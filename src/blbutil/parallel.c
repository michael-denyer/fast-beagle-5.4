/*
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "blbutil/parallel.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdlib.h>

#include "blbutil/utilities.h"

typedef struct {
    atomic_int next;
    int n_items;
    parallel_fn fn;
} shared;

typedef struct {
    shared *sh;
    void *worker;
} thread_arg;

static void *run(void *p) {
    thread_arg *a = p;
    for (int item = atomic_fetch_add(&a->sh->next, 1); item < a->sh->n_items; item = atomic_fetch_add(&a->sh->next, 1)) {
        a->sh->fn(a->worker, item);
    }
    return NULL;
}

static void *worker_at(void *workers, size_t worker_size, int t) {
    return (char *)workers + (size_t)t * worker_size;
}

void parallel_for(int n_threads, int n_items, void *workers, size_t worker_size, parallel_fn fn) {
    if (n_threads <= 1 || n_items <= 1) {
        for (int item = 0; item < n_items; ++item) fn(workers, item);
        return;
    }
    shared sh;
    atomic_init(&sh.next, 0);
    sh.n_items = n_items;
    sh.fn = fn;
    pthread_t *threads = util_malloc((size_t)n_threads * sizeof *threads);
    thread_arg *args = util_malloc((size_t)n_threads * sizeof *args);
    for (int t = 0; t < n_threads; ++t) {
        args[t] = (thread_arg){&sh, worker_at(workers, worker_size, t)};
        if (t > 0 && pthread_create(&threads[t], NULL, run, &args[t]) != 0) util_exit(PROGRAM ": cannot create thread");
    }
    run(&args[0]);
    for (int t = 1; t < n_threads; ++t) pthread_join(threads[t], NULL);
    free(threads);
    free(args);
}

typedef struct {
    pthread_mutex_t mu;
    pthread_cond_t built_cv;  /* an item finished building */
    pthread_cond_t room_cv;   /* an item was consumed, or every item is claimed */
    int n_items;
    int window;
    int next;       /* the next item to claim */
    int consumed;   /* items consumed, in item order */
    bool *built;    /* per slot: item % window is built and not yet consumed */
    parallel_fn build;
} ordered;

typedef struct {
    ordered *o;
    void *worker;
} ordered_arg;

static void *run_ordered(void *p) {
    ordered_arg *a = p;
    ordered *o = a->o;
    pthread_mutex_lock(&o->mu);
    for (;;) {
        while (o->next < o->n_items && o->next >= o->consumed + o->window) pthread_cond_wait(&o->room_cv, &o->mu);
        if (o->next >= o->n_items) break;
        int item = o->next++;
        pthread_mutex_unlock(&o->mu);
        o->build(a->worker, item);
        pthread_mutex_lock(&o->mu);
        o->built[item % o->window] = true;
        if (item == o->consumed) pthread_cond_signal(&o->built_cv);
    }
    pthread_mutex_unlock(&o->mu);
    return NULL;
}

void parallel_ordered(int n_threads, int n_items, int window, void *workers, size_t worker_size, parallel_fn build,
        void *ctx, parallel_fn consume) {
    ordered o = {.n_items = n_items, .window = window, .next = 0, .consumed = 0, .build = build};
    if (pthread_mutex_init(&o.mu, NULL) != 0 || pthread_cond_init(&o.built_cv, NULL) != 0
            || pthread_cond_init(&o.room_cv, NULL) != 0) {
        util_exit(PROGRAM ": cannot create thread");
    }
    o.built = util_malloc((size_t)window * sizeof *o.built);
    for (int s = 0; s < window; ++s) o.built[s] = false;
    pthread_t *threads = util_malloc((size_t)n_threads * sizeof *threads);
    ordered_arg *args = util_malloc((size_t)n_threads * sizeof *args);
    for (int t = 0; t < n_threads; ++t) {
        args[t] = (ordered_arg){&o, worker_at(workers, worker_size, t)};
        if (pthread_create(&threads[t], NULL, run_ordered, &args[t]) != 0) util_exit(PROGRAM ": cannot create thread");
    }
    for (int item = 0; item < n_items; ++item) {
        pthread_mutex_lock(&o.mu);
        while (!o.built[item % window]) pthread_cond_wait(&o.built_cv, &o.mu);
        pthread_mutex_unlock(&o.mu);
        consume(ctx, item);
        pthread_mutex_lock(&o.mu);
        o.built[item % window] = false;
        o.consumed = item + 1;
        pthread_cond_broadcast(&o.room_cv);
        pthread_mutex_unlock(&o.mu);
    }
    for (int t = 0; t < n_threads; ++t) pthread_join(threads[t], NULL);
    free(threads);
    free(args);
    free(o.built);
    pthread_cond_destroy(&o.built_cv);
    pthread_cond_destroy(&o.room_cv);
    pthread_mutex_destroy(&o.mu);
}
