/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/CompHapSegment.java, with
 * the priority queue operations Beagle uses on it, and the queue
 * state and updateHeadOfQ that phase/BasicPhaseStates.java,
 * phase/LowFreqPhaseStates.java and imp/ImpStates.java each repeat;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "beagleutil/comp_hap_queue.h"

#include <stdlib.h>

#include "blbutil/utilities.h"

/* CompHapSegment.compareTo */
static int compare(const comp_hap_segment *a, const comp_hap_segment *b) {
    if (a->last_ibs_step != b->last_ibs_step) return a->last_ibs_step < b->last_ibs_step ? -1 : 1;
    return 0;
}

void comp_hap_queue_init(comp_hap_queue *q, int capacity) {
    q->cap = capacity > 0 ? capacity : 1;
    q->es = util_malloc((size_t)q->cap * sizeof *q->es);
    q->size = 0;
}

void comp_hap_queue_free(comp_hap_queue *q) {
    free(q->es);
}

/* Binary min-heap, root at 0, children of i at 2i+1 and 2i+2. The tie rules
 * (stop at an equal parent; prefer the left child on equal children) fix the
 * order in which equal keys leave, which Beagle's output depends on. */
void comp_hap_queue_offer(comp_hap_queue *q, comp_hap_segment *key) {
    if (q->size == q->cap) {
        q->cap *= 2;
        q->es = util_realloc(q->es, (size_t)q->cap * sizeof *q->es);
    }
    int k = q->size++;
    while (k > 0) {
        int parent = (k - 1) >> 1;
        if (compare(key, q->es[parent]) >= 0) break;
        q->es[k] = q->es[parent];
        k = parent;
    }
    q->es[k] = key;
}

comp_hap_segment *comp_hap_queue_poll(comp_hap_queue *q) {
    if (q->size == 0) return NULL;
    comp_hap_segment *root = q->es[0];
    int n = --q->size;
    if (n > 0) {
        comp_hap_segment *key = q->es[n];
        int k = 0;
        for (int child = 1; child < n; child = 2 * k + 1) {
            if (child + 1 < n && compare(q->es[child + 1], q->es[child]) < 0) ++child;
            if (compare(key, q->es[child]) <= 0) break;
            q->es[k] = q->es[child];
            k = child;
        }
        q->es[k] = key;
    }
    return root;
}

static const int NIL = -103;

void comp_hap_tracker_init(comp_hap_tracker *t, int max_states) {
    comp_hap_queue_init(&t->q, max_states);
    t->segs = util_malloc((size_t)max_states * sizeof *t->segs);
    t->last_ibs_step = int_int_map_new();
    t->max_states = max_states;
}

void comp_hap_tracker_free(comp_hap_tracker *t) {
    int_int_map_free(t->last_ibs_step);
    free(t->segs);
    comp_hap_queue_free(&t->q);
}

void comp_hap_tracker_clear(comp_hap_tracker *t) {
    comp_hap_queue_clear(&t->q);
    int_int_map_clear(t->last_ibs_step);
}

/* Re-sorts segments whose haplotype was an IBS neighbour again after they
 * were queued, until the head's step is current. */
static void update_head_of_q(comp_hap_tracker *t) {
    comp_hap_segment *head = comp_hap_queue_peek(&t->q);
    if (head == NULL) return;
    int last_ibs_step = int_int_map_get(t->last_ibs_step, head->hap, NIL);
    while (head->last_ibs_step != last_ibs_step) {
        head = comp_hap_queue_poll(&t->q);
        head->last_ibs_step = last_ibs_step;
        comp_hap_queue_offer(&t->q, head);
        head = comp_hap_queue_peek(&t->q);
        last_ibs_step = int_int_map_get(t->last_ibs_step, head->hap, NIL);
    }
}

static int add_segment(comp_hap_tracker *t, int hap, int step) {
    int index = t->q.size;
    t->segs[index] = (comp_hap_segment){hap, 0, step, index};
    comp_hap_queue_offer(&t->q, &t->segs[index]);
    return index;
}

comp_hap_change comp_hap_tracker_observe(comp_hap_tracker *t, int hap, int step, int min_steps) {
    comp_hap_change change = {-1, -1, 0, 0};
    if (int_int_map_get(t->last_ibs_step, hap, NIL) == NIL) {
        update_head_of_q(t);
        comp_hap_segment *head = comp_hap_queue_peek(&t->q);
        bool full = t->q.size == t->max_states;
        bool stale = head != NULL && step - head->last_ibs_step >= min_steps;
        if (full || stale) {
            head = comp_hap_queue_poll(&t->q);
            int mid_step = (int)((unsigned)(head->last_ibs_step + step) >> 1);
            change = (comp_hap_change){head->comp_hap_index, head->hap, head->start_step, mid_step};
            int_int_map_remove(t->last_ibs_step, head->hap);
            head->hap = hap;
            head->start_step = mid_step;
            head->last_ibs_step = step;
            comp_hap_queue_offer(&t->q, head);
        } else {
            change.index = add_segment(t, hap, step);
        }
    }
    int_int_map_put(t->last_ibs_step, hap, step);
    return change;
}

int comp_hap_tracker_seed(comp_hap_tracker *t, int hap) {
    return add_segment(t, hap, 0);
}

int comp_hap_tracker_size(const comp_hap_tracker *t) {
    return t->q.size;
}

const comp_hap_segment *comp_hap_tracker_segment(const comp_hap_tracker *t, int index) {
    return &t->segs[index];
}
