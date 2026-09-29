/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/Utilities.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "blbutil/utilities.h"

#include <stdarg.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static FILE *exit_log;
static atomic_flag exiting = ATOMIC_FLAG_INIT;

void util_exit_log(FILE *log) {
    exit_log = log;
}

void util_exit(const char *fmt, ...) {
    if (atomic_flag_test_and_set(&exiting)) {
        for (;;) pause();
    }
    fflush(stdout);   /* progress lines first, as Java's autoflushed System.out */
    va_list ap;
    va_start(ap, fmt);
    va_list ap2;
    va_copy(ap2, ap);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    if (exit_log != NULL) {
        vfprintf(exit_log, fmt, ap2);
        fputc('\n', exit_log);
    }
    va_end(ap2);
    va_end(ap);
    exit(1);
}

void *util_malloc(size_t size) {
    void *p = malloc(size == 0 ? 1 : size);
    if (p == NULL) util_exit("ERROR: out of memory");
    return p;
}

void *util_realloc(void *p, size_t size) {
    p = realloc(p, size == 0 ? 1 : size);
    if (p == NULL) util_exit("ERROR: out of memory");
    return p;
}

char *util_strndup(const char *s, size_t n) {
    char *d = util_malloc(n + 1);
    memcpy(d, s, n);
    d[n] = '\0';
    return d;
}

int util_compare_ints(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

void util_shuffle(int *a, int n, int n_elements, jrandom *r) {
    for (int j = 0; j < n_elements; ++j) {
        int x = jrandom_next_int_bound(r, n - j);
        int tmp = a[j];
        a[j] = a[j + x];
        a[j + x] = tmp;
    }
}
