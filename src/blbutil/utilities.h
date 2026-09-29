/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/Utilities.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BLBUTIL_UTILITIES_H
#define BLBUTIL_UTILITIES_H

#include <stddef.h>
#include <stdio.h>

#include "jcompat/jrandom.h"

/* Main.PROGRAM: this program's name, which also starts the messages of the
 * errors that only fast-beagle raises. */
#define PROGRAM "fast-beagle"

/* Utilities.exit: prints the message to standard error and exits with status 1.
 * Only the first caller reports; a later or concurrent caller blocks until the
 * process exits. */
_Noreturn void util_exit(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
/* fast-beagle only: util_exit also writes its message to log, or stops when log
 * is NULL. Beagle's log never records why a run failed. */
void util_exit_log(FILE *log);

void *util_malloc(size_t size);
void *util_realloc(void *p, size_t size);
char *util_strndup(const char *s, size_t n);

/* A qsort comparator for ascending ints. */
int util_compare_ints(const void *a, const void *b);

/* Utilities.shuffle(ia, nElements, random): moves a random n_elements of the
 * n entries of a to its front. */
void util_shuffle(int *a, int n, int n_elements, jrandom *r);

#endif
