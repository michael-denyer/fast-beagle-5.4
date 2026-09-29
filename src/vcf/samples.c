/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/Samples.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/samples.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"

static int compare_ids(const void *a, const void *b) {
    return strcmp(*(char *const *)a, *(char *const *)b);
}

/* Samples.checkForNullsAndDuplicates. strcmp orders ASCII as String.compareTo does. */
void samples_init(samples *s, int n, char **ids, bool *is_diploid) {
    char **sorted = util_malloc((size_t)n * sizeof *sorted);
    memcpy(sorted, ids, (size_t)n * sizeof *sorted);
    qsort(sorted, (size_t)n, sizeof *sorted, compare_ids);
    for (int j = 0; j < n; ++j) {
        if (sorted[j][0] == '\0') util_exit("java.lang.IllegalArgumentException: Empty string identifier");
        if (j > 0 && strcmp(sorted[j], sorted[j - 1]) == 0) {
            fprintf(stderr, "Warning: duplicate sample identifier: %s\n", sorted[j]);
        }
    }
    free(sorted);
    s->n = n;
    s->ids = ids;
    s->is_diploid = is_diploid;
}

void samples_free(samples *s) {
    for (int j = 0; j < s->n; ++j) free(s->ids[j]);
    free(s->ids);
    free(s->is_diploid);
}
