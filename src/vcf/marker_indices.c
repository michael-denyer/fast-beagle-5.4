/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/MarkerIndices.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/marker_indices.h"

#include <stdlib.h>

#include "blbutil/utilities.h"
#include "jcompat/jarrays.h"

/* MarkerIndices.targIndex: the first target marker on or after the specified marker. */
static int targ_index(const int *targ_to_marker, int n, int marker) {
    int ins_pt = jarrays_search_int(targ_to_marker, 0, n, marker);
    return ins_pt < 0 ? -ins_pt - 1 : ins_pt;
}

static void check_bounds(int overlap_end, int overlap_start, int n_markers) {
    if (overlap_end < 0 || overlap_end > n_markers) util_exit("java.lang.IndexOutOfBoundsException: %d", overlap_end);
    if (overlap_start < 0 || overlap_start > n_markers) util_exit("java.lang.IndexOutOfBoundsException: %d", overlap_start);
}

void marker_indices_init(marker_indices *ix, const bool *in_targ, int n_markers, int overlap_end, int overlap_start) {
    check_bounds(overlap_end, overlap_start, n_markers);
    ix->n_markers = n_markers;
    ix->prev_splice = overlap_end >> 1;
    ix->overlap_end = overlap_end;
    ix->overlap_start = overlap_start;
    ix->next_splice = (int)((unsigned)(n_markers + overlap_start) >> 1);
    int n_targ = 0;
    for (int j = 0; j < n_markers; ++j) n_targ += in_targ[j];
    ix->n_targ_markers = n_targ;
    ix->targ_marker_to_marker = util_malloc((size_t)(n_targ > 0 ? n_targ : 1) * sizeof(int));
    ix->marker_to_targ_marker = util_malloc((size_t)(n_markers > 0 ? n_markers : 1) * sizeof(int));
    for (int j = 0, t = 0; j < n_markers; ++j) {
        ix->marker_to_targ_marker[j] = in_targ[j] ? t : -1;
        if (in_targ[j]) ix->targ_marker_to_marker[t++] = j;
    }
    ix->targ_prev_splice = targ_index(ix->targ_marker_to_marker, n_targ, ix->prev_splice);
    ix->targ_overlap_end = targ_index(ix->targ_marker_to_marker, n_targ, overlap_end);
    ix->targ_overlap_start = targ_index(ix->targ_marker_to_marker, n_targ, overlap_start);
    ix->targ_next_splice = targ_index(ix->targ_marker_to_marker, n_targ, ix->next_splice);
}

void marker_indices_init_targ(marker_indices *ix, int overlap_end, int overlap_start, int n_markers) {
    if (n_markers < 0) util_exit("java.lang.IllegalArgumentException: %d", n_markers);
    check_bounds(overlap_end, overlap_start, n_markers);
    ix->n_markers = ix->n_targ_markers = n_markers;
    ix->prev_splice = ix->targ_prev_splice = overlap_end >> 1;
    ix->overlap_end = ix->targ_overlap_end = overlap_end;
    ix->overlap_start = ix->targ_overlap_start = overlap_start;
    ix->next_splice = ix->targ_next_splice = (int)((unsigned)(n_markers + overlap_start) >> 1);
    ix->targ_marker_to_marker = util_malloc((size_t)(n_markers > 0 ? n_markers : 1) * sizeof(int));
    ix->marker_to_targ_marker = util_malloc((size_t)(n_markers > 0 ? n_markers : 1) * sizeof(int));
    for (int j = 0; j < n_markers; ++j) ix->targ_marker_to_marker[j] = ix->marker_to_targ_marker[j] = j;
}

void marker_indices_free(marker_indices *ix) {
    free(ix->targ_marker_to_marker);
    free(ix->marker_to_targ_marker);
}
