/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/MarkerIndices.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_MARKER_INDICES_H
#define VCF_MARKER_INDICES_H

#include <stdbool.h>

/* A window's overlap and splice points, for all markers and for target
 * markers, and the maps between the two indexings. */
typedef struct {
    int n_markers, n_targ_markers;
    int prev_splice, overlap_end, overlap_start, next_splice;
    int targ_prev_splice, targ_overlap_end, targ_overlap_start, targ_next_splice;
    int *targ_marker_to_marker;   /* n_targ_markers entries */
    int *marker_to_targ_marker;   /* n_markers entries, -1 where not a target marker */
} marker_indices;

/* new MarkerIndices(inTarg, overlapEnd, overlapStart) */
void marker_indices_init(marker_indices *ix, const bool *in_targ, int n_markers, int overlap_end, int overlap_start);
/* new MarkerIndices(overlapEnd, overlapStart, nMarkers): every marker is a target marker. */
void marker_indices_init_targ(marker_indices *ix, int overlap_end, int overlap_start, int n_markers);
void marker_indices_free(marker_indices *ix);

#endif
