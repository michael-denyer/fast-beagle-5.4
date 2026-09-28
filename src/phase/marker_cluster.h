/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/MarkerCluster.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_MARKER_CLUSTER_H
#define PHASE_MARKER_CLUSTER_H

#include <stdbool.h>

#include "phase/phase_data.h"
#include "phase/sample_phase.h"

/* A sample's clusters as HMM positions: each cluster's end marker, the
 * clusters holding an unphased heterozygote or a missing genotype, and the
 * probability of a switch since the previous cluster's end (0 for the first
 * cluster; the second also covers the first). */
typedef struct {
    sample_phase *sp;
    int n_clusters;
    int *ends;
    int n_unph;
    int *unph_clusters;
    bool *has_missing;
    int n_missing_clusters;
    float *p_recomb;
} marker_cluster;

/* new MarkerCluster(phaseData, sample) */
void marker_cluster_init(marker_cluster *mc, const phase_data *pd, int sample);
static inline int marker_cluster_start(const marker_cluster *mc, int c) {
    return c == 0 ? 0 : mc->ends[c - 1];
}
void marker_cluster_free(marker_cluster *mc);

#endif
