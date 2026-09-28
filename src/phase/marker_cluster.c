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
#include "phase/marker_cluster.h"

#include <limits.h>
#include <stdlib.h>

#include "blbutil/utilities.h"

/* MarkerCluster.unphHetClusters: the clusters holding an unphased marker. */
static void set_unph_clusters(marker_cluster *mc, const sample_phase *sp) {
    int n_unph = sp->n_unphased;
    mc->unph_clusters = util_malloc((size_t)(n_unph > 0 ? n_unph : 1) * sizeof *mc->unph_clusters);
    mc->n_unph = 0;
    int unph_index = 0;
    int next_unph = unph_index < n_unph ? sp->unphased[unph_index++] : INT_MAX;
    for (int j = 0; j < mc->n_clusters; ++j) {
        int clust_end = mc->ends[j];
        if (next_unph < clust_end) {
            mc->unph_clusters[mc->n_unph++] = j;
            next_unph = unph_index < n_unph ? sp->unphased[unph_index++] : INT_MAX;
            while (next_unph < clust_end) next_unph = unph_index < n_unph ? sp->unphased[unph_index++] : INT_MAX;
        }
    }
}

/* MarkerCluster.setClustHasMissingGT */
static void set_has_missing(marker_cluster *mc, const sample_phase *sp) {
    mc->has_missing = util_malloc((size_t)mc->n_clusters * sizeof *mc->has_missing);
    mc->n_missing_clusters = 0;
    int i = 0;
    for (int c = 0; c < mc->n_clusters; ++c) {
        int end = mc->ends[c];
        mc->has_missing[c] = false;
        for (; i < sp->n_missing && sp->missing[i] < end; ++i) mc->has_missing[c] = true;
        if (mc->has_missing[c]) ++mc->n_missing_clusters;
    }
}

void marker_cluster_init(marker_cluster *mc, const phase_data *pd, int sample) {
    sample_phase *sp = &pd->phase[sample];
    int n = sp->n_clusters;
    mc->sp = sp;
    mc->n_clusters = n;
    mc->ends = util_malloc((size_t)n * sizeof *mc->ends);
    sample_phase_clust_ends(sp, mc->ends);
    set_unph_clusters(mc, sp);
    set_has_missing(mc, sp);
    mc->p_recomb = util_malloc((size_t)n * sizeof *mc->p_recomb);
    mc->p_recomb[0] = 0.0f;
    /* Beagle 5.4 starts the first interval at marker 0, not at the first
     * cluster's end. */
    int start = 0;
    for (int c = 1; c < n; ++c) {
        int end = mc->ends[c];
        float p_no_recomb = 1.0f;
        for (int k = start; k < end; ++k) p_no_recomb *= (1.0f - pd->p_recomb[k]);
        mc->p_recomb[c] = 1.0f - p_no_recomb;
        start = end;
    }
}

void marker_cluster_free(marker_cluster *mc) {
    free(mc->ends);
    free(mc->unph_clusters);
    free(mc->has_missing);
    free(mc->p_recomb);
}
