/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/BasicPhaseStates.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_BASIC_PHASE_STATES_H
#define PHASE_BASIC_PHASE_STATES_H

#include <stdint.h>

#include "beagleutil/comp_hap_queue.h"
#include "phase/coded_steps.h"
#include "phase/marker_cluster.h"
#include "phase/pbwt_phase_ibs.h"

/* The HMM states for one target sample: up to max_states composite reference
 * haplotypes, each a sequence of segments of the IBS neighbours found for the
 * sample's two haplotypes. A segment is replaced when the queue is full or its
 * haplotype has not been an IBS neighbour for min_steps steps, starting half
 * way between the two. One instance per thread; each call rewrites it. */
typedef struct {
    const pbwt_phase_ibs *ibs;   /* ibs->cs holds the current and reference haplotypes */
    int n_markers;
    int max_states;
    int min_steps;
    comp_hap_tracker t;
    uint64_t **comp_haps;       /* max_states composite haplotypes */
    uint64_t *column;           /* one word of every composite haplotype */
} basic_phase_states;

void basic_phase_states_init(basic_phase_states *bps, const pbwt_phase_ibs *ibs, int max_states);
/* BasicPhaseStates.ibsStates(sample, nMismatches): fills mismatch[0][m][j]
 * and mismatch[1][m][j] with whether composite haplotype j differs from the
 * sample's first and second haplotype at stage-1 marker m. Returns the number
 * of states. */
int basic_phase_states_ibs_states(basic_phase_states *bps, int sample, uint8_t ***mismatch);
/* BasicPhaseStates.ibsStates(mc, refAtMissingGT, nMismatches): for each of
 * the sample's clusters c, mismatch[0][c][j] says whether composite haplotype
 * j differs from a homozygous cluster, and mismatch[1] and mismatch[2] whether
 * it differs from the sample's first and second haplotype there. For the k-th
 * missing-genotype cluster, ref_at_missing[k][j] is haplotype j's allele.
 * Java fills three rows per cluster. Here mismatch[i][c] points at the row:
 * rows that are always equal share one row in rows, which must hold two rows
 * of max_states per cluster, and rows that are always zero are zero_row. */
int basic_phase_states_cluster_states(basic_phase_states *bps, const marker_cluster *mc, int **ref_at_missing,
        const uint8_t **mismatch[3], uint8_t *rows, const uint8_t *zero_row);
void basic_phase_states_free(basic_phase_states *bps);

#endif
