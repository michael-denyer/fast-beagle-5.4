/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PhaseBaum1.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PHASE_BAUM1_H
#define PHASE_PHASE_BAUM1_H

#include <stdbool.h>
#include <stdint.h>

#include "phase/basic_phase_states.h"

/* PhaseBaum1's swap counters: switches of phase between successive unphased
 * heterozygotes, as a proportion of unphased heterozygotes, over one
 * iteration. */
typedef struct {
    int64_t n_swaps;
    int64_t n_unph_hets;
} swap_rate;

/* One sample at a time, runs the backward then the forward algorithm over the
 * sample's clusters with three HMMs (the homozygous clusters, and each
 * haplotype since the last unphased heterozygote), phases each unphased
 * heterozygote from the posterior odds, and imputes missing genotypes. After
 * burn-in, a heterozygote whose likelihood ratio reaches the sample's quantile
 * threshold is committed. The model parameters are fixed when it is created. */
typedef struct {
    phase_data *pd;
    bool burnin;
    bool trace;
    int max_states;
    basic_phase_states states;
    int n_states;
    uint8_t ***mismatch;   /* [3][cluster][state]; rows 1 and 2 swap with the haplotypes */
    float p_mismatch;
    float em_probs[2];
    float *fwd[3];
    float *bwd[3];
    float sum[3];
    int n_miss_cap;
    int **ref_alleles;
    float **miss_probs1;
    float **miss_probs2;
    int n_het_cap;
    float **bwd_het1;
    float **bwd_het2;
    float *lr;             /* each unphased heterozygote's likelihood ratio */
    float *lr_sorted;
    int n_lr;
    bool swap_haps;
    int miss_index;
    int swap_cnt;
    uint64_t em_digest;    /* trace seam T3e */
    float threshold;       /* trace seam T3e: NaN when none applies */
    char **trace_lines;    /* trace seam T3e: one line per sample, or NULL */
} phase_baum1;

/* new PhaseBaum1(phaseIbs). With trace_lines, each sample's T3e line is
 * stored there for the caller to write in sample order. */
void phase_baum1_init(phase_baum1 *pb, const pbwt_phase_ibs *ibs, char **trace_lines);
/* PhaseBaum1.phase(sample): updates the sample's SamplePhase and adds its
 * swaps to rate. */
void phase_baum1_phase(phase_baum1 *pb, int sample, swap_rate *rate);
void phase_baum1_free(phase_baum1 *pb);

#endif
