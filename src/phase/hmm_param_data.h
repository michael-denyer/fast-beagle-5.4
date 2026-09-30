/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/HmmParamData.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_HMM_PARAM_DATA_H
#define PHASE_HMM_PARAM_DATA_H

#include <stdint.h>

#include "phase/basic_phase_states.h"
#include "phase/param_estimates.h"

/* Runs the forward-backward HMM of each of a sample's haplotypes against its
 * composite reference haplotypes and sums the posterior mismatch and switch
 * probabilities that the parameter estimates are made from. The two
 * haplotypes' HMMs step through the markers together. */
typedef struct {
    const phase_data *pd;
    int n_markers;
    int max_states;
    basic_phase_states states;
    uint8_t ***al_match;       /* [2][marker][state] */
    float *fwd[2];             /* per haplotype */
    float *bwd[2];
    float *block_end_bwd[2];   /* [block][state]: the backward values at each block's last marker */
    float *block_bwd[2];       /* [marker in block][state]: the block the forward pass is in */
    float *joint[2];           /* one marker's per-state terms, before they are summed */
    float *state[2];
    float *discord[2];
    float *hap2_mismatch;      /* [marker]: the second haplotype's terms, added after all of the first's */
    double *hap2_switch;
    float em_probs[2];
    int mismatch_cnt;
    double sum_mismatch_prob;
    double sum_gen_dist;
    double sum_switch_prob;
} hmm_param_data;

/* new HmmParamData(phaseIbs): the emission probabilities use the current
 * mismatch probability. */
void hmm_param_data_init(hmm_param_data *hpd, const pbwt_phase_ibs *ibs);
void hmm_param_data_update(hmm_param_data *hpd, int sample);
/* Moves the sums into pe and resets them. */
void hmm_param_data_add_estimation_data(hmm_param_data *hpd, param_estimates *pe);
void hmm_param_data_free(hmm_param_data *hpd);

#endif
