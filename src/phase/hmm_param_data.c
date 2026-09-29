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
#include "phase/hmm_param_data.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"
#include "phase/hmm_updater.h"

void hmm_param_data_init(hmm_param_data *hpd, const pbwt_phase_ibs *ibs) {
    const phase_data *pd = ibs->pd;
    int max_states = pd->par->phase_states;
    int n = pd->fpd->n_stage1;
    hpd->pd = pd;
    hpd->n_markers = n;
    basic_phase_states_init(&hpd->states, ibs, max_states);
    hpd->al_match = util_malloc(2 * sizeof *hpd->al_match);
    for (int i = 0; i < 2; ++i) {
        hpd->al_match[i] = util_malloc((size_t)n * sizeof **hpd->al_match);
        for (int m = 0; m < n; ++m) hpd->al_match[i][m] = util_malloc((size_t)max_states);
    }
    hpd->fwd = util_malloc((size_t)max_states * sizeof *hpd->fwd);
    hpd->bwd = util_malloc((size_t)max_states * sizeof *hpd->bwd);
    hpd->saved_bwd = util_malloc((size_t)n * sizeof *hpd->saved_bwd);
    for (int m = 0; m < n; ++m) hpd->saved_bwd[m] = util_malloc((size_t)max_states * sizeof **hpd->saved_bwd);
    hpd->em_probs[0] = 1.0f - pd->p_mismatch;
    hpd->em_probs[1] = pd->p_mismatch;
    hpd->mismatch_cnt = 0;
    hpd->sum_mismatch_prob = 0.0;
    hpd->sum_gen_dist = 0.0;
    hpd->sum_switch_prob = 0.0;
}

void hmm_param_data_free(hmm_param_data *hpd) {
    for (int i = 0; i < 2; ++i) {
        for (int m = 0; m < hpd->n_markers; ++m) free(hpd->al_match[i][m]);
        free(hpd->al_match[i]);
    }
    free(hpd->al_match);
    for (int m = 0; m < hpd->n_markers; ++m) free(hpd->saved_bwd[m]);
    free(hpd->saved_bwd);
    free(hpd->fwd);
    free(hpd->bwd);
    basic_phase_states_free(&hpd->states);
}

void hmm_param_data_add_estimation_data(hmm_param_data *hpd, param_estimates *pe) {
    param_estimates_add_mismatch_data(pe, hpd->mismatch_cnt, hpd->sum_mismatch_prob);
    param_estimates_add_switch_data(pe, hpd->sum_gen_dist, hpd->sum_switch_prob);
    hpd->mismatch_cnt = 0;
    hpd->sum_mismatch_prob = 0.0;
    hpd->sum_gen_dist = 0.0;
    hpd->sum_switch_prob = 0.0;
}

static float fwd_update(hmm_param_data *hpd, int m, const uint8_t *al_discord, int n_states, float last_sum, float h_factor) {
    float p_switch = hpd->pd->p_recomb[m];
    float shift = p_switch / n_states;
    float scale = (1.0f - p_switch) / last_sum;
    float no_switch_scale = ((1.0f - p_switch) + shift) / last_sum;
    float joint_state_sum = 0.0f;
    float state_sum = 0.0f;
    const float *bwd_m = hpd->saved_bwd[m];
    float *fwd = hpd->fwd;
    float fwd_sum = 0.0f;
    float mismatch_sum = 0.0f;
    for (int k = 0; k < n_states; ++k) {
        float em = hpd->em_probs[al_discord[k]];
        joint_state_sum += bwd_m[k] * em * no_switch_scale * fwd[k];
        fwd[k] = em * (scale * fwd[k] + shift);
        fwd_sum += fwd[k];
        float state_prob = fwd[k] * bwd_m[k];
        state_sum += state_prob;
        if (al_discord[k] > 0) mismatch_sum += state_prob;
    }
    ++hpd->mismatch_cnt;
    hpd->sum_mismatch_prob += mismatch_sum / state_sum;
    double switch_prob = h_factor * (1.0f - joint_state_sum / state_sum);
    if (switch_prob > 0.0) {
        hpd->sum_gen_dist += hpd->pd->fpd->stage1_map.gen_dist[m];
        hpd->sum_switch_prob += switch_prob;
    }
    return fwd_sum;
}

static void get_param_data(hmm_param_data *hpd, uint8_t **al_match, int n_states) {
    int n = hpd->n_markers;
    const float *p_recomb = hpd->pd->p_recomb;
    for (int k = 0; k < n_states; ++k) hpd->bwd[k] = 1.0f;
    for (int k = 0; k < n_states; ++k) hpd->saved_bwd[n - 1][k] = 1.0f;
    for (int m = n - 2; m >= 0; --m) {
        hmm_bwd_update(hpd->bwd, p_recomb[m + 1], hpd->em_probs, al_match[m + 1], n_states);
        memcpy(hpd->saved_bwd[m], hpd->bwd, (size_t)n_states * sizeof *hpd->bwd);
    }
    float h_factor = n_states / (n_states - 1.0f);
    for (int k = 0; k < n_states; ++k) hpd->fwd[k] = 1.0f / n_states;
    float sum = 1.0f;
    for (int m = 0; m < n; ++m) sum = fwd_update(hpd, m, al_match[m], n_states, sum, h_factor);
}

void hmm_param_data_update(hmm_param_data *hpd, int sample) {
    int n_states = basic_phase_states_ibs_states(&hpd->states, sample, hpd->al_match);
    if (n_states > 1) {
        get_param_data(hpd, hpd->al_match[0], n_states);
        get_param_data(hpd, hpd->al_match[1], n_states);
    }
}
