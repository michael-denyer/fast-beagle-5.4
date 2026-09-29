/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/HmmUpdater.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/hmm_updater.h"

float hmm_fwd_update(const float *prev, float *fwd, float fwd_sum, float p_switch, const float p_mismatch[2], const uint8_t *mismatch, int n_states) {
    float shift = p_switch / n_states;
    float scale = (1.0f - p_switch) / fwd_sum;
    fwd_sum = 0.0f;
    for (int k = 0; k < n_states; ++k) {
        fwd[k] = p_mismatch[mismatch[k]] * (scale * prev[k] + shift);
        fwd_sum += fwd[k];
    }
    return fwd_sum;
}

void hmm_bwd_update(float *bwd, float p_switch, const float p_mismatch[2], const uint8_t *mismatch, int n_states) {
    float sum = 0.0f;
    for (int k = 0; k < n_states; ++k) {
        bwd[k] *= p_mismatch[mismatch[k]];
        sum += bwd[k];
    }
    float shift = p_switch / n_states;
    float scale = (1.0f - p_switch) / sum;
    for (int k = 0; k < n_states; ++k) bwd[k] = scale * bwd[k] + shift;
}

void hmm_fwd_update3(float *fwd[3], float fwd_sums[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states) {
    float *restrict f0 = fwd[0], *restrict f1 = fwd[1], *restrict f2 = fwd[2];
    float shift = p_switch / n_states;
    float scale0 = (1.0f - p_switch) / fwd_sums[0];
    float scale1 = (1.0f - p_switch) / fwd_sums[1];
    float scale2 = (1.0f - p_switch) / fwd_sums[2];
    float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f;
    for (int k = 0; k < n_states; ++k) {
        f0[k] = p_mismatch[m0[k]] * (scale0 * f0[k] + shift);
        f1[k] = p_mismatch[m1[k]] * (scale1 * f1[k] + shift);
        f2[k] = p_mismatch[m2[k]] * (scale2 * f2[k] + shift);
        s0 += f0[k];
        s1 += f1[k];
        s2 += f2[k];
    }
    fwd_sums[0] = s0;
    fwd_sums[1] = s1;
    fwd_sums[2] = s2;
}

void hmm_bwd_update3(float *bwd[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states) {
    float *restrict b0 = bwd[0], *restrict b1 = bwd[1], *restrict b2 = bwd[2];
    float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f;
    for (int k = 0; k < n_states; ++k) {
        b0[k] *= p_mismatch[m0[k]];
        b1[k] *= p_mismatch[m1[k]];
        b2[k] *= p_mismatch[m2[k]];
        s0 += b0[k];
        s1 += b1[k];
        s2 += b2[k];
    }
    float shift = p_switch / n_states;
    float scale0 = (1.0f - p_switch) / s0;
    float scale1 = (1.0f - p_switch) / s1;
    float scale2 = (1.0f - p_switch) / s2;
    for (int k = 0; k < n_states; ++k) {
        b0[k] = scale0 * b0[k] + shift;
        b1[k] = scale1 * b1[k] + shift;
        b2[k] = scale2 * b2[k] + shift;
    }
}
