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

/* The steps below compute each state's value in loops with no sum, which the
 * compiler can vectorise, and then add the values in state order. `omp simd`
 * states that a value loop's iterations are independent, so GCC at -O2
 * vectorises it without a runtime alias check. A mismatch byte is 0 or 1, so
 * the select equals p_mismatch[mismatch[k]]. */

static void sum2(const float *restrict a, int na, const float *restrict b, int nb, float sum[2]) {
    int n = na < nb ? na : nb;
    float s0 = 0.0f, s1 = 0.0f;
    for (int k = 0; k < n; ++k) {
        s0 += a[k];
        s1 += b[k];
    }
    for (int k = n; k < na; ++k) s0 += a[k];
    for (int k = n; k < nb; ++k) s1 += b[k];
    sum[0] = s0;
    sum[1] = s1;
}

void hmm_fwd_update2(const float *const prev[2], float *const fwd[2], float fwd_sum[2], float p_switch, const float p_mismatch[2],
        const uint8_t *const mismatch[2], const int n_states[2]) {
    float match = p_mismatch[0], differ = p_mismatch[1];
    for (int i = 0; i < 2; ++i) {
        const float *restrict p = prev[i];
        float *restrict f = fwd[i];
        const uint8_t *restrict d = mismatch[i];
        float shift = p_switch / n_states[i];
        float scale = (1.0f - p_switch) / fwd_sum[i];
        #pragma omp simd
        for (int k = 0; k < n_states[i]; ++k) f[k] = (d[k] ? differ : match) * (scale * p[k] + shift);
    }
    sum2(fwd[0], n_states[0], fwd[1], n_states[1], fwd_sum);
}

void hmm_bwd_update2(float *const bwd[2], float p_switch, const float p_mismatch[2], const uint8_t *const mismatch[2], const int n_states[2]) {
    float match = p_mismatch[0], differ = p_mismatch[1];
    for (int i = 0; i < 2; ++i) {
        float *restrict b = bwd[i];
        const uint8_t *restrict d = mismatch[i];
        #pragma omp simd
        for (int k = 0; k < n_states[i]; ++k) b[k] *= d[k] ? differ : match;
    }
    float sum[2];
    sum2(bwd[0], n_states[0], bwd[1], n_states[1], sum);
    for (int i = 0; i < 2; ++i) {
        float *restrict b = bwd[i];
        float shift = p_switch / n_states[i];
        float scale = (1.0f - p_switch) / sum[i];
        #pragma omp simd
        for (int k = 0; k < n_states[i]; ++k) b[k] = scale * b[k] + shift;
    }
}

static inline void fwd3_values(float *fwd[3], const float fwd_sums[3], float p_switch, const float p_mismatch[2], const uint8_t *m0,
        const uint8_t *m1, const uint8_t *m2, int n_states) {
    float *restrict f0 = fwd[0], *restrict f1 = fwd[1], *restrict f2 = fwd[2];
    const uint8_t *restrict d0 = m0, *restrict d1 = m1, *restrict d2 = m2;
    float match = p_mismatch[0], differ = p_mismatch[1];
    float shift = p_switch / n_states;
    float scale0 = (1.0f - p_switch) / fwd_sums[0];
    float scale1 = (1.0f - p_switch) / fwd_sums[1];
    float scale2 = (1.0f - p_switch) / fwd_sums[2];
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) f0[k] = (d0[k] ? differ : match) * (scale0 * f0[k] + shift);
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) f1[k] = (d1[k] ? differ : match) * (scale1 * f1[k] + shift);
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) f2[k] = (d2[k] ? differ : match) * (scale2 * f2[k] + shift);
}

static inline void bwd3_emit(float *bwd[3], const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states) {
    float *restrict b0 = bwd[0], *restrict b1 = bwd[1], *restrict b2 = bwd[2];
    const uint8_t *restrict d0 = m0, *restrict d1 = m1, *restrict d2 = m2;
    float match = p_mismatch[0], differ = p_mismatch[1];
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) b0[k] *= d0[k] ? differ : match;
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) b1[k] *= d1[k] ? differ : match;
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) b2[k] *= d2[k] ? differ : match;
}

static inline void bwd3_transition(float *bwd[3], const float sums[3], float p_switch, int n_states) {
    float *restrict b0 = bwd[0], *restrict b1 = bwd[1], *restrict b2 = bwd[2];
    float shift = p_switch / n_states;
    float scale0 = (1.0f - p_switch) / sums[0];
    float scale1 = (1.0f - p_switch) / sums[1];
    float scale2 = (1.0f - p_switch) / sums[2];
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) {
        b0[k] = scale0 * b0[k] + shift;
        b1[k] = scale1 * b1[k] + shift;
        b2[k] = scale2 * b2[k] + shift;
    }
}

/* Adds each HMM's states in [from, to) onto its sum, in state order. */
static inline void sum3(float *const v[3], int from, int to, float sums[3]) {
    const float *restrict v0 = v[0], *restrict v1 = v[1], *restrict v2 = v[2];
    float s0 = sums[0], s1 = sums[1], s2 = sums[2];
    for (int k = from; k < to; ++k) {
        s0 += v0[k];
        s1 += v1[k];
        s2 += v2[k];
    }
    sums[0] = s0;
    sums[1] = s1;
    sums[2] = s2;
}

/* The sums of two samples' three HMMs: together over the states both have,
 * then the rest of the longer one. */
static inline void sum3x2(float *const a[3], int na, float *const b[3], int nb, float sa[3], float sb[3]) {
    const float *restrict a0 = a[0], *restrict a1 = a[1], *restrict a2 = a[2];
    const float *restrict b0 = b[0], *restrict b1 = b[1], *restrict b2 = b[2];
    int n = na < nb ? na : nb;
    float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, t0 = 0.0f, t1 = 0.0f, t2 = 0.0f;
    for (int k = 0; k < n; ++k) {
        s0 += a0[k];
        s1 += a1[k];
        s2 += a2[k];
        t0 += b0[k];
        t1 += b1[k];
        t2 += b2[k];
    }
    sa[0] = s0;
    sa[1] = s1;
    sa[2] = s2;
    sb[0] = t0;
    sb[1] = t1;
    sb[2] = t2;
    sum3(a, n, na, sa);
    sum3(b, n, nb, sb);
}

void hmm_fwd_update3(float *fwd[3], float fwd_sums[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states) {
    fwd3_values(fwd, fwd_sums, p_switch, p_mismatch, m0, m1, m2, n_states);
    fwd_sums[0] = fwd_sums[1] = fwd_sums[2] = 0.0f;
    sum3(fwd, 0, n_states, fwd_sums);
}

void hmm_bwd_update3(float *bwd[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states) {
    bwd3_emit(bwd, p_mismatch, m0, m1, m2, n_states);
    float sums[3] = {0.0f, 0.0f, 0.0f};
    sum3(bwd, 0, n_states, sums);
    bwd3_transition(bwd, sums, p_switch, n_states);
}

void hmm_fwd_update3x2(const hmm3_step *a, const hmm3_step *b) {
    fwd3_values(a->val, a->sums, a->p_switch, a->p_mismatch, a->m[0], a->m[1], a->m[2], a->n_states);
    fwd3_values(b->val, b->sums, b->p_switch, b->p_mismatch, b->m[0], b->m[1], b->m[2], b->n_states);
    sum3x2(a->val, a->n_states, b->val, b->n_states, a->sums, b->sums);
}

void hmm_bwd_update3x2(const hmm3_step *a, const hmm3_step *b) {
    bwd3_emit(a->val, a->p_mismatch, a->m[0], a->m[1], a->m[2], a->n_states);
    bwd3_emit(b->val, b->p_mismatch, b->m[0], b->m[1], b->m[2], b->n_states);
    float sa[3], sb[3];
    sum3x2(a->val, a->n_states, b->val, b->n_states, sa, sb);
    bwd3_transition(a->val, sa, a->p_switch, a->n_states);
    bwd3_transition(b->val, sb, b->p_switch, b->n_states);
}
