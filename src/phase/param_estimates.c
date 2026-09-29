/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/ParamEstimates.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/param_estimates.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "blbutil/utilities.h"

void param_estimates_init(param_estimates *pe) {
    *pe = (param_estimates){0};
}

void param_estimates_free(param_estimates *pe) {
    free(pe->mismatch);
    free(pe->recomb);
}

void param_estimates_add_mismatch_data(param_estimates *pe, int marker_cnt, double p_mismatch_sum) {
    if (marker_cnt > 0 && p_mismatch_sum > 0 && isfinite(p_mismatch_sum)) {
        if (pe->n_mismatch == pe->mismatch_cap) {
            pe->mismatch_cap = pe->mismatch_cap == 0 ? 64 : 2 * pe->mismatch_cap;
            pe->mismatch = util_realloc(pe->mismatch, (size_t)pe->mismatch_cap * sizeof *pe->mismatch);
        }
        pe->mismatch[pe->n_mismatch++] = (mismatch_data){marker_cnt, p_mismatch_sum};
    }
}

void param_estimates_add_switch_data(param_estimates *pe, double gen_distances, double switch_probs) {
    if (gen_distances > 0 && switch_probs > 0 && isfinite(gen_distances) && isfinite(switch_probs)) {
        if (pe->n_recomb == pe->recomb_cap) {
            pe->recomb_cap = pe->recomb_cap == 0 ? 64 : 2 * pe->recomb_cap;
            pe->recomb = util_realloc(pe->recomb, (size_t)pe->recomb_cap * sizeof *pe->recomb);
        }
        pe->recomb[pe->n_recomb++] = (recomb_data){gen_distances, switch_probs};
    }
}

/* The stored values are positive and finite, where Double.compare is numeric
 * order. */
static int compare_doubles(double a, double b) {
    return (a > b) - (a < b);
}

static int compare_mismatch(const void *a, const void *b) {
    const mismatch_data *x = a, *y = b;
    int val = compare_doubles(x->p_mismatch_sum, y->p_mismatch_sum);
    return val != 0 ? val : (x->marker_cnt > y->marker_cnt) - (x->marker_cnt < y->marker_cnt);
}

static int compare_recomb(const void *a, const void *b) {
    const recomb_data *x = a, *y = b;
    int val = compare_doubles(x->gen_distance, y->gen_distance);
    return val != 0 ? val : compare_doubles(x->switch_prob, y->switch_prob);
}

float param_estimates_p_mismatch(param_estimates *pe) {
    qsort(pe->mismatch, (size_t)pe->n_mismatch, sizeof *pe->mismatch, compare_mismatch);
    int64_t sum_markers = 0;
    double sum_p_mismatch = 0.0;
    for (int j = 0; j < pe->n_mismatch; ++j) {
        sum_markers += pe->mismatch[j].marker_cnt;
        sum_p_mismatch += pe->mismatch[j].p_mismatch_sum;
    }
    return sum_markers == 0 ? NAN : (float)(sum_p_mismatch / sum_markers);
}

float param_estimates_recomb_intensity(param_estimates *pe) {
    qsort(pe->recomb, (size_t)pe->n_recomb, sizeof *pe->recomb, compare_recomb);
    double sum_switches = 0.0;
    double sum_distances = 0.0;
    for (int j = 0; j < pe->n_recomb; ++j) {
        sum_switches += pe->recomb[j].switch_prob;
        sum_distances += pe->recomb[j].gen_distance;
    }
    return sum_distances == 0.0 ? NAN : (float)(sum_switches / sum_distances);
}
