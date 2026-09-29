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
#ifndef PHASE_PARAM_ESTIMATES_H
#define PHASE_PARAM_ESTIMATES_H

typedef struct {
    int marker_cnt;
    double p_mismatch_sum;
} mismatch_data;

typedef struct {
    double gen_distance;
    double switch_prob;
} recomb_data;

/* Per-sample sums from HmmParamData. The estimates sort them before summing,
 * as Java does, so the order samples are added in does not matter. */
typedef struct {
    mismatch_data *mismatch;
    int n_mismatch, mismatch_cap;
    recomb_data *recomb;
    int n_recomb, recomb_cap;
} param_estimates;

void param_estimates_init(param_estimates *pe);
void param_estimates_add_mismatch_data(param_estimates *pe, int marker_cnt, double p_mismatch_sum);
void param_estimates_add_switch_data(param_estimates *pe, double gen_distances, double switch_probs);
/* NaN without data. Sorts the stored data. */
float param_estimates_p_mismatch(param_estimates *pe);
float param_estimates_recomb_intensity(param_estimates *pe);
void param_estimates_free(param_estimates *pe);

#endif
