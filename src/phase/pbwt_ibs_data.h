/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtIbsData.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PBWT_IBS_DATA_H
#define PHASE_PBWT_IBS_DATA_H

#include <stdbool.h>

#include "beagleutil/pbwt_div_updater.h"
#include "phase/coded_steps.h"
#include "phase/phase_data.h"

/* The parameters of an IBS neighbour search over n_steps steps: the number
 * of candidates by iteration (10 after the stage-1 iterations), the buffer
 * and the stage-2 backoff in steps, and the batches of steps, which depend on
 * nthreads. */
typedef struct {
    int n_steps;
    int n_candidates;
    int n_overlap_steps;
    int max_backoff_steps;
    int steps_per_batch;
    int n_batches;
} pbwt_ibs_data;

void pbwt_ibs_data_init(pbwt_ibs_data *d, const phase_data *pd, int n_steps);

/* One batch's PBWT over the coded steps. pbwt_batch_init restarts the PBWT
 * and runs it over the batch's buffer steps; each pbwt_batch_next adds the
 * batch's next step (backward if bwd) to a and d and sets step, until the
 * batch is done. */
typedef struct {
    const coded_steps *cs;
    bool bwd;
    int step;
    int stop;
    int *a;   /* prefix array */
    int *d;   /* divergence array, n_haps + 1 entries */
    pbwt_div_updater u;
} pbwt_batch;

void pbwt_batch_init(pbwt_batch *b, const pbwt_ibs_data *data, const coded_steps *cs, int batch, bool bwd);
bool pbwt_batch_next(pbwt_batch *b);
void pbwt_batch_free(pbwt_batch *b);

#endif
