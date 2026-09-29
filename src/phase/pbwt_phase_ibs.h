/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtPhaseIbs.java and
 * phase/PbwtIbsData.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PBWT_PHASE_IBS_H
#define PHASE_PBWT_PHASE_IBS_H

#include <stdbool.h>

#include "phase/coded_steps.h"
#include "phase/phase_data.h"

/* For each stage-1 step and target haplotype, one haplotype chosen at random
 * among its nearest PBWT neighbours over the coded steps, excluding its own
 * sample and IBS2 samples, or -1. The PBWT runs backward if use_bwd (Beagle:
 * even iterations), else forward, restarting for each batch of steps with a
 * buffer of buffer= cM. Beagle's output depends on nthreads through the batches. */
typedef struct {
    phase_data *pd;            /* borrowed */
    const coded_steps *cs;     /* borrowed */
    int n_steps;
    int n_targ_haps;
    int n_candidates;
    int steps_per_batch;
    int n_overlap_steps;
    int **ibs_haps;   /* [step][target hap] */
} pbwt_phase_ibs;

/* new PbwtPhaseIbs(phaseData, codedSteps, useBwd). Writes trace seam T3d
 * (T3d-<it>-bwd or T3d-<it>-fwd). */
void pbwt_phase_ibs_init(pbwt_phase_ibs *pi, phase_data *pd, const coded_steps *cs, bool use_bwd);
void pbwt_phase_ibs_free(pbwt_phase_ibs *pi);

#endif
