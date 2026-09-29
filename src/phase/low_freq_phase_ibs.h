/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/LowFreqPhaseIbs.java and
 * phase/LowFreqPbwtPhaseIbs.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_LOW_FREQ_PHASE_IBS_H
#define PHASE_LOW_FREQ_PHASE_IBS_H

#include "phase/coded_steps.h"
#include "phase/phase_data.h"

/* The stage-2 IBS neighbours of each target haplotype at each stage-1 step,
 * from forward and backward PBWT searches over the final stage-1 haplotypes.
 * A neighbour that also carries one of the rare alleles between this step
 * and the next is preferred when its match reaches back (or forward) far
 * enough; otherwise one of the 10 nearest is chosen at random. -1 if none. */
typedef struct {
    const phase_data *pd;   /* borrowed */
    coded_steps cs;   /* the final stage-1 haplotypes (LowFreqPhaseIbs.allHaps) */
    int n_steps;
    int n_targ_haps;
    int **fwd;   /* [step][target hap] */
    int **bwd;
} low_freq_phase_ibs;

/* new LowFreqPhaseIbs(phaseData): builds its own CodedSteps. Writes trace
 * seam T4a. */
void low_freq_phase_ibs_init(low_freq_phase_ibs *lf, const phase_data *pd);
void low_freq_phase_ibs_free(low_freq_phase_ibs *lf);

#endif
