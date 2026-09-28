/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PhaseData.java,
 * phase/EstPhase.java and vcf/MarkerMap.java pRecomb; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PHASE_DATA_H
#define PHASE_PHASE_DATA_H

#include <stdint.h>

#include "jcompat/jrandom.h"
#include "main/par.h"
#include "phase/fixed_phase_data.h"
#include "phase/sample_phase.h"

/* The phasing state of a window: each target sample's current phase
 * (EstPhase) and the model parameters, which the iterations update. */
typedef struct {
    const fixed_phase_data *fpd;
    const par *par;
    int n_samples;
    sample_phase *phase;
    int64_t seed;
    int it;
    float *leave_unph_prop;   /* per sample, from its initial unphased count */
    float recomb_intensity;
    float *p_recomb;       /* per stage-1 marker */
    float p_mismatch;
} phase_data;

/* new PhaseData(fpd, seed). Writes trace seam T3b1. */
void phase_data_init(phase_data *pd, const fixed_phase_data *fpd, const par *p, int64_t seed);
void phase_data_free(phase_data *pd);
/* Trace seam T3b: the iteration count, the swap rate (raw bits), and each
 * sample's haplotypes, clusters, and unphased and missing markers. */
void phase_data_trace_iteration(const phase_data *pd, double swap_rate);
/* PhaseData.incrementIt and advanceToFirstPhasingIt. */
void phase_data_increment_it(phase_data *pd);
void phase_data_advance_to_first_phasing_it(phase_data *pd);
/* PhaseData.updatePMismatch: p_mismatch must be in [0, 1]. */
void phase_data_update_p_mismatch(phase_data *pd, float p_mismatch);
/* PhaseData.updateRecombIntensity: recomputes pRecomb. */
void phase_data_update_recomb_intensity(phase_data *pd, float recomb_intensity);

/* PhaseData.seed(): the window seed plus the iteration */
static inline int64_t phase_data_seed(const phase_data *pd) {
    return jrandom_seed_plus(pd->seed, pd->it);
}

#endif
