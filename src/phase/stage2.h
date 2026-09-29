/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/HmmStateProbs.java,
 * phase/Stage2Baum.java, phase/Stage2Haps.java and PhaseLS.runStage2;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_STAGE2_H
#define PHASE_STAGE2_H

#include <pthread.h>

#include "blbutil/int_list.h"
#include "phase/phase_data.h"

/* The stage-2 phase of the markers left out of stage 1, kept as Java's
 * Stage2Haps keeps it: for each such marker and each allele that is not
 * high-frequency, the target haplotypes carrying it, in increasing order.
 * Every other haplotype carries the marker's one high-frequency allele.
 * Stage-1 markers keep the final stage-1 phase. */
typedef struct {
    const phase_data *pd;
    int_list **carriers;   /* [marker][allele]: NULL rows for stage-1 markers; unused for alleles FixedPhaseData marks high-frequency */
    pthread_mutex_t *locks; /* per marker, guarding its carrier lists while samples are phased in parallel */
    int *major;            /* per stage-2 marker, the allele of haplotypes in no list */
} stage2_haps;

/* PhaseLS.runStage2(pd). Writes trace seams T4a, T4b, T4c and T4d. */
void phase_ls_run_stage2(stage2_haps *s2, const phase_data *pd);
/* Stage2Haps.toGTRecs: target haplotype hap's allele at target marker m, from
 * the final stage-1 phase at stage-1 markers. */
int stage2_haps_allele(const void *s2, int m, int hap);
/* Frees s2, or does nothing for a zero-initialised s2. */
void stage2_haps_free(stage2_haps *s2);

#endif
