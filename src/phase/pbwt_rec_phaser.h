/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtRecPhaser.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PBWT_REC_PHASER_H
#define PHASE_PBWT_REC_PHASER_H

#include <stdbool.h>

#include "beagleutil/pbwt_updater.h"
#include "phase/fixed_phase_data.h"

/* Phases one stage-1 marker at a time from the haplotypes adjacent in a PBWT
 * of the markers already visited, in either direction. Heterozygotes whose
 * neighbours agree on a phase are phased first; the agreement required drops
 * from 2 to 0 as passes stop changing anything. Missing alleles are copied
 * from agreeing neighbours. */
typedef struct {
    const fixed_phase_data *fpd;
    int n_targ_haps;
    int n_targ_samples;
    int n_haps;
    int *a;          /* PBWT prefix array */
    int *inv_a;
    pbwt_updater pbwt;
    int *cdf;
    int cdf_cap;
} pbwt_rec_phaser;

void pbwt_rec_phaser_init(pbwt_rec_phaser *p, const fixed_phase_data *fpd);
/* PbwtRecPhaser.phase(currentMkr, alleles, nextMkr, missing, unphHet): adds
 * the alleles of current_mkr (-1 for none) to the PBWT, loads next_mkr's
 * alleles into alleles[], and phases them at or after the stage-1 overlap.
 * Returns the cumulative allele counts at next_mkr, valid until the next call. */
const int *pbwt_rec_phaser_phase(pbwt_rec_phaser *p, int current_mkr, int *alleles, int next_mkr, bool *missing, bool *unph_het);
void pbwt_rec_phaser_free(pbwt_rec_phaser *p);

#endif
