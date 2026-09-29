/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/RevPbwtPhaser.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_REV_PBWT_PHASER_H
#define PHASE_REV_PBWT_PHASER_H

#include <stdint.h>

#include "phase/fixed_phase_data.h"

/* The target haplotypes at stage-1 markers [start, end) phased by a PBWT run
 * from the last marker to the first. Heterozygotes left unphased get a random
 * phase, and missing alleles left unimputed are drawn from the allele
 * frequencies, from new Random(seed). */
typedef struct {
    int start;
    int end;
    int *bits_per_allele;   /* per marker */
    uint64_t **bits;        /* per marker, bits_per_allele bits per target haplotype */
} rev_pbwt_phaser;

/* new RevPbwtPhaser(fpd, start, end, seed). Writes trace seam T3b0. */
void rev_pbwt_phaser_init(rev_pbwt_phaser *rp, const fixed_phase_data *fpd, int start, int end, int64_t seed);
int rev_pbwt_phaser_allele(const rev_pbwt_phaser *rp, int m, int hap);
void rev_pbwt_phaser_free(rev_pbwt_phaser *rp);

#endif
