/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/FwdPbwtPhaser.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_FWD_PBWT_PHASER_H
#define PHASE_FWD_PBWT_PHASER_H

#include <stdint.h>

#include "phase/fixed_phase_data.h"

/* The target haplotypes at stage-1 markers [start, end) phased by a PBWT run
 * from the first marker to the last. Heterozygotes the PBWT leaves unphased
 * take their phase relative to the sample's previous heterozygote from a
 * reverse PBWT phase, which also supplies missing alleles left unimputed. */
typedef struct {
    int start;
    int end;
    int *bits_per_allele;   /* per marker */
    uint64_t **bits;        /* per marker, bits_per_allele bits per target haplotype */
} fwd_pbwt_phaser;

/* new FwdPbwtPhaser(fpd, start, end, seed) */
void fwd_pbwt_phaser_init(fwd_pbwt_phaser *fp, const fixed_phase_data *fpd, int start, int end, int64_t seed);
int fwd_pbwt_phaser_allele(const fwd_pbwt_phaser *fp, int m, int hap);
void fwd_pbwt_phaser_free(fwd_pbwt_phaser *fp);

#endif
