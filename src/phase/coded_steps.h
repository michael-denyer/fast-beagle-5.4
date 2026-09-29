/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/CodedSteps.java and the parts
 * of vcf/XRefGT.java it uses; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_CODED_STEPS_H
#define PHASE_CODED_STEPS_H

#include <stdint.h>

#include "phase/phase_data.h"

/* Each haplotype's sequence index within each stage-1 step, for the current
 * target haplotypes followed by the reference haplotypes. Sequences are
 * numbered in order of first appearance and told apart by BitArray.hash of
 * the step's bits, so haplotypes whose bits collide share an index. */
typedef struct {
    int n_steps;
    int n_haps;
    const uint64_t **haps;   /* per haplotype: rows of targ_bits, then the reference haplotypes of FixedPhaseData */
    uint64_t *targ_bits;     /* owned copies of the target haplotypes, one row each */
    int **hap_to_seq;        /* [step][hap] */
    int *n_seq;              /* per step */
} coded_steps;

/* new CodedSteps(estPhase) */
void coded_steps_init(coded_steps *cs, const phase_data *pd);
void coded_steps_free(coded_steps *cs);

#endif
