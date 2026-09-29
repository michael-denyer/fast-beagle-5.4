/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) imp/CodedSteps.java and
 * imp/ImpIbs.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef IMP_IMP_IBS_H
#define IMP_IMP_IBS_H

#include "imp/imp_data.h"

/* A shared list of reference haplotypes. */
typedef struct {
    int n;
    int *haps;
} ibs_set;

/* For each imputation step (imp-step= cM of clusters, the first half as
 * long) and target haplotype, up to imp-states / (imp-segment / imp-step)
 * reference haplotypes that share its allele sequence over the next
 * imp-nsteps steps, found by partitioning the haplotypes step by step. Sets
 * are shared by the target haplotypes that reach them together. */
typedef struct {
    int n_steps;
    int *step_starts;       /* first cluster of each step */
    int **step_seq;         /* [step][hap]: CodedSteps, over reference then target haplotypes */
    int *step_n_seq;
    int n_targ_haps;
    const ibs_set ***ibs;   /* [step][target hap] */
    ibs_set **owned;        /* [step]: the distinct sets allocated for the step */
    int *n_owned;
} imp_ibs;

/* new ImpIbs(impData). Writes trace seam T5b. */
void imp_ibs_init(imp_ibs *ib, const imp_data *id, const par *p);
void imp_ibs_free(imp_ibs *ib);

#endif
