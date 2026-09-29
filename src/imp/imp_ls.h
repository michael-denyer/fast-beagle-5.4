/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) imp/ImpLS.java, imp/ImpStates.java,
 * imp/ImpLSBaum.java, imp/StateProbsFactory.java and imp/StateProbs.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef IMP_IMP_LS_H
#define IMP_IMP_LS_H

#include "imp/imp_data.h"

/* One state kept for a target haplotype at a cluster: its reference
 * haplotype and its probability there and at the next cluster (the same
 * cluster for the last). The writer reads the three together. */
typedef struct {
    int ref_hap;
    float prob;
    float prob_p1;
} kept_state;

/* One target haplotype's state probabilities kept at one cluster: the
 * states whose probability there or at the next cluster exceeds
 * min(0.005, 0.9999 / nStates). */
typedef struct {
    int n;
    kept_state *s;
} cluster_probs;

/* StateProbs of every target haplotype: a haplotype's kept states live in
 * one block, and cluster (hap, c) is at cluster[hap * n_clusters + c]. */
typedef struct {
    int n_targ_haps;
    int n_clusters;
    cluster_probs *cluster;
} state_probs;

static inline const cluster_probs *state_probs_at(const state_probs *sp, int hap, int cluster) {
    return &sp->cluster[(size_t)hap * (size_t)sp->n_clusters + (size_t)cluster];
}

/* ImpLS.stateProbs(impData): each target haplotype's Li and Stephens HMM
 * over its composite reference haplotypes (ImpStates), from the IBS sets it
 * builds and frees. Writes trace seams T5b and T5c. */
state_probs *imp_ls_state_probs(const imp_data *id, const par *p);
void state_probs_free(state_probs *sp);

#endif
