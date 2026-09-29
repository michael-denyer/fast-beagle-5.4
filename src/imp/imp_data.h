/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) imp/ImpData.java and
 * imp/HaplotypeCoder.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef IMP_IMP_DATA_H
#define IMP_IMP_DATA_H

#include "main/par.h"
#include "vcf/genetic_map.h"
#include "vcf/sliding_window.h"

/* One cluster's haplotype-to-sequence map over the reference haplotypes then
 * the target haplotypes. Sequences are numbered from 1 in order of first
 * appearance among the target haplotypes; a reference haplotype whose
 * sequence no target haplotype carries is 0. When the cluster lies within one
 * SeqCoder3 group, the reference part is the group's map composed with a
 * per-sequence table, as in Java; otherwise it is stored per haplotype. */
typedef struct {
    const seq_group *group;   /* NULL when stored per haplotype */
    int *ref_seq;             /* group ? per group sequence : per reference haplotype */
    int *targ_seq;            /* per target haplotype */
    int n_seq;                /* IndexArray.valueSize */
} hap_seq_map;

/* The data for imputing one window: target markers grouped into clusters
 * of at most cluster= cM that do not cross a change of reference sequence
 * group, each cluster's mean cM position, mismatch and switch
 * probabilities, and the interpolation weights of the reference markers
 * between clusters. */
typedef struct {
    const window *w;
    phased_allele_fn targ_allele;   /* the phased target */
    const void *targ_ctx;
    int n_clusters;
    int *targ_clust_start_end;   /* n_clusters + 1 target marker boundaries */
    int *ref_cluster_start;
    int *ref_cluster_end;
    hap_seq_map *hap_to_seq;
    float *err_prob;
    double *pos;
    float *p_recomb;
    float *weight;               /* per reference marker; NaN in and outside clusters */
    int n_ref_haps;
    int n_targ_haps;
    int n_input_targ_haps;       /* haploid samples count once */
    int n_haps;
} imp_data;

/* new ImpData(par, window, phasedTarg, map). Writes trace seam T5a. */
void imp_data_init(imp_data *id, const par *p, const window *w, const samples *targ_samples,
        phased_allele_fn targ_allele, const void *targ_ctx, const genetic_map *gm);
/* ImpData.allele(cluster, hap): hap's sequence in the cluster. */
static inline int imp_data_seq(const imp_data *id, int cluster, int hap) {
    const hap_seq_map *hs = &id->hap_to_seq[cluster];
    if (hap >= id->n_ref_haps) return hs->targ_seq[hap - id->n_ref_haps];
    return hs->group != NULL ? hs->ref_seq[hs->group->hap_to_seq[hap]] : hs->ref_seq[hap];
}
void imp_data_free(imp_data *id);

#endif
