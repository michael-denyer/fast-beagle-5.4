/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/FixedPhaseData.java and
 * vcf/Window.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_FIXED_PHASE_DATA_H
#define PHASE_FIXED_PHASE_DATA_H

#include <stdint.h>

#include "main/par.h"
#include "phase/ibs2.h"
#include "vcf/marker_map.h"
#include "vcf/sliding_window.h"

/* The carriers of one allele. Java tells the two empty kinds apart by
 * identity (Window.ZERO_FREQ_ARRAY, Window.HIGH_FREQ_ARRAY). */
typedef enum {
    CARRIERS_LOW_FREQ,    /* at most max_carriers samples, listed */
    CARRIERS_ZERO_FREQ,   /* no carriers */
    CARRIERS_HIGH_FREQ,   /* more than max_carriers samples, or all alleles when rare carriers are ignored */
} carriers_kind;

typedef struct {
    carriers_kind kind;
    int n;
    int *samples;   /* increasing; reference samples follow the target samples */
} carriers;

typedef struct {
    int n_alleles;
    carriers *allele;
} marker_carriers;

/* FixedPhaseData for one window. Stage-1 markers
 * are those with at least two high-frequency alleles, unless there are fewer
 * than two such markers or they exceed 75% of the markers, in which case
 * every marker is a stage-1 marker. */
typedef struct fixed_phase_data {
    const window *win;        /* owned by the caller; outlives this */
    int window;
    int n_markers;            /* target markers */
    int n_haps;               /* target plus reference */
    int overlap;              /* markers phased by the previous window */
    const int *const *overlap_alleles; /* [hap][marker]: that phase, owned by the caller */
    marker_map map;
    marker_carriers *carriers;
    int n_stage1;
    int *stage1_to2;          /* stage-1 marker to marker */
    int *stage1_hap_bits;     /* Markers.sumHapBits: first bit of each stage-1 marker in a haplotype, n_stage1 + 1 entries */
    int n_ref_haps;
    uint64_t **stage1_ref_haps; /* stage1XRefGT: each reference haplotype's stage-1 alleles at stage1_hap_bits */
    marker_map stage1_map;
    float ibs_step;
    steps stage1_steps;
    int stage1_overlap;
    int *prev_stage1_marker;  /* per marker, the last stage-1 marker at or before it */
    float *prev_stage1_wt;    /* per marker, the interpolation weight of that stage-1 marker */
    float *stage1_maf;        /* per stage-1 marker, from at most 10,000 haplotypes */
    ibs2 stage1_ibs2;
} fixed_phase_data;

/* new FixedPhaseData(par, ped, window, phasedOverlap): overlap_alleles[hap]
 * holds the previous window's phased alleles at this window's first overlap
 * target markers (not read when overlap is 0). */
void fixed_phase_data_init(fixed_phase_data *fpd, const par *p, const genetic_map *gm, const window *w, int overlap,
        const int *const *overlap_alleles);
void fixed_phase_data_free(fixed_phase_data *fpd);

/* The stage-1 genotypes (stage1TargGT, stage1RefGT) by stage-1 marker. */
static inline int fpd_n_targ_haps(const fixed_phase_data *fpd) {
    return fpd->win->targ[0]->n_haps;
}
static inline const marker *fpd_marker(const fixed_phase_data *fpd, int m) {
    return &fpd->win->targ[fpd->stage1_to2[m]]->marker;
}
static inline int fpd_n_alleles(const fixed_phase_data *fpd, int m) {
    return marker_n_alleles(fpd_marker(fpd, m));
}
/* SplicedGT.allele: FixedPhaseData.targGT at target marker m, the previous
 * window's phase in the overlap and the input after it. */
static inline int fpd_spliced_allele(const void *ctx, int m, int hap) {
    const fixed_phase_data *fpd = ctx;
    return m < fpd->overlap ? fpd->overlap_alleles[hap][m] : gt_rec_get(fpd->win->targ[m], hap);
}
static inline int fpd_targ_allele(const fixed_phase_data *fpd, int m, int hap) {
    return fpd_spliced_allele(fpd, fpd->stage1_to2[m], hap);
}
static inline int fpd_ref_allele(const fixed_phase_data *fpd, int m, int hap) {
    const window *w = fpd->win;
    return ref_gt_rec_get(w->ref[w->indices.targ_marker_to_marker[fpd->stage1_to2[m]]], hap);
}

/* Trace seam T2b. */
void fixed_phase_data_trace(const fixed_phase_data *fpd);

#endif
