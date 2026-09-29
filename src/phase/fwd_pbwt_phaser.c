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
#include "phase/fwd_pbwt_phaser.h"

#include <stdbool.h>
#include <stdlib.h>

#include "blbutil/bit_array.h"
#include "blbutil/utilities.h"
#include "phase/pbwt_rec_phaser.h"
#include "phase/rev_pbwt_phaser.h"

static int impute_allele(const fwd_pbwt_phaser *fp, const rev_pbwt_phaser *rev, int last_het, int m, int hap) {
    if (last_het < 0) return rev_pbwt_phaser_allele(rev, m, hap);
    int comp_hap = hap ^ 1;
    int a1 = rev_pbwt_phaser_allele(rev, last_het, hap);
    int a2 = rev_pbwt_phaser_allele(rev, last_het, comp_hap);
    int b1 = fwd_pbwt_phaser_allele(fp, last_het, hap);
    int b2 = fwd_pbwt_phaser_allele(fp, last_het, comp_hap);
    return (a1 < a2) == (b1 < b2) ? rev_pbwt_phaser_allele(rev, m, hap) : rev_pbwt_phaser_allele(rev, m, comp_hap);
}

static void finish_phasing(const fwd_pbwt_phaser *fp, const rev_pbwt_phaser *rev, int m, int *alleles,
        const int *last_het, bool *unph_het, int n_samples) {
    for (int s = 0; s < n_samples; ++s) {
        int h1 = s << 1;
        int h2 = h1 | 1;
        if (unph_het[s]) {
            int prev_het = last_het[s];
            if (prev_het >= 0) {
                int a1 = rev_pbwt_phaser_allele(rev, prev_het, h1);
                int a2 = rev_pbwt_phaser_allele(rev, prev_het, h2);
                int b1 = rev_pbwt_phaser_allele(rev, m, h1);
                int b2 = rev_pbwt_phaser_allele(rev, m, h2);
                bool rev_same_phase = (a1 < a2) == (b1 < b2);
                int c1 = fwd_pbwt_phaser_allele(fp, prev_het, h1);
                int c2 = fwd_pbwt_phaser_allele(fp, prev_het, h2);
                bool fwd_same_phase = (c1 < c2) == (alleles[h1] < alleles[h2]);
                if (rev_same_phase != fwd_same_phase) {
                    int tmp = alleles[h1];
                    alleles[h1] = alleles[h2];
                    alleles[h2] = tmp;
                }
            }
            unph_het[s] = false;
        } else {
            if (alleles[h1] == -1) alleles[h1] = impute_allele(fp, rev, last_het[s], m, h1);
            if (alleles[h2] == -1) alleles[h2] = impute_allele(fp, rev, last_het[s], m, h2);
        }
    }
}

void fwd_pbwt_phaser_init(fwd_pbwt_phaser *fp, const fixed_phase_data *fpd, int start, int end, int64_t seed) {
    if (start < 0 || end > fpd->n_markers || start >= end) util_exit("java.lang.IllegalArgumentException: %d", start);
    int n_targ_haps = fpd_n_targ_haps(fpd);
    int n_samples = n_targ_haps >> 1;
    fp->start = start;
    fp->end = end;
    fp->bits_per_allele = util_malloc((size_t)(end - start) * sizeof *fp->bits_per_allele);
    fp->bits = util_malloc((size_t)(end - start) * sizeof *fp->bits);

    pbwt_rec_phaser rec;
    pbwt_rec_phaser_init(&rec, fpd);
    rev_pbwt_phaser rev;
    rev_pbwt_phaser_init(&rev, fpd, start, end, seed);
    bool *missing = util_malloc((size_t)n_samples * sizeof *missing);
    bool *unph_het = util_malloc((size_t)n_samples * sizeof *unph_het);
    int *last_het = util_malloc((size_t)n_samples * sizeof *last_het);
    for (int s = 0; s < n_samples; ++s) last_het[s] = -1;
    int *alleles = util_malloc((size_t)fpd->n_haps * sizeof *alleles);
    int last_m = -1;
    for (int m = start; m < end; ++m) {
        pbwt_rec_phaser_phase(&rec, last_m, alleles, m, missing, unph_het);
        if (m >= fpd->stage1_overlap) finish_phasing(fp, &rev, m, alleles, last_het, unph_het, n_samples);
        int bpa = marker_bits_per_allele(fpd_marker(fpd, m));
        fp->bits_per_allele[m - start] = bpa;
        fp->bits[m - start] = allele_bits_store(alleles, n_targ_haps, bpa);
        for (int s = 0; s < n_samples; ++s) {
            if (!missing[s] && alleles[s << 1] != alleles[(s << 1) | 1]) last_het[s] = m;
        }
        last_m = m;
    }
    free(alleles);
    free(last_het);
    free(unph_het);
    free(missing);
    rev_pbwt_phaser_free(&rev);
    pbwt_rec_phaser_free(&rec);
}

int fwd_pbwt_phaser_allele(const fwd_pbwt_phaser *fp, int m, int hap) {
    return allele_bits_get(fp->bits[m - fp->start], hap, fp->bits_per_allele[m - fp->start]);
}

void fwd_pbwt_phaser_free(fwd_pbwt_phaser *fp) {
    for (int m = 0; m < fp->end - fp->start; ++m) free(fp->bits[m]);
    free(fp->bits);
    free(fp->bits_per_allele);
}
