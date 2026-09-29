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
#include "phase/rev_pbwt_phaser.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include <htslib/kstring.h>

#include "blbutil/bit_array.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jrandom.h"
#include "phase/pbwt_rec_phaser.h"

static int impute_allele(const int *cdf, int n_alleles, jrandom *r) {
    int bound = cdf[n_alleles - 1];
    if (bound == 0) return 0;
    int x = jrandom_next_int_bound(r, bound);
    int allele = 0;
    while (x >= cdf[allele]) ++allele;
    return allele;
}

static void finish_phasing(int *alleles, bool *unph_het, int n_samples, const int *cdf, int n_alleles, jrandom *r) {
    for (int s = 0; s < n_samples; ++s) {
        int h1 = s << 1;
        int h2 = h1 | 1;
        if (unph_het[s]) {
            if (jrandom_next_boolean(r)) {
                int tmp = alleles[h1];
                alleles[h1] = alleles[h2];
                alleles[h2] = tmp;
            }
            unph_het[s] = false;
        } else {
            if (alleles[h1] == -1) alleles[h1] = impute_allele(cdf, n_alleles, r);
            if (alleles[h2] == -1) alleles[h2] = impute_allele(cdf, n_alleles, r);
        }
    }
}

static void trace(const rev_pbwt_phaser *rp, int n_targ_haps, int64_t seed) {
    char seam[32];
    snprintf(seam, sizeof seam, "T3b0-%d", rp->start);
    trace_line(seam, "rev\t%d\t%d\t%" PRId64, rp->start, rp->end, seed);
    kstring_t s = {0, 0, NULL};
    for (int m = rp->start; m < rp->end; ++m) {
        s.l = 0;
        kputw(m, &s);
        kputc('\t', &s);
        for (int h = 0; h < n_targ_haps; ++h) {
            if (h > 0) kputc(',', &s);
            kputw(rev_pbwt_phaser_allele(rp, m, h), &s);
        }
        trace_line(seam, "%s", s.s);
    }
    free(s.s);
}

void rev_pbwt_phaser_init(rev_pbwt_phaser *rp, const fixed_phase_data *fpd, int start, int end, int64_t seed) {
    if (start < 0 || end > fpd->n_stage1 || start >= end) util_exit("java.lang.IllegalArgumentException: %d", start);
    int n_targ_haps = fpd_n_targ_haps(fpd);
    int n_samples = n_targ_haps >> 1;
    rp->start = start;
    rp->end = end;
    rp->bits_per_allele = util_malloc((size_t)(end - start) * sizeof *rp->bits_per_allele);
    rp->bits = util_malloc((size_t)(end - start) * sizeof *rp->bits);

    jrandom r;
    jrandom_init(&r, seed);
    pbwt_rec_phaser rec;
    pbwt_rec_phaser_init(&rec, fpd);
    bool *missing = util_malloc((size_t)n_samples * sizeof *missing);
    bool *unph_het = util_malloc((size_t)n_samples * sizeof *unph_het);
    int *alleles = util_malloc((size_t)fpd->n_haps * sizeof *alleles);
    int last_m = -1;
    for (int m = end - 1; m >= start; --m) {
        const int *cdf = pbwt_rec_phaser_phase(&rec, last_m, alleles, m, missing, unph_het);
        if (m >= fpd->stage1_overlap) finish_phasing(alleles, unph_het, n_samples, cdf, fpd_n_alleles(fpd, m), &r);
        int bpa = marker_bits_per_allele(fpd_marker(fpd, m));
        rp->bits_per_allele[m - start] = bpa;
        rp->bits[m - start] = allele_bits_store(alleles, n_targ_haps, bpa);
        last_m = m;
    }
    free(alleles);
    free(unph_het);
    free(missing);
    pbwt_rec_phaser_free(&rec);
    if (trace_on()) trace(rp, n_targ_haps, seed);
}

int rev_pbwt_phaser_allele(const rev_pbwt_phaser *rp, int m, int hap) {
    return allele_bits_get(rp->bits[m - rp->start], hap, rp->bits_per_allele[m - rp->start]);
}

void rev_pbwt_phaser_free(rev_pbwt_phaser *rp) {
    for (int m = 0; m < rp->end - rp->start; ++m) free(rp->bits[m]);
    free(rp->bits);
    free(rp->bits_per_allele);
}
