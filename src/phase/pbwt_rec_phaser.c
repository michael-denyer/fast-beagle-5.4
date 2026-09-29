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
#include "phase/pbwt_rec_phaser.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"

void pbwt_rec_phaser_init(pbwt_rec_phaser *p, const fixed_phase_data *fpd) {
    p->fpd = fpd;
    p->n_targ_haps = fpd_n_targ_haps(fpd);
    p->n_targ_samples = p->n_targ_haps >> 1;
    p->n_haps = fpd->n_haps;
    p->a = util_malloc((size_t)p->n_haps * sizeof *p->a);
    p->inv_a = util_malloc((size_t)p->n_haps * sizeof *p->inv_a);
    for (int h = 0; h < p->n_haps; ++h) p->a[h] = p->inv_a[h] = h;
    pbwt_updater_init(&p->pbwt, p->n_haps);
    p->cdf = NULL;
    p->cdf_cap = 0;
}

void pbwt_rec_phaser_free(pbwt_rec_phaser *p) {
    free(p->a);
    free(p->inv_a);
    pbwt_updater_free(&p->pbwt);
    free(p->cdf);
}

static bool is_phased_neighbour(const pbwt_rec_phaser *p, const bool *unph_het, int h) {
    int s = h >> 1;
    return s >= p->n_targ_samples || !unph_het[s];
}

static int phase_cnt(int adjacent, int a1, int a2) {
    if (adjacent == a1) return 1;
    if (adjacent == a2) return -1;
    return 0;
}

static int neighbour_phase_cnt(const pbwt_rec_phaser *p, const int *alleles, const bool *unph_het, int ai, int a1, int a2) {
    int cnt = 0;
    if (ai > 0) {
        int h = p->a[ai - 1];
        if (is_phased_neighbour(p, unph_het, h)) cnt += phase_cnt(alleles[h], a1, a2);
    }
    if (ai + 1 < p->n_haps) {
        int h = p->a[ai + 1];
        if (is_phased_neighbour(p, unph_het, h)) cnt += phase_cnt(alleles[h], a1, a2);
    }
    return cnt;
}

static bool phase_het(const pbwt_rec_phaser *p, int s, int threshold, int *alleles, bool *unph_het) {
    int h1 = s << 1;
    int h2 = h1 | 1;
    int a1 = alleles[h1];
    int a2 = alleles[h2];
    int cnt = neighbour_phase_cnt(p, alleles, unph_het, p->inv_a[h1], a1, a2)
            + neighbour_phase_cnt(p, alleles, unph_het, p->inv_a[h2], a2, a1);
    if (cnt >= threshold) {
        unph_het[s] = false;
        return true;
    }
    if (cnt <= -threshold) {
        alleles[h1] = a2;
        alleles[h2] = a1;
        unph_het[s] = false;
        return true;
    }
    return false;
}

static int impute(const pbwt_rec_phaser *p, const int *alleles, const bool *unph_het, int ai) {
    int prev = -1;
    int next = -1;
    if (ai > 0) {
        int h = p->a[ai - 1];
        if (is_phased_neighbour(p, unph_het, h)) prev = alleles[h];
    }
    if (ai + 1 < p->n_haps) {
        int h = p->a[ai + 1];
        if (is_phased_neighbour(p, unph_het, h)) next = alleles[h];
    }
    if (prev >= 0 && (prev == next || next < 0)) return prev;
    if (prev < 0 && next >= 0) return next;
    return -1;
}

static void phase(pbwt_rec_phaser *p, int *alleles, bool *unph_het) {
    for (int j = 0; j < p->n_haps; ++j) p->inv_a[p->a[j]] = j;
    int threshold = 2;
    bool change_made = true;
    while (threshold > 0 || change_made) {
        change_made = false;
        for (int s = 0; s < p->n_targ_samples; ++s) {
            if (unph_het[s]) {
                change_made |= phase_het(p, s, threshold, alleles, unph_het);
            } else {
                int h1 = s << 1;
                int h2 = h1 | 1;
                if (alleles[h1] == -1) {
                    alleles[h1] = impute(p, alleles, unph_het, p->inv_a[h1]);
                    change_made |= alleles[h1] >= 0;
                }
                if (alleles[h2] == -1) {
                    alleles[h2] = impute(p, alleles, unph_het, p->inv_a[h2]);
                    change_made |= alleles[h2] >= 0;
                }
            }
        }
        if (!change_made) --threshold;
    }
}

static const int *set_alleles(pbwt_rec_phaser *p, int m, int *alleles, bool *unph_het, bool *missing) {
    const fixed_phase_data *fpd = p->fpd;
    int n_alleles = fpd_n_alleles(fpd, m);
    if (n_alleles > p->cdf_cap) {
        p->cdf_cap = n_alleles;
        p->cdf = util_realloc(p->cdf, (size_t)p->cdf_cap * sizeof *p->cdf);
    }
    memset(p->cdf, 0, (size_t)n_alleles * sizeof *p->cdf);
    for (int s = 0; s < p->n_targ_samples; ++s) {
        int h1 = s << 1;
        int h2 = h1 | 1;
        int a1 = fpd_targ_allele(fpd, m, h1);
        int a2 = fpd_targ_allele(fpd, m, h2);
        alleles[h1] = a1;
        alleles[h2] = a2;
        unph_het[s] = m >= fpd->stage1_overlap && a1 >= 0 && a2 >= 0 && a1 != a2;
        missing[s] = a1 < 0 || a2 < 0;
        if (a1 >= 0) ++p->cdf[a1];
        if (a2 >= 0) ++p->cdf[a2];
    }
    for (int h = p->n_targ_haps; h < p->n_haps; ++h) {
        int a = fpd_ref_allele(fpd, m, h - p->n_targ_haps);
        alleles[h] = a;
        ++p->cdf[a];
    }
    for (int j = 1; j < n_alleles; ++j) p->cdf[j] += p->cdf[j - 1];
    return p->cdf;
}

const int *pbwt_rec_phaser_phase(pbwt_rec_phaser *p, int current_mkr, int *alleles, int next_mkr, bool *missing, bool *unph_het) {
    if (current_mkr != -1) pbwt_updater_update(&p->pbwt, alleles, fpd_n_alleles(p->fpd, current_mkr), p->a);
    const int *cdf = set_alleles(p, next_mkr, alleles, unph_het, missing);
    if (next_mkr >= p->fpd->stage1_overlap) phase(p, alleles, unph_het);
    return cdf;
}
