/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/PbwtDivUpdater.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "beagleutil/pbwt_div_updater.h"

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"

void pbwt_div_updater_init(pbwt_div_updater *u, int n_haps) {
    size_t n = (size_t)(n_haps > 0 ? n_haps : 1);
    u->n_haps = n_haps;
    u->p_cap = 4;
    u->p = util_malloc((size_t)u->p_cap * sizeof *u->p);
    u->counts = util_malloc(((size_t)u->p_cap + 1) * sizeof *u->counts);
    u->allele = util_malloc(n * sizeof *u->allele);
    u->div = util_malloc(n * sizeof *u->div);
    u->prefix = util_malloc(n * sizeof *u->prefix);
}

void pbwt_div_updater_free(pbwt_div_updater *u) {
    free(u->p);
    free(u->counts);
    free(u->allele);
    free(u->div);
    free(u->prefix);
}

static void init_p(pbwt_div_updater *u, int n_alleles, int init_value) {
    if (n_alleles < 1) util_exit("java.lang.IllegalArgumentException: %d", n_alleles);
    if (n_alleles > u->p_cap) {
        u->p_cap = n_alleles;
        u->p = util_realloc(u->p, (size_t)u->p_cap * sizeof *u->p);
        u->counts = util_realloc(u->counts, ((size_t)u->p_cap + 1) * sizeof *u->counts);
    }
    for (int j = 0; j < n_alleles; ++j) u->p[j] = init_value;
}

/* Java's per-allele lists: each haplotype and its new divergence go to its
 * allele's list in prefix order, then the lists are concatenated. */
static void regroup(pbwt_div_updater *u, int n_alleles, int *prefix, int *div) {
    memset(u->counts, 0, ((size_t)n_alleles + 1) * sizeof *u->counts);
    for (int i = 0; i < u->n_haps; ++i) ++u->counts[u->allele[i] + 1];
    for (int a = 1; a <= n_alleles; ++a) u->counts[a] += u->counts[a - 1];
    for (int i = 0; i < u->n_haps; ++i) {
        int k = u->counts[u->allele[i]]++;
        u->prefix[k] = prefix[i];
        div[k] = u->div[i];
    }
    memcpy(prefix, u->prefix, (size_t)u->n_haps * sizeof *prefix);
}

static void update(pbwt_div_updater *u, const int *rec, int n_alleles, int init_value, bool fwd, int *prefix, int *div) {
    init_p(u, n_alleles, init_value);
    for (int i = 0; i < u->n_haps; ++i) {
        int allele = rec[prefix[i]];
        if (allele >= n_alleles) util_exit("java.lang.IndexOutOfBoundsException: %d", n_alleles);
        for (int j = 0; j < n_alleles; ++j) {
            if (fwd ? div[i] > u->p[j] : div[i] < u->p[j]) u->p[j] = div[i];
        }
        u->allele[i] = allele;
        u->div[i] = u->p[allele];
        u->p[allele] = fwd ? INT_MIN : INT_MAX;
    }
    regroup(u, n_alleles, prefix, div);
}

void pbwt_div_updater_fwd(pbwt_div_updater *u, const int *rec, int n_alleles, int marker, int *prefix, int *div) {
    update(u, rec, n_alleles, marker + 1, true, prefix, div);
}

void pbwt_div_updater_bwd(pbwt_div_updater *u, const int *rec, int n_alleles, int marker, int *prefix, int *div) {
    update(u, rec, n_alleles, marker - 1, false, prefix, div);
}
