/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/PbwtUpdater.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "beagleutil/pbwt_updater.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"

void pbwt_updater_init(pbwt_updater *u, int n_haps) {
    if (n_haps < 0) util_exit("java.lang.IllegalArgumentException: %d", n_haps);
    u->n_haps = n_haps;
    u->scratch = util_malloc((size_t)(n_haps > 0 ? n_haps : 1) * sizeof *u->scratch);
    u->counts_cap = 4;
    u->counts = util_malloc((size_t)u->counts_cap * sizeof *u->counts);
}

void pbwt_updater_update(pbwt_updater *u, const int *alleles, int n_alleles, int *prefix) {
    if (n_alleles < 1) util_exit("java.lang.IllegalArgumentException: %d", n_alleles);
    if (n_alleles + 1 > u->counts_cap) {
        u->counts_cap = n_alleles + 1;
        u->counts = util_realloc(u->counts, (size_t)u->counts_cap * sizeof *u->counts);
    }
    memset(u->counts, 0, (size_t)(n_alleles + 1) * sizeof *u->counts);
    for (int j = 0; j < u->n_haps; ++j) {
        int allele = alleles[prefix[j]];
        if (allele >= n_alleles) util_exit("java.lang.IndexOutOfBoundsException: %d", allele);
        ++u->counts[allele + 1];
    }
    for (int a = 1; a <= n_alleles; ++a) u->counts[a] += u->counts[a - 1];
    for (int j = 0; j < u->n_haps; ++j) u->scratch[u->counts[alleles[prefix[j]]]++] = prefix[j];
    memcpy(prefix, u->scratch, (size_t)u->n_haps * sizeof *prefix);
}

void pbwt_updater_free(pbwt_updater *u) {
    free(u->scratch);
    free(u->counts);
}
