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
#ifndef BEAGLEUTIL_PBWT_UPDATER_H
#define BEAGLEUTIL_PBWT_UPDATER_H

/* Advances a positional Burrows-Wheeler transform prefix array by one marker:
 * haplotypes are regrouped by allele, keeping their prefix order within each
 * allele. */
typedef struct {
    int n_haps;
    int *scratch;
    int *counts;
    int counts_cap;
} pbwt_updater;

void pbwt_updater_init(pbwt_updater *u, int n_haps);
/* PbwtUpdater.update(alleles, nAlleles, prefix): alleles[h] is haplotype h's
 * allele, in [0, n_alleles). */
void pbwt_updater_update(pbwt_updater *u, const int *alleles, int n_alleles, int *prefix);
void pbwt_updater_free(pbwt_updater *u);

#endif
