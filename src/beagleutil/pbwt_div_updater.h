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
#ifndef BEAGLEUTIL_PBWT_DIV_UPDATER_H
#define BEAGLEUTIL_PBWT_DIV_UPDATER_H

/* Advances a PBWT prefix array and its divergence array by one step, forward
 * (divergence is the start of each haplotype's match with its predecessor) or
 * backward (the end of the match). */
typedef struct {
    int n_haps;
    int *p;          /* per allele, the running divergence bound */
    int p_cap;
    int *allele;     /* per prefix position, scratch */
    int *div;        /* per prefix position, scratch */
    int *prefix;     /* scratch */
    int *counts;     /* per allele + 1, scratch */
} pbwt_div_updater;

void pbwt_div_updater_init(pbwt_div_updater *u, int n_haps);
/* PbwtDivUpdater.fwdUpdate(rec, nAlleles, marker, prefix, div): rec[h] is
 * haplotype h's allele in [0, n_alleles). div has n_haps + 1 entries; the
 * last is not changed. */
void pbwt_div_updater_fwd(pbwt_div_updater *u, const int *rec, int n_alleles, int marker, int *prefix, int *div);
/* PbwtDivUpdater.bwdUpdate(rec, nAlleles, marker, prefix, div) */
void pbwt_div_updater_bwd(pbwt_div_updater *u, const int *rec, int n_alleles, int marker, int *prefix, int *div);
void pbwt_div_updater_free(pbwt_div_updater *u);

#endif
