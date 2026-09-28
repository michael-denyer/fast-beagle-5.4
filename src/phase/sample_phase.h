/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/SamplePhase.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_SAMPLE_PHASE_H
#define PHASE_SAMPLE_PHASE_H

#include <stdint.h>

#include "phase/fixed_phase_data.h"

/* One target sample's current haplotypes at the stage-1 markers, its missing
 * genotypes, its unphased heterozygotes (those whose phase relative to the
 * previous heterozygote is still open), and the partition of the markers into
 * clusters. A missing genotype or a heterozygote is a cluster of its own; runs
 * of homozygotes are grouped into clusters of at most 255 markers spanning at
 * most 0.005 cM. The marker lists are increasing. */
typedef struct {
    int sample;
    const fixed_phase_data *fpd;
    uint64_t *hap1;
    uint64_t *hap2;
    int n_unphased;
    int *unphased;
    int n_missing;
    int *missing;
    int n_clusters;
    uint8_t *clust_size;
} sample_phase;

/* new SamplePhase(markers, genPos, hap1, hap2, unphased, missing). */
void sample_phase_init(sample_phase *sp, int sample, const fixed_phase_data *fpd, const int *hap1, const int *hap2,
        const int *unphased, int n_unphased, const int *missing, int n_missing);
int sample_phase_allele1(const sample_phase *sp, int m);
int sample_phase_allele2(const sample_phase *sp, int m);
void sample_phase_set_allele1(sample_phase *sp, int m, int allele);
void sample_phase_set_allele2(sample_phase *sp, int m, int allele);
/* Swaps the two haplotypes' alleles at stage-1 markers [start, end). */
void sample_phase_swap_haps(sample_phase *sp, int start, int end);
/* Each cluster's exclusive end marker, n_clusters entries. */
void sample_phase_clust_ends(const sample_phase *sp, int *ends);
void sample_phase_free(sample_phase *sp);

#endif
