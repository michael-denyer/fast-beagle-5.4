/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/SamplePhase.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/sample_phase.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/bit_array.h"
#include "blbutil/int_list.h"
#include "blbutil/utilities.h"

static const float MAX_CLUSTER_CM = 0.005f;

/* Markers.allelesToBits */
static uint64_t *alleles_to_bits(const fixed_phase_data *fpd, const int *alleles) {
    const int *hap_bits = fpd->stage1_hap_bits;
    uint64_t *bits = bit_array_new((size_t)hap_bits[fpd->n_stage1]);
    for (int m = 0; m < fpd->n_stage1; ++m) bit_array_set_allele(bits, hap_bits, m, alleles[m]);
    return bits;
}

static int *copy_ints(const int *a, int n) {
    int *copy = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *copy);
    if (n > 0) memcpy(copy, a, (size_t)n * sizeof *copy);
    return copy;
}

/* SamplePhase.clustSize */
static void set_clusters(sample_phase *sp, const int *hap1, const int *hap2) {
    const double *gen_pos = sp->fpd->stage1_map.gen_pos;
    int n_markers = sp->fpd->n_stage1;
    int_list sizes = {0};
    double max_clust_end = gen_pos[0] + (double)MAX_CLUSTER_CM;
    bool prev_is_missing_or_het = false;
    int last_end = 0;
    int miss_index = 0;
    int next_miss = miss_index < sp->n_missing ? sp->missing[miss_index++] : -1;
    for (int m = 0; m < n_markers; ++m) {
        int size = m - last_end;
        bool is_missing = m == next_miss;
        if (is_missing) next_miss = miss_index < sp->n_missing ? sp->missing[miss_index++] : -1;
        bool is_missing_or_het = is_missing || hap1[m] != hap2[m];
        if (prev_is_missing_or_het || is_missing_or_het || gen_pos[m] > max_clust_end || size == 255) {
            if (m > 0) {
                int_list_add(&sizes, size);
                max_clust_end = gen_pos[m] + (double)MAX_CLUSTER_CM;
                last_end = m;
            }
        }
        prev_is_missing_or_het = is_missing_or_het;
    }
    int_list_add(&sizes, n_markers - last_end);

    sp->n_clusters = sizes.n;
    sp->clust_size = util_malloc((size_t)sizes.n * sizeof *sp->clust_size);
    for (int c = 0; c < sizes.n; ++c) sp->clust_size[c] = (uint8_t)sizes.v[c];
    free(sizes.v);
}

void sample_phase_init(sample_phase *sp, int sample, const fixed_phase_data *fpd, const int *hap1, const int *hap2,
        const int *unphased, int n_unphased, const int *missing, int n_missing) {
    sp->sample = sample;
    sp->fpd = fpd;
    sp->hap1 = alleles_to_bits(fpd, hap1);
    sp->hap2 = alleles_to_bits(fpd, hap2);
    sp->n_unphased = n_unphased;
    sp->unphased = copy_ints(unphased, n_unphased);
    sp->n_missing = n_missing;
    sp->missing = copy_ints(missing, n_missing);
    set_clusters(sp, hap1, hap2);
}

int sample_phase_allele1(const sample_phase *sp, int m) {
    return bit_array_allele(sp->hap1, sp->fpd->stage1_hap_bits, m);
}

int sample_phase_allele2(const sample_phase *sp, int m) {
    return bit_array_allele(sp->hap2, sp->fpd->stage1_hap_bits, m);
}

void sample_phase_free(sample_phase *sp) {
    free(sp->hap1);
    free(sp->hap2);
    free(sp->unphased);
    free(sp->missing);
    free(sp->clust_size);
}

void sample_phase_set_allele1(sample_phase *sp, int m, int allele) {
    bit_array_set_allele(sp->hap1, sp->fpd->stage1_hap_bits, m, allele);
}

void sample_phase_set_allele2(sample_phase *sp, int m, int allele) {
    bit_array_set_allele(sp->hap2, sp->fpd->stage1_hap_bits, m, allele);
}

void sample_phase_swap_haps(sample_phase *sp, int start, int end) {
    const int *hap_bits = sp->fpd->stage1_hap_bits;
    for (int b = hap_bits[start]; b < hap_bits[end]; ++b) {
        uint64_t mask = (uint64_t)1 << (b & 63);
        uint64_t diff = (sp->hap1[b >> 6] ^ sp->hap2[b >> 6]) & mask;
        sp->hap1[b >> 6] ^= diff;
        sp->hap2[b >> 6] ^= diff;
    }
}

void sample_phase_clust_ends(const sample_phase *sp, int *ends) {
    int cum_sum = 0;
    for (int c = 0; c < sp->n_clusters; ++c) {
        cum_sum += sp->clust_size[c];
        ends[c] = cum_sum;
    }
}
