/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PhaseBaum1.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/phase_baum1.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "phase/hmm_updater.h"
#include "phase/marker_cluster.h"

void phase_baum1_init(phase_baum1 *pb, const pbwt_phase_ibs *ibs, char **trace_lines) {
    phase_data *pd = ibs->pd;
    const par *p = pd->par;
    int n_markers = pd->fpd->n_stage1;
    pb->pd = pd;
    pb->burnin = pd->it < p->burnin;
    pb->trace = trace_lines != NULL;
    pb->trace_lines = trace_lines;
    pb->max_states = p->phase_states;
    basic_phase_states_init(&pb->states, ibs, pb->max_states);
    pb->n_states = 0;
    pb->mismatch = util_malloc(3 * sizeof *pb->mismatch);
    for (int i = 0; i < 3; ++i) {
        pb->mismatch[i] = util_malloc((size_t)n_markers * sizeof **pb->mismatch);
        for (int m = 0; m < n_markers; ++m) pb->mismatch[i][m] = util_malloc((size_t)pb->max_states);
        pb->fwd[i] = util_malloc((size_t)pb->max_states * sizeof *pb->fwd[i]);
        pb->bwd[i] = util_malloc((size_t)pb->max_states * sizeof *pb->bwd[i]);
    }
    pb->p_mismatch = pd->p_mismatch;
    pb->em_probs[0] = 1.0f - pb->p_mismatch;
    pb->em_probs[1] = pb->p_mismatch;
    pb->n_miss_cap = pb->n_het_cap = 0;
    pb->ref_alleles = NULL;
    pb->miss_probs1 = pb->miss_probs2 = pb->bwd_het1 = pb->bwd_het2 = NULL;
    pb->lr = pb->lr_sorted = NULL;
    pb->n_lr = 0;
    pb->swap_haps = false;
    pb->miss_index = -1;
    pb->swap_cnt = 0;
}

void phase_baum1_free(phase_baum1 *pb) {
    int n_markers = pb->pd->fpd->n_stage1;
    for (int i = 0; i < 3; ++i) {
        for (int m = 0; m < n_markers; ++m) free(pb->mismatch[i][m]);
        free(pb->mismatch[i]);
        free(pb->fwd[i]);
        free(pb->bwd[i]);
    }
    free(pb->mismatch);
    for (int j = 0; j < pb->n_miss_cap; ++j) {
        free(pb->ref_alleles[j]);
        free(pb->miss_probs1[j]);
        free(pb->miss_probs2[j]);
    }
    for (int j = 0; j < pb->n_het_cap; ++j) {
        free(pb->bwd_het1[j]);
        free(pb->bwd_het2[j]);
    }
    free(pb->ref_alleles);
    free(pb->miss_probs1);
    free(pb->miss_probs2);
    free(pb->bwd_het1);
    free(pb->bwd_het2);
    free(pb->lr);
    free(pb->lr_sorted);
    basic_phase_states_free(&pb->states);
}

static void ensure_capacity(phase_baum1 *pb, int n_unph, int n_miss) {
    size_t row = (size_t)pb->max_states;
    if (pb->n_miss_cap < n_miss) {
        pb->ref_alleles = util_realloc(pb->ref_alleles, (size_t)n_miss * sizeof *pb->ref_alleles);
        pb->miss_probs1 = util_realloc(pb->miss_probs1, (size_t)n_miss * sizeof *pb->miss_probs1);
        pb->miss_probs2 = util_realloc(pb->miss_probs2, (size_t)n_miss * sizeof *pb->miss_probs2);
        for (int j = pb->n_miss_cap; j < n_miss; ++j) {
            pb->ref_alleles[j] = util_malloc(row * sizeof **pb->ref_alleles);
            pb->miss_probs1[j] = util_malloc(row * sizeof **pb->miss_probs1);
            pb->miss_probs2[j] = util_malloc(row * sizeof **pb->miss_probs2);
        }
        pb->n_miss_cap = n_miss;
    }
    if (pb->n_het_cap < n_unph) {
        pb->bwd_het1 = util_realloc(pb->bwd_het1, (size_t)n_unph * sizeof *pb->bwd_het1);
        pb->bwd_het2 = util_realloc(pb->bwd_het2, (size_t)n_unph * sizeof *pb->bwd_het2);
        for (int j = pb->n_het_cap; j < n_unph; ++j) {
            pb->bwd_het1[j] = util_malloc(row * sizeof **pb->bwd_het1);
            pb->bwd_het2[j] = util_malloc(row * sizeof **pb->bwd_het2);
        }
        pb->lr = util_realloc(pb->lr, (size_t)n_unph * sizeof *pb->lr);
        pb->lr_sorted = util_realloc(pb->lr_sorted, (size_t)n_unph * sizeof *pb->lr_sorted);
        pb->n_het_cap = n_unph;
    }
}

static void copy_states(const phase_baum1 *pb, float *dst, const float *src) {
    memcpy(dst, src, (size_t)pb->n_states * sizeof *dst);
}

/* The emission probabilities for a cluster of n markers. Beagle 5.4 does not
 * cap them, so a long cluster can have a negative match probability. */
static void set_em_probs(phase_baum1 *pb, int n_markers_in_cluster) {
    pb->em_probs[1] = n_markers_in_cluster * pb->p_mismatch;
    pb->em_probs[0] = 1.0f - pb->em_probs[1];
    if (pb->trace) {
        pb->em_digest = trace_fold(pb->em_digest, jnum_float_bits(pb->em_probs[0]));
        pb->em_digest = trace_fold(pb->em_digest, jnum_float_bits(pb->em_probs[1]));
    }
}

static int cluster_size(const marker_cluster *mc, int c) {
    return mc->ends[c] - marker_cluster_start(mc, c);
}

/* Backward over clusters [start_clust, end_clust), the haplotype HMMs
 * restarting from the homozygous one. The step from cluster c + 1 to c takes
 * the size of cluster c with the mismatches of cluster c + 1, as Beagle 5.4
 * does. */
static void bwd_segment(phase_baum1 *pb, const marker_cluster *mc, int start_clust, int end_clust) {
    copy_states(pb, pb->bwd[1], pb->bwd[0]);
    copy_states(pb, pb->bwd[2], pb->bwd[0]);
    for (int c = end_clust - 1; c >= start_clust; --c) {
        int c_p1 = c + 1;
        float p_rec = mc->p_recomb[c_p1];
        set_em_probs(pb, cluster_size(mc, c));
        hmm_bwd_update3(pb->bwd, p_rec, pb->em_probs, pb->mismatch[0][c_p1], pb->mismatch[1][c_p1], pb->mismatch[2][c_p1], pb->n_states);
        if (mc->has_missing[c]) {
            --pb->miss_index;
            copy_states(pb, pb->miss_probs1[pb->miss_index], pb->bwd[1]);
            copy_states(pb, pb->miss_probs2[pb->miss_index], pb->bwd[2]);
        }
    }
}

/* Stores the haplotype HMMs' backward values after each unphased
 * heterozygote's segment, and at each missing genotype. */
static void bwd_alg(phase_baum1 *pb, const marker_cluster *mc) {
    for (int k = 0; k < pb->n_states; ++k) pb->bwd[0][k] = 1.0f / pb->n_states;
    int end_clust = mc->n_clusters - 1;
    if (mc->has_missing[end_clust]) {
        --pb->miss_index;
        copy_states(pb, pb->miss_probs1[pb->miss_index], pb->bwd[0]);
        copy_states(pb, pb->miss_probs2[pb->miss_index], pb->bwd[0]);
    }
    for (int j = mc->n_unph - 1; j >= 0; --j) {
        int start_clust = mc->unph_clusters[j] - 1;
        bwd_segment(pb, mc, start_clust, end_clust);
        copy_states(pb, pb->bwd_het1[j], pb->bwd[1]);
        copy_states(pb, pb->bwd_het2[j], pb->bwd[2]);
        end_clust = start_clust;
    }
    bwd_segment(pb, mc, 0, end_clust);
}

static void swap_haps(phase_baum1 *pb, sample_phase *sp, const marker_cluster *mc, int start_clust, int end_clust) {
    for (int c = start_clust; c < end_clust; ++c) {
        uint8_t *tmp = pb->mismatch[1][c];
        pb->mismatch[1][c] = pb->mismatch[2][c];
        pb->mismatch[2][c] = tmp;
    }
    sample_phase_swap_haps(sp, marker_cluster_start(mc, start_clust), mc->ends[end_clust - 1]);
}

static void impute_alleles(phase_baum1 *pb, sample_phase *sp, const marker_cluster *mc, int cluster) {
    int mi = pb->miss_index;
    if (pb->swap_haps) {
        float *tmp = pb->miss_probs1[mi];
        pb->miss_probs1[mi] = pb->miss_probs2[mi];
        pb->miss_probs2[mi] = tmp;
    }
    float *state_probs1 = pb->miss_probs1[mi];
    float *state_probs2 = pb->miss_probs2[mi];
    const int *ref_al = pb->ref_alleles[mi];
    for (int k = 0; k < pb->n_states; ++k) {
        state_probs1[k] *= pb->fwd[1][k];
        state_probs2[k] *= pb->fwd[2][k];
    }
    int marker = marker_cluster_start(mc, cluster);
    int n_alleles = fpd_n_alleles(pb->pd->fpd, marker);
    float *al_freq1 = util_malloc((size_t)n_alleles * sizeof *al_freq1);
    float *al_freq2 = util_malloc((size_t)n_alleles * sizeof *al_freq2);
    for (int a = 0; a < n_alleles; ++a) al_freq1[a] = al_freq2[a] = 0.0f;
    for (int k = 0; k < pb->n_states; ++k) {
        al_freq1[ref_al[k]] += state_probs1[k];
        al_freq2[ref_al[k]] += state_probs2[k];
    }
    int a1 = 0;
    int a2 = 0;
    for (int j = 1; j < n_alleles; ++j) {
        if (al_freq1[j] > al_freq1[a1]) a1 = j;
        if (al_freq2[j] > al_freq2[a2]) a2 = j;
    }
    sample_phase_set_allele1(sp, marker, a1);
    sample_phase_set_allele2(sp, marker, a2);
    free(al_freq1);
    free(al_freq2);
    ++pb->miss_index;
}

/* Forward over clusters [start_clust, end_clust), first applying the phase
 * switch decided at the previous heterozygote. */
static void fwd_segment(phase_baum1 *pb, sample_phase *sp, const marker_cluster *mc, int start_clust, int end_clust) {
    if (pb->swap_haps) swap_haps(pb, sp, mc, start_clust, end_clust);
    copy_states(pb, pb->fwd[1], pb->fwd[0]);
    copy_states(pb, pb->fwd[2], pb->fwd[0]);
    pb->sum[1] = pb->sum[2] = pb->sum[0];
    for (int c = start_clust; c < end_clust; ++c) {
        float p_rec = mc->p_recomb[c];
        set_em_probs(pb, cluster_size(mc, c));
        hmm_fwd_update3(pb->fwd, pb->sum, p_rec, pb->em_probs, pb->mismatch[0][c], pb->mismatch[1][c], pb->mismatch[2][c], pb->n_states);
        if (mc->has_missing[c]) impute_alleles(pb, sp, mc, c);
    }
}

static void phase_het(phase_baum1 *pb, int het_index) {
    const float *b1 = pb->bwd_het1[het_index];
    const float *b2 = pb->bwd_het2[het_index];
    float p11 = 0.0f, p12 = 0.0f, p21 = 0.0f, p22 = 0.0f;
    for (int k = 0; k < pb->n_states; ++k) {
        p11 += pb->fwd[1][k] * b1[k];
        p12 += pb->fwd[1][k] * b2[k];
        p21 += pb->fwd[2][k] * b1[k];
        p22 += pb->fwd[2][k] * b2[k];
    }
    bool last_swap_haps = pb->swap_haps;
    float num = p11 * p22;
    float den = p12 * p21;
    pb->swap_haps = num < den;
    if (pb->swap_haps != last_swap_haps) ++pb->swap_cnt;
    pb->lr[pb->n_lr++] = pb->swap_haps ? den / num : num / den;
}

static void fwd_alg(phase_baum1 *pb, sample_phase *sp, const marker_cluster *mc) {
    for (int k = 0; k < pb->n_states; ++k) pb->fwd[0][k] = 1.0f / pb->n_states;
    pb->sum[0] = 1.0f;
    int start_clust = 0;
    for (int j = 0; j < mc->n_unph; ++j) {
        int end_clust = mc->unph_clusters[j];
        fwd_segment(pb, sp, mc, start_clust, end_clust);
        phase_het(pb, j);
        start_clust = end_clust;
    }
    if (start_clust < mc->n_clusters) fwd_segment(pb, sp, mc, start_clust, mc->n_clusters);
}

/* Float.compare order, as Arrays.sort(float[]) leaves it: -0.0 before 0.0,
 * NaN last. */
static int float_compare(const void *pa, const void *pb) {
    float a = *(const float *)pa, b = *(const float *)pb;
    if (a < b) return -1;
    if (a > b) return 1;
    int32_t a_bits = isnan(a) ? 0x7fc00000 : (int32_t)jnum_float_bits(a);
    int32_t b_bits = isnan(b) ? 0x7fc00000 : (int32_t)jnum_float_bits(b);
    return a_bits == b_bits ? 0 : (a_bits < b_bits ? -1 : 1);
}

/* PhaseBaum1.threshold: the likelihood ratio at rank floor(p*n + 0.5),
 * or the largest. */
static float threshold(phase_baum1 *pb, float leave_unphased_prop) {
    int n = pb->n_lr;
    memcpy(pb->lr_sorted, pb->lr, (size_t)n * sizeof *pb->lr_sorted);
    qsort(pb->lr_sorted, (size_t)n, sizeof *pb->lr_sorted, float_compare);
    int rank = jnum_d2i(floor((double)(leave_unphased_prop * n + 0.5f)));
    return pb->lr_sorted[rank < n ? rank : n - 1];
}

/* After burn-in, keeps unphased only the heterozygotes whose likelihood
 * ratio is below the sample's threshold. */
static void update_phase(phase_baum1 *pb, int sample, sample_phase *sp) {
    if (sp->n_unphased > 0 && !pb->burnin) {
        float t = threshold(pb, pb->pd->leave_unph_prop[sample]);
        pb->threshold = t;
        int n_next = 0;
        for (int j = 0; j < sp->n_unphased; ++j) {
            if (pb->lr[j] < t) sp->unphased[n_next++] = sp->unphased[j];
        }
        sp->n_unphased = n_next;
    }
}

static void put_ints(kstring_t *s, const int *a, int n) {
    if (n == 0) kputc('-', s);
    for (int j = 0; j < n; ++j) {
        if (j > 0) kputc(',', s);
        kputw(a[j], s);
    }
}

/* T3e: the cluster count, state count, a digest of the state mismatches
 * (rows 0 to 2, each cluster, each state) then of the state alleles at each
 * missing-genotype cluster, and a digest of the cluster recombination
 * probabilities. */
static void trace_states(const phase_baum1 *pb, const marker_cluster *mc, kstring_t *s) {
    uint64_t h = TRACE_FNV_BASIS;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < mc->n_clusters; ++c) {
            for (int k = 0; k < pb->n_states; ++k) h = trace_fold(h, pb->mismatch[r][c][k]);
        }
    }
    for (int j = 0; j < mc->n_missing_clusters; ++j) {
        for (int k = 0; k < pb->n_states; ++k) h = trace_fold(h, (uint32_t)pb->ref_alleles[j][k]);
    }
    uint64_t r = TRACE_FNV_BASIS;
    for (int c = 0; c < mc->n_clusters; ++c) r = trace_fold(r, jnum_float_bits(mc->p_recomb[c]));
    ksprintf(s, "\t%d\t%d\t%" PRIx64 "\t%" PRIx64, mc->n_clusters, pb->n_states, h, r);
}

static void trace_update(const phase_baum1 *pb, const sample_phase *sp, kstring_t *s) {
    ksprintf(s, "\t%" PRIx64 "\t", pb->em_digest);
    if (pb->n_lr == 0) kputc('-', s);
    for (int j = 0; j < pb->n_lr; ++j) {
        if (j > 0) kputc(',', s);
        ksprintf(s, "%" PRIx32, jnum_float_bits(pb->lr[j]));
    }
    ksprintf(s, "\t%d\t", pb->swap_cnt);
    if (isnan(pb->threshold)) kputc('-', s);
    else ksprintf(s, "%" PRIx32, jnum_float_bits(pb->threshold));
    kputc('\t', s);
    put_ints(s, sp->unphased, sp->n_unphased);
}

void phase_baum1_phase(phase_baum1 *pb, int sample, swap_rate *rate) {
    sample_phase *sp = &pb->pd->phase[sample];
    pb->swap_haps = false;
    pb->swap_cnt = 0;
    int n_unph = sp->n_unphased;
    int n_miss = sp->n_missing;
    kstring_t s = {0, 0, NULL};
    if (pb->trace) ksprintf(&s, "sample\t%d\t%d\t%d", sample, n_unph, n_miss);
    if (n_miss > 0 || n_unph > 0) {
        pb->n_lr = 0;
        ensure_capacity(pb, n_unph, n_miss);
        marker_cluster mc;
        marker_cluster_init(&mc, pb->pd, sample);
        pb->miss_index = mc.n_missing_clusters;
        pb->n_states = basic_phase_states_cluster_states(&pb->states, &mc, pb->ref_alleles, pb->mismatch);
        if (pb->trace) trace_states(pb, &mc, &s);
        pb->em_digest = TRACE_FNV_BASIS;
        pb->threshold = NAN;
        bwd_alg(pb, &mc);
        fwd_alg(pb, sp, &mc);
        update_phase(pb, sample, sp);
        if (pb->trace) trace_update(pb, sp, &s);
        marker_cluster_free(&mc);
    }
    if (pb->trace) pb->trace_lines[sample] = s.s;
    rate->n_unph_hets += n_unph;
    rate->n_swaps += pb->swap_cnt;
}
