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
    pb->rows = NULL;
    pb->rows_cap = 0;
    pb->zero_row = util_malloc((size_t)pb->max_states);
    memset(pb->zero_row, 0, (size_t)pb->max_states);
    for (int i = 0; i < 3; ++i) {
        pb->mismatch[i] = util_malloc((size_t)n_markers * sizeof *pb->mismatch[i]);
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
    for (int i = 0; i < 3; ++i) {
        free(pb->mismatch[i]);
        free(pb->fwd[i]);
        free(pb->bwd[i]);
    }
    free(pb->rows);
    free(pb->zero_row);
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

static void ensure_capacity(phase_baum1 *pb, int n_clusters, int n_unph, int n_miss) {
    size_t row = (size_t)pb->max_states;
    size_t n_rows_bytes = 2 * (size_t)n_clusters * row;
    if (pb->rows_cap < n_rows_bytes) {
        free(pb->rows);
        pb->rows = util_malloc(n_rows_bytes);
        pb->rows_cap = n_rows_bytes;
    }
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

/* The backward pass over a sample's clusters, held as a position so that
 * two samples can take their steps together. It stores the haplotype HMMs'
 * backward values at each missing genotype and at the cluster just before
 * each unphased heterozygote, where the haplotype HMMs restart from the
 * homozygous one. */
typedef struct {
    phase_baum1 *pb;
    const marker_cluster *mc;
    int c;   /* the cluster the next step reaches, from n_clusters - 2 down to 0 */
    int unph_index;
} bwd_pass;

static void bwd_begin(bwd_pass *bp, phase_baum1 *pb, const marker_cluster *mc) {
    *bp = (bwd_pass){pb, mc, mc->n_clusters - 2, mc->n_unph - 1};
    for (int k = 0; k < pb->n_states; ++k) pb->bwd[0][k] = 1.0f / pb->n_states;
    if (mc->has_missing[mc->n_clusters - 1]) {
        --pb->miss_index;
        copy_states(pb, pb->miss_probs1[pb->miss_index], pb->bwd[0]);
        copy_states(pb, pb->miss_probs2[pb->miss_index], pb->bwd[0]);
    }
    copy_states(pb, pb->bwd[1], pb->bwd[0]);
    copy_states(pb, pb->bwd[2], pb->bwd[0]);
}

/* The update that takes the backward values from cluster c + 1 to cluster c.
 * It takes the size of cluster c with the mismatches of cluster c + 1, as
 * Beagle 5.4 does. */
static hmm3_step bwd_step(bwd_pass *bp) {
    phase_baum1 *pb = bp->pb;
    int c_p1 = bp->c + 1;
    set_em_probs(pb, cluster_size(bp->mc, bp->c));
    return (hmm3_step){pb->bwd, NULL, bp->mc->p_recomb[c_p1], pb->em_probs,
            {pb->mismatch[0][c_p1], pb->mismatch[1][c_p1], pb->mismatch[2][c_p1]}, pb->n_states};
}

static void bwd_end_step(bwd_pass *bp) {
    phase_baum1 *pb = bp->pb;
    const marker_cluster *mc = bp->mc;
    int c = bp->c;
    if (mc->has_missing[c]) {
        --pb->miss_index;
        copy_states(pb, pb->miss_probs1[pb->miss_index], pb->bwd[1]);
        copy_states(pb, pb->miss_probs2[pb->miss_index], pb->bwd[2]);
    }
    int j = bp->unph_index;
    if (j >= 0 && c + 1 == mc->unph_clusters[j]) {
        copy_states(pb, pb->bwd_het1[j], pb->bwd[1]);
        copy_states(pb, pb->bwd_het2[j], pb->bwd[2]);
        copy_states(pb, pb->bwd[1], pb->bwd[0]);
        copy_states(pb, pb->bwd[2], pb->bwd[0]);
        --bp->unph_index;
    }
    --bp->c;
}

static void bwd_finish(bwd_pass *bp) {
    while (bp->c >= 0) {
        hmm3_step s = bwd_step(bp);
        hmm_bwd_update3(s.val, s.p_switch, s.p_mismatch, s.m[0], s.m[1], s.m[2], s.n_states);
        bwd_end_step(bp);
    }
}

static void swap_haps(phase_baum1 *pb, const marker_cluster *mc, int start_clust, int end_clust) {
    for (int c = start_clust; c < end_clust; ++c) {
        const uint8_t *tmp = pb->mismatch[1][c];
        pb->mismatch[1][c] = pb->mismatch[2][c];
        pb->mismatch[2][c] = tmp;
    }
    sample_phase_swap_haps(mc->sp, marker_cluster_start(mc, start_clust), mc->ends[end_clust - 1]);
}

static void impute_alleles(phase_baum1 *pb, const marker_cluster *mc, int cluster) {
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
    sample_phase_set_allele1(mc->sp, marker, a1);
    sample_phase_set_allele2(mc->sp, marker, a2);
    free(al_freq1);
    free(al_freq2);
    ++pb->miss_index;
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

/* The forward pass, held as a position in the same way. */
typedef struct {
    phase_baum1 *pb;
    const marker_cluster *mc;
    int c;   /* the cluster the next step reaches */
    int unph_index;
} fwd_pass;

static void fwd_begin(fwd_pass *fp, phase_baum1 *pb, const marker_cluster *mc) {
    *fp = (fwd_pass){pb, mc, 0, 0};
    for (int k = 0; k < pb->n_states; ++k) pb->fwd[0][k] = 1.0f / pb->n_states;
    pb->sum[0] = 1.0f;
    copy_states(pb, pb->fwd[1], pb->fwd[0]);
    copy_states(pb, pb->fwd[2], pb->fwd[0]);
    pb->sum[1] = pb->sum[2] = pb->sum[0];
}

/* At an unphased heterozygote's cluster, phases it, applies the phase switch
 * to the clusters up to the next one, and restarts the haplotype HMMs from the
 * homozygous one. Then gives the update that takes the forward values to
 * cluster c. */
static hmm3_step fwd_step(fwd_pass *fp) {
    phase_baum1 *pb = fp->pb;
    const marker_cluster *mc = fp->mc;
    int c = fp->c;
    if (fp->unph_index < mc->n_unph && c == mc->unph_clusters[fp->unph_index]) {
        phase_het(pb, fp->unph_index);
        ++fp->unph_index;
        if (pb->swap_haps) {
            int swap_end = fp->unph_index < mc->n_unph ? mc->unph_clusters[fp->unph_index] : mc->n_clusters;
            swap_haps(pb, mc, c, swap_end);
        }
        copy_states(pb, pb->fwd[1], pb->fwd[0]);
        copy_states(pb, pb->fwd[2], pb->fwd[0]);
        pb->sum[1] = pb->sum[2] = pb->sum[0];
    }
    set_em_probs(pb, cluster_size(mc, c));
    return (hmm3_step){pb->fwd, pb->sum, mc->p_recomb[c], pb->em_probs,
            {pb->mismatch[0][c], pb->mismatch[1][c], pb->mismatch[2][c]}, pb->n_states};
}

static void fwd_end_step(fwd_pass *fp) {
    if (fp->mc->has_missing[fp->c]) impute_alleles(fp->pb, fp->mc, fp->c);
    ++fp->c;
}

static void fwd_finish(fwd_pass *fp) {
    while (fp->c < fp->mc->n_clusters) {
        hmm3_step s = fwd_step(fp);
        hmm_fwd_update3(s.val, s.sums, s.p_switch, s.p_mismatch, s.m[0], s.m[1], s.m[2], s.n_states);
        fwd_end_step(fp);
    }
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

/* Clusters a sample and builds its states. False if the sample has nothing
 * to phase or impute. */
static bool begin_sample(phase_baum1 *pb, int sample, marker_cluster *mc, kstring_t *s) {
    sample_phase *sp = &pb->pd->phase[sample];
    pb->swap_haps = false;
    pb->swap_cnt = 0;
    int n_unph = sp->n_unphased;
    int n_miss = sp->n_missing;
    if (pb->trace) ksprintf(s, "sample\t%d\t%d\t%d", sample, n_unph, n_miss);
    if (n_miss == 0 && n_unph == 0) return false;
    pb->n_lr = 0;
    marker_cluster_init(mc, pb->pd, sample);
    ensure_capacity(pb, mc->n_clusters, n_unph, n_miss);
    pb->miss_index = mc->n_missing_clusters;
    pb->n_states = basic_phase_states_cluster_states(&pb->states, mc, pb->ref_alleles, pb->mismatch, pb->rows, pb->zero_row);
    if (pb->trace) trace_states(pb, mc, s);
    pb->em_digest = TRACE_FNV_BASIS;
    pb->threshold = NAN;
    return true;
}

static void end_sample(phase_baum1 *pb, int sample, marker_cluster *mc, bool ran, kstring_t *s, swap_rate *rate) {
    sample_phase *sp = &pb->pd->phase[sample];
    int n_unph = sp->n_unphased;   /* before update_phase drops the committed heterozygotes */
    if (ran) {
        update_phase(pb, sample, sp);
        if (pb->trace) trace_update(pb, sp, s);
        marker_cluster_free(mc);
    }
    if (pb->trace) pb->trace_lines[sample] = s->s;
    rate->n_unph_hets += n_unph;
    rate->n_swaps += pb->swap_cnt;
}

/* Each sample is phased as PhaseBaum1.phase phases it. The samples are
 * independent, so their backward passes, and then their forward passes, take
 * their steps together until the shorter one ends. */
void phase_baum1_phase_pair(phase_baum1 *pb0, int sample0, phase_baum1 *pb1, int sample1, swap_rate *rate) {
    marker_cluster mc0 = {0}, mc1 = {0};
    kstring_t s0 = {0, 0, NULL}, s1 = {0, 0, NULL};
    bool run0 = begin_sample(pb0, sample0, &mc0, &s0);
    bool run1 = sample1 >= 0 && begin_sample(pb1, sample1, &mc1, &s1);
    bwd_pass b0 = {0}, b1 = {0};
    if (run0) bwd_begin(&b0, pb0, &mc0);
    if (run1) bwd_begin(&b1, pb1, &mc1);
    while (run0 && run1 && b0.c >= 0 && b1.c >= 0) {
        hmm3_step t0 = bwd_step(&b0), t1 = bwd_step(&b1);
        hmm_bwd_update3x2(&t0, &t1);
        bwd_end_step(&b0);
        bwd_end_step(&b1);
    }
    if (run0) bwd_finish(&b0);
    if (run1) bwd_finish(&b1);
    fwd_pass f0 = {0}, f1 = {0};
    if (run0) fwd_begin(&f0, pb0, &mc0);
    if (run1) fwd_begin(&f1, pb1, &mc1);
    while (run0 && run1 && f0.c < mc0.n_clusters && f1.c < mc1.n_clusters) {
        hmm3_step t0 = fwd_step(&f0), t1 = fwd_step(&f1);
        hmm_fwd_update3x2(&t0, &t1);
        fwd_end_step(&f0);
        fwd_end_step(&f1);
    }
    if (run0) fwd_finish(&f0);
    if (run1) fwd_finish(&f1);
    end_sample(pb0, sample0, &mc0, run0, &s0, rate);
    if (sample1 >= 0) end_sample(pb1, sample1, &mc1, run1, &s1, rate);
}
