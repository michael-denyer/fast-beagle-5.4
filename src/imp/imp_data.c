/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) imp/ImpData.java and
 * imp/HaplotypeCoder.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "imp/imp_data.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/int_list.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jmath.h"
#include "jcompat/jnum.h"

static const double MIN_CM_DIST = 1e-7;

/* ImpData.cumPos: cM positions from 0, each step at least MIN_CM_DIST. */
static double *cum_pos(const genetic_map *gm, const marker *const *markers, int n) {
    double *cum = util_malloc((size_t)n * sizeof *cum);
    double last_gen_pos = genetic_map_gen_pos(gm, markers[0]->chrom_index, markers[0]->pos);
    cum[0] = 0.0;
    for (int j = 1; j < n; ++j) {
        double gen_pos = genetic_map_gen_pos(gm, markers[j]->chrom_index, markers[j]->pos);
        double gen_dist = fabs(gen_pos - last_gen_pos);
        if (!(gen_dist > MIN_CM_DIST)) gen_dist = MIN_CM_DIST;
        cum[j] = cum[j - 1] + gen_dist;
        last_gen_pos = gen_pos;
    }
    return cum;
}

static const ref_gt_rec *restricted_ref(const window *w, int targ_marker) {
    return w->ref[w->indices.targ_marker_to_marker[targ_marker]];
}

/* ImpData.targBlockEnd: target marker indices where the reference records'
 * SeqCoder3 group changes, then the marker count. */
static int_list targ_block_end(const window *w) {
    int_list ends = {0};
    const seq_group *last = NULL;
    for (int j = 0; j < w->n_targ; ++j) {
        const ref_gt_rec *rec = restricted_ref(w, j);
        if (rec->kind == REF_HAP && rec->group != last) {
            if (last != NULL) int_list_add(&ends, j);
            last = rec->group;
        }
    }
    int_list_add(&ends, w->n_targ);
    return ends;
}

static int_list targ_clust_start_end(const double *pos, const int_list *block_end, float cluster_dist) {
    int_list se = {0};
    int_list_add(&se, 0);
    for (int j = 0; j < block_end->n; ++j) {
        int clust_start = se.v[se.n - 1];
        int end = block_end->v[j];
        double start_pos = pos[clust_start];
        for (int m = clust_start + 1; m < end; ++m) {
            if (pos[m] - start_pos > (double)cluster_dist) {
                int_list_add(&se, m);
                start_pos = pos[m];
            }
        }
        int_list_add(&se, end);
    }
    return se;
}

/* HaplotypeCoder.codeTarg: numbers the target haplotypes' allele sequences
 * over [start, end) from 1 in order of first appearance. seq_map[j] maps a
 * (sequence before marker start + j, allele) pair to the next sequence. */
static int code_targ(const imp_data *id, int start, int end, int **seq_map, int *targ_seq) {
    const window *w = id->w;
    for (int h = 0; h < id->n_targ_haps; ++h) targ_seq[h] = 1;
    int seq_cnt = 2;
    for (int j = 0; j < end - start; ++j) {
        int m = start + j;
        int n_alleles = marker_n_alleles(&w->targ[m]->marker);
        size_t size = (size_t)seq_cnt * (size_t)n_alleles;
        seq_map[j] = util_malloc(size * sizeof **seq_map);
        memset(seq_map[j], 0, size * sizeof **seq_map);
        seq_cnt = 1;
        for (int h = 0; h < id->n_targ_haps; ++h) {
            int index = n_alleles * targ_seq[h] + id->targ_allele(id->targ_ctx, m, h);
            if (seq_map[j][index] == 0) seq_map[j][index] = seq_cnt++;
            targ_seq[h] = seq_map[j][index];
        }
    }
    return seq_cnt;
}

static bool is_hap_coded(const window *w, int start, int end) {
    const ref_gt_rec *first = restricted_ref(w, start);
    if (first->kind != REF_HAP) return false;
    for (int m = start + 1; m < end; ++m) {
        const ref_gt_rec *rec = restricted_ref(w, m);
        if (rec->kind != REF_HAP || rec->group != first->group) return false;
    }
    return true;
}

/* HaplotypeCoder.run: the reference haplotypes follow the target sequences'
 * transitions; one that leaves them is 0 from then on. */
static void code_cluster(const imp_data *id, int start, int end, hap_seq_map *hs) {
    const window *w = id->w;
    int n = end - start;
    int **seq_map = util_malloc((size_t)n * sizeof *seq_map);
    hs->targ_seq = util_malloc((size_t)id->n_targ_haps * sizeof *hs->targ_seq);
    hs->n_seq = code_targ(id, start, end, seq_map, hs->targ_seq);
    if (is_hap_coded(w, start, end)) {
        const ref_gt_rec *first = restricted_ref(w, start);
        hs->group = first->group;
        int n_seq1 = first->group->n_seq;
        hs->ref_seq = util_malloc((size_t)n_seq1 * sizeof *hs->ref_seq);
        for (int s = 0; s < n_seq1; ++s) hs->ref_seq[s] = 1;
        for (int j = 0; j < n; ++j) {
            const ref_gt_rec *rec = restricted_ref(w, start + j);
            int n_alleles = marker_n_alleles(&rec->marker);
            for (int s = 0; s < n_seq1; ++s) {
                if (hs->ref_seq[s] > 0) hs->ref_seq[s] = seq_map[j][hs->ref_seq[s] * n_alleles + rec->seq_to_allele[s]];
            }
        }
    } else {
        hs->group = NULL;
        hs->ref_seq = util_malloc((size_t)id->n_ref_haps * sizeof *hs->ref_seq);
        for (int h = 0; h < id->n_ref_haps; ++h) hs->ref_seq[h] = 1;
        for (int j = 0; j < n; ++j) {
            const ref_gt_rec *rec = restricted_ref(w, start + j);
            int n_alleles = marker_n_alleles(&rec->marker);
            for (int h = 0; h < id->n_ref_haps; ++h) {
                if (hs->ref_seq[h] > 0) hs->ref_seq[h] = seq_map[j][hs->ref_seq[h] * n_alleles + ref_gt_rec_get(rec, h)];
            }
        }
    }
    for (int j = 0; j < n; ++j) free(seq_map[j]);
    free(seq_map);
}

/* ImpData.wts: between clusters, a reference marker's weight on the
 * preceding cluster falls linearly with cM from 1 to 0. */
static float *weights(const imp_data *id, const genetic_map *gm) {
    const window *w = id->w;
    const marker **markers = util_malloc((size_t)w->n_ref * sizeof *markers);
    for (int m = 0; m < w->n_ref; ++m) markers[m] = &w->ref[m]->marker;
    double *cum = cum_pos(gm, markers, w->n_ref);
    free(markers);
    float *wts = util_malloc((size_t)w->n_ref * sizeof *wts);
    int last = id->n_clusters - 1;
    for (int m = 0; m < id->ref_cluster_start[0]; ++m) wts[m] = NAN;
    for (int j = 0; j < last; ++j) {
        int start = id->ref_cluster_start[j];
        int end = id->ref_cluster_end[j];
        int next_start = id->ref_cluster_start[j + 1];
        double total_length = cum[next_start] - cum[end - 1];
        for (int m = start; m < end; ++m) wts[m] = NAN;
        for (int m = end; m < next_start; ++m) wts[m] = (float)((cum[next_start] - cum[m]) / total_length);
    }
    for (int m = id->ref_cluster_start[last]; m < w->n_ref; ++m) wts[m] = NAN;
    free(cum);
    return wts;
}

static void trace(const imp_data *id, float err_rate) {
    trace_line("T5a", "imp\t%d\t%d\t%d\t%d\t%d\t%" PRIx32, id->w->index, id->n_clusters, id->n_ref_haps, id->n_targ_haps,
            id->n_input_targ_haps, jnum_float_bits(err_rate));
    for (int c = 0; c < id->n_clusters; ++c) {
        uint64_t h = TRACE_FNV_BASIS;
        for (int hap = 0; hap < id->n_haps; ++hap) h = trace_fold(h, (uint32_t)imp_data_seq(id, c, hap));
        uint64_t pos_bits;
        memcpy(&pos_bits, &id->pos[c], sizeof pos_bits);
        trace_line("T5a", "c\t%d\t%d\t%d\t%d\t%d\t%" PRIx64 "\t%" PRIx32 "\t%" PRIx32 "\t%d\t%" PRIx64, c,
                id->targ_clust_start_end[c], id->targ_clust_start_end[c + 1], id->ref_cluster_start[c],
                id->ref_cluster_end[c], pos_bits, jnum_float_bits(id->err_prob[c]), jnum_float_bits(id->p_recomb[c]),
                id->hap_to_seq[c].n_seq, h);
    }
    uint64_t h = TRACE_FNV_BASIS;
    for (int m = 0; m < id->w->n_ref; ++m) h = trace_fold(h, jnum_float_bits(id->weight[m]));
    trace_line("T5a", "weight\t%" PRIx64, h);
}

void imp_data_init(imp_data *id, const par *p, const window *w, const samples *targ_samples,
        phased_allele_fn targ_allele, const void *targ_ctx, const genetic_map *gm) {
    id->w = w;
    id->targ_allele = targ_allele;
    id->targ_ctx = targ_ctx;
    const marker **markers = util_malloc((size_t)w->n_targ * sizeof *markers);
    for (int m = 0; m < w->n_targ; ++m) markers[m] = &w->targ[m]->marker;
    double *targ_pos = cum_pos(gm, markers, w->n_targ);
    free(markers);
    int_list block_end = targ_block_end(w);
    int_list se = targ_clust_start_end(targ_pos, &block_end, p->cluster);
    free(block_end.v);
    id->targ_clust_start_end = se.v;
    id->n_clusters = se.n - 1;
    id->pos = util_malloc((size_t)id->n_clusters * sizeof *id->pos);
    for (int c = 0; c < id->n_clusters; ++c) id->pos[c] = (targ_pos[se.v[c]] + targ_pos[se.v[c + 1] - 1]) / 2;
    free(targ_pos);
    id->n_ref_haps = w->ref[0]->n_haps;
    id->n_targ_haps = w->targ[0]->n_haps;
    id->n_input_targ_haps = 0;
    for (int s = 0; s < targ_samples->n; ++s) id->n_input_targ_haps += targ_samples->is_diploid[s] ? 2 : 1;
    id->n_haps = id->n_ref_haps + id->n_targ_haps;
    id->hap_to_seq = util_malloc((size_t)id->n_clusters * sizeof *id->hap_to_seq);
    for (int c = 0; c < id->n_clusters; ++c) code_cluster(id, se.v[c], se.v[c + 1], &id->hap_to_seq[c]);
    const int *targ_to_ref = w->indices.targ_marker_to_marker;
    id->ref_cluster_start = util_malloc((size_t)id->n_clusters * sizeof *id->ref_cluster_start);
    id->ref_cluster_end = util_malloc((size_t)id->n_clusters * sizeof *id->ref_cluster_end);
    for (int c = 0; c < id->n_clusters; ++c) {
        id->ref_cluster_start[c] = targ_to_ref[se.v[c]];
        id->ref_cluster_end[c] = targ_to_ref[se.v[c + 1] - 1] + 1;
    }
    float err_rate = par_err(p, id->n_haps);
    id->err_prob = util_malloc((size_t)id->n_clusters * sizeof *id->err_prob);
    for (int c = 0; c < id->n_clusters; ++c) {
        id->err_prob[c] = err_rate * (se.v[c + 1] - se.v[c]);
        if (id->err_prob[c] > 0.5f) id->err_prob[c] = 0.5f;
    }
    double c_rec = -(0.04 * p->ne / id->n_ref_haps);
    id->p_recomb = util_malloc((size_t)id->n_clusters * sizeof *id->p_recomb);
    id->p_recomb[0] = 0.0f;
    for (int c = 1; c < id->n_clusters; ++c) id->p_recomb[c] = (float)-jmath_expm1(c_rec * (id->pos[c] - id->pos[c - 1]));
    id->weight = weights(id, gm);
    if (trace_on()) trace(id, err_rate);
}

void imp_data_free(imp_data *id) {
    for (int c = 0; c < id->n_clusters; ++c) {
        free(id->hap_to_seq[c].ref_seq);
        free(id->hap_to_seq[c].targ_seq);
    }
    free(id->hap_to_seq);
    free(id->targ_clust_start_end);
    free(id->ref_cluster_start);
    free(id->ref_cluster_end);
    free(id->pos);
    free(id->err_prob);
    free(id->p_recomb);
    free(id->weight);
}
