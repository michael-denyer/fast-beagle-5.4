/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/FixedPhaseData.java and
 * vcf/Window.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/fixed_phase_data.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/bit_array.h"
#include "blbutil/int_list.h"
#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"

static const float MAX_HIFREQ_PROP = 0.75f;

static void add_carrier(int_list *list, int sample, int max_carriers) {
    if (list->n <= max_carriers) int_list_add(list, sample);
}

/* Window.carriers(m, maxCarriers): lists stop growing one past max_carriers. */
static marker_carriers carriers_at(const window *w, int m, int max_carriers) {
    const gt_rec *targ = w->targ[m];
    int n_alleles = marker_n_alleles(&targ->marker);
    int_list *lists = util_malloc((size_t)n_alleles * sizeof *lists);
    for (int a = 0; a < n_alleles; ++a) lists[a] = (int_list){0};
    int n_targ_samples = targ->n_haps >> 1;
    for (int s = 0; s < n_targ_samples; ++s) {
        int a1 = gt_rec_get(targ, s << 1);
        int a2 = gt_rec_get(targ, (s << 1) | 1);
        if (a1 >= 0) add_carrier(&lists[a1], s, max_carriers);
        if (a2 >= 0 && a2 != a1) add_carrier(&lists[a2], s, max_carriers);
    }
    if (w->n_ref > 0) {
        const ref_gt_rec *ref = w->ref[w->indices.targ_marker_to_marker[m]];
        for (int s = 0, n = ref->n_haps >> 1; s < n; ++s) {
            int a1 = ref_gt_rec_get(ref, s << 1);
            int a2 = ref_gt_rec_get(ref, (s << 1) | 1);
            add_carrier(&lists[a1], n_targ_samples + s, max_carriers);
            if (a2 != a1) add_carrier(&lists[a2], n_targ_samples + s, max_carriers);
        }
    }
    carriers *c = util_malloc((size_t)n_alleles * sizeof *c);
    for (int a = 0; a < n_alleles; ++a) {
        if (lists[a].n == 0) {
            c[a] = (carriers){CARRIERS_ZERO_FREQ, 0, NULL};
        } else if (lists[a].n <= max_carriers) {
            c[a] = (carriers){CARRIERS_LOW_FREQ, lists[a].n, lists[a].v};
        } else {
            c[a] = (carriers){CARRIERS_HIGH_FREQ, 0, NULL};
            free(lists[a].v);
        }
    }
    free(lists);
    return (marker_carriers){n_alleles, c};
}

typedef struct {
    fixed_phase_data *fpd;
    int max_carriers;
} carriers_ctx;

static void carriers_task(void *worker, int m) {
    const carriers_ctx *c = worker;
    c->fpd->carriers[m] = carriers_at(c->fpd->win, m, c->max_carriers);
}

static int max_carriers(const par *p, const window *w) {
    int n_samples = (w->targ[0]->n_haps >> 1) + (w->n_ref > 0 ? w->ref[0]->n_haps >> 1 : 0);
    int max = jnum_d2i(floor((double)((float)n_samples * p->rare)));
    return max > 3 ? max : 3;
}

static bool is_hi_freq(marker_carriers mc) {
    int n_high = 0;
    for (int a = 0; a < mc.n_alleles; ++a) n_high += mc.allele[a].kind == CARRIERS_HIGH_FREQ;
    return n_high > 1;
}

static void ignore_low_freq_carriers(fixed_phase_data *fpd) {
    for (int m = 0; m < fpd->n_markers; ++m) {
        for (int a = 0; a < fpd->carriers[m].n_alleles; ++a) {
            free(fpd->carriers[m].allele[a].samples);
            fpd->carriers[m].allele[a] = (carriers){CARRIERS_HIGH_FREQ, 0, NULL};
        }
    }
}

static int compare_doubles(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* The distances are positive and finite, so Java's sort order is numeric order. */
static float median_diff(const marker_map *map) {
    int n = map->n - 1;
    double *diffs = util_malloc((size_t)n * sizeof *diffs);
    for (int j = 1; j < map->n; ++j) diffs[j - 1] = map->gen_pos[j] - map->gen_pos[j - 1];
    qsort(diffs, (size_t)n, sizeof *diffs, compare_doubles);
    float median = 0.5f * (float)(diffs[(n - 1) >> 1] + diffs[n >> 1]);
    free(diffs);
    return median;
}

static int *prev_stage1_marker(int n_markers, const int *stage1, int n_stage1) {
    int *mkr = util_malloc((size_t)n_markers * sizeof *mkr);
    for (int m = 0; m < stage1[1]; ++m) mkr[m] = 0;
    int start = stage1[1];
    for (int j = 2; j < n_stage1; ++j) {
        for (int m = start; m < stage1[j]; ++m) mkr[m] = j - 1;
        start = stage1[j];
    }
    for (int m = start; m < n_markers; ++m) mkr[m] = n_stage1 - 1;
    return mkr;
}

static float *prev_wt(const marker_map *map, const int *stage1, int n_stage1) {
    float *wt = util_malloc((size_t)map->n * sizeof *wt);
    for (int m = 0; m < stage1[0]; ++m) wt[m] = 1.0f;
    int start = stage1[0];
    for (int j = 1; j < n_stage1; ++j) {
        int end = stage1[j];
        double pos_b = map->gen_pos[end];
        double d = pos_b - map->gen_pos[start];
        wt[start] = 1.0f;
        for (int m = start + 1; m < end; ++m) wt[m] = (float)((pos_b - map->gen_pos[m]) / d);
        start = end;
    }
    for (int m = start; m < map->n; ++m) wt[m] = 1.0f;
    return wt;
}

/* FixedPhaseData.randHaps: all n_haps haplotypes, or a sorted random
 * max_haps of them. */
static int *rand_haps(int n_haps, int max_haps, jrandom *r, int *n_out) {
    int *haps = util_malloc((size_t)n_haps * sizeof *haps);
    for (int h = 0; h < n_haps; ++h) haps[h] = h;
    *n_out = n_haps;
    if (n_haps > max_haps) {
        util_shuffle(haps, n_haps, max_haps, r);
        *n_out = max_haps;
        qsort(haps, (size_t)max_haps, sizeof *haps, util_compare_ints);
    }
    return haps;
}

typedef struct {
    const fixed_phase_data *fpd;
    const int *targ_haps, *ref_haps;
    int n_targ, n_ref;
    float *maf;
} maf_ctx;

/* The second most frequent allele's share of the non-missing sampled alleles. */
static void maf_task(void *worker, int j) {
    const maf_ctx *c = worker;
    const window *w = c->fpd->win;
    int m = c->fpd->stage1_to2[j];
    int n_cnts = marker_n_alleles(&w->targ[m]->marker) + 1;
    int *cnts = util_malloc((size_t)n_cnts * sizeof *cnts);
    memset(cnts, 0, (size_t)n_cnts * sizeof *cnts);
    for (int k = 0; k < c->n_targ; ++k) ++cnts[fpd_spliced_allele(c->fpd, m, c->targ_haps[k]) + 1];
    if (c->n_ref > 0) {
        const ref_gt_rec *ref = w->ref[w->indices.targ_marker_to_marker[m]];
        for (int k = 0; k < c->n_ref; ++k) ++cnts[ref_gt_rec_get(ref, c->ref_haps[k]) + 1];
    }
    cnts[0] = 0;
    qsort(cnts, (size_t)n_cnts, sizeof *cnts, util_compare_ints);
    int den = 0;
    for (int a = 1; a < n_cnts; ++a) den += cnts[a];
    c->maf[j] = (float)(den == 0 ? 0.0 : (double)cnts[n_cnts - 2] / den);
    free(cnts);
}

/* The haplotypes are drawn serially, as in Java; only the per-marker counts
 * run in parallel. */
static float *stage1_maf(const fixed_phase_data *fpd, int max_haps, int64_t seed, int nthreads) {
    const window *w = fpd->win;
    jrandom r;
    jrandom_init(&r, seed);
    int n_targ, n_ref = 0;
    int *targ_haps = rand_haps(w->targ[0]->n_haps, max_haps, &r, &n_targ);
    int *ref_haps = NULL;
    if (n_targ < max_haps && w->n_ref > 0) ref_haps = rand_haps(w->ref[0]->n_haps, max_haps - n_targ, &r, &n_ref);
    float *maf = util_malloc((size_t)fpd->n_stage1 * sizeof *maf);
    maf_ctx c = {fpd, targ_haps, ref_haps, n_targ, n_ref, maf};
    parallel_for(parallel_threads(nthreads, fpd->n_stage1), fpd->n_stage1, &c, 0, maf_task);
    free(targ_haps);
    free(ref_haps);
    return maf;
}

/* XRefGT.fromPhasedGT: reference haplotype h's stage-1 alleles. */
static void ref_hap_task(void *worker, int h) {
    fixed_phase_data *fpd = worker;
    fpd->stage1_ref_haps[h] = bit_array_new((size_t)fpd->stage1_hap_bits[fpd->n_stage1]);
    for (int j = 0; j < fpd->n_stage1; ++j) {
        bit_array_set_allele(fpd->stage1_ref_haps[h], fpd->stage1_hap_bits, j, fpd_ref_allele(fpd, j, h));
    }
}

/* The number of stage-1 markers before marker index `overlap`. */
static int stage1_overlap(const int *stage1, int n_stage1, int overlap) {
    int low = 0;
    while (low < n_stage1 && stage1[low] < overlap) ++low;
    return low;
}

void fixed_phase_data_init(fixed_phase_data *fpd, const par *p, const genetic_map *gm, const window *w, int overlap,
        const int *const *overlap_alleles) {
    int n = w->n_targ;
    fpd->win = w;
    fpd->window = w->index;
    fpd->n_markers = n;
    fpd->n_haps = w->targ[0]->n_haps + (w->n_ref > 0 ? w->ref[0]->n_haps : 0);
    fpd->overlap = overlap;
    fpd->overlap_alleles = overlap_alleles;

    const marker **markers = util_malloc((size_t)n * sizeof *markers);
    for (int m = 0; m < n; ++m) markers[m] = &w->targ[m]->marker;
    marker_map_init(&fpd->map, gm, markers, n);
    free(markers);

    int max = max_carriers(p, w);
    fpd->carriers = util_malloc((size_t)n * sizeof *fpd->carriers);
    carriers_ctx cc = {fpd, max};
    parallel_for(parallel_threads(p->nthreads, n), n, &cc, 0, carriers_task);
    int_list hi_freq = {0};
    for (int m = 0; m < n; ++m) {
        if (is_hi_freq(fpd->carriers[m])) int_list_add(&hi_freq, m);
    }

    if (hi_freq.n < 2 || (float)hi_freq.n > MAX_HIFREQ_PROP * (float)n) {
        ignore_low_freq_carriers(fpd);
        fpd->n_stage1 = n;
        fpd->stage1_to2 = util_malloc((size_t)n * sizeof *fpd->stage1_to2);
        for (int m = 0; m < n; ++m) fpd->stage1_to2[m] = m;
        marker_map_restrict(&fpd->stage1_map, &fpd->map, fpd->stage1_to2, n);
        fpd->stage1_overlap = overlap;
        fpd->prev_stage1_marker = util_malloc((size_t)n * sizeof *fpd->prev_stage1_marker);
        fpd->prev_stage1_wt = util_malloc((size_t)n * sizeof *fpd->prev_stage1_wt);
        for (int m = 0; m < n; ++m) {
            fpd->prev_stage1_marker[m] = m;
            fpd->prev_stage1_wt[m] = 1.0f;
        }
        free(hi_freq.v);
    } else {
        fpd->n_stage1 = hi_freq.n;
        fpd->stage1_to2 = hi_freq.v;
        marker_map_restrict(&fpd->stage1_map, &fpd->map, fpd->stage1_to2, fpd->n_stage1);
        fpd->stage1_overlap = stage1_overlap(fpd->stage1_to2, fpd->n_stage1, overlap);
        fpd->prev_stage1_marker = prev_stage1_marker(n, fpd->stage1_to2, fpd->n_stage1);
        fpd->prev_stage1_wt = prev_wt(&fpd->map, fpd->stage1_to2, fpd->n_stage1);
    }
    fpd->stage1_hap_bits = util_malloc((size_t)(fpd->n_stage1 + 1) * sizeof *fpd->stage1_hap_bits);
    fpd->stage1_hap_bits[0] = 0;
    for (int j = 0; j < fpd->n_stage1; ++j) {
        fpd->stage1_hap_bits[j + 1] = fpd->stage1_hap_bits[j] + marker_bits_per_allele(fpd_marker(fpd, j));
    }
    fpd->n_ref_haps = w->n_ref > 0 ? w->ref[0]->n_haps : 0;
    fpd->stage1_ref_haps = util_malloc((size_t)(fpd->n_ref_haps > 0 ? fpd->n_ref_haps : 1) * sizeof *fpd->stage1_ref_haps);
    parallel_for(parallel_threads(p->nthreads, fpd->n_ref_haps), fpd->n_ref_haps, fpd, 0, ref_hap_task);
    fpd->ibs_step = p->step_scale * median_diff(&fpd->stage1_map);
    steps_init(&fpd->stage1_steps, &fpd->stage1_map, fpd->ibs_step);
    int max_maf_haps = 10000;
    fpd->stage1_maf = stage1_maf(fpd, max_maf_haps, p->seed, p->nthreads);
    ibs2_init(&fpd->stage1_ibs2, fpd, p->nthreads);
}

void fixed_phase_data_free(fixed_phase_data *fpd) {
    for (int m = 0; m < fpd->n_markers; ++m) {
        for (int a = 0; a < fpd->carriers[m].n_alleles; ++a) free(fpd->carriers[m].allele[a].samples);
        free(fpd->carriers[m].allele);
    }
    free(fpd->carriers);
    marker_map_free(&fpd->map);
    marker_map_free(&fpd->stage1_map);
    steps_free(&fpd->stage1_steps);
    free(fpd->stage1_to2);
    free(fpd->stage1_hap_bits);
    for (int h = 0; h < fpd->n_ref_haps; ++h) free(fpd->stage1_ref_haps[h]);
    free(fpd->stage1_ref_haps);
    free(fpd->prev_stage1_marker);
    free(fpd->prev_stage1_wt);
    free(fpd->stage1_maf);
    ibs2_free(&fpd->stage1_ibs2);
}

void fixed_phase_data_trace(const fixed_phase_data *fpd) {
    trace_line("T2b", "fpd\t%d\t%d\t%d\t%d\t%d\t%" PRIx32, fpd->window, fpd->n_markers, fpd->n_haps,
            fpd->overlap, fpd->stage1_overlap, jnum_float_bits(fpd->ibs_step));
    kstring_t s = {0, 0, NULL};
    for (int m = 0; m < fpd->n_markers; ++m) {
        uint64_t pos_bits;
        memcpy(&pos_bits, &fpd->map.gen_pos[m], sizeof pos_bits);
        s.l = 0;
        ksprintf(&s, "marker\t%d\t%" PRIx64 "\t%" PRIx32 "\t%d\t%" PRIx32 "\t", m, pos_bits,
                jnum_float_bits(fpd->map.gen_dist[m]), fpd->prev_stage1_marker[m], jnum_float_bits(fpd->prev_stage1_wt[m]));
        for (int a = 0; a < fpd->carriers[m].n_alleles; ++a) {
            if (a > 0) kputc(';', &s);
            const carriers *c = &fpd->carriers[m].allele[a];
            if (c->kind == CARRIERS_HIGH_FREQ) {
                kputc('H', &s);
            } else if (c->kind == CARRIERS_ZERO_FREQ) {
                kputc('Z', &s);
            } else {
                for (int j = 0; j < c->n; ++j) {
                    if (j > 0) kputc(',', &s);
                    kputw(c->samples[j], &s);
                }
            }
        }
        trace_line("T2b", "%s", s.s);
    }
    s.l = 0;
    kputs("stage1\t", &s);
    for (int j = 0; j < fpd->n_stage1; ++j) {
        if (j > 0) kputc(',', &s);
        kputw(fpd->stage1_to2[j], &s);
    }
    trace_line("T2b", "%s", s.s);
    s.l = 0;
    kputs("stage1Dist\t", &s);
    for (int j = 0; j < fpd->n_stage1; ++j) {
        if (j > 0) kputc(',', &s);
        ksprintf(&s, "%" PRIx32, jnum_float_bits(fpd->stage1_map.gen_dist[j]));
    }
    trace_line("T2b", "%s", s.s);
    s.l = 0;
    kputs("steps\t", &s);
    for (int j = 0; j < fpd->stage1_steps.n; ++j) {
        if (j > 0) kputc(',', &s);
        kputw(fpd->stage1_steps.ends[j], &s);
    }
    trace_line("T2b", "%s", s.s);
    free(s.s);
}
