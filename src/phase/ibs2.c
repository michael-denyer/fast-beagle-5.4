/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/Ibs2.java,
 * phase/Ibs2Markers.java, phase/Ibs2Sets.java and phase/SampleSeg.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/ibs2.h"

#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/int_list.h"
#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "phase/fixed_phase_data.h"

static const float MIN_IBS2_CM = 2.0f;
static const float MAX_IBD_GAP_CM = 4.0f;
static const float MAX_MISS_FREQ = 0.1f;
static const float MIN_MINOR_FREQ = 0.1f;
static const int MIN_MARKER_CNT = 50;
static const float MIN_INTERMARKER_CM = 0.02f;
static const float MAX_MISS_STEP_FREQ = 0.1f;

typedef struct {
    const fixed_phase_data *fpd;
    int n_markers;
    int n_samples;
} stage1_gt;

static int allele(const stage1_gt *gt, int m, int hap) {
    return fpd_targ_allele(gt->fpd, m, hap);
}

static int n_alleles(const stage1_gt *gt, int m) {
    return fpd_n_alleles(gt->fpd, m);
}

/* Ibs2Markers: the markers used for IBS2 discovery and the step starts. */

static bool use_marker(const stage1_gt *gt, int m, const float *maf, int max_miss) {
    if (!(maf[m] >= MIN_MINOR_FREQ)) return false;
    int miss = 0;
    for (int h = 0, n = gt->n_samples << 1; h < n; ++h) miss += allele(gt, m, h) < 0;
    return miss <= max_miss;
}

/* Clears use[] for markers closer than MIN_INTERMARKER_CM to the last used
 * marker, as Java's nextStart does. */
static int next_start(const marker_map *map, int start, bool *use) {
    double cm_pos = map->gen_pos[start];
    double min_cm_pos = cm_pos + (double)MIN_INTERMARKER_CM;
    int next = start + 1;
    int mkr_cnt = 0;
    while (next < map->n && mkr_cnt < MIN_MARKER_CNT) {
        if (use[next]) {
            cm_pos = map->gen_pos[next];
            if (cm_pos < min_cm_pos) {
                use[next] = false;
            } else {
                ++mkr_cnt;
                min_cm_pos = cm_pos + (double)MIN_INTERMARKER_CM;
            }
        }
        ++next;
    }
    return next;
}

/* The last two steps are combined. */
static int_list step_starts(const marker_map *map, bool *use) {
    int_list starts = {0};
    int last_start = 0;
    int next = next_start(map, last_start, use);
    while (next < map->n) {
        int_list_add(&starts, last_start);
        last_start = next;
        next = next_start(map, next, use);
    }
    return starts;
}

/* Ibs2Sets: per step, the partition of target samples by their genotypes at
 * every used marker. */

typedef struct {
    int *samples;   /* increasing */
    int n;
    bool is_hom;    /* every genotype so far is homozygous */
} samp_clust;

typedef struct {
    samp_clust *v;
    int n, cap;
} clust_list;

static void clust_list_add(clust_list *l, samp_clust c) {
    if (l->n == l->cap) {
        l->cap = l->cap == 0 ? 8 : 2 * l->cap;
        l->v = util_realloc(l->v, (size_t)l->cap * sizeof *l->v);
    }
    l->v[l->n++] = c;
}

static int genotype_index(const stage1_gt *gt, int m, int s) {
    int a1 = allele(gt, m, s << 1);
    int a2 = allele(gt, m, (s << 1) | 1);
    if (a1 < 0 || a2 < 0) return -1;
    return a1 <= a2 ? ((a2 * (a2 + 1)) >> 1) + a1 : ((a1 * (a1 + 1)) >> 1) + a2;
}

/* Samples missing a genotype join every list: those that exist and those
 * created later. Children keep genotype-index order and need two samples. */
static void partition(const stage1_gt *gt, const samp_clust *parent, int m, clust_list *out) {
    int n_al = n_alleles(gt, m);
    int n_gt = (n_al * (n_al + 1)) >> 1;
    int_list *lists = util_malloc((size_t)n_gt * sizeof *lists);
    bool *created = util_malloc((size_t)n_gt * sizeof *created);
    for (int k = 0; k < n_gt; ++k) {
        lists[k] = (int_list){0};
        created[k] = false;
    }
    int_list missing = {0};
    for (int j = 0; j < parent->n; ++j) {
        int s = parent->samples[j];
        int g = genotype_index(gt, m, s);
        if (g < 0) {
            int_list_add(&missing, s);
            for (int k = 0; k < n_gt; ++k) {
                if (created[k]) int_list_add(&lists[k], s);
            }
        } else {
            if (!created[g]) {
                created[g] = true;
                for (int i = 0; i < missing.n; ++i) int_list_add(&lists[g], missing.v[i]);
            }
            int_list_add(&lists[g], s);
        }
    }
    for (int k = 0; k < n_gt; ++k) {
        if (created[k] && lists[k].n > 1) {
            bool is_hom_gt = false;
            for (int a = 0; a < n_al; ++a) is_hom_gt |= k == ((a * (a + 1)) >> 1) + a;
            clust_list_add(out, (samp_clust){lists[k].v, lists[k].n, parent->is_hom && is_hom_gt});
        } else {
            free(lists[k].v);
        }
    }
    free(missing.v);
    free(created);
    free(lists);
}

/* Excludes samples missing a genotype at more than MAX_MISS_STEP_FREQ of the
 * step's markers. */
static samp_clust init_cluster(const stage1_gt *gt, const int_list *step_markers) {
    int *miss = util_malloc((size_t)gt->n_samples * sizeof *miss);
    memset(miss, 0, (size_t)gt->n_samples * sizeof *miss);
    for (int j = 0; j < step_markers->n; ++j) {
        int m = step_markers->v[j];
        for (int s = 0; s < gt->n_samples; ++s) {
            miss[s] += allele(gt, m, s << 1) == -1 || allele(gt, m, (s << 1) | 1) == -1;
        }
    }
    int max_miss = jnum_d2i(floor((double)(MAX_MISS_STEP_FREQ * (float)step_markers->n)));
    int_list samples = {0};
    for (int s = 0; s < gt->n_samples; ++s) {
        if (miss[s] <= max_miss) int_list_add(&samples, s);
    }
    free(miss);
    return (samp_clust){samples.v, samples.n, true};
}

/* Ibs2Sets.results: a sample's set is its first heterozygous cluster, or the
 * sorted distinct union when missing genotypes put it in several. */
static void add_to_set(int_list *set, const samp_clust *c) {
    bool was_empty = set->n == 0;
    for (int j = 0; j < c->n; ++j) int_list_add(set, c->samples[j]);
    if (was_empty) return;
    qsort(set->v, (size_t)set->n, sizeof *set->v, util_compare_ints);
    int n = 0;
    for (int j = 0; j < set->n; ++j) {
        if (n == 0 || set->v[n - 1] != set->v[j]) set->v[n++] = set->v[j];
    }
    set->n = n;
}

typedef struct {
    sample_seg *v;
    int n, cap;
} seg_list;

static void seg_list_add(seg_list *l, sample_seg ss) {
    if (l->n == l->cap) {
        l->cap = l->cap == 0 ? 8 : 2 * l->cap;
        l->v = util_realloc(l->v, (size_t)l->cap * sizeof *l->v);
    }
    l->v[l->n++] = ss;
}

/* Ibs2Sets.ibs2Sets for one step: each sample's IBS2 set, other than the
 * sample itself, as (sample, other) pairs in sample order. */
static int_list step_pairs(const stage1_gt *gt, const bool *use, int start, int end) {
    int_list step_markers = {0};
    for (int m = start; m < end; ++m) {
        if (use[m]) int_list_add(&step_markers, m);
    }
    clust_list part = {0};
    clust_list_add(&part, init_cluster(gt, &step_markers));
    for (int j = 0; j < step_markers.n; ++j) {
        clust_list next = {0};
        for (int c = 0; c < part.n; ++c) {
            partition(gt, &part.v[c], step_markers.v[j], &next);
            free(part.v[c].samples);
        }
        free(part.v);
        part = next;
    }
    int_list *sets = util_malloc((size_t)gt->n_samples * sizeof *sets);
    for (int s = 0; s < gt->n_samples; ++s) sets[s] = (int_list){0};
    for (int c = 0; c < part.n; ++c) {
        if (!part.v[c].is_hom) {
            for (int j = 0; j < part.v[c].n; ++j) {
                int s = part.v[c].samples[j];
                add_to_set(&sets[s], &part.v[c]);
            }
        }
        free(part.v[c].samples);
    }
    int_list pairs = {0};
    for (int s = 0; s < gt->n_samples; ++s) {
        for (int j = 0; j < sets[s].n; ++j) {
            if (sets[s].v[j] != s) {
                int_list_add(&pairs, s);
                int_list_add(&pairs, sets[s].v[j]);
            }
        }
        free(sets[s].v);
    }
    free(sets);
    free(part.v);
    free(step_markers.v);
    return pairs;
}

/* Ibs2: sort, merge, extend, merge again, then drop short segments. */

static int compare_segs(const void *a, const void *b) {
    const sample_seg *x = a, *y = b;
    if (x->sample != y->sample) return x->sample < y->sample ? -1 : 1;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    if (x->incl_end != y->incl_end) return x->incl_end < y->incl_end ? -1 : 1;
    return 0;
}

static void merge_segments(seg_list *l, const marker_map *map) {
    if (l->n < 2) return;
    int n = 0;
    sample_seg prev = l->v[0];
    for (int j = 1; j < l->n; ++j) {
        sample_seg next = l->v[j];
        double gap = map->gen_pos[next.start] - map->gen_pos[prev.incl_end];
        if (prev.sample == next.sample && gap <= (double)MAX_IBD_GAP_CM) {
            prev.incl_end = next.incl_end;
        } else {
            l->v[n++] = prev;
            prev = next;
        }
    }
    l->v[n++] = prev;
    l->n = n;
}

static bool phase_consistent(int a1, int a2, int b1, int b2) {
    return (a1 < 0 || b1 < 0 || a1 == b1) && (a2 < 0 || b2 < 0 || a2 == b2);
}

static bool are_ibs2(const stage1_gt *gt, int m, int s1, int s2) {
    int a1 = allele(gt, m, s1 << 1);
    int a2 = allele(gt, m, (s1 << 1) | 1);
    int b1 = allele(gt, m, s2 << 1);
    int b2 = allele(gt, m, (s2 << 1) | 1);
    return phase_consistent(a1, a2, b1, b2) || phase_consistent(a1, a2, b2, b1);
}

static void extend_segments(const stage1_gt *gt, int sample, seg_list *l) {
    for (int j = 0; j < l->n; ++j) {
        sample_seg *ss = &l->v[j];
        int incl_start = ss->start;
        int excl_end = ss->incl_end + 1;
        while (incl_start > 0 && are_ibs2(gt, incl_start - 1, sample, ss->sample)) --incl_start;
        while (excl_end < gt->n_markers && are_ibs2(gt, excl_end, sample, ss->sample)) ++excl_end;
        ss->start = incl_start;
        ss->incl_end = excl_end - 1;
    }
}

static void apply_length_filter(seg_list *l, const marker_map *map) {
    int n = 0;
    for (int j = 0; j < l->n; ++j) {
        if (map->gen_pos[l->v[j].incl_end] - map->gen_pos[l->v[j].start] >= (double)MIN_IBS2_CM) l->v[n++] = l->v[j];
    }
    l->n = n;
}

static void trace(const ibs2 *ib, const float *maf, const bool *use, const int_list *starts) {
    kstring_t s = {0, 0, NULL};
    kputs("maf\t", &s);
    for (int m = 0; m < ib->n_markers; ++m) {
        if (m > 0) kputc(',', &s);
        ksprintf(&s, "%" PRIx32, jnum_float_bits(maf[m]));
    }
    trace_line("T3a", "%s", s.s);
    s.l = 0;
    kputs("markers\t", &s);
    bool first = true;
    for (int m = 0; m < ib->n_markers; ++m) {
        if (!use[m]) continue;
        if (!first) kputc(',', &s);
        kputw(m, &s);
        first = false;
    }
    trace_line("T3a", "%s", s.s);
    s.l = 0;
    kputs("stepStarts\t", &s);
    for (int j = 0; j < starts->n; ++j) {
        if (j > 0) kputc(',', &s);
        kputw(starts->v[j], &s);
    }
    trace_line("T3a", "%s", s.s);
    for (int t = 0; t < ib->n_samples; ++t) {
        s.l = 0;
        ksprintf(&s, "segs\t%d\t", t);
        for (int j = 0; j < ib->n_segs[t]; ++j) {
            const sample_seg *ss = &ib->segs[t][j];
            if (j > 0) kputc(',', &s);
            ksprintf(&s, "%d:%d-%d", ss->sample, ss->start, ss->incl_end);
        }
        trace_line("T3a", "%s", s.s);
    }
    free(s.s);
}

/* The Java streams over markers (Ibs2Markers), steps (Ibs2Sets) and samples
 * (Ibs2) run in parallel; each writes only its own index. */
typedef struct {
    const stage1_gt *gt;
    const marker_map *map;
    const float *maf;
    int max_miss;
    bool *use;
    const int_list *starts;
    int_list *pairs;    /* per step */
    seg_list *segs;     /* per sample */
} ibs2_ctx;

static void use_task(void *worker, int m) {
    const ibs2_ctx *c = worker;
    c->use[m] = use_marker(c->gt, m, c->maf, c->max_miss);
}

static void step_task(void *worker, int j) {
    const ibs2_ctx *c = worker;
    int end = j + 1 < c->starts->n ? c->starts->v[j + 1] : c->gt->n_markers;
    c->pairs[j] = step_pairs(c->gt, c->use, c->starts->v[j], end);
}

static void sample_task(void *worker, int s) {
    const ibs2_ctx *c = worker;
    seg_list *l = &c->segs[s];
    if (l->n > 1) qsort(l->v, (size_t)l->n, sizeof *l->v, compare_segs);
    merge_segments(l, c->map);
    extend_segments(c->gt, s, l);
    merge_segments(l, c->map);
    apply_length_filter(l, c->map);
}

void ibs2_init(ibs2 *ib, const fixed_phase_data *fpd, int nthreads) {
    const marker_map *map = &fpd->stage1_map;
    const float *maf = fpd->stage1_maf;
    stage1_gt gt = {fpd, map->n, fpd_n_targ_haps(fpd) >> 1};
    int max_miss = jnum_d2i(ceil((double)(MAX_MISS_FREQ * (float)(gt.n_samples << 1))));
    bool *use = util_malloc((size_t)gt.n_markers * sizeof *use);
    ibs2_ctx c = {&gt, map, maf, max_miss, use, NULL, NULL, NULL};
    parallel_for(parallel_threads(nthreads, gt.n_markers), gt.n_markers, &c, 0, use_task);
    int_list starts = step_starts(map, use);
    int_list *pairs = util_malloc((size_t)(starts.n > 0 ? starts.n : 1) * sizeof *pairs);
    c.starts = &starts;
    c.pairs = pairs;
    parallel_for(parallel_threads(nthreads, starts.n), starts.n, &c, 0, step_task);

    /* Ibs2Sets.segList: each sample's segments in step order. */
    seg_list *segs = util_malloc((size_t)gt.n_samples * sizeof *segs);
    for (int s = 0; s < gt.n_samples; ++s) segs[s] = (seg_list){0};
    for (int j = 0; j < starts.n; ++j) {
        int end = j + 1 < starts.n ? starts.v[j + 1] : gt.n_markers;
        for (int k = 0; k < pairs[j].n; k += 2) {
            seg_list_add(&segs[pairs[j].v[k]], (sample_seg){pairs[j].v[k + 1], starts.v[j], end - 1});
        }
        free(pairs[j].v);
    }
    free(pairs);
    c.segs = segs;
    parallel_for(parallel_threads(nthreads, gt.n_samples), gt.n_samples, &c, 0, sample_task);

    ib->n_markers = gt.n_markers;
    ib->n_samples = gt.n_samples;
    ib->n_segs = util_malloc((size_t)gt.n_samples * sizeof *ib->n_segs);
    ib->segs = util_malloc((size_t)gt.n_samples * sizeof *ib->segs);
    for (int s = 0; s < gt.n_samples; ++s) {
        seg_list *l = &segs[s];
        ib->n_segs[s] = l->n;
        ib->segs[s] = l->v;
    }
    if (trace_on()) trace(ib, maf, use, &starts);
    free(segs);
    free(starts.v);
    free(use);
}

void ibs2_free(ibs2 *ib) {
    for (int s = 0; s < ib->n_samples; ++s) free(ib->segs[s]);
    free(ib->segs);
    free(ib->n_segs);
}

bool ibs2_are_ibs2(const ibs2 *ib, int targ_sample, int other_sample, int start, int incl_end) {
    if (start > incl_end) util_exit("java.lang.IndexOutOfBoundsException: %d", start);
    if (targ_sample == other_sample) return true;
    for (int j = 0; j < ib->n_segs[targ_sample]; ++j) {
        const sample_seg *ss = &ib->segs[targ_sample][j];
        if (ss->sample == other_sample && start <= ss->incl_end && ss->start <= incl_end) return true;
    }
    return false;
}
