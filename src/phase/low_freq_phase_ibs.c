/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/LowFreqPhaseIbs.java and
 * phase/LowFreqPbwtPhaseIbs.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/low_freq_phase_ibs.h"

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>

#include <htslib/kstring.h>

#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jrandom.h"
#include "phase/pbwt_ibs_data.h"

typedef struct {
    const phase_data *pd;
    const coded_steps *cs;
    const pbwt_ibs_data *data;
    int n_haps;
    int n_targ_haps;
    int *a_inv;
    int *i_to_prev_i;
    int *i_to_next_i;
    int *scratch;       /* sorted prefix positions of one carrier list */
} search;

static int min_int(int a, int b) {
    return a < b ? a : b;
}

static int max_int(int a, int b) {
    return a > b ? a : b;
}

/* setIToPrevNextI: for each haplotype position in the PBWT, the nearest
 * position above and below whose haplotype carries a rare allele it carries,
 * at a marker from this stage-1 step up to the next. */
static void set_i_to_prev_next_i(search *s, int step, const int *a) {
    const fixed_phase_data *fpd = s->pd->fpd;
    for (int j = 0; j < s->n_haps; ++j) s->a_inv[a[j]] = j;
    for (int j = 0; j < s->n_haps; ++j) {
        s->i_to_prev_i[j] = INT_MIN;
        s->i_to_next_i[j] = INT_MAX;
    }
    int start = step == 0 ? 0 : fpd->stage1_to2[steps_start(&fpd->stage1_steps, step)];
    int end = step + 1 < fpd->stage1_steps.n ? fpd->stage1_to2[steps_start(&fpd->stage1_steps, step + 1)] : fpd->n_markers;
    for (int m = start; m < end; ++m) {
        const marker_carriers *mc = &fpd->carriers[m];
        for (int al = 0; al < mc->n_alleles; ++al) {
            const carriers *c = &mc->allele[al];
            if (c->n <= 1) continue;
            int n = 0;
            for (int j = 0; j < c->n; ++j) {
                int h1 = c->samples[j] << 1;
                s->scratch[n++] = s->a_inv[h1];
                s->scratch[n++] = s->a_inv[h1 | 1];
            }
            qsort(s->scratch, (size_t)n, sizeof *s->scratch, util_compare_ints);
            for (int k = 1; k < n; ++k) {
                int i0 = s->scratch[k - 1];
                int i1 = s->scratch[k];
                if (i0 > s->i_to_prev_i[i1]) s->i_to_prev_i[i1] = i0;
                if (i1 < s->i_to_next_i[i0]) s->i_to_next_i[i0] = i1;
            }
        }
    }
}

static bool are_ibs2(const search *s, int step, int h1, int h2) {
    const fixed_phase_data *fpd = s->pd->fpd;
    return ibs2_are_ibs2(&fpd->stage1_ibs2, h1 >> 1, h2 >> 1, steps_start(&fpd->stage1_steps, step), fpd->stage1_steps.ends[step] - 1);
}

/* The nearest carrier-sharing neighbour above or below position i, when not
 * IBS2 with it; the one whose match starts later loses unless its partner is
 * absent. -1 if neither qualifies. */
static int best_fwd_stage2_index(const search *s, int step, int i, const int *a, const int *d) {
    int best_prev_match = -1;
    int best_next_match = -1;
    int prev_match_start = 0;
    int next_match_start = 0;
    int min_match_start = i + 1 < s->n_haps ? min_int(d[i], d[i + 1]) : d[i];
    int d_max = min_int(min_match_start + s->data->max_backoff_steps, step);
    int prev_i = s->i_to_prev_i[i];
    while (prev_i > INT_MIN && are_ibs2(s, step, a[i], a[prev_i])) prev_i = s->i_to_prev_i[prev_i];
    if (prev_i > INT_MIN) {
        int u = i;
        while (u - 1 != prev_i && d[u] <= d_max) prev_match_start = max_int(prev_match_start, d[u--]);
        if (u - 1 == prev_i && d[u] <= d_max) {
            prev_match_start = max_int(prev_match_start, d[u]);
            best_prev_match = prev_i;
        }
    }
    int next_i = s->i_to_next_i[i];
    while (next_i < INT_MAX && are_ibs2(s, step, a[i], a[next_i])) next_i = s->i_to_next_i[next_i];
    if (next_i < INT_MAX) {
        int v = i;
        while (v + 1 != next_i && d[v + 1] <= d_max) next_match_start = max_int(next_match_start, d[++v]);
        if (v + 1 == next_i && d[v + 1] <= d_max) {
            next_match_start = max_int(next_match_start, d[++v]);
            best_next_match = next_i;
        }
    }
    return prev_match_start < next_match_start && best_prev_match != -1 ? best_prev_match : best_next_match;
}

static int best_bwd_stage2_index(const search *s, int step, int i, const int *a, const int *d) {
    int n_steps_m1 = s->pd->fpd->stage1_steps.n - 1;
    int best_prev_match = -1;
    int best_next_match = -1;
    int prev_match_incl_end = n_steps_m1;
    int next_match_incl_end = n_steps_m1;
    int max_match_start = i + 1 < s->n_haps ? max_int(d[i], d[i + 1]) : d[i];
    int d_min = max_int(max_match_start - s->data->max_backoff_steps, step);
    int prev_i = s->i_to_prev_i[i];
    while (prev_i > INT_MIN && are_ibs2(s, step, a[i], a[prev_i])) prev_i = s->i_to_prev_i[prev_i];
    if (prev_i > INT_MIN) {
        int u = i;
        while (u - 1 != prev_i && d[u] >= d_min) prev_match_incl_end = min_int(prev_match_incl_end, d[u--]);
        if (u - 1 == prev_i && d[u] >= d_min) {
            prev_match_incl_end = min_int(prev_match_incl_end, d[u]);
            best_prev_match = prev_i;
        }
    }
    int next_i = s->i_to_next_i[i];
    while (next_i < INT_MAX && are_ibs2(s, step, a[i], a[next_i])) next_i = s->i_to_next_i[next_i];
    if (next_i < INT_MAX) {
        int v = i;
        while (v + 1 != next_i && d[v + 1] >= d_min) next_match_incl_end = min_int(next_match_incl_end, d[++v]);
        if (v + 1 == next_i && d[v + 1] >= d_min) {
            next_match_incl_end = min_int(next_match_incl_end, d[++v]);
            best_next_match = next_i;
        }
    }
    return prev_match_incl_end > next_match_incl_end && best_prev_match != -1 ? best_prev_match : best_next_match;
}

/* A random haplotype in positions [i_start, i_end) not IBS2 with a[i]. */
static int get_match(const search *s, int step, int i, int i_start, int i_end, const int *a, jrandom *r) {
    int i_length = i_end - i_start;
    if (i_length == 1) return -1;
    int match = -1;
    int index = i_start + jrandom_next_int_bound(r, i_length);
    for (int j = 0; j < i_length && match == -1; ++j) {
        if (!are_ibs2(s, step, a[i], a[index])) match = a[index];
        if (++index == i_end) index = i_start;
    }
    return match;
}

static int *bwd_ibs_haps(search *s, int step, const int *a, int *d) {
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(s->pd), step));
    int *selected = util_malloc((size_t)s->n_targ_haps * sizeof *selected);
    d[s->n_haps] = step - 1;
    for (int i = 0; i < s->n_haps; ++i) {
        if (a[i] >= s->n_targ_haps) continue;
        int best_i = best_bwd_stage2_index(s, step, i, a, d);
        if (best_i >= 0) {
            selected[a[i]] = a[best_i];
            continue;
        }
        int u = i;
        int v = i + 1;
        int u_next_match_end = d[u];
        int v_next_match_end = d[v];
        while (v - u < s->data->n_candidates && (step <= u_next_match_end || step <= v_next_match_end)) {
            if (u_next_match_end <= v_next_match_end) {
                v_next_match_end = min_int(d[++v], v_next_match_end);
            } else {
                u_next_match_end = min_int(d[--u], u_next_match_end);
            }
        }
        selected[a[i]] = get_match(s, step, i, u, v, a, &r);
    }
    return selected;
}

static int *fwd_ibs_haps(search *s, int step, const int *a, int *d) {
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(s->pd), step));
    int *selected = util_malloc((size_t)s->n_targ_haps * sizeof *selected);
    d[s->n_haps] = step + 1;
    for (int i = 0; i < s->n_haps; ++i) {
        if (a[i] >= s->n_targ_haps) continue;
        int best_i = best_fwd_stage2_index(s, step, i, a, d);
        if (best_i >= 0) {
            selected[a[i]] = a[best_i];
            continue;
        }
        int u = i;
        int v = i + 1;
        int u_next_match_start = d[u];
        int v_next_match_start = d[v];
        while (v - u < s->data->n_candidates && (u_next_match_start <= step || v_next_match_start <= step)) {
            if (v_next_match_start <= u_next_match_start) {
                v_next_match_start = max_int(d[++v], v_next_match_start);
            } else {
                u_next_match_start = max_int(d[--u], u_next_match_start);
            }
        }
        selected[a[i]] = get_match(s, step, i, u, v, a, &r);
    }
    return selected;
}

typedef struct {
    search s;
    int **haps;
    bool use_bwd;
} batch_worker;

static void batch_task(void *worker, int b) {
    batch_worker *bw = worker;
    pbwt_batch pb;
    pbwt_batch_init(&pb, bw->s.data, bw->s.cs, b, bw->use_bwd);
    while (pbwt_batch_next(&pb)) {
        set_i_to_prev_next_i(&bw->s, pb.step, pb.a);
        bw->haps[pb.step] = bw->use_bwd ? bwd_ibs_haps(&bw->s, pb.step, pb.a, pb.d)
                                        : fwd_ibs_haps(&bw->s, pb.step, pb.a, pb.d);
    }
    pbwt_batch_free(&pb);
}

static int **ibs_haps(const phase_data *pd, const coded_steps *cs, bool use_bwd) {
    pbwt_ibs_data data;
    pbwt_ibs_data_init(&data, pd, cs->n_steps);
    size_t n = (size_t)pd->fpd->n_haps;
    int **haps = util_malloc((size_t)(cs->n_steps > 0 ? cs->n_steps : 1) * sizeof *haps);
    int n_threads = parallel_threads(pd->par->nthreads, data.n_batches);
    /* Each batch restarts the PBWT and writes only its own steps; each
     * thread has its own scratch arrays. */
    batch_worker *workers = util_malloc((size_t)n_threads * sizeof *workers);
    for (int t = 0; t < n_threads; ++t) {
        workers[t] = (batch_worker){{pd, cs, &data, pd->fpd->n_haps, fpd_n_targ_haps(pd->fpd), util_malloc(n * sizeof(int)),
                util_malloc(n * sizeof(int)), util_malloc(n * sizeof(int)), util_malloc(n * sizeof(int))}, haps, use_bwd};
    }
    parallel_for(n_threads, data.n_batches, workers, sizeof *workers, batch_task);
    for (int t = 0; t < n_threads; ++t) {
        free(workers[t].s.a_inv);
        free(workers[t].s.i_to_prev_i);
        free(workers[t].s.i_to_next_i);
        free(workers[t].s.scratch);
    }
    free(workers);
    return haps;
}

static void trace_haps(const char *label, int **haps, int n_steps, int n_targ_haps) {
    kstring_t s = {0, 0, NULL};
    for (int j = 0; j < n_steps; ++j) {
        s.l = 0;
        ksprintf(&s, "%s\t%d\t", label, j);
        for (int h = 0; h < n_targ_haps; ++h) {
            if (h > 0) kputc(',', &s);
            kputw(haps[j][h], &s);
        }
        trace_line("T4a", "%s", s.s);
    }
    free(s.s);
}

void low_freq_phase_ibs_init(low_freq_phase_ibs *lf, const phase_data *pd) {
    lf->pd = pd;
    coded_steps_init(&lf->cs, pd);
    lf->n_steps = lf->cs.n_steps;
    lf->n_targ_haps = fpd_n_targ_haps(pd->fpd);
    lf->fwd = ibs_haps(pd, &lf->cs, false);
    lf->bwd = ibs_haps(pd, &lf->cs, true);
    if (trace_on()) {
        trace_line("T4a", "stage2\t%d", pd->it);
        trace_haps("fwd", lf->fwd, lf->n_steps, lf->n_targ_haps);
        trace_haps("bwd", lf->bwd, lf->n_steps, lf->n_targ_haps);
    }
}

void low_freq_phase_ibs_free(low_freq_phase_ibs *lf) {
    for (int j = 0; j < lf->n_steps; ++j) {
        free(lf->fwd[j]);
        free(lf->bwd[j]);
    }
    free(lf->fwd);
    free(lf->bwd);
    coded_steps_free(&lf->cs);
}
