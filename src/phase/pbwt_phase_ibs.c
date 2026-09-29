/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtPhaseIbs.java and
 * phase/PbwtIbsData.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/pbwt_phase_ibs.h"

#include <stdio.h>
#include <stdlib.h>

#include <htslib/kstring.h>

#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jrandom.h"
#include "phase/pbwt_ibs_data.h"

static int min_int(int a, int b) {
    return a < b ? a : b;
}

static int max_int(int a, int b) {
    return a > b ? a : b;
}

/* Scanning a[u..v) cyclically from a random start, the first haplotype
 * outside a[i]'s sample and not IBS2 with it, or -1. */
static int select_hap(const pbwt_phase_ibs *pi, int step, const int *a, int i, int u, int v, jrandom *r) {
    int n = v - u;
    if (n <= 1) return -1;
    const fixed_phase_data *fpd = pi->pd->fpd;
    int m_start = steps_start(&fpd->stage1_steps, step);
    int m_incl_end = fpd->stage1_steps.ends[step] - 1;
    int index = u + jrandom_next_int_bound(r, n);
    for (int j = 0; j < n; ++j, ++index) {
        if (index == v) index = u;
        if (index != i && !ibs2_are_ibs2(&fpd->stage1_ibs2, a[i] >> 1, a[index] >> 1, m_start, m_incl_end)) {
            return a[index];
        }
    }
    return -1;
}

static int *bwd_ibs_haps(const pbwt_phase_ibs *pi, int step, const int *a, int *d) {
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(pi->pd), step));
    int n_haps = pi->pd->fpd->n_haps;
    int *selected = util_malloc((size_t)pi->n_targ_haps * sizeof *selected);
    d[0] = d[n_haps] = step - 2;
    for (int i = 0; i < n_haps; ++i) {
        if (a[i] >= pi->n_targ_haps) continue;
        int u = i;
        int v = i + 1;
        int u_next_match_end = d[u];
        int v_next_match_end = d[v];
        while (v - u < pi->n_candidates && (step <= u_next_match_end || step <= v_next_match_end)) {
            if (u_next_match_end <= v_next_match_end) {
                v_next_match_end = min_int(d[++v], v_next_match_end);
            } else {
                u_next_match_end = min_int(d[--u], u_next_match_end);
            }
        }
        selected[a[i]] = select_hap(pi, step, a, i, u, v, &r);
    }
    return selected;
}

static int *fwd_ibs_haps(const pbwt_phase_ibs *pi, int step, const int *a, int *d) {
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(pi->pd), step));
    int n_haps = pi->pd->fpd->n_haps;
    int *selected = util_malloc((size_t)pi->n_targ_haps * sizeof *selected);
    d[0] = d[n_haps] = step + 2;
    for (int i = 0; i < n_haps; ++i) {
        if (a[i] >= pi->n_targ_haps) continue;
        int u = i;
        int v = i + 1;
        int u_next_match_start = d[u];
        int v_next_match_start = d[v];
        while (v - u < pi->n_candidates && (u_next_match_start <= step || v_next_match_start <= step)) {
            if (v_next_match_start <= u_next_match_start) {
                v_next_match_start = max_int(d[++v], v_next_match_start);
            } else {
                u_next_match_start = max_int(d[--u], u_next_match_start);
            }
        }
        selected[a[i]] = select_hap(pi, step, a, i, u, v, &r);
    }
    return selected;
}

static void trace(const pbwt_phase_ibs *pi, bool use_bwd) {
    const phase_data *pd = pi->pd;
    const coded_steps *cs = pi->cs;
    char seam[32];
    snprintf(seam, sizeof seam, "T3d-%d-%s", pd->it, use_bwd ? "bwd" : "fwd");
    trace_line(seam, "ibs\t%d\t%s\t%d\t%d\t%d", pd->it, use_bwd ? "true" : "false", pi->n_candidates,
            pi->steps_per_batch, pi->n_overlap_steps);
    kstring_t s = {0, 0, NULL};
    for (int j = 0; j < cs->n_steps; ++j) {
        s.l = 0;
        ksprintf(&s, "code\t%d\t%d\t", j, cs->n_seq[j]);
        for (int h = 0; h < cs->n_haps; ++h) {
            if (h > 0) kputc(',', &s);
            kputw(cs->hap_to_seq[j][h], &s);
        }
        trace_line(seam, "%s", s.s);
    }
    for (int j = 0; j < pi->n_steps; ++j) {
        s.l = 0;
        ksprintf(&s, "ibsHaps\t%d\t", j);
        for (int h = 0; h < pi->n_targ_haps; ++h) {
            if (h > 0) kputc(',', &s);
            kputw(pi->ibs_haps[j][h], &s);
        }
        trace_line(seam, "%s", s.s);
    }
    free(s.s);
}

typedef struct {
    const pbwt_ibs_data *data;
    pbwt_phase_ibs *pi;
    bool use_bwd;
} batch_job;

static void batch_task(void *worker, int b) {
    const batch_job *job = worker;
    pbwt_batch pb;
    pbwt_batch_init(&pb, job->data, job->pi->cs, b, job->use_bwd);
    while (pbwt_batch_next(&pb)) {
        job->pi->ibs_haps[pb.step] = job->use_bwd ? bwd_ibs_haps(job->pi, pb.step, pb.a, pb.d)
                                                  : fwd_ibs_haps(job->pi, pb.step, pb.a, pb.d);
    }
    pbwt_batch_free(&pb);
}

void pbwt_phase_ibs_init(pbwt_phase_ibs *pi, phase_data *pd, const coded_steps *cs, bool use_bwd) {
    const fixed_phase_data *fpd = pd->fpd;
    pbwt_ibs_data data;
    pbwt_ibs_data_init(&data, pd, cs->n_steps);
    pi->pd = pd;
    pi->cs = cs;
    pi->n_steps = cs->n_steps;
    pi->n_targ_haps = fpd_n_targ_haps(fpd);
    pi->n_candidates = data.n_candidates;
    pi->n_overlap_steps = data.n_overlap_steps;
    pi->steps_per_batch = data.steps_per_batch;
    int n_batches = data.n_batches;
    pi->ibs_haps = util_malloc((size_t)pi->n_steps * sizeof *pi->ibs_haps);
    /* Each batch restarts the PBWT and writes only its own steps. */
    batch_job job = {&data, pi, use_bwd};
    parallel_for(parallel_threads(pd->par->nthreads, n_batches), n_batches, &job, 0, batch_task);
    if (trace_on()) trace(pi, use_bwd);
}

void pbwt_phase_ibs_free(pbwt_phase_ibs *pi) {
    for (int j = 0; j < pi->n_steps; ++j) free(pi->ibs_haps[j]);
    free(pi->ibs_haps);
}
