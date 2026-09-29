/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtIbsData.java, with the
 * batch PBWT loop that PbwtPhaseIbs and LowFreqPbwtPhaseIbs share; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/pbwt_ibs_data.h"

#include <math.h>
#include <stdlib.h>

#include "blbutil/utilities.h"
#include "jcompat/jnum.h"

static const int BURNIN_CANDIDATES = 100;
static const int MAX_PHASE_CANDIDATES = 90;
static const int MIN_PHASE_CANDIDATES = 5;
static const float MAX_BACKOFF_CM = 0.3f;

static int min_int(int a, int b) {
    return a < b ? a : b;
}

static int max_int(int a, int b) {
    return a > b ? a : b;
}

/* PbwtIbsData.nCandidates1 */
static int n_candidates1(const phase_data *pd) {
    const par *p = pd->par;
    int n = BURNIN_CANDIDATES;
    if (pd->it >= p->burnin) {
        double n_its_remaining = p->burnin + p->iterations - pd->it;
        double prop = n_its_remaining / p->iterations;
        n = (int)jnum_round_d(prop * MAX_PHASE_CANDIDATES);
        if (n < MIN_PHASE_CANDIDATES) n = MIN_PHASE_CANDIDATES;
    }
    return min_int(n, pd->fpd->n_haps);
}

/* PbwtIbsData.nCandidates2: stage 2 scales with rare= and the haplotype
 * count. */
static int n_candidates2(const phase_data *pd) {
    int n_haps = pd->fpd->n_haps;
    float scale_factor = 0.5f;
    int n = jnum_d2i(floor((double)(scale_factor * pd->par->rare * n_haps)));
    return min_int(max_int(n, MIN_PHASE_CANDIDATES), n_haps);
}

void pbwt_ibs_data_init(pbwt_ibs_data *d, const phase_data *pd, int n_steps) {
    const fixed_phase_data *fpd = pd->fpd;
    const par *p = pd->par;
    int n_its = p->burnin + p->iterations;
    d->n_steps = n_steps;
    d->n_candidates = pd->it < n_its ? n_candidates1(pd) : n_candidates2(pd);
    d->n_overlap_steps = jnum_d2i(rint((double)(p->buffer / fpd->ibs_step)));
    d->max_backoff_steps = jnum_d2i(rint((double)(MAX_BACKOFF_CM / fpd->ibs_step)));
    d->steps_per_batch = (n_steps + p->nthreads - 1) / p->nthreads;
    d->n_batches = (n_steps + d->steps_per_batch - 1) / d->steps_per_batch;
}

static void advance(pbwt_batch *b, int j) {
    if (b->bwd) pbwt_div_updater_bwd(&b->u, b->cs->hap_to_seq[j], b->cs->n_seq[j], j, b->a, b->d);
    else pbwt_div_updater_fwd(&b->u, b->cs->hap_to_seq[j], b->cs->n_seq[j], j, b->a, b->d);
}

/* The batch and buffer bounds are PbwtIbsData.startStep, endStep,
 * bufferStartStep and bufferEndStep. */
void pbwt_batch_init(pbwt_batch *b, const pbwt_ibs_data *data, const coded_steps *cs, int batch, bool bwd) {
    int n = cs->n_haps;
    int start = batch * data->steps_per_batch;
    int end = min_int((batch + 1) * data->steps_per_batch, data->n_steps);
    b->cs = cs;
    b->bwd = bwd;
    b->a = util_malloc((size_t)n * sizeof *b->a);
    b->d = util_malloc(((size_t)n + 1) * sizeof *b->d);
    pbwt_div_updater_init(&b->u, n);
    for (int h = 0; h < n; ++h) b->a[h] = h;
    if (bwd) {
        int buffer_end = min_int(end + data->n_overlap_steps, data->n_steps);
        for (int h = 0; h <= n; ++h) b->d[h] = buffer_end - 1;
        for (int j = buffer_end - 1; j >= end; --j) advance(b, j);
        b->step = end;
        b->stop = start - 1;
    } else {
        int buffer_start = max_int(0, start - data->n_overlap_steps);
        for (int h = 0; h <= n; ++h) b->d[h] = buffer_start;
        for (int j = buffer_start; j < start; ++j) advance(b, j);
        b->step = start - 1;
        b->stop = end;
    }
}

bool pbwt_batch_next(pbwt_batch *b) {
    int next = b->bwd ? b->step - 1 : b->step + 1;
    if (next == b->stop) return false;
    b->step = next;
    advance(b, next);
    return true;
}

void pbwt_batch_free(pbwt_batch *b) {
    pbwt_div_updater_free(&b->u);
    free(b->d);
    free(b->a);
}
