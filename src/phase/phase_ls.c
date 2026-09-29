/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PhaseLS.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/phase_ls.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"
#include "phase/coded_steps.h"
#include "phase/hmm_param_data.h"
#include "phase/param_estimates.h"
#include "phase/pbwt_phase_ibs.h"

/* PhaseLS.samplesToAnalyze: all target samples, or a random 500. */
static int *samples_to_analyze(const phase_data *pd, jrandom *r, int *n_out) {
    int max_samples = 500;
    int *ia = util_malloc((size_t)pd->n_samples * sizeof *ia);
    for (int s = 0; s < pd->n_samples; ++s) ia[s] = s;
    *n_out = pd->n_samples;
    if (pd->n_samples > max_samples) {
        util_shuffle(ia, pd->n_samples, max_samples, r);
        *n_out = max_samples;
    }
    return ia;
}

static uint64_t double_bits(double d) {
    uint64_t bits;
    memcpy(&bits, &d, sizeof bits);
    return bits;
}

static void trace_estimates(const phase_data *pd, const param_estimates *pe, float p_mismatch, float recomb_intensity) {
    kstring_t s = {0, 0, NULL};
    kputs("mismatch\t", &s);
    for (int j = 0; j < pe->n_mismatch; ++j) {
        if (j > 0) kputc(',', &s);
        ksprintf(&s, "%d:%" PRIx64, pe->mismatch[j].marker_cnt, double_bits(pe->mismatch[j].p_mismatch_sum));
    }
    trace_line("T3c", "%s", s.s);
    s.l = 0;
    kputs("switch\t", &s);
    for (int j = 0; j < pe->n_recomb; ++j) {
        if (j > 0) kputc(',', &s);
        ksprintf(&s, "%" PRIx64 ":%" PRIx64, double_bits(pe->recomb[j].gen_distance), double_bits(pe->recomb[j].switch_prob));
    }
    trace_line("T3c", "%s", s.s);
    free(s.s);
    trace_line("T3c", "params\t%d\t%" PRIx32 "\t%" PRIx32 "\t%" PRIx32 "\t%" PRIx32, pd->it, jnum_float_bits(p_mismatch),
            jnum_float_bits(recomb_intensity), jnum_float_bits(pd->p_mismatch), jnum_float_bits(pd->recomb_intensity));
}

typedef struct {
    hmm_param_data hpd;
    param_estimates pe;
    const int *samples;
} param_worker;

static void param_task(void *worker, int j) {
    param_worker *pw = worker;
    hmm_param_data_update(&pw->hpd, pw->samples[j]);
    hmm_param_data_add_estimation_data(&pw->hpd, &pw->pe);
}

static void add_all(param_estimates *pe, const param_estimates *from) {
    for (int j = 0; j < from->n_mismatch; ++j) {
        param_estimates_add_mismatch_data(pe, from->mismatch[j].marker_cnt, from->mismatch[j].p_mismatch_sum);
    }
    for (int j = 0; j < from->n_recomb; ++j) {
        param_estimates_add_switch_data(pe, from->recomb[j].gen_distance, from->recomb[j].switch_prob);
    }
}

/* PhaseLS.getParamEst. Each sample's sums are independent of the others and
 * the estimates sort them before adding, so the threads can take samples in
 * any order. Java's stopping rule on the switch probability sum never stops
 * early: the sum is reset after every sample. */
static void get_param_est(param_estimates *pe, const pbwt_phase_ibs *ibs, jrandom *r) {
    const phase_data *pd = ibs->pd;
    int n_samples;
    int *samples = samples_to_analyze(pd, r, &n_samples);
    int n_threads = parallel_threads(pd->par->nthreads, n_samples);
    param_worker *workers = util_malloc((size_t)n_threads * sizeof *workers);
    for (int t = 0; t < n_threads; ++t) {
        hmm_param_data_init(&workers[t].hpd, ibs);
        param_estimates_init(&workers[t].pe);
        workers[t].samples = samples;
    }
    parallel_for(n_threads, n_samples, workers, sizeof *workers, param_task);
    for (int t = 0; t < n_threads; ++t) {
        add_all(pe, &workers[t].pe);
        param_estimates_free(&workers[t].pe);
        hmm_param_data_free(&workers[t].hpd);
    }
    free(workers);
    free(samples);
}

static void update_parameters(const pbwt_phase_ibs *ibs, jrandom *r) {
    phase_data *pd = ibs->pd;
    param_estimates pe;
    param_estimates_init(&pe);
    get_param_est(&pe, ibs, r);
    float prev_p_mismatch = pd->p_mismatch;
    float p_mismatch = param_estimates_p_mismatch(&pe);
    float recomb_intensity = param_estimates_recomb_intensity(&pe);
    if (isfinite(p_mismatch) && p_mismatch > prev_p_mismatch) phase_data_update_p_mismatch(pd, p_mismatch);
    if (isfinite(recomb_intensity) && recomb_intensity > 0.0f) phase_data_update_recomb_intensity(pd, recomb_intensity);
    if (trace_on()) trace_estimates(pd, &pe, p_mismatch, recomb_intensity);
    param_estimates_free(&pe);
}

static void initialize_parameters(const pbwt_phase_ibs *ibs, jrandom *r) {
    phase_data *pd = ibs->pd;
    float prev_rec_int = pd->recomb_intensity;
    int max_initial_its = 15;
    for (int j = 0; j < max_initial_its; ++j) {
        update_parameters(ibs, r);
        float rec_int = pd->recomb_intensity;
        if (fabsf(rec_int - prev_rec_int) <= 0.1 * prev_rec_int) break;
        prev_rec_int = rec_int;
    }
}

typedef struct {
    phase_baum1 baum;
    swap_rate rate;
} baum_worker;

static void baum_task(void *worker, int sample) {
    baum_worker *bw = worker;
    phase_baum1_phase(&bw->baum, sample, &bw->rate);
}

void phase_ls_run_stage1(phase_data *pd, swap_rate *rate) {
    coded_steps cs;
    coded_steps_init(&cs, pd);
    pbwt_phase_ibs ibs;
    pbwt_phase_ibs_init(&ibs, pd, &cs, (pd->it & 1) == 0);
    if (trace_on() && pd->it == 0) {
        /* The forward search, traced before Beagle first runs it at iteration 1. */
        pbwt_phase_ibs fwd;
        pbwt_phase_ibs_init(&fwd, pd, &cs, false);
        pbwt_phase_ibs_free(&fwd);
    }
    if (pd->par->em) {
        jrandom r;
        jrandom_init(&r, phase_data_seed(pd));
        if (pd->it == 0) initialize_parameters(&ibs, &r);
        else if (pd->it < pd->par->burnin) update_parameters(&ibs, &r);
    }
    /* Each sample is phased only against the haplotypes copied at the start
     * of the iteration, so the threads can take samples in any order. */
    int n_threads = parallel_threads(pd->par->nthreads, pd->n_samples);
    char **trace_lines = trace_on() ? util_malloc((size_t)pd->n_samples * sizeof *trace_lines) : NULL;
    baum_worker *workers = util_malloc((size_t)n_threads * sizeof *workers);
    for (int t = 0; t < n_threads; ++t) {
        phase_baum1_init(&workers[t].baum, &ibs, trace_lines);
        workers[t].rate = (swap_rate){0, 0};
    }
    parallel_for(n_threads, pd->n_samples, workers, sizeof *workers, baum_task);
    for (int t = 0; t < n_threads; ++t) {
        rate->n_swaps += workers[t].rate.n_swaps;
        rate->n_unph_hets += workers[t].rate.n_unph_hets;
        phase_baum1_free(&workers[t].baum);
    }
    free(workers);
    if (trace_lines != NULL) {
        trace_line("T3e", "it\t%d", pd->it);
        for (int s = 0; s < pd->n_samples; ++s) {
            trace_line("T3e", "%s", trace_lines[s]);
            free(trace_lines[s]);
        }
        free(trace_lines);
    }
    pbwt_phase_ibs_free(&ibs);
    coded_steps_free(&cs);
}
