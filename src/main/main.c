/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) main/Main.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "imp/imp_data.h"
#include "imp/imp_ls.h"
#include "imp/imputed_writer.h"
#include "jcompat/jrandom.h"
#include "main/par.h"
#include "main/run_stats.h"
#include "main/window_writer.h"
#include "phase/fixed_phase_data.h"
#include "phase/phase_data.h"
#include "phase/phase_ls.h"
#include "phase/stage2.h"
#include "vcf/sliding_window.h"

/* BasicGT.isPhased: every target genotype in the window is phased, so
 * Beagle writes the input phase without running the phasing iterations. */
static bool targ_is_phased(const window *w) {
    for (int j = 0; j < w->n_targ; ++j) {
        if (!w->targ[j]->is_phased) return false;
    }
    return true;
}

/* Main.phaseStage1Variants: burn-in and phasing iterations, leaving burn-in
 * early once at most 1% of unphased heterozygotes switch phase. */
static void phase_stage1_variants(phase_data *pd, const par *p, run_stats *rs) {
    int n_its = p->burnin + p->iterations;
    double max_burnin_swap_rate = 0.01;
    while (pd->it < n_its) {
        swap_rate rate = {0, 0};
        int64_t t0 = run_stats_nanos();
        phase_ls_run_stage1(pd, &rate);
        run_stats_stage1(rs, pd, run_stats_nanos() - t0);
        phase_data_increment_it(pd);
        double swap = (double)rate.n_swaps / rate.n_unph_hets;
        if (trace_on()) phase_data_trace_iteration(pd, swap);
        if (pd->it < p->burnin && swap <= max_burnin_swap_rate) phase_data_advance_to_first_phasing_it(pd);
    }
}

static int stage1_allele(const void *ctx, int m, int hap) {
    const phase_data *pd = ctx;
    const sample_phase *sp = &pd->phase[hap >> 1];
    return (hap & 1) == 0 ? sample_phase_allele1(sp, m) : sample_phase_allele2(sp, m);
}

/* The next window's phased overlap (Main.phasedOverlap): each target
 * haplotype's alleles at target markers [targOverlapStart, nextTargSplice). */
typedef struct {
    int n_markers;
    int n_haps;
    int **alleles;   /* [hap][marker] */
} phased_overlap;

static void phased_overlap_free(phased_overlap *po) {
    for (int h = 0; h < po->n_haps; ++h) free(po->alleles[h]);
    free(po->alleles);
    *po = (phased_overlap){0, 0, NULL};
}

static void set_phased_overlap(phased_overlap *po, const window *w, phased_allele_fn allele, const void *ctx) {
    int start = w->indices.targ_overlap_start;
    int end = w->indices.targ_next_splice;
    po->n_markers = end - start;
    po->n_haps = w->targ[0]->n_haps;
    po->alleles = util_malloc((size_t)po->n_haps * sizeof *po->alleles);
    for (int h = 0; h < po->n_haps; ++h) {
        po->alleles[h] = util_malloc((size_t)po->n_markers * sizeof **po->alleles);
        for (int m = start; m < end; ++m) po->alleles[h][m - start] = allele(ctx, m, h);
    }
}

/* Main.phaseAndImpute for one window: phases the target genotypes, then
 * writes target markers [prevTargSplice, nextTargSplice), or when the window
 * has markers to impute, the imputed reference markers [prevSplice,
 * nextSplice). Returns the phased overlap in next. */
static void phase_window(const fixed_phase_data *fpd, const par *p, int64_t seed, window_writer *ww, bool impute,
        const samples *targ_samples, const genetic_map *gm, run_stats *rs, phased_overlap *next) {
    const window *w = fpd->win;
    phase_data pd;
    phase_data_init(&pd, fpd, p, seed);
    stage2_haps s2 = {0};
    phased_allele_fn allele;
    const void *ctx;
    if (targ_is_phased(w)) {
        allele = fpd_spliced_allele;
        ctx = fpd;
    } else {
        phase_stage1_variants(&pd, p, rs);
        if (fpd->n_stage1 == fpd->n_markers) {
            allele = stage1_allele;
            ctx = &pd;
        } else {
            int64_t t0 = run_stats_nanos();
            phase_ls_run_stage2(&s2, &pd);
            run_stats_stage2(rs, run_stats_nanos() - t0);
            allele = stage2_haps_allele;
            ctx = &s2;
        }
    }
    if (impute) {
        int64_t t0 = run_stats_nanos();
        imp_data id;
        imp_data_init(&id, p, w, targ_samples, allele, ctx, gm);
        state_probs *sp = imp_ls_state_probs(&id, p);
        imputed_writer_print(ww, &id, sp, p, w->indices.prev_splice, w->indices.next_splice);
        state_probs_free(sp);
        imp_data_free(&id);
        run_stats_imputation(rs, run_stats_nanos() - t0);
    } else {
        window_writer_print_phased(ww, w, w->indices.targ_prev_splice, w->indices.targ_next_splice, allele, ctx);
    }
    set_phased_overlap(next, w, allele, ctx);
    stage2_haps_free(&s2);
    phase_data_free(&pd);
}

int main(int argc, char **argv) {
    par p;
    par_parse(&p, argc - 1, argv + 1);
    if (p.trace != NULL) trace_init(p.trace);
    run_stats rs;
    run_stats_open(&rs, &p, argv[0]);
    sliding_window *sw = sliding_window_open(&p);
    window_writer ww;
    window_writer_open(&ww, &p, sliding_window_targ_samples(sw));
    jrandom rand;
    jrandom_init(&rand, p.seed);
    phased_overlap overlap = {0, 0, NULL};
    int64_t cum_targ_markers = 0, cum_markers = 0;
    window *w;
    while ((w = sliding_window_next(sw)) != NULL) {
        if (w->index == 1) {
            run_stats_sample_summary(&rs, sliding_window_n_ref_samples(sw), sliding_window_targ_samples(sw)->n);
        }
        cum_targ_markers += w->indices.n_targ_markers - w->indices.targ_overlap_end;
        cum_markers += w->indices.n_markers - w->indices.overlap_end;
        window_writer_begin_window(&ww, w);
        run_stats_window_update(&rs, w);
        fixed_phase_data fpd;
        fixed_phase_data_init(&fpd, &p, sliding_window_gen_map(sw), w, overlap.n_markers,
                (const int *const *)overlap.alleles);
        if (trace_on()) fixed_phase_data_trace(&fpd);
        int64_t seed = jrandom_next_long(&rand);
        bool impute = w->indices.n_markers != w->indices.n_targ_markers;
        phased_overlap next;
        phase_window(&fpd, &p, seed, &ww, impute, sliding_window_targ_samples(sw), sliding_window_gen_map(sw), &rs, &next);
        phased_overlap_free(&overlap);
        overlap = next;
        fixed_phase_data_free(&fpd);
        window_free(w);
    }
    phased_overlap_free(&overlap);
    window_writer_close(&ww);   /* before the samples it holds are freed */
    run_stats_close(&rs, cum_targ_markers, cum_markers);
    sliding_window_close(sw);
    trace_close();
    return 0;
}
