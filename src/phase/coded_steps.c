/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/CodedSteps.java and the parts
 * of vcf/XRefGT.java it uses; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/coded_steps.h"

#include <stdlib.h>
#include <string.h>

#include <htslib/khash.h>

#include "blbutil/bit_array.h"
#include "blbutil/parallel.h"
#include "blbutil/utilities.h"

KHASH_MAP_INIT_INT(seq_index, int)

typedef struct {
    coded_steps *cs;
    const fixed_phase_data *fpd;
    khash_t(seq_index) *map;
} step_worker;

/* Numbers a step's distinct haplotype sequences in order of first
 * appearance. Each step is independent, so the threads can take steps in any
 * order. */
static void step_task(void *worker, int j) {
    step_worker *w = worker;
    coded_steps *cs = w->cs;
    const steps *st = &w->fpd->stage1_steps;
    int from = w->fpd->stage1_hap_bits[steps_start(st, j)];
    int to = w->fpd->stage1_hap_bits[st->ends[j]];
    int *seq = util_malloc((size_t)cs->n_haps * sizeof *seq);
    int n_seq = 0;
    kh_clear(seq_index, w->map);
    for (int h = 0; h < cs->n_haps; ++h) {
        int absent;
        khiter_t k = kh_put(seq_index, w->map, (khint32_t)bit_array_hash(cs->haps[h], from, to), &absent);
        if (absent < 0) util_exit("fast-beagle: out of memory");
        if (absent) kh_value(w->map, k) = n_seq++;
        seq[h] = kh_value(w->map, k);
    }
    cs->hap_to_seq[j] = seq;
    cs->n_seq[j] = n_seq;
}

void coded_steps_init(coded_steps *cs, const phase_data *pd) {
    const fixed_phase_data *fpd = pd->fpd;
    const steps *st = &fpd->stage1_steps;
    int n_targ_haps = fpd_n_targ_haps(fpd);
    cs->n_steps = st->n;
    cs->n_haps = n_targ_haps + fpd->n_ref_haps;
    cs->haps = util_malloc((size_t)cs->n_haps * sizeof *cs->haps);
    /* EstPhase.phasedHaps copies the target haplotypes: phasing a sample
     * during the iteration must not change what the others see. */
    size_t n_words = ((size_t)fpd->stage1_hap_bits[fpd->n_stage1] + 63) >> 6;
    size_t n_bytes = n_words * sizeof(uint64_t);
    cs->targ_bits = util_malloc((size_t)n_targ_haps * n_bytes);
    for (int s = 0; s < pd->n_samples; ++s) {
        uint64_t *h1 = cs->targ_bits + (size_t)(s << 1) * n_words;
        uint64_t *h2 = h1 + n_words;
        memcpy(h1, pd->phase[s].hap1, n_bytes);
        memcpy(h2, pd->phase[s].hap2, n_bytes);
        cs->haps[s << 1] = h1;
        cs->haps[(s << 1) | 1] = h2;
    }
    for (int h = 0; h < fpd->n_ref_haps; ++h) cs->haps[n_targ_haps + h] = fpd->stage1_ref_haps[h];

    cs->hap_to_seq = util_malloc((size_t)(cs->n_steps > 0 ? cs->n_steps : 1) * sizeof *cs->hap_to_seq);
    cs->n_seq = util_malloc((size_t)(cs->n_steps > 0 ? cs->n_steps : 1) * sizeof *cs->n_seq);
    int n_threads = parallel_threads(pd->par->nthreads, cs->n_steps);
    step_worker *workers = util_malloc((size_t)n_threads * sizeof *workers);
    for (int t = 0; t < n_threads; ++t) workers[t] = (step_worker){cs, fpd, kh_init(seq_index)};
    parallel_for(n_threads, cs->n_steps, workers, sizeof *workers, step_task);
    for (int t = 0; t < n_threads; ++t) kh_destroy(seq_index, workers[t].map);
    free(workers);
}

void coded_steps_free(coded_steps *cs) {
    for (int j = 0; j < cs->n_steps; ++j) free(cs->hap_to_seq[j]);
    free(cs->hap_to_seq);
    free(cs->n_seq);
    free(cs->haps);
    free(cs->targ_bits);
}
