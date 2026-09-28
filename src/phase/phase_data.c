/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PhaseData.java,
 * phase/EstPhase.java and vcf/MarkerMap.java pRecomb; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/phase_data.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jmath.h"
#include "jcompat/jnum.h"
#include "phase/pbwt_phaser.h"

/* PhaseData.leaveUnphasedProp: per sample, the proportion of unphased
 * heterozygotes to leave unphased at each phasing iteration, n^(-1/iterations)
 * for its initial count n. */
static float *leave_unphased_prop(const phase_data *pd) {
    float *prop = util_malloc((size_t)pd->n_samples * sizeof *prop);
    for (int s = 0; s < pd->n_samples; ++s) {
        prop[s] = (float)jmath_pow((double)pd->phase[s].n_unphased, -1.0 / pd->par->iterations);
    }
    return prop;
}

/* MarkerMap.pRecomb(recombIntensity) */
static float *p_recomb(const marker_map *map, float recomb_intensity) {
    double c = -recomb_intensity;
    float *p = util_malloc((size_t)map->n * sizeof *p);
    for (int m = 0; m < map->n; ++m) p[m] = (float)-jmath_expm1(c * map->gen_dist[m]);
    return p;
}

static void put_ints(kstring_t *s, const int *a, int n) {
    if (n == 0) kputc('-', s);
    for (int j = 0; j < n; ++j) {
        if (j > 0) kputc(',', s);
        kputw(a[j], s);
    }
}

/* Each target sample's haplotypes, cluster sizes, and unphased and missing
 * markers. */
static void trace_samples(const phase_data *pd, const char *seam) {
    const fixed_phase_data *fpd = pd->fpd;
    kstring_t s = {0, 0, NULL};
    for (int j = 0; j < pd->n_samples; ++j) {
        const sample_phase *sp = &pd->phase[j];
        s.l = 0;
        ksprintf(&s, "sample\t%d\t", sp->sample);
        for (int m = 0; m < fpd->n_stage1; ++m) {
            kputw(sample_phase_allele1(sp, m), &s);
            kputc(m + 1 < fpd->n_stage1 ? ',' : '\t', &s);
        }
        for (int m = 0; m < fpd->n_stage1; ++m) {
            kputw(sample_phase_allele2(sp, m), &s);
            kputc(m + 1 < fpd->n_stage1 ? ',' : '\t', &s);
        }
        for (int c = 0; c < sp->n_clusters; ++c) {
            if (c > 0) kputc(',', &s);
            kputw(sp->clust_size[c], &s);
        }
        kputc('\t', &s);
        put_ints(&s, sp->unphased, sp->n_unphased);
        kputc('\t', &s);
        put_ints(&s, sp->missing, sp->n_missing);
        trace_line(seam, "%s", s.s);
    }
    free(s.s);
}

static void trace(const phase_data *pd) {
    const fixed_phase_data *fpd = pd->fpd;
    trace_line("T3b1", "pd\t%" PRId64 "\t%" PRIx32 "\t%" PRIx32, pd->seed, jnum_float_bits(pd->recomb_intensity),
            jnum_float_bits(pd->p_mismatch));
    kstring_t s = {0, 0, NULL};
    kputs("pRecomb\t", &s);
    for (int m = 0; m < fpd->n_stage1; ++m) {
        if (m > 0) kputc(',', &s);
        ksprintf(&s, "%" PRIx32, jnum_float_bits(pd->p_recomb[m]));
    }
    trace_line("T3b1", "%s", s.s);
    s.l = 0;
    kputs("leaveUnph\t", &s);
    for (int j = 0; j < pd->n_samples; ++j) {
        if (j > 0) kputc(',', &s);
        ksprintf(&s, "%" PRIx32, jnum_float_bits(pd->leave_unph_prop[j]));
    }
    trace_line("T3b1", "%s", s.s);
    free(s.s);
    trace_samples(pd, "T3b1");
}

void phase_data_trace_iteration(const phase_data *pd, double swap_rate) {
    uint64_t bits;
    memcpy(&bits, &swap_rate, sizeof bits);
    trace_line("T3b", "it\t%d\t%" PRIx64, pd->it, bits);
    trace_samples(pd, "T3b");
}

void phase_data_init(phase_data *pd, const fixed_phase_data *fpd, const par *p, int64_t seed) {
    pd->fpd = fpd;
    pd->par = p;
    pd->n_samples = fpd_n_targ_haps(fpd) >> 1;
    pd->phase = pbwt_phaser_init_phase(fpd, p->nthreads, seed);
    pd->seed = seed;
    pd->it = 0;
    pd->leave_unph_prop = leave_unphased_prop(pd);
    pd->recomb_intensity = 0.04f * p->ne / (float)fpd->n_haps;
    pd->p_recomb = p_recomb(&fpd->stage1_map, pd->recomb_intensity);
    pd->p_mismatch = par_li_stephens_p_mismatch(fpd->n_haps);
    if (trace_on()) trace(pd);
}

void phase_data_free(phase_data *pd) {
    for (int s = 0; s < pd->n_samples; ++s) sample_phase_free(&pd->phase[s]);
    free(pd->phase);
    free(pd->leave_unph_prop);
    free(pd->p_recomb);
}

void phase_data_update_p_mismatch(phase_data *pd, float p_mismatch) {
    if (p_mismatch < 0.0 || p_mismatch > 1.0 || !isfinite(p_mismatch)) {
        util_exit("fast-beagle: mismatch probability %g outside [0, 1]", (double)p_mismatch);
    }
    pd->p_mismatch = p_mismatch;
}

void phase_data_update_recomb_intensity(phase_data *pd, float recomb_intensity) {
    free(pd->p_recomb);
    pd->recomb_intensity = recomb_intensity;
    pd->p_recomb = p_recomb(&pd->fpd->stage1_map, recomb_intensity);
}

void phase_data_increment_it(phase_data *pd) {
    ++pd->it;
}

void phase_data_advance_to_first_phasing_it(phase_data *pd) {
    if (pd->it < pd->par->burnin) pd->it = pd->par->burnin;
}
