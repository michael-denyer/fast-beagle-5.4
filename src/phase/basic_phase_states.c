/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/BasicPhaseStates.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/basic_phase_states.h"

#include <math.h>
#include <stdlib.h>

#include "blbutil/bit_array.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"

void basic_phase_states_init(basic_phase_states *bps, const pbwt_phase_ibs *ibs, int max_states) {
    if (max_states < 1) util_exit("java.lang.IllegalArgumentException: %d", max_states);
    const fixed_phase_data *fpd = ibs->pd->fpd;
    bps->ibs = ibs;
    bps->n_markers = fpd->n_stage1;
    bps->max_states = max_states;
    int ceil_steps = jnum_d2i(ceil((double)(1.0f / fpd->ibs_step)));
    bps->min_steps = ceil_steps > 200 ? ceil_steps : 200;
    comp_hap_tracker_init(&bps->t, max_states, ibs->cs->n_haps);
    bps->comp_haps = util_malloc((size_t)max_states * sizeof *bps->comp_haps);
    for (int j = 0; j < max_states; ++j) bps->comp_haps[j] = bit_array_new((size_t)fpd->stage1_hap_bits[fpd->n_stage1]);
    bps->column = util_malloc((size_t)max_states * sizeof *bps->column);
}

void basic_phase_states_free(basic_phase_states *bps) {
    for (int j = 0; j < bps->max_states; ++j) free(bps->comp_haps[j]);
    free(bps->comp_haps);
    free(bps->column);
    comp_hap_tracker_free(&bps->t);
}

/* XRefGT.copyTo: haplotype hap's bits for markers [start, end) into
 * composite haplotype index. */
static void copy_to(basic_phase_states *bps, int hap, int start, int end, int index) {
    const int *hap_bits = bps->ibs->pd->fpd->stage1_hap_bits;
    bit_array_copy_range(bps->comp_haps[index], bps->ibs->cs->haps[hap], hap_bits[start], hap_bits[end]);
}

static void add_ibs_hap(basic_phase_states *bps, int ibs_hap, int step) {
    comp_hap_change c = comp_hap_tracker_observe(&bps->t, ibs_hap, step, bps->min_steps);
    if (c.old_hap >= 0) {
        const steps *st = &bps->ibs->pd->fpd->stage1_steps;
        copy_to(bps, c.old_hap, steps_start(st, c.start_step), steps_start(st, c.end_step), c.index);
    }
}

/* Used when no IBS neighbour was found: random haplotypes outside the sample,
 * a haplotype drawn twice used twice. Beagle 5.4 queues them at the last step,
 * not step 0; nothing reads the queue order afterwards. */
static void fill_q_with_random_haps(basic_phase_states *bps, int sample) {
    int n_haps = bps->ibs->cs->n_haps;
    int n_states = n_haps - 2 < bps->max_states ? n_haps - 2 : bps->max_states;
    if (n_states <= 0) util_exit("ERROR: there is only one sample");
    const phase_data *pd = bps->ibs->pd;
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(pd), sample));
    for (int j = 0; j < n_states; ++j) {
        int h = jrandom_next_int_bound(&r, n_haps);
        while ((h >> 1) == sample) h = jrandom_next_int_bound(&r, n_haps);
        comp_hap_tracker_seed(&bps->t, h);
    }
}

static int copy_final_ref_segs(basic_phase_states *bps) {
    int n_comp_haps = comp_hap_tracker_size(&bps->t);
    for (int j = 0; j < n_comp_haps; ++j) {
        const comp_hap_segment *seg = comp_hap_tracker_segment(&bps->t, j);
        int start = steps_start(&bps->ibs->pd->fpd->stage1_steps, seg->start_step);
        copy_to(bps, seg->hap, start, bps->n_markers, j);
    }
    return n_comp_haps;
}

static int set_comp_ref_haps(basic_phase_states *bps, int sample) {
    int h1 = sample << 1;
    int h2 = h1 | 1;
    comp_hap_tracker_clear(&bps->t);
    for (int step = 0; step < bps->ibs->n_steps; ++step) {
        int ibs_hap1 = bps->ibs->ibs_haps[step][h1];
        if (ibs_hap1 >= 0) add_ibs_hap(bps, ibs_hap1, step);
        int ibs_hap2 = bps->ibs->ibs_haps[step][h2];
        if (ibs_hap2 >= 0) add_ibs_hap(bps, ibs_hap2, step);
    }
    if (comp_hap_tracker_size(&bps->t) == 0) fill_q_with_random_haps(bps, sample);
    return copy_final_ref_segs(bps);
}

/* One word of every composite haplotype, gathered once for all the bit
 * ranges that lie within it. Start a new cursor whenever the composite
 * haplotypes change. */
typedef struct {
    const basic_phase_states *bps;
    int n_comp_haps;
    int word;   /* the word held in bps->column, or -1 */
} column_cursor;

static inline const uint64_t *column_word(column_cursor *cc, int word) {
    const basic_phase_states *bps = cc->bps;
    if (word != cc->word) {
        for (int j = 0; j < cc->n_comp_haps; ++j) bps->column[j] = bps->comp_haps[j][word];
        cc->word = word;
    }
    return bps->column;
}

/* out1[j] and out2[j]: whether composite haplotype j differs from hap1 and
 * from hap2 over bits [from, to). A range within one word compares that word
 * of every composite haplotype at once. */
static inline __attribute__((always_inline)) void set_mismatches(column_cursor *cc, const uint64_t *hap1, const uint64_t *hap2, int from, int to,
        uint8_t *restrict out1, uint8_t *restrict out2) {
    const basic_phase_states *bps = cc->bps;
    int n = cc->n_comp_haps;
    int word = from >> 6;
    if (from < to && word == (to - 1) >> 6) {
        const uint64_t *column = column_word(cc, word);
        uint64_t mask = (~(uint64_t)0 << (from & 63)) & (~(uint64_t)0 >> ((unsigned)-to & 63));
        uint64_t t1 = hap1[word] & mask;
        uint64_t t2 = hap2[word] & mask;
        for (int j = 0; j < n; ++j) {
            uint64_t v = column[j] & mask;
            out1[j] = v != t1;
            out2[j] = v != t2;
        }
    } else {
        for (int j = 0; j < n; ++j) {
            out1[j] = !bit_array_equal_range(hap1, bps->comp_haps[j], from, to);
            out2[j] = !bit_array_equal_range(hap2, bps->comp_haps[j], from, to);
        }
    }
}

static inline __attribute__((always_inline)) void set_mismatch(column_cursor *cc, const uint64_t *hap, int from, int to, uint8_t *restrict out) {
    const basic_phase_states *bps = cc->bps;
    int n = cc->n_comp_haps;
    int word = from >> 6;
    if (from < to && word == (to - 1) >> 6) {
        const uint64_t *column = column_word(cc, word);
        uint64_t mask = (~(uint64_t)0 << (from & 63)) & (~(uint64_t)0 >> ((unsigned)-to & 63));
        uint64_t t = hap[word] & mask;
        for (int j = 0; j < n; ++j) out[j] = (column[j] & mask) != t;
    } else {
        for (int j = 0; j < n; ++j) out[j] = !bit_array_equal_range(hap, bps->comp_haps[j], from, to);
    }
}

int basic_phase_states_ibs_states(basic_phase_states *bps, int sample, uint8_t ***mismatch) {
    int n_comp_haps = set_comp_ref_haps(bps, sample);
    const int *hap_bits = bps->ibs->pd->fpd->stage1_hap_bits;
    const uint64_t *hap1 = bps->ibs->cs->haps[sample << 1];
    const uint64_t *hap2 = bps->ibs->cs->haps[(sample << 1) | 1];
    column_cursor cc = {bps, n_comp_haps, -1};
    for (int m = 0; m < bps->n_markers; ++m) {
        set_mismatches(&cc, hap1, hap2, hap_bits[m], hap_bits[m + 1], mismatch[0][m], mismatch[1][m]);
    }
    return n_comp_haps;
}

/* Row 0 is the homozygous-cluster HMM: it matches rows 1 and 2 at a
 * homozygous cluster and has no mismatches at a heterozygous one. So a
 * homozygous cluster has one distinct row, a heterozygous cluster two, and a
 * missing-genotype cluster none. */
int basic_phase_states_cluster_states(basic_phase_states *bps, const marker_cluster *mc, int **ref_at_missing,
        const uint8_t **mismatch[3], uint8_t *rows, const uint8_t *zero_row) {
    int n_comp_haps = set_comp_ref_haps(bps, mc->sp->sample);
    size_t n = (size_t)n_comp_haps;
    const int *hap_bits = bps->ibs->pd->fpd->stage1_hap_bits;
    const uint64_t *hap1 = mc->sp->hap1;
    const uint64_t *hap2 = mc->sp->hap2;
    column_cursor cc = {bps, n_comp_haps, -1};
    int miss_index = 0;
    for (int c = 0; c < mc->n_clusters; ++c) {
        int m_start = marker_cluster_start(mc, c);
        int b_start = hap_bits[m_start];
        int b_end = hap_bits[mc->ends[c]];
        if (mc->has_missing[c]) {
            mismatch[0][c] = mismatch[1][c] = mismatch[2][c] = zero_row;
            int *ref_alleles = ref_at_missing[miss_index++];
            for (int j = 0; j < n_comp_haps; ++j) ref_alleles[j] = bit_array_allele(bps->comp_haps[j], hap_bits, m_start);
        } else if (bit_array_equal_range(hap1, hap2, b_start, b_end)) {
            set_mismatch(&cc, hap1, b_start, b_end, rows);
            mismatch[0][c] = mismatch[1][c] = mismatch[2][c] = rows;
            rows += n;
        } else {
            set_mismatches(&cc, hap1, hap2, b_start, b_end, rows, rows + n);
            mismatch[0][c] = zero_row;
            mismatch[1][c] = rows;
            mismatch[2][c] = rows + n;
            rows += 2 * n;
        }
    }
    return n_comp_haps;
}
