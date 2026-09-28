/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/LowFreqPhaseStates.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/low_freq_phase_states.h"

#include <math.h>
#include <stdlib.h>

#include "blbutil/bit_array.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"

void low_freq_phase_states_init(low_freq_phase_states *st, const low_freq_phase_ibs *ibs, int max_states) {
    const phase_data *pd = ibs->pd;
    if (max_states < 1) util_exit("java.lang.IllegalArgumentException: %d", max_states);
    st->ibs = ibs;
    st->n_markers = pd->fpd->n_stage1;
    st->max_states = max_states;
    int ceil_steps = jnum_d2i(ceil((double)(1.0f / pd->fpd->ibs_step)));
    st->min_steps = ceil_steps > 200 ? ceil_steps : 200;
    comp_hap_tracker_init(&st->t, max_states);
    size_t n = (size_t)max_states;
    st->comp_hap_hap = util_malloc(n * sizeof *st->comp_hap_hap);
    st->comp_hap_end = util_malloc(n * sizeof *st->comp_hap_end);
    for (int j = 0; j < max_states; ++j) st->comp_hap_hap[j] = st->comp_hap_end[j] = (int_list){0};
    st->segment_index = util_malloc(n * sizeof *st->segment_index);
    st->comp_hap_to_hap = util_malloc(n * sizeof *st->comp_hap_to_hap);
    st->comp_hap_to_end = util_malloc(n * sizeof *st->comp_hap_to_end);
}

void low_freq_phase_states_free(low_freq_phase_states *st) {
    for (int j = 0; j < st->max_states; ++j) {
        free(st->comp_hap_hap[j].v);
        free(st->comp_hap_end[j].v);
    }
    free(st->comp_hap_hap);
    free(st->comp_hap_end);
    free(st->segment_index);
    free(st->comp_hap_to_hap);
    free(st->comp_hap_to_end);
    comp_hap_tracker_free(&st->t);
}

static void add_ibs_hap(low_freq_phase_states *st, int ibs_hap, int step) {
    if (ibs_hap < 0) return;
    comp_hap_change c = comp_hap_tracker_observe(&st->t, ibs_hap, step, st->min_steps);
    if (c.index >= 0) {
        int_list_add(&st->comp_hap_hap[c.index], ibs_hap);
        if (c.old_hap >= 0) {
            int next_start = steps_start(&st->ibs->pd->fpd->stage1_steps, c.end_step);
            int_list_add(&st->comp_hap_end[c.index], next_start);
        }
    }
}

/* A haplotype drawn twice is used twice. Beagle 5.4 queues these segments at
 * the last step, not step 0; nothing reads the queue order afterwards. */
static void fill_q_with_random_haps(low_freq_phase_states *st, int hap) {
    int n_haps = st->ibs->cs.n_haps;
    int n_states = n_haps - 2 < st->max_states ? n_haps - 2 : st->max_states;
    if (n_states <= 0) util_exit("ERROR: there is only one sample");
    const phase_data *pd = st->ibs->pd;
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(pd), hap));
    int sample = hap >> 1;
    for (int j = 0; j < n_states; ++j) {
        int h = jrandom_next_int_bound(&r, n_haps);
        while ((h >> 1) == sample) h = jrandom_next_int_bound(&r, n_haps);
        int_list_add(&st->comp_hap_hap[comp_hap_tracker_seed(&st->t, h)], h);
    }
}

static int set_final_ref_segs(low_freq_phase_states *st) {
    int n_comp_haps = comp_hap_tracker_size(&st->t);
    for (int c = 0; c < n_comp_haps; ++c) {
        int_list_add(&st->comp_hap_end[c], st->n_markers);
        st->segment_index[c] = 0;
        st->comp_hap_to_hap[c] = st->comp_hap_hap[c].v[0];
        st->comp_hap_to_end[c] = st->comp_hap_end[c].v[0];
    }
    return n_comp_haps;
}

int low_freq_phase_states_ibs_states(low_freq_phase_states *st, int targ_hap, int **haps, uint8_t **mismatch) {
    comp_hap_tracker_clear(&st->t);
    for (int j = 0; j < st->max_states; ++j) st->comp_hap_hap[j].n = st->comp_hap_end[j].n = 0;
    for (int step = 0; step < st->ibs->n_steps; ++step) {
        add_ibs_hap(st, st->ibs->fwd[step][targ_hap], step);
        add_ibs_hap(st, st->ibs->bwd[step][targ_hap], step);
    }
    if (comp_hap_tracker_size(&st->t) == 0) fill_q_with_random_haps(st, targ_hap);
    int n_comp_haps = set_final_ref_segs(st);

    const int *hap_bits = st->ibs->pd->fpd->stage1_hap_bits;
    const uint64_t *const *all_haps = st->ibs->cs.haps;
    for (int m = 0; m < st->n_markers; ++m) {
        int obs_allele = bit_array_allele(all_haps[targ_hap], hap_bits, m);
        for (int j = 0; j < n_comp_haps; ++j) {
            if (m == st->comp_hap_to_end[j]) {
                ++st->segment_index[j];
                st->comp_hap_to_hap[j] = st->comp_hap_hap[j].v[st->segment_index[j]];
                st->comp_hap_to_end[j] = st->comp_hap_end[j].v[st->segment_index[j]];
            }
            int ref_hap = st->comp_hap_to_hap[j];
            haps[m][j] = ref_hap;
            mismatch[m][j] = bit_array_allele(all_haps[ref_hap], hap_bits, m) == obs_allele ? 0 : 1;
        }
    }
    return n_comp_haps;
}
