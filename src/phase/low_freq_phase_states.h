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
#ifndef PHASE_LOW_FREQ_PHASE_STATES_H
#define PHASE_LOW_FREQ_PHASE_STATES_H

#include <stdint.h>

#include "beagleutil/comp_hap_queue.h"
#include "blbutil/int_list.h"
#include "phase/low_freq_phase_ibs.h"

/* The stage-2 HMM states for one target haplotype: up to max_states composite
 * haplotypes, each a list of reference haplotype segments built from the
 * forward and backward stage-2 IBS neighbours, as in BasicPhaseStates. */
typedef struct {
    const low_freq_phase_ibs *ibs;
    int n_markers;
    int max_states;
    int min_steps;
    comp_hap_tracker t;
    int_list *comp_hap_hap;   /* per state, the haplotype of each segment */
    int_list *comp_hap_end;   /* per state, the end marker of each segment */
    int *segment_index;
    int *comp_hap_to_hap;
    int *comp_hap_to_end;
} low_freq_phase_states;

void low_freq_phase_states_init(low_freq_phase_states *st, const low_freq_phase_ibs *ibs, int max_states);
/* LowFreqPhaseStates.ibsStates(targHap, haps, nMismatches). Sets
 * mismatch[m * max_states + j], whether state j's allele at stage-1 marker m
 * differs from targ_hap's, for every marker m, and copies the states'
 * haplotypes at each marker with slot[m] >= 0 to row slot[m] of haps. Rows
 * are max_states long. Returns the number of states. */
int low_freq_phase_states_ibs_states(low_freq_phase_states *st, int targ_hap, const int *slot, int *haps, uint8_t *mismatch);
void low_freq_phase_states_free(low_freq_phase_states *st);

#endif
