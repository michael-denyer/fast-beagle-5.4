/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtPhaser.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PBWT_PHASER_H
#define PHASE_PBWT_PHASER_H

#include "blbutil/int_list.h"
#include "phase/fixed_phase_data.h"
#include "phase/sample_phase.h"

/* PbwtPhaser.hiFreqWindows: overlapping [start, end) ranges of stage-1
 * markers, each advancing max(3 cM, (total cM - 1.5) / nthreads) and
 * overlapping the previous by 1.5 cM, with Beagle 5.4's index rounding.
 * Beagle's output depends on nthreads through this split. Returns start, end
 * pairs. */
int_list pbwt_phaser_windows(const fixed_phase_data *fpd, int nthreads);

/* PbwtPhaser.initPhase(fpd, seed): each target sample's initial phase from
 * forward PBWT phasers over those windows, window j seeded with seed + j and
 * aligned to the previous window at a heterozygote in their overlap. Returns
 * one SamplePhase per target sample. */
sample_phase *pbwt_phaser_init_phase(const fixed_phase_data *fpd, int nthreads, int64_t seed);

#endif
