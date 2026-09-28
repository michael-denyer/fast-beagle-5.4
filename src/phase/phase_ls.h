/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PhaseLS.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_PHASE_LS_H
#define PHASE_PHASE_LS_H

#include "phase/phase_baum1.h"
#include "phase/phase_data.h"

/* PhaseLS.runStage1: the IBS neighbour search, then with em=true the
 * parameter estimates (repeated up to 15 times at iteration 0 until the
 * recombination intensity settles, once per later burn-in iteration), then
 * each sample's phasing. Adds the iteration's swaps to rate. Writes trace
 * seam T3e. */
void phase_ls_run_stage1(phase_data *pd, swap_rate *rate);

#endif
