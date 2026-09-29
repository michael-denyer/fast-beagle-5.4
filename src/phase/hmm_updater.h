/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/HmmUpdater.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_HMM_UPDATER_H
#define PHASE_HMM_UPDATER_H

#include <stdint.h>

/* One step of the scaled forward and backward Li and Stephens HMM
 * recursions, in float. p_mismatch[mismatch[k]] is state k's emission
 * probability. hmm_fwd_update writes the step from prev into fwd. */
float hmm_fwd_update(const float *prev, float *fwd, float fwd_sum, float p_switch, const float p_mismatch[2], const uint8_t *mismatch, int n_states);
void hmm_bwd_update(float *bwd, float p_switch, const float p_mismatch[2], const uint8_t *mismatch, int n_states);

/* The same steps for three HMMs sharing p_switch and n_states, in one pass
 * over the states. Each HMM's sum accumulates in state order, so results
 * equal three separate calls. */
void hmm_fwd_update3(float *fwd[3], float fwd_sums[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states);
void hmm_bwd_update3(float *bwd[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states);

#endif
