/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/Ibs2.java,
 * phase/Ibs2Markers.java, phase/Ibs2Sets.java and phase/SampleSeg.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_IBS2_H
#define PHASE_IBS2_H

#include <stdbool.h>

#include "vcf/marker_map.h"

/* Markers start..incl_end over which a target sample and `sample` share both
 * alleles. */
typedef struct {
    int sample;
    int start;
    int incl_end;
} sample_seg;

/* The IBS2 segments of each target sample at the stage-1 markers, at least
 * 2 cM long. Phasing uses them only to veto IBS neighbours. */
typedef struct {
    int n_markers;
    int n_samples;
    int *n_segs;
    sample_seg **segs;   /* per target sample, sorted by sample, start, end */
} ibs2;

struct fixed_phase_data;

/* new Ibs2(stage1TargGT, stage1Map, stage1Maf) from FixedPhaseData's
 * stage-1 genotypes, map and frequencies. Writes trace seam T3a. */
void ibs2_init(ibs2 *ib, const struct fixed_phase_data *fpd);
/* Ibs2.areIbs2(targSample, otherSample, start, inclEnd): whether the samples
 * share an IBS2 segment overlapping markers [start, incl_end]. A sample is
 * IBS2 with itself. */
bool ibs2_are_ibs2(const ibs2 *ib, int targ_sample, int other_sample, int start, int incl_end);
void ibs2_free(ibs2 *ib);

#endif
