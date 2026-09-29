/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/VcfIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_VCF_IT_H
#define VCF_VCF_IT_H

#include "blbutil/sample_file_it.h"
#include "blbutil/str_set.h"
#include "vcf/marker.h"

/* VcfIt: the target VCF records (gt_rec) that pass the marker filter, and the
 * samples of the header. Blank lines between records are skipped. With tracing
 * on, writes seams T1a-target, T1b-target and T1c-target. */
sample_file_it vcf_it_open(const char *path, const str_set *exclude_samples, const str_set *exclude_markers,
        int n_threads);

/* Trace seam T1a: a marker as an iterator returns it. Shared with ref_gt_rec. */
void vcf_it_trace_marker(const char *seam, const marker *m);

#endif
