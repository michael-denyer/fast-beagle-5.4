/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/RefIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_REF_IT_H
#define VCF_REF_IT_H

#include "blbutil/sample_file_it.h"
#include "blbutil/str_set.h"

/* RefIt: the reference VCF records (ref_gt_rec) that pass the marker filter,
 * and the samples of the header. Records whose major
 * allele is rare enough are grouped by SeqCoder3 into sequence-coded records;
 * a group is released on a chromosome change, when it is full, and at the end
 * of the file, so records are read ahead to the end of their group. Blank
 * lines are not skipped. Lines are parsed on n_threads threads, a batch at a
 * time (vcf/block_reader.h). With tracing on, writes seams T1a-ref, T1b-ref and
 * T1d-ref. */
sample_file_it ref_it_open(const char *path, const str_set *exclude_samples, const str_set *exclude_markers,
        int n_threads);

#endif
