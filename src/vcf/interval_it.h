/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/IntervalVcfIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_INTERVAL_IT_H
#define VCF_INTERVAL_IT_H

#include "beagleutil/chrom_interval.h"
#include "blbutil/sample_file_it.h"

/* The records of src from its first record in interval up to the first record
 * after that outside interval, which ends the stream. One record is read ahead:
 * the first in-interval record on open, then each next record as the previous
 * is returned. Exits with Java's message when src has no record in interval.
 * Takes ownership of src and closes it on close. With interval NULL (no chrom=
 * parameter), returns src itself. */
sample_file_it interval_it_open(sample_file_it src, const chrom_interval *interval);

#endif
