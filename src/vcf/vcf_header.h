/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/VcfHeader.java and
 * vcf/VcfMetaInfo.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_VCF_HEADER_H
#define VCF_VCF_HEADER_H

#include <stddef.h>

#include "blbutil/line_reader.h"
#include "blbutil/str_set.h"
#include "vcf/samples.h"

typedef struct {
    const char *src;      /* the file name */
    int n_header_fields;
    int *included;        /* unfiltered column index of each kept sample */
    samples samples;
} vcf_header;

/* new VcfHeader(src, lines, VcfHeader.isDiploid(firstDataLine), sampleFilter).
 * lines are the '#' lines; the last must be the #CHROM header line. Samples in
 * exclude (which may be NULL) are dropped. Exits on a format error. */
void vcf_header_init(vcf_header *h, const char *src, char **lines, int n_lines,
        const char *first_data_line, size_t first_data_len, const str_set *exclude);
/* VcfIt.head followed by new VcfHeader: reads the '#' lines and the first data
 * line, which is left in first_line. Exits if the file has no data lines. */
void vcf_header_read(vcf_header *h, line_reader *r, kstring_t *first_line, const str_set *exclude);
/* Trace seam T1b: one line per sample with its ploidy. */
void vcf_header_trace(const vcf_header *h, const char *seam);
void vcf_header_free(vcf_header *h);

#endif
