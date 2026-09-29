/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/VcfRecGTParser.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_GT_PARSE_H
#define VCF_GT_PARSE_H

#include <stdbool.h>
#include <stddef.h>

#include "vcf/marker.h"
#include "vcf/vcf_header.h"

/* The walk over a record's sample columns that VcfRecGTParser's methods share:
 * columns of excluded samples are skipped without being parsed. */
typedef struct {
    const char *rec;
    size_t len;
    const vcf_header *h;
    const marker *m;
    bool is_ref;   /* selects Java's "reference sample" error wording */
    long tab;      /* the tab before the next column, or -1 */
    int unfilt;
} gt_fields;

/* The GT subfield of one sample: the first allele is rec[start, end1) and a
 * second allele, if any, is rec[end1 + 1, end2). */
typedef struct {
    size_t start, end1, end2;
    bool diploid;
} gt_field;

void gt_fields_init(gt_fields *f, const char *rec, size_t len, const vcf_header *h, const marker *m, bool is_ref);
/* Moves to sample s, exiting on too few columns or an empty GT field. */
void gt_fields_next(gt_fields *f, int s, gt_field *gt);
/* VcfRecGTParser.parseAllele: -1 for '.'. */
int gt_parse_allele(const gt_fields *f, size_t start, size_t end);
_Noreturn void gt_ploidy_error(const gt_fields *f, int s, bool diploid);
_Noreturn void gt_sample_error(const gt_fields *f, const char *msg, int s);

#endif
