/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/VcfRecGTParser.java (hapListRep),
 * vcf/VcfIt.java (TO_LOWMEM_GT_REC), vcf/LowMafDiallelicGTRec.java,
 * vcf/LowMafGTRec.java and vcf/BitArrayGTRec.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_GT_REC_H
#define VCF_GT_REC_H

#include <stdbool.h>
#include <stdint.h>

#include "vcf/marker.h"
#include "vcf/vcf_header.h"

/* The Java class TO_LOWMEM_GT_REC chooses: a low-frequency record keeps the
 * haplotypes carrying each non-major allele, any other keeps every allele in
 * a bit array. The choice is part of the trace because later stages branch on
 * the record class. */
typedef enum { GT_LOW_MAF_DIALLELIC, GT_LOW_MAF, GT_BIT_ARRAY } gt_rec_kind;

/* A target VCF record. A haploid sample's allele fills both of its haplotypes.
 * Phase is one flag per record: false if any genotype is unphased or missing. */
typedef struct {
    int refs;               /* owners; see gt_rec_release */
    marker marker;
    gt_rec_kind kind;
    bool is_phased;
    int n_haps;
    int major_allele;
    int n_missing;
    int *missing;           /* increasing sample indices with a missing allele */
    int *hap_list_len;      /* low-MAF kinds: per allele, 0 for the major allele */
    int **hap_lists;        /* low-MAF kinds: increasing haplotypes carrying each non-major allele */
    int bits_per_allele;    /* bit-array kind */
    uint64_t *allele_bits;  /* bit-array kind */
    uint64_t *missing_bits; /* bit-array kind: one bit per sample */
} gt_rec;

/* Parses a target record as TO_LOWMEM_GT_REC does, exiting with Java's message
 * on a format error. */
void gt_rec_parse(gt_rec *rec, const char *line, size_t len, const vcf_header *h);
void gt_rec_free(gt_rec *rec);
/* For heap records shared by windows: add an owner, or drop one and free the
 * record with its last owner. */
gt_rec *gt_rec_retain(gt_rec *rec);
void gt_rec_release(gt_rec *rec);

/* GTRec.get: the allele on a haplotype, or -1 if the sample has a missing allele. */
int gt_rec_get(const gt_rec *rec, int hap);
const char *gt_rec_class_name(const gt_rec *rec);

#endif
