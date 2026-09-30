/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/RefGTRec.java,
 * vcf/LowMafRefDiallelicGTRec.java, vcf/LowMafRefGTRec.java,
 * vcf/SeqCodedRefGTRec.java and vcf/VcfRecGTParser.java (phasedAlleles,
 * nonMajRefIndices); modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_REF_GT_REC_H
#define VCF_REF_GT_REC_H

#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>

#include "vcf/marker.h"
#include "vcf/vcf_header.h"

/* The haplotype-to-sequence map that SeqCoder3 shares between the
 * SeqCodedRefGTRec records of one group. Later stages compare groups by
 * identity. */
typedef struct {
    atomic_int refs;
    int n_haps;
    int *hap_to_seq;
    int n_seq;
    int trace_id;   /* order of first appearance in trace seam T1d, or -1 */
} seq_group;

/* A group owned by the caller, its hap_to_seq allocated but not filled. Each
 * record that shares the group retains it (ref_gt_rec_set_seq_coded); the
 * creator releases its own reference once the records are built. */
seq_group *seq_group_new(int n_haps, int n_seq);
seq_group *seq_group_retain(seq_group *g);
void seq_group_release(seq_group *g);

typedef enum { REF_TWO_ALLELE, REF_ALLELE, REF_HAP } ref_gt_rec_kind;

/* A phased, non-missing reference record. The allele-coded kinds keep the
 * increasing haplotype lists of each non-major allele; the sequence-coded kind
 * keeps its group and the allele of each sequence. */
typedef struct {
    atomic_int refs;        /* owners; see ref_gt_rec_release */
    marker marker;
    ref_gt_rec_kind kind;
    int n_haps;
    int major_allele;       /* allele-coded kinds */
    int *hap_list_len;      /* allele-coded kinds: 0 for the major allele */
    int **hap_lists;        /* allele-coded kinds */
    seq_group *group;       /* sequence-coded kind */
    uint8_t *seq_to_allele; /* sequence-coded kind, group->n_seq entries */
} ref_gt_rec;

/* RefGTRec.alleleRefGTRec(new VcfRecGTParser(header, line, parser)). Exits with
 * Java's message on an unphased or missing genotype or a format error. */
void ref_gt_rec_parse(ref_gt_rec *rec, const char *line, size_t len, const vcf_header *h);
/* SeqCoder3.getCompressedList and Bref3Reader.readHapRecord: replace an
 * allele-coded record's haplotype lists with a sequence coding. Retains g and
 * takes ownership of seq_to_allele, which has g->n_seq entries. */
void ref_gt_rec_set_seq_coded(ref_gt_rec *rec, seq_group *g, uint8_t *seq_to_allele);
/* For heap records shared by windows, whose owners can be on the window reader
 * and caller threads: add an owner, or drop one and free the record with its
 * last owner. */
ref_gt_rec *ref_gt_rec_retain(ref_gt_rec *rec);
void ref_gt_rec_release(ref_gt_rec *rec);

int ref_gt_rec_get(const ref_gt_rec *rec, int hap);
/* RefGTRec.majorAllele: for a sequence-coded record, the most frequent allele
 * with the lowest index. */
int ref_gt_rec_major_allele(const ref_gt_rec *rec);
const char *ref_gt_rec_class_name(const ref_gt_rec *rec);
/* Trace seams T1a-ref and T1d-ref: a record as a reference reader returns it.
 * A group's line precedes the first record of the group, and groups are
 * numbered in that order over the whole run. */
void ref_gt_rec_trace(const ref_gt_rec *rec);

/* The sample_file_it_ops record accessors of the reference readers. */
const marker *ref_gt_rec_it_marker(const void *rec);
void ref_gt_rec_it_release(void *rec);

#endif
