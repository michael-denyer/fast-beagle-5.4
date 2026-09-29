/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) bref/SeqCoder3.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BREF_SEQ_CODER3_H
#define BREF_SEQ_CODER3_H

#include <stdbool.h>

#include "vcf/ref_gt_rec.h"

#define SEQ_CODER3_MAX_NALLELES 255
#define SEQ_CODER3_COMPRESS_FREQ_THRESHOLD 0.995f

/* Groups consecutive allele-coded reference records by the distinct allele
 * sequences their haplotypes carry, up to maxNSeq sequences per group. */
typedef struct seq_coder3 seq_coder3;

/* SeqCoder3.defaultMaxNSeq, using StrictMath as Java's Math does at this site. */
int seq_coder3_default_max_nseq(int n_samples);
seq_coder3 *seq_coder3_new(int n_samples, int max_nseq);
int seq_coder3_max_nseq(const seq_coder3 *c);
/* SeqCoder3.add: takes ownership of rec and returns true if the group can hold
 * it; on false the coder is unchanged and the caller keeps rec. */
bool seq_coder3_add(seq_coder3 *c, ref_gt_rec *rec);
/* SeqCoder3.getCompressedList: converts the added records, in order, to
 * sequence-coded records sharing one group, and empties the coder. Returns the
 * number of records written to out, which must have room for all of them. */
int seq_coder3_flush(seq_coder3 *c, ref_gt_rec **out);
int seq_coder3_n_recs(const seq_coder3 *c);
void seq_coder3_free(seq_coder3 *c);

#endif
