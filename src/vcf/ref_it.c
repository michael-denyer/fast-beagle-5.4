/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) vcf/RefIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/ref_it.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/line_reader.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "bref/seq_coder3.h"
#include "jcompat/jnum.h"
#include "vcf/block_reader.h"
#include "vcf/filter_util.h"
#include "vcf/ref_gt_rec.h"
#include "vcf/vcf_header.h"

/* RefIt.DEFAULT_BUFFER_SIZE: the lines in each BlockLineReader block. */
#define REF_IT_BLOCK_LINES 1024

/* A FIFO of owned record pointers; NULL entries are placeholders. */
typedef struct {
    ref_gt_rec **v;
    int head, n, cap;
} rec_queue;

static void queue_push(rec_queue *q, ref_gt_rec *rec) {
    if (q->head + q->n == q->cap) {
        if (q->head > 0) {
            memmove(q->v, q->v + q->head, (size_t)q->n * sizeof *q->v);
            q->head = 0;
        }
        if (q->n == q->cap) {
            q->cap = q->cap == 0 ? 64 : 2 * q->cap;
            q->v = util_realloc(q->v, (size_t)q->cap * sizeof *q->v);
        }
    }
    q->v[q->head + q->n++] = rec;
}

static ref_gt_rec *queue_pop(rec_queue *q) {
    --q->n;
    return q->v[q->head++];
}

typedef struct {
    line_reader *reader;
    vcf_header header;
    block_reader *blocks;
    const str_set *exclude;
    seq_coder3 *coder;
    int max_seq_coded_alleles;
    int max_seq_coding_major_cnt;
    int last_chrom;
    int n_blocks;         /* BlockLineReader blocks coded so far */
    rec_queue low_freq;   /* records in file order, NULL where SeqCoder3 holds one */
    rec_queue ready;      /* records that next() can return */
} ref_it;

/* RefIt.flushCompressedRecords: fill the placeholders with the group's
 * sequence-coded records and release everything buffered, in order. */
static void flush_compressed(ref_it *it) {
    ref_gt_rec **coded = util_malloc((size_t)(seq_coder3_n_recs(it->coder) + 1) * sizeof *coded);
    int n_coded = seq_coder3_flush(it->coder, coded);
    int index = 0;
    while (it->low_freq.n > 0) {
        ref_gt_rec *rec = queue_pop(&it->low_freq);
        if (rec == NULL) {
            /* A record that SeqCoder3 refused left a slot with no coded
             * record: List.get's exception, or Collections.emptyList's. */
            if (n_coded == 0) util_exit("java.lang.IndexOutOfBoundsException: Index: 0");
            if (index == n_coded) {
                util_exit("java.lang.IndexOutOfBoundsException: Index %d out of bounds for length %d", index, n_coded);
            }
            rec = coded[index++];
        }
        queue_push(&it->ready, rec);
    }
    it->low_freq.head = 0;
    free(coded);
}

/* RefIt.applySeqCoding */
static bool apply_seq_coding(const ref_it *it, const ref_gt_rec *rec) {
    int n_alleles = marker_n_alleles(&rec->marker);
    if (n_alleles > it->max_seq_coded_alleles) return false;
    int maj_cnt = rec->n_haps;
    for (int a = 0; a < n_alleles; ++a) {
        if (a != rec->major_allele) maj_cnt -= rec->hap_list_len[a];
    }
    return maj_cnt <= it->max_seq_coding_major_cnt;
}

/* The loop body of RefIt.fillRecBuffer for one parsed record. */
static void add_rec(ref_it *it, ref_gt_rec *rec) {
    if (!filter_accept_marker(it->exclude, &rec->marker)) {
        ref_gt_rec_release(rec);
        return;
    }
    int chrom = rec->marker.chrom_index;
    if (it->last_chrom == -1) it->last_chrom = chrom;
    if (chrom != it->last_chrom) {
        flush_compressed(it);
        it->last_chrom = chrom;
    }
    if (!apply_seq_coding(it, rec)) {
        queue_push(&it->low_freq, rec);
    } else {
        if (!seq_coder3_add(it->coder, rec)) {
            flush_compressed(it);
            /* RefIt only asserts that the retry succeeds: without -ea the
             * record is dropped and the next flush throws. */
            if (!seq_coder3_add(it->coder, rec)) ref_gt_rec_release(rec);
        }
        queue_push(&it->low_freq, NULL);
    }
}

/* RefIt.fillRecBuffer: codes whole blocks of BlockLineReader lines, the first
 * one line longer (RefIt.combine), until a flush has made records ready. A
 * flush that throws anywhere in a block therefore throws when Java's does,
 * before the records ahead of it are used. */
static void fill(ref_it *it) {
    while (it->ready.n == 0) {
        int block = it->n_blocks++ == 0 ? REF_IT_BLOCK_LINES + 1 : REF_IT_BLOCK_LINES;
        int n = 0;
        for (ref_gt_rec *rec; n < block && (rec = block_reader_next(it->blocks)) != NULL; ++n) add_rec(it, rec);
        if (n == 0) {   /* BlockLineReader.SENTINAL */
            flush_compressed(it);
            return;
        }
    }
}

static const sample_file_it_ops ref_it_ops;

/* RefIt.create(...): reading, parsing and sequence coding the VCF records. */
sample_file_it ref_it_open(const char *path, const str_set *exclude_samples, const str_set *exclude_markers,
        int n_threads) {
    ref_it *it = util_malloc(sizeof *it);
    *it = (ref_it){0};
    it->reader = line_reader_open(path, n_threads);
    it->exclude = exclude_markers;
    it->last_chrom = -1;
    kstring_t line = {0, 0, NULL};
    vcf_header_read(&it->header, it->reader, &line, exclude_samples);
    if (trace_on()) vcf_header_trace(&it->header, "T1b-ref");

    int n_samples = it->header.samples.n;
    it->coder = seq_coder3_new(n_samples, seq_coder3_default_max_nseq(n_samples));
    int max_nseq = seq_coder3_max_nseq(it->coder);
    it->max_seq_coded_alleles = max_nseq < SEQ_CODER3_MAX_NALLELES ? max_nseq : SEQ_CODER3_MAX_NALLELES;
    /* (int) Math.floor(nHaps*0.995f - 1): the product and difference are float. */
    float threshold = (float)(n_samples << 1) * SEQ_CODER3_COMPRESS_FREQ_THRESHOLD - 1.0f;
    it->max_seq_coding_major_cnt = jnum_d2i(floor((double)threshold));
    it->blocks = block_reader_open(it->reader, line, &it->header, n_threads);
    fill(it);
    return (sample_file_it){&ref_it_ops, it};
}

static const samples *ref_it_samples(const void *self) {
    const ref_it *it = self;
    return &it->header.samples;
}

/* RefIt.next */
static void *ref_it_next(void *self) {
    ref_it *it = self;
    if (it->ready.n == 0) return NULL;
    ref_gt_rec *rec = queue_pop(&it->ready);
    if (trace_on()) ref_gt_rec_trace(rec);
    if (it->ready.n == 0) fill(it);
    return rec;
}

static void ref_it_close(void *self) {
    ref_it *it = self;
    block_reader_close(it->blocks);
    while (it->ready.n > 0) ref_gt_rec_release(queue_pop(&it->ready));
    while (it->low_freq.n > 0) ref_gt_rec_release(queue_pop(&it->low_freq));
    free(it->ready.v);
    free(it->low_freq.v);
    seq_coder3_free(it->coder);
    vcf_header_free(&it->header);
    line_reader_close(it->reader);
    free(it);
}

static const sample_file_it_ops ref_it_ops = {.samples = ref_it_samples, .next = ref_it_next, .close = ref_it_close,
        .marker = ref_gt_rec_it_marker, .release = ref_gt_rec_it_release};
