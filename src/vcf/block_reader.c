/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/BlockLineReader.java and
 * the parseLines step of vcf/RefIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/block_reader.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/parallel.h"
#include "blbutil/utilities.h"
#include "jcompat/jutf8.h"

/* Lines parsed together on the worker threads, and the batches in flight:
 * one being read, one being parsed, one being consumed and one buffered
 * between stages. */
#define BLOCK_READER_BATCH 1024
#define BLOCK_READER_SLOTS 4

typedef struct {
    kstring_t lines[BLOCK_READER_BATCH];
    ref_gt_rec *recs[BLOCK_READER_BATCH];
    int n;   /* 0 marks the end of the file, as BlockLineReader.SENTINAL does */
} batch;

struct block_reader {
    line_reader *reader;
    const vcf_header *header;
    kstring_t line;       /* the reader thread's lookahead line */
    bool have_line;
    int n_threads;
    pthread_t reader_thread;
    pthread_t parser_thread;
    pthread_mutex_t mutex;
    pthread_cond_t changed;
    batch slots[BLOCK_READER_SLOTS];
    batch *free_slots[BLOCK_READER_SLOTS];
    int n_free;
    batch *read[BLOCK_READER_SLOTS];   /* batches of lines in file order, not yet parsed */
    int read_head, n_read;
    batch *full[BLOCK_READER_SLOTS];   /* parsed batches in file order */
    int full_head, n_full;
    bool stop;
    batch *cur;           /* the batch being consumed */
    int cur_next;
};

typedef struct {
    const vcf_header *header;
    batch *b;
} parse_ctx;

static void parse_task(void *worker, int j) {
    const parse_ctx *c = worker;
    ref_gt_rec *rec = util_malloc(sizeof *rec);
    jutf8_sanitize(&c->b->lines[j]);
    ref_gt_rec_parse(rec, c->b->lines[j].s, c->b->lines[j].l, c->header);
    rec->refs = 1;
    c->b->recs[j] = rec;
}

/* Indexes the line's chromosome id, sanitised as the parser will see it. The
 * lines are sanitised on the parse workers; only the id before the first tab
 * is needed here, and it is ASCII for every real chromosome name. */
static void index_chrom(const char *line, size_t len) {
    const char *tab = memchr(line, '\t', len);
    if (tab == NULL || tab == line) return;
    size_t id_len = (size_t)(tab - line);
    bool ascii = true;
    for (size_t i = 0; i < id_len && ascii; ++i) ascii = (unsigned char)line[i] < 0x80;
    if (ascii) {
        chrom_ids_index(line, id_len);
        return;
    }
    kstring_t id = {0, 0, NULL};
    kputsn(line, id_len, &id);
    jutf8_sanitize(&id);
    chrom_ids_index(id.s, id.l);
    free(id.s);
}

/* BlockLineReader.fillBuffer: reads the next lines. Chromosome indices are
 * assigned here, in file order, so the parsers only look them up. */
static void read_lines(block_reader *r, batch *b) {
    int n = 0;
    while (n < BLOCK_READER_BATCH && r->have_line) {
        kstring_t t = b->lines[n];
        b->lines[n++] = r->line;
        r->line = t;
        r->line.l = 0;
        r->have_line = line_reader_next_raw(r->reader, &r->line);
    }
    for (int j = 0; j < n; ++j) index_chrom(b->lines[j].s, b->lines[j].l);
    b->n = n;
}

/* VcfIt.fillEmissionBuffer: parses a batch's lines in parallel. */
static void parse_batch(block_reader *r, batch *b) {
    parse_ctx ctx = {r->header, b};
    parallel_for(parallel_threads(r->n_threads, b->n), b->n, &ctx, 0, parse_task);
}

/* The reader thread reads lines into free batches and the parser thread
 * parses them, so reading overlaps parsing. Both run until the batch with no
 * lines that marks the end of the file, or until block_reader_close stops
 * them. */
static void *read_batches(void *arg) {
    block_reader *r = arg;
    for (;;) {
        pthread_mutex_lock(&r->mutex);
        while (r->n_free == 0 && !r->stop) pthread_cond_wait(&r->changed, &r->mutex);
        if (r->stop) {
            pthread_mutex_unlock(&r->mutex);
            return NULL;
        }
        batch *b = r->free_slots[--r->n_free];
        pthread_mutex_unlock(&r->mutex);
        read_lines(r, b);
        pthread_mutex_lock(&r->mutex);
        r->read[(r->read_head + r->n_read++) % BLOCK_READER_SLOTS] = b;
        pthread_cond_broadcast(&r->changed);
        pthread_mutex_unlock(&r->mutex);
        if (b->n == 0) return NULL;
    }
}

static void *parse_batches(void *arg) {
    block_reader *r = arg;
    for (;;) {
        pthread_mutex_lock(&r->mutex);
        while (r->n_read == 0 && !r->stop) pthread_cond_wait(&r->changed, &r->mutex);
        if (r->stop) {
            pthread_mutex_unlock(&r->mutex);
            return NULL;
        }
        batch *b = r->read[r->read_head];
        r->read_head = (r->read_head + 1) % BLOCK_READER_SLOTS;
        --r->n_read;
        pthread_mutex_unlock(&r->mutex);
        parse_batch(r, b);
        /* The consumer may recycle b as soon as it is published. */
        bool at_end = b->n == 0;
        pthread_mutex_lock(&r->mutex);
        r->full[(r->full_head + r->n_full++) % BLOCK_READER_SLOTS] = b;
        pthread_cond_broadcast(&r->changed);
        pthread_mutex_unlock(&r->mutex);
        if (at_end) return NULL;
    }
}

block_reader *block_reader_open(line_reader *reader, kstring_t line, const vcf_header *header, int n_threads) {
    block_reader *r = util_malloc(sizeof *r);
    *r = (block_reader){0};
    r->reader = reader;
    r->header = header;
    r->line = line;
    r->have_line = true;
    r->n_threads = n_threads;
    pthread_mutex_init(&r->mutex, NULL);
    pthread_cond_init(&r->changed, NULL);
    for (int j = 0; j < BLOCK_READER_SLOTS; ++j) r->free_slots[r->n_free++] = &r->slots[j];
    if (pthread_create(&r->reader_thread, NULL, read_batches, r) != 0) util_exit("fast-beagle: cannot create thread");
    if (pthread_create(&r->parser_thread, NULL, parse_batches, r) != 0) util_exit("fast-beagle: cannot create thread");
    return r;
}

ref_gt_rec *block_reader_next(block_reader *r) {
    if (r->cur != NULL && r->cur->n == 0) return NULL;
    if (r->cur == NULL || r->cur_next == r->cur->n) {
        pthread_mutex_lock(&r->mutex);
        if (r->cur != NULL) r->free_slots[r->n_free++] = r->cur;
        while (r->n_full == 0) pthread_cond_wait(&r->changed, &r->mutex);
        r->cur = r->full[r->full_head];
        r->full_head = (r->full_head + 1) % BLOCK_READER_SLOTS;
        --r->n_full;
        pthread_cond_broadcast(&r->changed);
        pthread_mutex_unlock(&r->mutex);
        r->cur_next = 0;
        if (r->cur->n == 0) return NULL;
    }
    return r->cur->recs[r->cur_next++];
}

void block_reader_close(block_reader *r) {
    pthread_mutex_lock(&r->mutex);
    r->stop = true;
    pthread_cond_broadcast(&r->changed);
    pthread_mutex_unlock(&r->mutex);
    pthread_join(r->reader_thread, NULL);
    pthread_join(r->parser_thread, NULL);
    if (r->cur != NULL) {
        for (int j = r->cur_next; j < r->cur->n; ++j) ref_gt_rec_release(r->cur->recs[j]);
    }
    for (int k = 0; k < r->n_full; ++k) {
        batch *b = r->full[(r->full_head + k) % BLOCK_READER_SLOTS];
        for (int j = 0; j < b->n; ++j) ref_gt_rec_release(b->recs[j]);
    }
    for (int k = 0; k < BLOCK_READER_SLOTS; ++k) {
        for (int j = 0; j < BLOCK_READER_BATCH; ++j) free(r->slots[k].lines[j].s);
    }
    pthread_mutex_destroy(&r->mutex);
    pthread_cond_destroy(&r->changed);
    free(r->line.s);
    free(r);
}
