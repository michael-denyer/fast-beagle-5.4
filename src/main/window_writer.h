/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) main/WindowWriter.java,
 * vcf/VcfWriter.java, vcf/VcfRecBuilder.java and imp/ImputedRecBuilder.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef MAIN_WINDOW_WRITER_H
#define MAIN_WINDOW_WRITER_H

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "bgen/bgen_writer.h"
#include "main/par.h"
#include "main/vcf_index.h"
#include "vcf/samples.h"
#include "vcf/sliding_window.h"

/* The output VCF, <out>.vcf.gz, compressed with BGZF, its index when
 * tbi=true, and the BGEN output when bgen= is set. Only window_writer.c reads the fields. */
typedef struct {
    BGZF *out;
    char *path;
    const samples *samples;
    kstring_t line;
    bgen_writer *bgen;   /* NULL without bgen= */
    vcf_index *index;    /* NULL without tbi=true */
    uint64_t n_written;  /* uncompressed bytes */
    bool ap, gp;
    int n_haps;
    kstring_t hom_ref[5];
} window_writer;

/* new WindowWriter(par, samples): writes the meta-information and header
 * lines, including the DS and INFO lines Beagle writes even when it only
 * phases. Opens the BGEN output when bgen= is set. */
void window_writer_open(window_writer *ww, const par *p, const samples *s);
/* Before a window's work: fails now if the outputs cannot take its first
 * target marker. */
void window_writer_begin_window(window_writer *ww, const window *w);

/* A phased-only record has GT fields. Genotyped and imputed records have
 * dosage fields; only the latter carries the IMP flag. */
typedef enum { OUT_PHASED, OUT_GENOTYPED, OUT_IMPUTED } out_rec_kind;

/* One record under construction, then encoded and put in output order.
 * Only window_writer.c reads its fields. */
typedef struct {
    const window_writer *writer;
    const marker *mk;
    out_rec_kind kind;
    kstring_t info;
    kstring_t fields;
    bgen_rec *bgen;   /* NULL without bgen= */
    int n_alleles, n_samples, hap_cnt;
    float *sum_al_probs, *sum_al_probs2;
    uint64_t digest;
} out_rec;

/* Readies a zero-initialised record, or one already put, keeping its buffers.
 * ww and mk must outlive the record. All construction can run on workers. */
void window_writer_rec_begin(const window_writer *ww, out_rec *r, const marker *mk, out_rec_kind kind);
void window_writer_rec_free(out_rec *r);
/* Add samples in the writer's sample order. _gt is for OUT_PHASED; a2 is
 * ignored for haploid samples. _probs is for dosage records, scales its
 * arrays in place and captures both outputs before returning. p2 is ignored
 * for haploid samples. Arrays may be reused immediately after the call. */
void window_writer_rec_gt(out_rec *r, int a1, int a2);
void window_writer_rec_probs(out_rec *r, float *p1, float *p2);
/* T5d, emitted in record order after worker construction. */
void window_writer_rec_trace(const out_rec *r, int cluster, int marker_index);

/* One thread's working buffers for window_writer_encode. */
typedef struct out_worker out_worker;
out_worker *window_writer_worker_new(const window_writer *ww);
void window_writer_worker_free(out_worker *wk);
/* Completes INFO, checks the sample count and encodes BGEN on any thread,
 * once every sample has been added. */
void window_writer_encode(out_worker *wk, out_rec *r);
/* Writes an encoded record to every output. Records are put in VCF order. */
void window_writer_put(window_writer *ww, out_rec *r);

/* WindowWriter.printPhased: records for target markers [start, end) of w,
 * haploid samples with one allele. */
void window_writer_print_phased(window_writer *ww, const window *w, int start, int end, phased_allele_fn allele, const void *ctx);
void window_writer_close(window_writer *ww);

#endif
