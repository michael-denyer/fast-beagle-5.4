/*
 * Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang
 * Ported to C from PLINK 2.0 v2.0.0-a.7.8 plink2_import.cc (VCF dosage and
 * hardcall import), plink2_data.cc (ApplyHardCallThreshPhased),
 * plink2_pvar.cc (--extract-if-info), plink2_filter.cc (--maf) and
 * plink2_export.cc (--export bgen-1.2, .sample); modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BGEN_BGEN_WRITER_H
#define BGEN_BGEN_WRITER_H

#include <stddef.h>
#include <stdint.h>

#include "main/par.h"
#include "main/run_outputs.h"
#include "vcf/marker.h"
#include "vcf/samples.h"

/* The BGEN writer, storing each probability in N = bgen-bits bits:
 * - bgen=plink2: <out>.bgen and <out>.sample, byte-identical to
 *   plink2 --vcf <out>.vcf.gz dosage=DS --export bgen-1.2 bits=N
 *   [--extract-if-info "DR2 >= x"] [--maf x], with multiallelic records
 *   skipped as --import-max-alleles 2 skips them.
 * - bgen=phased: layout 2 with one probability per haplotype and allele,
 *   every chromosome and allele count, alleles in VCF order (REF first),
 *   each value in units of 1/(2^N - 1). */
typedef struct bgen_writer bgen_writer;

/* The n values v[0, n), each below 2^bits, as BGEN v1.2 packs them: value k
 * in bits k*bits to (k+1)*bits - 1 counting from the low bit of out[0], and
 * zero bits after the last value. Writes (n*bits + 7)/8 bytes. */
void bgen_pack_bits(const uint32_t *v, size_t n, int bits, uint8_t *out);

/* Probabilities p[0, n) summing to 1 as integers summing to max: each
 * p * max rounded down, then one more to the largest remainders, lower
 * alleles first on ties (the BGEN v1.2 specification's suggested method).
 * Writes the values for alleles 0 to n-2; BGEN leaves the last implied. */
void bgen_quantise(const float *p, int n, uint32_t max, uint32_t *out);

/* In a VCF INFO field (semicolon-separated entries, as printed): the value
 * of entry key=value, or "." if there is none; and whether a flag entry
 * that is exactly key is present. Keys match whole entries only. */
span bgen_info_value(span info, const char *key);
bool bgen_info_flag(span info, const char *key);

/* Working buffers for bgen_writer_encode; one per thread. */
typedef struct bgen_scratch bgen_scratch;

/* One VCF record's BGEN output: its sample calls until bgen_writer_encode,
 * then its .bgen bytes and .info row until bgen_writer_put. A record may be
 * begun again once put, keeping its output buffers. */
typedef struct bgen_rec bgen_rec;

/* s must outlive the writer. */
bgen_writer *bgen_writer_open(const par *p, const run_outputs *out, const samples *s);
/* For bgen=plink2, exits, removing the partial .bgen, unless chrom is an
 * autosome: plink2 needs sex information for the others. */
void bgen_writer_check_chrom(bgen_writer *bw, const char *chrom);
bgen_scratch *bgen_scratch_new(const bgen_writer *bw);
void bgen_scratch_free(bgen_scratch *sc);
/* Readies rec, or a new record when rec is NULL, for marker mk, which must
 * outlive bgen_writer_put. Set each sample's call before encoding. */
bgen_rec *bgen_rec_begin(const bgen_writer *bw, bgen_rec *rec, const marker *mk);
/* Sample s's fields as the VCF record prints them: its GT alleles (a2
 * unused for a haploid sample), its DS string (NULL when the record prints
 * no DS) and each haplotype's allele probabilities, summing to 1 (p2 NULL
 * for a haploid sample, both NULL when the GT alleles are certain). */
void bgen_rec_set_call(bgen_rec *rec, int s, int a1, int a2, const char *ds, const float *p1, const float *p2);
void bgen_rec_free(bgen_rec *rec);
/* Encodes a record whose calls are all set, given its INFO field as
 * printed, and releases the calls. Reads only the writer's settings, so
 * threads with their own scratch may encode records at once. */
void bgen_writer_encode(const bgen_writer *bw, bgen_scratch *sc, span info, bgen_rec *rec);
/* Writes an encoded record, in VCF order. */
void bgen_writer_put(bgen_writer *bw, bgen_rec *rec);
/* Writes the variant count and the .sample file. */
void bgen_writer_close(bgen_writer *bw);

#endif
