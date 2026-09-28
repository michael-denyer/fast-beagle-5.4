/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) bref/Bref3It.java and
 * bref/Bref3Reader.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "bref/bref3_it.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "vcf/filter_util.h"
#include "vcf/ref_gt_rec.h"

#define MAGIC_NUMBER_V3 2055763188
#define SEQ_CODED 0
#define ALLELE_CODED 1

typedef struct {
    FILE *f;
    char *path;
    const str_set *exclude_markers;
    int n_haps;
    samples samples;
    uint8_t *bytes;          /* 2 bytes per haplotype */
    ref_gt_rec **buf;        /* the current block's records not yet returned */
    int buf_head, buf_n, buf_cap;
    bool at_end;             /* the end-of-data block has been read */
    kstring_t str, prefix;
} bref3_it;

static void read_fully(bref3_it *it, void *p, size_t n) {
    if (fread(p, 1, n, it->f) != n) util_exit("\nError reading file\n\njava.io.EOFException: %s", it->path);
}

static int32_t read_int(bref3_it *it) {
    uint8_t b[4];
    read_fully(it, b, 4);
    return (int32_t)((uint32_t)b[0] << 24 | (uint32_t)b[1] << 16 | (uint32_t)b[2] << 8 | b[3]);
}

static int read_unsigned_short(bref3_it *it) {
    uint8_t b[2];
    read_fully(it, b, 2);
    return b[0] << 8 | b[1];
}

static int read_unsigned_byte(bref3_it *it) {
    uint8_t b;
    read_fully(it, &b, 1);
    return b;
}

static void put_code_point(kstring_t *s, unsigned c) {
    if (c < 0x80) {
        kputc((int)c, s);
    } else if (c < 0x800) {
        kputc((int)(0xc0 | c >> 6), s);
        kputc((int)(0x80 | (c & 0x3f)), s);
    } else if (c < 0x10000) {
        kputc((int)(0xe0 | c >> 12), s);
        kputc((int)(0x80 | (c >> 6 & 0x3f)), s);
        kputc((int)(0x80 | (c & 0x3f)), s);
    } else {
        kputc((int)(0xf0 | c >> 18), s);
        kputc((int)(0x80 | (c >> 12 & 0x3f)), s);
        kputc((int)(0x80 | (c >> 6 & 0x3f)), s);
        kputc((int)(0x80 | (c & 0x3f)), s);
    }
}

/* DataInput.readUTF, converted to the UTF-8 that Java's output writers print:
 * a surrogate pair becomes one 4-byte sequence and a lone surrogate becomes
 * '?'. A NUL character, which the C strings cannot hold, is an error. */
static void read_utf(bref3_it *it, kstring_t *out) {
    int n = read_unsigned_short(it);
    uint8_t *b = util_malloc((size_t)n);
    read_fully(it, b, (size_t)n);
    unsigned *units = util_malloc((size_t)n * sizeof *units);
    int n_units = 0;
    for (int i = 0; i < n;) {
        unsigned c = b[i];
        if (c < 0x80) {
            units[n_units++] = c;
            i += 1;
        } else if ((c & 0xe0) == 0xc0 && i + 1 < n && (b[i + 1] & 0xc0) == 0x80) {
            units[n_units++] = (c & 0x1f) << 6 | (b[i + 1] & 0x3f);
            i += 2;
        } else if ((c & 0xf0) == 0xe0 && i + 2 < n && (b[i + 1] & 0xc0) == 0x80 && (b[i + 2] & 0xc0) == 0x80) {
            units[n_units++] = (c & 0x0f) << 12 | (b[i + 1] & 0x3fu) << 6 | (b[i + 2] & 0x3f);
            i += 3;
        } else {
            util_exit("java.io.UTFDataFormatException: malformed input in %s", it->path);
        }
    }
    out->l = 0;
    for (int j = 0; j < n_units; ++j) {
        unsigned u = units[j];
        if (u == 0) util_exit("fast-beagle: a string in %s contains a NUL character", it->path);
        if (u >= 0xd800 && u <= 0xdbff && j + 1 < n_units && units[j + 1] >= 0xdc00 && units[j + 1] <= 0xdfff) {
            put_code_point(out, 0x10000 + ((u - 0xd800) << 10) + (units[j + 1] - 0xdc00));
            ++j;
        } else if (u >= 0xd800 && u <= 0xdfff) {
            kputc('?', out);
        } else {
            put_code_point(out, u);
        }
    }
    if (out->s == NULL) kputs("", out);
    free(units);
    free(b);
}

/* The Bref3Reader constructor */
static void read_header(bref3_it *it) {
    if (read_int(it) != MAGIC_NUMBER_V3) {
        util_exit("\nERROR: Unrecognized input file.  Was input file created \n"
                "with a different version of the bref program?\n\nTerminating program.");
    }
    read_utf(it, &it->str);   /* the program string */
    int n_ids = read_int(it);
    if (n_ids < 0) util_exit("java.lang.NullPointerException: %s has no sample list", it->path);
    if (n_ids > INT_MAX / 4) {
        /* Java reads every ID before new byte[2*nHaps] (nHaps = 2*nSamples)
         * throws on the overflowed size. */
        for (int j = 0; j < n_ids; ++j) read_utf(it, &it->str);
        util_exit("java.lang.NegativeArraySizeException: %d", (int)((unsigned)n_ids << 2));
    }
    char **ids = util_malloc((size_t)(n_ids > 0 ? n_ids : 1) * sizeof *ids);
    int n = n_ids;
    for (int j = 0; j < n; ++j) {
        read_utf(it, &it->str);
        ids[j] = util_strndup(it->str.s, it->str.l);
    }
    bool *is_diploid = util_malloc((size_t)n * sizeof *is_diploid);
    for (int j = 0; j < n; ++j) is_diploid[j] = true;
    samples_init(&it->samples, n, ids, is_diploid);
    it->n_haps = n << 1;
}

static void buf_add(bref3_it *it, ref_gt_rec *rec) {
    if (it->buf_n == it->buf_cap) {
        it->buf_cap = it->buf_cap == 0 ? 512 : 2 * it->buf_cap;
        it->buf = util_realloc(it->buf, (size_t)it->buf_cap * sizeof *it->buf);
    }
    it->buf[it->buf_n++] = rec;
}

/* Bref3Reader.readByteLengthStringArray joined with ';', leaving out the
 * empty and "." IDs that the BasicMarker constructor removes. */
static void read_ids(bref3_it *it, kstring_t *out) {
    int n = read_unsigned_byte(it);
    size_t start = out->l;
    for (int j = 0; j < n; ++j) {
        read_utf(it, &it->str);
        if (it->str.l == 0 || (it->str.l == 1 && it->str.s[0] == '.')) continue;
        if (out->l > start) kputc(';', out);
        kputsn(it->str.s, it->str.l, out);
    }
    if (out->l == start) kputc('.', out);
}

/* Bref3Reader.readMarker: the marker is parsed from a VCF record prefix
 * holding the fields that Java passes to the BasicMarker constructor. */
static void read_marker(bref3_it *it, marker *m, const char *chrom) {
    static const char BASES[] = "ACGT";
    kstring_t *s = &it->prefix;
    s->l = 0;
    int32_t pos = read_int(it);
    ksprintf(s, "%s\t%d\t", chrom, pos);
    read_ids(it, s);
    kputc('\t', s);
    int allele_code = (int8_t)read_unsigned_byte(it);
    int32_t end = -1;
    int n_alleles;
    if (allele_code == -1) {
        n_alleles = read_int(it);
        if (n_alleles <= 0) util_exit("java.lang.NullPointerException: a marker in %s has no alleles", it->path);
        for (int j = 0; j < n_alleles; ++j) {
            read_utf(it, &it->str);
            if (j > 0) kputc(j == 1 ? '\t' : ',', s);
            kputsn(it->str.s, it->str.l, s);
        }
        end = read_int(it);
    } else {
        /* The first nAlleles bases of permutation permIndex of ACGT, the 24
         * permutations in lexicographic order. */
        n_alleles = 1 + (allele_code & 0x3);
        int perm = allele_code >> 2;
        if (perm < 0 || perm >= 24) util_exit("java.lang.ArrayIndexOutOfBoundsException: SNV code %d in %s", allele_code, it->path);
        char avail[4];
        memcpy(avail, BASES, 4);
        int n_avail = 4, fact = 6;
        for (int j = 0; j < n_alleles; ++j) {
            int k = perm / fact;
            perm %= fact;
            if (n_avail > 1) fact /= n_avail - 1;
            if (j > 0) kputc(j == 1 ? '\t' : ',', s);
            kputc(avail[k], s);
            memmove(avail + k, avail + k + 1, (size_t)(n_avail - k - 1));
            --n_avail;
        }
    }
    if (n_alleles == 1) kputs("\t.", s);
    kputs("\t.\t.\t", s);
    if (end >= 0) ksprintf(s, "END=%d", end);
    else kputc('.', s);
    kputs("\tGT", s);
    marker_parse(m, s->s, s->l);
}

static seq_group *read_hap_to_seq(bref3_it *it, int n_seq) {
    read_fully(it, it->bytes, 2 * (size_t)it->n_haps);
    seq_group *g = seq_group_new(it->n_haps, n_seq);
    for (int k = 0; k < it->n_haps; ++k) {
        g->hap_to_seq[k] = it->bytes[2 * k] << 8 | it->bytes[2 * k + 1];
        if (n_seq > 0 && g->hap_to_seq[k] >= n_seq) util_exit("fast-beagle: inconsistent data in %s", it->path);
    }
    return g;
}

/* Bref3Reader.readHapRecord */
static void read_hap_record(bref3_it *it, ref_gt_rec *rec, seq_group *g) {
    if (g->n_seq == 0) util_exit("fast-beagle: inconsistent data in %s", it->path);
    int n_alleles = marker_n_alleles(&rec->marker);
    uint8_t *seq_to_allele = util_malloc((size_t)g->n_seq);
    read_fully(it, seq_to_allele, (size_t)g->n_seq);
    for (int s = 0; s < g->n_seq; ++s) {
        if (seq_to_allele[s] >= n_alleles) util_exit("fast-beagle: inconsistent data in %s", it->path);
    }
    ref_gt_rec_set_seq_coded(rec, g, seq_to_allele);
}

/* Bref3Reader.readHapCodedRec and RefGTRec.hapCodedInstance: the major allele
 * is the one whose list is absent. The lists are checked after all are read. */
static void read_allele_record(bref3_it *it, ref_gt_rec *rec) {
    int n_alleles = marker_n_alleles(&rec->marker);
    rec->kind = n_alleles == 2 ? REF_TWO_ALLELE : REF_ALLELE;
    rec->hap_list_len = util_malloc((size_t)n_alleles * sizeof *rec->hap_list_len);
    rec->hap_lists = util_malloc((size_t)n_alleles * sizeof *rec->hap_lists);
    rec->major_allele = -1;
    for (int a = 0; a < n_alleles; ++a) {
        int len = read_int(it);
        rec->hap_list_len[a] = 0;
        rec->hap_lists[a] = NULL;
        if (len == -1) continue;
        if (len < 0) util_exit("java.lang.NegativeArraySizeException: %d", len);
        int *list = util_malloc((size_t)(len > 0 ? len : 1) * sizeof *list);
        for (int c = 0; c < len; ++c) list[c] = read_int(it);
        rec->hap_list_len[a] = len;
        rec->hap_lists[a] = list;
    }
    /* LowMafRefGTRec.checkIndicesAndReturnNullIndex */
    for (int a = 0; a < n_alleles; ++a) {
        if (rec->hap_lists[a] == NULL) {
            if (rec->major_allele != -1) util_exit("java.lang.IllegalArgumentException: invalid array");
            rec->major_allele = a;
            continue;
        }
        const int *list = rec->hap_lists[a];
        int len = rec->hap_list_len[a];
        if (len > 0 && (list[0] < 0 || list[len - 1] >= it->n_haps)) {
            util_exit("java.lang.IllegalArgumentException: invalid array");
        }
        for (int k = 1; k < len; ++k) {
            if (list[k - 1] >= list[k]) util_exit("java.lang.IllegalArgumentException: invalid array");
        }
    }
    if (rec->major_allele == -1) util_exit("java.lang.IllegalArgumentException: invalid array");
}

/* Bref3Reader.readBlock(dataIn, buffer, nRecs) */
static void read_block(bref3_it *it, int n_recs) {
    read_utf(it, &it->str);
    int chrom_index = chrom_ids_index(it->str.s, it->str.l);
    const char *chrom = chrom_ids_id(chrom_index);
    int n_seq = read_unsigned_short(it);
    seq_group *g = read_hap_to_seq(it, n_seq);
    for (int j = 0; j < n_recs; ++j) {
        ref_gt_rec *rec = util_malloc(sizeof *rec);
        memset(rec, 0, sizeof *rec);
        rec->refs = 1;
        rec->n_haps = it->n_haps;
        read_marker(it, &rec->marker, chrom);
        int flag = (int8_t)read_unsigned_byte(it);
        if (flag == SEQ_CODED) read_hap_record(it, rec, g);
        else if (flag == ALLELE_CODED) read_allele_record(it, rec);
        else util_exit("\nError reading file\n\nTerminating program.");
        if (filter_accept_marker(it->exclude_markers, &rec->marker)) buf_add(it, rec);
        else ref_gt_rec_release(rec);
    }
    seq_group_release(g);
}

/* Bref3Reader.readBlock(bref, buffer): reads blocks until one has a record
 * that passes the filter or the end-of-data block is read. */
static void fill(bref3_it *it) {
    it->buf_head = it->buf_n = 0;
    while (it->buf_n == 0 && !it->at_end) {
        int n_recs = read_int(it);
        if (n_recs == 0) it->at_end = true;
        else read_block(it, n_recs);
    }
}

static const sample_file_it_ops bref3_it_ops;

sample_file_it bref3_it_open(const char *path, const str_set *exclude_markers) {
    bref3_it *it = util_malloc(sizeof *it);
    *it = (bref3_it){0};
    it->path = util_strndup(path, strlen(path));
    it->f = fopen(path, "rb");
    if (it->f == NULL) util_exit("\nError: file not found [%s]\n\njava.io.FileNotFoundException: %s", path, path);
    setvbuf(it->f, NULL, _IOFBF, 1 << 20);
    it->exclude_markers = exclude_markers;
    read_header(it);
    it->bytes = util_malloc(2 * (size_t)(it->n_haps > 0 ? it->n_haps : 1));
    fill(it);
    if (trace_on()) {
        for (int j = 0; j < it->samples.n; ++j) trace_line("T1b-ref", "%s\tdiploid", it->samples.ids[j]);
    }
    return (sample_file_it){&bref3_it_ops, it};
}

static const samples *bref3_it_samples(const void *self) {
    const bref3_it *it = self;
    return &it->samples;
}

static void *bref3_it_next(void *self) {
    bref3_it *it = self;
    if (it->buf_head == it->buf_n) return NULL;
    ref_gt_rec *rec = it->buf[it->buf_head++];
    if (trace_on()) ref_gt_rec_trace(rec);
    if (it->buf_head == it->buf_n) fill(it);
    return rec;
}

static void bref3_it_close(void *self) {
    bref3_it *it = self;
    for (int j = it->buf_head; j < it->buf_n; ++j) ref_gt_rec_release(it->buf[j]);
    free(it->buf);
    fclose(it->f);
    samples_free(&it->samples);
    free(it->bytes);
    free(it->str.s);
    free(it->prefix.s);
    free(it->path);
    free(it);
}

static const sample_file_it_ops bref3_it_ops = {.samples = bref3_it_samples, .next = bref3_it_next,
        .close = bref3_it_close, .marker = ref_gt_rec_it_marker, .release = ref_gt_rec_it_release};
