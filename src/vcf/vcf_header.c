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
#include "vcf/vcf_header.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/trace.h"
#include "blbutil/utilities.h"

#define SAMPLE_OFFSET 9
static const char HEADER_PREFIX[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

/* VcfMetaInfo: after trimming, "##" and then an '=' that is not the last char. */
static void check_meta_info(const char *line) {
    size_t start = 0, end = strlen(line);
    while (start < end && (unsigned char)line[start] <= ' ') ++start;
    while (end > start && (unsigned char)line[end - 1] <= ' ') --end;
    const char *s = line + start;
    size_t n = end - start;
    if (n < 2 || s[0] != '#' || s[1] != '#') {
        util_exit("VCF meta-information line: missing starting \"##\": %.*s", (int)n, s);
    }
    const char *eq = memchr(s, '=', n);
    if (eq == NULL || eq == s || (size_t)(eq - s) == n - 1) {
        util_exit("VCF meta-information line: missing \"=\"");
    }
}

/* VcfHeader.isDiploid: a sample column is diploid when it contains '/' or '|'
 * anywhere, not only in its GT subfield. */
static int is_diploid_columns(const char *rec, size_t len, bool **out) {
    const char *p = rec, *end = rec + len;
    for (int j = 0; j < 9; ++j) {
        p = memchr(p, '\t', (size_t)(end - p));
        if (p == NULL) util_exit("VCF record format error: %.*s", (int)len, rec);
        ++p;
    }
    int n = 1;
    for (const char *c = p; c < end; ++c) n += *c == '\t';
    bool *is_diploid = util_malloc((size_t)n * sizeof *is_diploid);
    int k = 0;
    bool sep = false;
    for (const char *c = p; c < end; ++c) {
        if (*c == '\t') {
            is_diploid[k++] = sep;
            sep = false;
        } else if (*c == '/' || *c == '|') {
            sep = true;
        }
    }
    is_diploid[k] = sep;
    *out = is_diploid;
    return n;
}

void vcf_header_init(vcf_header *h, const char *src, char **lines, int n_lines,
        const char *first_data_line, size_t first_data_len, const str_set *exclude) {
    bool *col_diploid;
    int n_cols = is_diploid_columns(first_data_line, first_data_len, &col_diploid);

    if (n_lines == 0) {
        util_exit("\n\nERROR: Missing the VCF meta information lines and the VCF header line\nVCF source: %s\n", src);
    }
    const char *header = lines[n_lines - 1];
    if (strncmp(header, HEADER_PREFIX, sizeof HEADER_PREFIX - 1) != 0) {
        util_exit("\n\nERROR: Missing the VCF header line.\nVCF source: %s\n"
                "The VCF header line must immediately follow the meta-information lines.\n"
                "The fields of the VCF header line must be tab-delimited and begin with:\n%s\n",
                src, HEADER_PREFIX);
    }
    for (int j = 0; j < n_lines - 1; ++j) check_meta_info(lines[j]);

    /* StringUtil.getFields(header, '\t') keeps empty fields. */
    int n_fields = 1;
    for (const char *c = header; *c != '\0'; ++c) n_fields += *c == '\t';
    const char **fields = util_malloc((size_t)n_fields * sizeof *fields);
    size_t *field_len = util_malloc((size_t)n_fields * sizeof *field_len);
    const char *f = header;
    for (int j = 0; j < n_fields; ++j) {
        const char *tab = strchr(f, '\t');
        fields[j] = f;
        field_len[j] = tab == NULL ? strlen(f) : (size_t)(tab - f);
        f = tab == NULL ? f + field_len[j] : tab + 1;
    }

    h->src = src;
    h->n_header_fields = n_fields;
    int n_unfiltered = n_fields > SAMPLE_OFFSET ? n_fields - SAMPLE_OFFSET : 0;
    h->included = util_malloc((size_t)n_unfiltered * sizeof *h->included);
    int n = 0;
    for (int j = 0; j < n_unfiltered; ++j) {
        const char *id = fields[SAMPLE_OFFSET + j];
        size_t len = field_len[SAMPLE_OFFSET + j];
        if (exclude == NULL || str_set_find(exclude, id, len) < 0) h->included[n++] = j;
    }
    if (n == 0) {
        util_exit("\nError      :  All samples in the VCF file are excluded\nFile       :  %s", src);
    }
    char **ids = util_malloc((size_t)n * sizeof *ids);
    bool *is_diploid = util_malloc((size_t)n * sizeof *is_diploid);
    for (int j = 0; j < n; ++j) {
        int col = h->included[j];
        ids[j] = util_strndup(fields[SAMPLE_OFFSET + col], field_len[SAMPLE_OFFSET + col]);
        if (col >= n_cols) util_exit("java.lang.ArrayIndexOutOfBoundsException: Index %d out of bounds for length %d", col, n_cols);
        is_diploid[j] = col_diploid[col];
    }
    samples_init(&h->samples, n, ids, is_diploid);
    free(fields);
    free(field_len);
    free(col_diploid);
}

void vcf_header_read(vcf_header *h, line_reader *r, kstring_t *first_line, const str_set *exclude) {
    char **lines = NULL;
    int n_lines = 0;
    bool have = line_reader_next(r, first_line);
    while (have && first_line->l > 0 && first_line->s[0] == '#') {
        lines = util_realloc(lines, (size_t)(n_lines + 1) * sizeof *lines);
        lines[n_lines++] = util_strndup(first_line->s, first_line->l);
        have = line_reader_next(r, first_line);
    }
    if (!have) util_exit("ERROR: missing VCF data lines (%s)", line_reader_name(r));
    vcf_header_init(h, line_reader_name(r), lines, n_lines, first_line->s, first_line->l, exclude);
    for (int j = 0; j < n_lines; ++j) free(lines[j]);
    free(lines);
}

void vcf_header_trace(const vcf_header *h, const char *seam) {
    const samples *s = &h->samples;
    for (int j = 0; j < s->n; ++j) {
        trace_line(seam, "%s\t%s", s->ids[j], s->is_diploid[j] ? "diploid" : "haploid");
    }
}

void vcf_header_free(vcf_header *h) {
    free(h->included);
    samples_free(&h->samples);
}
