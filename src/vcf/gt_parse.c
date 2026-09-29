/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/VcfRecGTParser.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/gt_parse.h"

#include <stdint.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jutf8.h"

static long index_of_tab(const gt_fields *f, long from) {
    if (from < 0 || (size_t)from >= f->len) return -1;
    const char *p = memchr(f->rec + from, '\t', f->len - (size_t)from);
    return p == NULL ? -1 : p - f->rec;
}

static _Noreturn void field_count_error(const gt_fields *f) {
    int n_fields = 1;
    for (size_t j = 0; j < f->len; ++j) n_fields += f->rec[j] == '\t';
    util_exit("ERROR: CF header line has %d fields, but data line has %d fields\nFile source: %s\n",
            f->h->n_header_fields, n_fields, f->h->src);
}

/* Marker.toString: CHROM, POS, ID, REF and ALT, tab-separated. */
static char *marker_string(const marker *m) {
    kstring_t s = {0, 0, NULL};
    span id = marker_id(m), al = marker_alleles(m);
    ksprintf(&s, "%s\t%d\t%.*s\t%.*s", marker_chrom(m), m->pos, id.n, id.s, al.n, al.s);
    return s.s;
}

void gt_fields_init(gt_fields *f, const char *rec, size_t len, const vcf_header *h, const marker *m, bool is_ref) {
    *f = (gt_fields){rec, len, h, m, is_ref, -1, -1};
    for (int j = 0; j < 9; ++j) {
        f->tab = index_of_tab(f, f->tab + 1);
        if (f->tab == -1) util_exit("VCF record format error: %.*s", (int)len, rec);
    }
}

void gt_sample_error(const gt_fields *f, const char *msg, int s) {
    util_exit("ERROR: %s%s%s at marker [%s]", msg, f->is_ref ? " for reference sample " : " for sample ",
            f->h->samples.ids[s], marker_string(f->m));
}

void gt_ploidy_error(const gt_fields *f, int s, bool diploid) {
    const samples *smp = &f->h->samples;
    util_exit("%s%s has an inconsistent number of alleles. The first genotype is %s, but the genotype at position %s:%d is %s",
            f->is_ref ? "Reference sample " : "Sample ", smp->ids[s], smp->is_diploid[s] ? "diploid" : "haploid",
            marker_chrom(f->m), f->m->pos, diploid ? "diploid" : "haploid");
}

void gt_fields_next(gt_fields *f, int s, gt_field *gt) {
    if (f->tab == -1) field_count_error(f);
    int next_unfiltered = f->h->included[s];
    while (++f->unfilt < next_unfiltered) {
        f->tab = index_of_tab(f, f->tab + 1);
        if (f->tab == -1) field_count_error(f);
    }
    /* VcfRecGTParser.alEnd1 and alEnd2 (exclusive ends). */
    size_t start = (size_t)f->tab + 1;
    if (start == f->len) {
        util_exit("ERROR: genotype is missing allele separator:\n%.*s\nExiting Program\n", (int)f->len, f->rec);
    }
    size_t end1 = start;
    while (end1 < f->len && f->rec[end1] != '/' && f->rec[end1] != '|' && f->rec[end1] != '\t' && f->rec[end1] != ':') {
        ++end1;
    }
    if (start == end1) gt_sample_error(f, "missing data", s);
    size_t end2 = end1;
    while (end2 < f->len && f->rec[end2] != ':' && f->rec[end2] != '\t') ++end2;
    *gt = (gt_field){start, end1, end2, end1 != end2};
    f->tab = index_of_tab(f, (long)end2);
}

/* One Java char is c - '0' with no digit check; longer text goes through
 * Integer.parseInt. */
int gt_parse_allele(const gt_fields *f, size_t start, size_t end) {
    if (start == end) util_exit("ERROR: Missing sample allele: %.*s", (int)f->len, f->rec);
    const char *s = f->rec + start;
    size_t n = end - start, one = 0;
    int32_t c = jutf8_next_bmp(s, n, &one);
    int32_t al;
    if (c >= 0 && one == n) {
        if (c == '.') return -1;
        al = c - '0';
    } else if (!jnum_parse_int(s, n, &al)) {
        util_exit("java.lang.NumberFormatException: For input string: \"%.*s\"", (int)n, s);
    }
    if (al < 0 || al >= marker_n_alleles(f->m)) {
        util_exit("ERROR: Invalid allele [%.*s] at character %zu in record \"%s\t...\"", (int)(end - start),
                f->rec + start, start, marker_string(f->m));
    }
    return al;
}
