/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/BasicMarker.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/marker.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jutf8.h"

#include <htslib/kstring.h>

#define ID_STORED ((uint16_t)(1 << 15))
#define ALLELES_STORED ((uint16_t)(1 << 14))
#define SNV_INDEX_MASK 0x7f

#define N_SNV_PERMS 100
#define MAX_SNV_ALLELES 5

/* A compact allele store taken from Beagle 5.5 MarkerUtils.snvPerms(), which
 * has no 5.4 counterpart: REF '\t' ALT for every ordering of {*, A, C, G, T}
 * whose REF is not '*', plus the four REF-only records, sorted. */
static char snv_perms[N_SNV_PERMS][2 * MAX_SNV_ALLELES];
static int n_snv_perms;

static void permute(char *start, int n_start, const char *end, int n_end) {
    if (n_end == 0 && start[0] != '*') {
        char *s = snv_perms[n_snv_perms++];
        int k = 0;
        s[k++] = start[0];
        for (int j = 1; j < n_start; ++j) {
            s[k++] = j == 1 ? '\t' : ',';
            s[k++] = start[j];
        }
        s[k] = '\0';
        return;
    }
    for (int j = 0; j < n_end; ++j) {
        start[n_start] = end[j];
        char rest[MAX_SNV_ALLELES];
        int n = 0;
        for (int k = 0; k < n_end; ++k) {
            if (k != j) rest[n++] = end[k];
        }
        permute(start, n_start + 1, rest, n);
    }
}

static int compare_str(const void *a, const void *b) {
    return strcmp((const char *)a, (const char *)b);
}

static void init_snv_perms(void) {
    char start[MAX_SNV_ALLELES];
    permute(start, 0, "*ACGT", 5);
    for (const char *b = "ACGT"; *b != '\0'; ++b) {
        char *s = snv_perms[n_snv_perms++];
        s[0] = *b;
        s[1] = '\t';
        s[2] = '.';
        s[3] = '\0';
    }
    qsort(snv_perms, N_SNV_PERMS, sizeof snv_perms[0], compare_str);
}

/* The text an SNV permutation index stands for: the REF '\t' ALT prefix of the
 * permutation that has n_alleles alleles. */
static span snv_alleles(int snv, int n_alleles) {
    const char *perm = snv_perms[snv];
    return (span){perm, n_alleles == 1 ? (int)strlen(perm) : 2 * n_alleles - 1};
}

/* The index of the permutation that stands for REF '\t' ALT s, or -1. The
 * binary search finds the first entry >= s, which starts with s if any does. */
static int snv_index(const char *s, size_t len, int n_alleles) {
    if (n_alleles > MAX_SNV_ALLELES) return -1;
    int lo = 0, hi = N_SNV_PERMS;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        int c = strncmp(snv_perms[mid], s, len);
        if (c == 0 && snv_perms[mid][len] != '\0') c = 1;  /* longer entry sorts after s */
        if (c < 0) lo = mid + 1;
        else hi = mid;
    }
    if (lo == N_SNV_PERMS) return -1;
    span sp = snv_alleles(lo, n_alleles);
    return (size_t)sp.n == len && memcmp(sp.s, s, len) == 0 ? lo : -1;
}

static const char *index_of(const char *s, const char *end, char c) {
    const void *p = memchr(s, c, (size_t)(end - s));
    return p;
}

/* The Java char at s[*i] of UTF-8 that the line reader sanitized, advancing *i
 * past it: its code point, or -1 for a supplementary character, which is a
 * surrogate pair (the only 4-byte sequence), two Java chars that none of the
 * checks below accepts. */
static int32_t next_char(const char *s, size_t len, size_t *i) {
    int32_t c = jutf8_next_bmp(s, len, i);
    if (c < 0) *i += (unsigned char)s[*i] >= 0xf0 && len - *i >= 4 ? 4 : 1;
    return c;
}

/* Character.isWhitespace */
static bool is_java_whitespace(int32_t c) {
    return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f) || c == 0x1680
            || (c >= 0x2000 && c <= 0x200a && c != 0x2007) || c == 0x2028 || c == 0x2029 || c == 0x205f
            || c == 0x3000;
}

/* Character.isDigit on each char: Integer.parseInt reads one char as a number
 * exactly when it is a decimal digit. */
static bool all_digits(const char *s, size_t len) {
    int32_t digit;
    for (size_t i = 0; i < len;) {
        size_t start = i;
        next_char(s, len, &i);
        if (!jnum_parse_int(s + start, i - start, &digit)) return false;
    }
    return true;
}

/* Integer.parseInt on digits, which throws for "" and on overflow. */
static int32_t parse_int(const char *s, size_t len) {
    int32_t v;
    if (!jnum_parse_int(s, len, &v)) {
        util_exit("java.lang.NumberFormatException: For input string: \"%.*s\"", (int)len, s);
    }
    return v;
}

/* msg, a newline and vcfRecord.substring(0, 80), which throws when the record
 * is shorter; or with clamp, substring(0, min(80, length)). Java prints the
 * high half of a surrogate pair cut at 80 as '?'. */
static _Noreturn void record_error(kstring_t *msg, const char *rec, size_t len, bool clamp) {
    kputc('\n', msg);
    int n = 0;
    for (size_t i = 0; i < len && n < 80;) {
        size_t start = i;
        next_char(rec, len, &i);
        int units = i - start == 4 ? 2 : 1;
        if (n + units > 80) {
            kputc('?', msg);
            n = 80;
        } else {
            kputsn(rec + start, i - start, msg);
            n += units;
        }
    }
    if (n < 80 && !clamp) util_exit("java.lang.StringIndexOutOfBoundsException: Range [0, 80) out of bounds for length %d", n);
    util_exit("%s", msg->s);
}

/* extractChrom */
static int extract_chrom(const char *rec, size_t len, size_t chrom_len) {
    kstring_t msg = {0, 0, NULL};
    if (chrom_len == 0 || (chrom_len == 1 && rec[0] == '.')) {
        kputs("ERROR: missing CHROM field: ", &msg);
        record_error(&msg, rec, len, false);
    }
    for (size_t i = 0; i < chrom_len;) {
        size_t start = i;
        int32_t c = next_char(rec, chrom_len, &i);
        if (c == ':' || is_java_whitespace(c)) {
            kputs("invalid character in CHROM field ['", &msg);
            kputsn(rec + start, i - start, &msg);
            kputs("']: ", &msg);
            record_error(&msg, rec, len, false);
        }
    }
    return chrom_ids_index(rec, chrom_len);
}

/* extractPos */
static int32_t extract_pos(const char *pos, size_t pos_len, const char *rec, size_t len) {
    if (!all_digits(pos, pos_len)) {
        kstring_t msg = {0, 0, NULL};
        ksprintf(&msg, "ERROR: invalid POS field: \"%.*s\"", (int)pos_len, pos);
        record_error(&msg, rec, len, true);
    }
    return parse_int(pos, pos_len);
}

/* extractIds: the ';'-separated IDs after removeMissingIds, whose result
 * extractIds discards. That method moves the IDs that are neither empty nor
 * "." to the front in place, so the stored list is those IDs followed by the
 * original entries from that count on: "rs1;.;rs2" becomes "rs1;rs2;rs2". */
static void put_ids(kstring_t *sb, const char *id, size_t len, const char *coord) {
    int n = 1;
    for (size_t j = 0; j < len; ++j) n += id[j] == ';';
    span *ids = util_malloc((size_t)n * sizeof *ids);
    const char *f = id, *end = id + len;
    int index = 0;
    for (int j = 0; j < n; ++j) {
        const char *semi = index_of(f, end, ';');
        const char *f_end = semi == NULL ? end : semi;
        span s = {f, (int)(f_end - f)};
        f = f_end + 1;
        for (size_t i = 0; i < (size_t)s.n;) {
            if (is_java_whitespace(next_char(s.s, (size_t)s.n, &i))) {
                util_exit("ERROR: ID field contains white-space at %s [%.*s]", coord, s.n, s.s);
            }
        }
        ids[j] = s;
        if (s.n > 0 && !(s.n == 1 && s.s[0] == '.')) ids[index++] = s;
    }
    for (int j = 0; j < n; ++j) {
        if (j > 0) kputc(';', sb);
        kputsn(ids[j].s, (size_t)ids[j].n, sb);
    }
    free(ids);
}

/* extractAlleles and checkAlleles: returns the allele count. */
static int check_alleles(const char *ref, const char *alt, const char *alt_end, const char *coord) {
    size_t ref_len = (size_t)(alt - 1 - ref), alt_len = (size_t)(alt_end - alt);
    if (ref_len == 0) util_exit("ERROR: missing REF field at %s", coord);
    if (alt_len == 0) util_exit("ERROR: missing ALT field: at %s", coord);
    int n = 1;
    if (!(alt_len == 1 && alt[0] == '.')) {
        n = 2;
        for (const char *c = alt; c < alt_end; ++c) n += *c == ',';
    }
    span *alleles = util_malloc((size_t)n * sizeof *alleles);
    alleles[0] = (span){ref, (int)ref_len};
    const char *f = alt;
    for (int j = 1; j < n; ++j) {
        const char *comma = index_of(f, alt_end, ',');
        const char *f_end = comma == NULL ? alt_end : comma;
        alleles[j] = (span){f, (int)(f_end - f)};
        f = f_end + 1;
    }
    for (int j = 1; j < n; ++j) {
        for (int k = 0; k < j; ++k) {
            if (alleles[j].n == alleles[k].n && memcmp(alleles[j].s, alleles[k].s, (size_t)alleles[j].n) == 0) {
                kstring_t list = {0, 0, NULL};
                for (int a = 0; a < n; ++a) {
                    kputs(a == 0 ? "[" : ", ", &list);
                    kputsn(alleles[a].s, (size_t)alleles[a].n, &list);
                }
                util_exit("ERROR: duplicate allele at %s %s]", coord, list.s);
            }
        }
    }
    free(alleles);
    /* checkREF: Character.toUpperCase maps no other char to one of these. */
    for (size_t j = 0; j < ref_len; ++j) {
        if (memchr("ACGTNacgtn", ref[j], 10) == NULL) {
            util_exit("ERROR: REF field is not a sequence of A, C, T, G, or N characters at %s [%.*s]", coord,
                    (int)ref_len, ref);
        }
    }
    return n;
}

/* extractEnd: the last END= subfield of INFO, or -1. */
static int32_t extract_end(const char *info, const char *info_end, int32_t pos, const char *coord) {
    int32_t end = -1;
    for (const char *f = info; f <= info_end;) {
        const char *semi = index_of(f, info_end, ';');
        const char *f_end = semi == NULL ? info_end : semi;
        if (f_end - f >= 4 && memcmp(f, "END=", 4) == 0) {
            const char *value = f + 4;
            size_t value_len = (size_t)(f_end - value);
            if (!all_digits(value, value_len)) {
                util_exit("ERROR: invalid INFO:END field at %s [END=%.*s]", coord, (int)value_len, value);
            }
            end = parse_int(value, value_len);
            if (end < pos) util_exit("ERROR: invalid INFO:END field at %s [%d]", coord, end);
        }
        f = f_end + 1;
    }
    return end;
}

void marker_parse(marker *m, const char *rec, size_t len) {
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    pthread_once(&once, init_snv_perms);
    /* StringUtil.getFields(vcfRecord, '\t', 9): INFO ends at the eighth tab,
     * or at the end of a record that has only seven. */
    const char *rec_end = rec + len;
    const char *tabs[8];
    const char *p = rec;
    int n_tabs = 0;
    while (n_tabs < 8 && (p = index_of(p, rec_end, '\t')) != NULL) tabs[n_tabs++] = p++;
    if (n_tabs < 7) util_exit("VCF record does not contain at least 8 tab-delimited fields: %.*s", (int)len, rec);

    m->chrom_index = extract_chrom(rec, len, (size_t)(tabs[0] - rec));
    m->pos = extract_pos(tabs[0] + 1, (size_t)(tabs[1] - tabs[0] - 1), rec, len);
    kstring_t coord = {0, 0, NULL};
    ksprintf(&coord, "%.*s:%d", (int)(tabs[0] - rec), rec, m->pos);

    kstring_t sb = {0, 0, NULL};
    uint16_t info = 0;
    const char *id = tabs[1] + 1;
    size_t id_len = (size_t)(tabs[2] - id);
    if (id_len == 0) util_exit("ERROR: missing ID field at %s", coord.s);
    if (!(id_len == 1 && id[0] == '.')) {
        put_ids(&sb, id, id_len, coord.s);
        info |= ID_STORED;
    }

    const char *alleles = tabs[2] + 1;
    size_t alleles_len = (size_t)(tabs[4] - alleles);
    m->n_alleles = check_alleles(alleles, tabs[3] + 1, tabs[4], coord.s);
    int snv = snv_index(alleles, alleles_len, m->n_alleles);
    if (snv >= 0) {
        info |= (uint16_t)snv;
    } else {
        if (sb.l > 0) kputc('\t', &sb);
        kputsn(alleles, alleles_len, &sb);
        info |= ALLELES_STORED;
    }
    m->end = extract_end(tabs[6] + 1, n_tabs == 8 ? tabs[7] : rec_end, m->pos, coord.s);
    free(coord.s);

    m->field_info = info;
    m->fields_len = sb.l;
    m->fields = sb.s;
}

void marker_free(marker *m) {
    free(m->fields);
    m->fields = NULL;
}

const char *marker_chrom(const marker *m) {
    return chrom_ids_id(m->chrom_index);
}

int marker_n_alleles(const marker *m) {
    return m->n_alleles;
}

bool marker_has_id(const marker *m) {
    return (m->field_info & ID_STORED) != 0;
}

/* The length of the stored ID list. */
static size_t id_len(const marker *m) {
    const char *tab = memchr(m->fields, '\t', m->fields_len);
    return tab == NULL ? m->fields_len : (size_t)(tab - m->fields);
}

span marker_id(const marker *m) {
    if (!marker_has_id(m)) return (span){".", 1};
    return (span){m->fields, (int)id_len(m)};
}

span marker_alleles(const marker *m) {
    if (m->field_info & ALLELES_STORED) {
        size_t start = marker_has_id(m) ? id_len(m) + 1 : 0;
        return (span){m->fields + start, (int)(m->fields_len - start)};
    }
    return snv_alleles(m->field_info & SNV_INDEX_MASK, m->n_alleles);
}

int marker_bits_per_allele(const marker *m) {
    int n_alleles = marker_n_alleles(m);
    int bits = 0;
    while ((1 << bits) < n_alleles) ++bits;
    return bits;
}
