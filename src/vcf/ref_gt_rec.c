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
#include "vcf/ref_gt_rec.h"

#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jarrays.h"
#include "vcf/gt_parse.h"
#include "vcf/vcf_it.h"

void ref_gt_rec_parse(ref_gt_rec *rec, const char *line, size_t len, const vcf_header *h) {
    memset(rec, 0, sizeof *rec);
    marker_parse(&rec->marker, line, len);
    const samples *smp = &h->samples;
    int n_alleles = marker_n_alleles(&rec->marker);
    gt_fields f;
    gt_fields_init(&f, line, len, h, &rec->marker, true);

    /* VcfRecGTParser.phasedAlleles: here the alleles are parsed before the ploidy check. */
    int n_haps = smp->n << 1;
    int *alleles = util_malloc((size_t)n_haps * sizeof *alleles);
    for (int s = 0; s < smp->n; ++s) {
        /* A phased diploid "d|d" column with no skipped column before it
         * parses directly; anything else, errors included, takes the general
         * path below. */
        size_t p = (size_t)f.tab + 1;
        if (f.tab != -1 && h->included[s] == f.unfilt + 1 && smp->is_diploid[s] && p + 3 <= len
                && (p + 3 == len || line[p + 3] == '\t') && line[p + 1] == '|'
                && line[p] >= '0' && line[p] <= '9' && line[p + 2] >= '0' && line[p + 2] <= '9'
                && line[p] - '0' < n_alleles && line[p + 2] - '0' < n_alleles) {
            alleles[s << 1] = line[p] - '0';
            alleles[(s << 1) | 1] = line[p + 2] - '0';
            f.unfilt = h->included[s];
            f.tab = p + 3 == len ? -1 : (long)(p + 3);
            continue;
        }
        gt_field gt;
        gt_fields_next(&f, s, &gt);
        int a1 = gt_parse_allele(&f, gt.start, gt.end1);
        int a2 = gt.diploid ? gt_parse_allele(&f, gt.end1 + 1, gt.end2) : a1;
        if (gt.diploid != smp->is_diploid[s]) gt_ploidy_error(&f, s, gt.diploid);
        if ((gt.diploid && line[gt.end1] != '|') || a1 == -1 || a2 == -1) {
            gt_sample_error(&f, "unphased or missing genotype", s);
        }
        alleles[s << 1] = a1;
        alleles[(s << 1) | 1] = a2;
    }

    /* VcfRecGTParser.nonMajRefIndices */
    int *counts = util_malloc((size_t)n_alleles * sizeof *counts);
    memset(counts, 0, (size_t)n_alleles * sizeof *counts);
    for (int j = 0; j < n_haps; ++j) ++counts[alleles[j]];
    int major = 0;
    for (int a = 1; a < n_alleles; ++a) {
        if (counts[a] > counts[major]) major = a;
    }
    rec->hap_list_len = util_malloc((size_t)n_alleles * sizeof *rec->hap_list_len);
    rec->hap_lists = util_malloc((size_t)n_alleles * sizeof *rec->hap_lists);
    for (int a = 0; a < n_alleles; ++a) {
        rec->hap_list_len[a] = a == major ? 0 : counts[a];
        rec->hap_lists[a] = a == major ? NULL : util_malloc((size_t)counts[a] * sizeof(int));
    }
    memset(counts, 0, (size_t)n_alleles * sizeof *counts);
    for (int j = 0; j < n_haps; ++j) {
        int a = alleles[j];
        if (a != major) rec->hap_lists[a][counts[a]++] = j;
    }
    free(counts);
    free(alleles);

    rec->kind = n_alleles == 2 ? REF_TWO_ALLELE : REF_ALLELE;
    rec->n_haps = n_haps;
    rec->major_allele = major;
}

seq_group *seq_group_new(int n_haps, int n_seq) {
    seq_group *g = util_malloc(sizeof *g);
    g->refs = 1;
    g->n_haps = n_haps;
    g->n_seq = n_seq;
    g->trace_id = -1;
    g->hap_to_seq = util_malloc((size_t)n_haps * sizeof *g->hap_to_seq);
    return g;
}

seq_group *seq_group_retain(seq_group *g) {
    ++g->refs;
    return g;
}

void seq_group_release(seq_group *g) {
    if (g != NULL && --g->refs == 0) {
        free(g->hap_to_seq);
        free(g);
    }
}

static void free_hap_lists(ref_gt_rec *rec) {
    if (rec->hap_lists != NULL) {
        for (int a = 0; a < marker_n_alleles(&rec->marker); ++a) free(rec->hap_lists[a]);
    }
    free(rec->hap_lists);
    free(rec->hap_list_len);
    rec->hap_lists = NULL;
    rec->hap_list_len = NULL;
}

void ref_gt_rec_set_seq_coded(ref_gt_rec *rec, seq_group *g, uint8_t *seq_to_allele) {
    free_hap_lists(rec);
    rec->kind = REF_HAP;
    rec->group = seq_group_retain(g);
    rec->seq_to_allele = seq_to_allele;
}

static void ref_gt_rec_free(ref_gt_rec *rec) {
    marker_free(&rec->marker);
    free_hap_lists(rec);
    seq_group_release(rec->group);
    free(rec->seq_to_allele);
    memset(rec, 0, sizeof *rec);
}

ref_gt_rec *ref_gt_rec_retain(ref_gt_rec *rec) {
    ++rec->refs;
    return rec;
}

void ref_gt_rec_release(ref_gt_rec *rec) {
    if (rec != NULL && --rec->refs == 0) {
        ref_gt_rec_free(rec);
        free(rec);
    }
}

int ref_gt_rec_get(const ref_gt_rec *rec, int hap) {
    if (hap < 0 || hap >= rec->n_haps) util_exit("java.lang.IndexOutOfBoundsException: %d", hap);
    if (rec->kind == REF_HAP) return rec->seq_to_allele[rec->group->hap_to_seq[hap]];
    for (int a = 0; a < marker_n_alleles(&rec->marker); ++a) {
        if (a != rec->major_allele && jarrays_search_int(rec->hap_lists[a], 0, rec->hap_list_len[a], hap) >= 0) return a;
    }
    return rec->major_allele;
}

int ref_gt_rec_major_allele(const ref_gt_rec *rec) {
    if (rec->kind != REF_HAP) return rec->major_allele;
    int n_alleles = marker_n_alleles(&rec->marker);
    int *counts = util_malloc((size_t)n_alleles * sizeof *counts);
    memset(counts, 0, (size_t)n_alleles * sizeof *counts);
    for (int h = 0; h < rec->n_haps; ++h) ++counts[ref_gt_rec_get(rec, h)];
    int major = 0;
    for (int a = 1; a < n_alleles; ++a) {
        if (counts[a] > counts[major]) major = a;
    }
    free(counts);
    return major;
}

const char *ref_gt_rec_class_name(const ref_gt_rec *rec) {
    switch (rec->kind) {
        case REF_TWO_ALLELE: return "LowMafRefDiallelicGTRec";
        case REF_ALLELE: return "LowMafRefGTRec";
        default: return "SeqCodedRefGTRec";
    }
}

void ref_gt_rec_trace(const ref_gt_rec *rec) {
    vcf_it_trace_marker("T1a-ref", &rec->marker);
    kstring_t s = {0, 0, NULL};
    if (rec->kind == REF_HAP) {
        seq_group *g = rec->group;
        if (g->trace_id < 0) {
            static int n_groups;
            g->trace_id = n_groups++;
            ksprintf(&s, "group\t%d\t", g->trace_id);
            for (int h = 0; h < g->n_haps; ++h) {
                if (h > 0) kputc(',', &s);
                kputw(g->hap_to_seq[h], &s);
            }
            trace_line("T1d-ref", "%s", s.s);
            s.l = 0;
        }
    }
    ksprintf(&s, "rec\t%s\t%d\t", ref_gt_rec_class_name(rec), ref_gt_rec_major_allele(rec));
    if (rec->kind == REF_HAP) {
        kputw(rec->group->trace_id, &s);
        kputc('\t', &s);
        for (int j = 0; j < rec->group->n_seq; ++j) {
            if (j > 0) kputc(',', &s);
            kputw(rec->seq_to_allele[j], &s);
        }
    } else {
        kputs("-\t-", &s);
    }
    kputc('\t', &s);
    for (int h = 0; h < rec->n_haps; ++h) {
        if (h > 0) kputc(',', &s);
        kputw(ref_gt_rec_get(rec, h), &s);
    }
    trace_line("T1d-ref", "%s", s.s);
    free(s.s);
}

const marker *ref_gt_rec_it_marker(const void *rec) {
    return &((const ref_gt_rec *)rec)->marker;
}

void ref_gt_rec_it_release(void *rec) {
    ref_gt_rec_release(rec);
}
