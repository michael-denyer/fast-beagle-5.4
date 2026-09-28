/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) bref/SeqCoder3.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "bref/seq_coder3.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/int_list.h"
#include "blbutil/utilities.h"
#include "jcompat/jmath.h"
#include "jcompat/jnum.h"

#define CHARACTER_MAX_VALUE 65535

struct seq_coder3 {
    int n_haps;
    int max_nseq;
    ref_gt_rec **recs;
    int n_recs, cap_recs;
    int *hap_to_seq;
    int_list seq_to_cnt;
    int_list *seq_map;   /* per sequence: (allele, new sequence) pairs */
    int n_seq_map, cap_seq_map;
};

int seq_coder3_default_max_nseq(int n_samples) {
    if (n_samples < 1) util_exit("java.lang.IllegalArgumentException: %d", n_samples);
    if (n_samples == 1) return 3;
    double exponent = 2 * jmath_log10(n_samples) + 1;
    int64_t max_nseq = jnum_d2l(floor(jmath_pow(2.0, exponent)));
    return max_nseq > CHARACTER_MAX_VALUE ? CHARACTER_MAX_VALUE : (int)max_nseq;
}

static void seq_map_push(seq_coder3 *c) {
    if (c->n_seq_map == c->cap_seq_map) {
        c->cap_seq_map = c->cap_seq_map == 0 ? 16 : 2 * c->cap_seq_map;
        c->seq_map = util_realloc(c->seq_map, (size_t)c->cap_seq_map * sizeof *c->seq_map);
        memset(c->seq_map + c->n_seq_map, 0, (size_t)(c->cap_seq_map - c->n_seq_map) * sizeof *c->seq_map);
    }
    c->seq_map[c->n_seq_map++].n = 0;
}

/* SeqCoder3.initialize: every haplotype carries the empty sequence 0. */
static void initialize(seq_coder3 *c) {
    c->n_recs = 0;
    memset(c->hap_to_seq, 0, (size_t)c->n_haps * sizeof *c->hap_to_seq);
    c->seq_to_cnt.n = 0;
    int_list_add(&c->seq_to_cnt, c->n_haps);
    c->n_seq_map = 0;
    seq_map_push(c);
}

seq_coder3 *seq_coder3_new(int n_samples, int max_nseq) {
    if (max_nseq < 0 || max_nseq >= CHARACTER_MAX_VALUE) util_exit("java.lang.IllegalArgumentException: %d", max_nseq);
    seq_coder3 *c = util_malloc(sizeof *c);
    *c = (seq_coder3){0};
    c->n_haps = 2 * n_samples;
    c->max_nseq = max_nseq;
    c->hap_to_seq = util_malloc((size_t)c->n_haps * sizeof *c->hap_to_seq);
    initialize(c);
    return c;
}

int seq_coder3_max_nseq(const seq_coder3 *c) {
    return c->max_nseq;
}

int seq_coder3_n_recs(const seq_coder3 *c) {
    return c->n_recs;
}

/* The index in a sequence's pair list of the pair for allele a, or the list length. */
static int find_allele(const int_list *list, int a) {
    int index = 0;
    while (index < list->n && list->v[index] != a) index += 2;
    return index;
}

/* SeqCoder3.addMajorAllele: a sequence some of whose haplotypes keep the major
 * allele gives its index to the major allele, and the allele that first
 * claimed it moves to a new sequence. */
static void add_major_allele(seq_coder3 *c, const int *seq_non_major_cnt, int n_start_seq, int major) {
    for (int seq = 0; seq < n_start_seq; ++seq) {
        if (seq_non_major_cnt[seq] < c->seq_to_cnt.v[seq]) {
            int_list *list = &c->seq_map[seq];
            if (list->n == 0) {
                int_list_add(list, major);
                int_list_add(list, seq);
            } else {
                int_list_add(list, list->v[0]);
                int_list_add(list, c->n_seq_map);
                list->v[0] = major;
                seq_map_push(c);
            }
        }
    }
}

/* SeqCoder3.setAlleleMap */
static bool set_allele_map(seq_coder3 *c, const ref_gt_rec *rec) {
    int n_start_seq = c->seq_to_cnt.n;
    int *seq_non_major_cnt = util_malloc((size_t)(n_start_seq > 0 ? n_start_seq : 1) * sizeof(int));
    memset(seq_non_major_cnt, 0, (size_t)n_start_seq * sizeof(int));
    for (int j = 0; j < c->n_seq_map; ++j) c->seq_map[j].n = 0;
    int n_alleles = marker_n_alleles(&rec->marker);
    for (int a = 0; a < n_alleles; ++a) {
        if (a == rec->major_allele) continue;
        for (int k = 0; k < rec->hap_list_len[a]; ++k) {
            int seq = c->hap_to_seq[rec->hap_lists[a][k]];
            ++seq_non_major_cnt[seq];
            if (c->seq_map[seq].n == 0) {
                int_list_add(&c->seq_map[seq], a);
                int_list_add(&c->seq_map[seq], seq);
            } else if (find_allele(&c->seq_map[seq], a) == c->seq_map[seq].n) {
                int_list_add(&c->seq_map[seq], a);
                int_list_add(&c->seq_map[seq], c->n_seq_map);
                seq_map_push(c);
            }
        }
    }
    add_major_allele(c, seq_non_major_cnt, n_start_seq, rec->major_allele);
    free(seq_non_major_cnt);
    if (c->n_seq_map >= c->max_nseq) {
        c->n_seq_map = n_start_seq;
        return false;
    }
    return true;
}

bool seq_coder3_add(seq_coder3 *c, ref_gt_rec *rec) {
    if (rec->kind == REF_HAP) util_exit("java.lang.IllegalArgumentException: class vcf.SeqCodedRefGTRec");
    if (!set_allele_map(c, rec)) return false;
    if (c->n_recs == c->cap_recs) {
        c->cap_recs = c->cap_recs == 0 ? 64 : 2 * c->cap_recs;
        c->recs = util_realloc(c->recs, (size_t)c->cap_recs * sizeof *c->recs);
    }
    c->recs[c->n_recs++] = rec;
    int n_alleles = marker_n_alleles(&rec->marker);
    for (int a = 0; a < n_alleles; ++a) {
        if (a == rec->major_allele) continue;
        for (int k = 0; k < rec->hap_list_len[a]; ++k) {
            int h = rec->hap_lists[a][k];
            int old_seq = c->hap_to_seq[h];
            const int_list *list = &c->seq_map[old_seq];
            int new_seq = list->v[find_allele(list, a) + 1];
            if (new_seq != old_seq) {
                while (new_seq >= c->seq_to_cnt.n) int_list_add(&c->seq_to_cnt, 0);
                c->hap_to_seq[h] = new_seq;
                --c->seq_to_cnt.v[old_seq];
                ++c->seq_to_cnt.v[new_seq];
            }
        }
    }
    return true;
}

int seq_coder3_flush(seq_coder3 *c, ref_gt_rec **out) {
    int n = c->n_recs;
    if (n == 0) return 0;
    int n_seq = c->n_seq_map;
    seq_group *group = seq_group_new(c->n_haps, n_seq);
    memcpy(group->hap_to_seq, c->hap_to_seq, (size_t)c->n_haps * sizeof *group->hap_to_seq);
    for (int j = 0; j < n; ++j) {
        ref_gt_rec *rec = c->recs[j];
        /* Every haplotype in a sequence carries the same allele, so the
         * non-major haplotype lists give each sequence's allele directly. */
        uint8_t *seq_to_allele = util_malloc((size_t)n_seq);
        memset(seq_to_allele, rec->major_allele, (size_t)n_seq);
        int n_alleles = marker_n_alleles(&rec->marker);
        for (int a = 0; a < n_alleles; ++a) {
            if (a == rec->major_allele) continue;
            for (int k = 0; k < rec->hap_list_len[a]; ++k) seq_to_allele[c->hap_to_seq[rec->hap_lists[a][k]]] = (uint8_t)a;
        }
        ref_gt_rec_set_seq_coded(rec, group, seq_to_allele);
        out[j] = rec;
    }
    seq_group_release(group);
    initialize(c);
    return n;
}

void seq_coder3_free(seq_coder3 *c) {
    for (int j = 0; j < c->n_recs; ++j) ref_gt_rec_release(c->recs[j]);
    free(c->recs);
    free(c->hap_to_seq);
    free(c->seq_to_cnt.v);
    for (int j = 0; j < c->cap_seq_map; ++j) free(c->seq_map[j].v);
    free(c->seq_map);
    free(c);
}
