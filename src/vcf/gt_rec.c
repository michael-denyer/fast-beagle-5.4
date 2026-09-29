/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/VcfRecGTParser.java (hapListRep),
 * vcf/VcfIt.java (TO_LOWMEM_GT_REC), vcf/LowMafDiallelicGTRec.java,
 * vcf/LowMafGTRec.java and vcf/BitArrayGTRec.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/gt_rec.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/int_list.h"
#include "blbutil/utilities.h"
#include "jcompat/jarrays.h"
#include "vcf/gt_parse.h"

void gt_rec_parse(gt_rec *rec, const char *line, size_t len, const vcf_header *h) {
    memset(rec, 0, sizeof *rec);
    marker_parse(&rec->marker, line, len);
    const samples *smp = &h->samples;
    int n_alleles = marker_n_alleles(&rec->marker);
    gt_fields f;
    gt_fields_init(&f, line, len, h, &rec->marker, false);

    /* VcfRecGTParser.hapListRep */
    int_list *lists = util_malloc((size_t)n_alleles * sizeof *lists);
    memset(lists, 0, (size_t)n_alleles * sizeof *lists);
    int_list missing = {0};
    bool is_phased = true;
    for (int s = 0; s < smp->n; ++s) {
        gt_field gt;
        gt_fields_next(&f, s, &gt);
        is_phased &= !gt.diploid || line[gt.end1] == '|';
        if (gt.diploid != smp->is_diploid[s]) gt_ploidy_error(&f, s, gt.diploid);
        int a1 = gt_parse_allele(&f, gt.start, gt.end1);
        int a2 = gt.diploid ? gt_parse_allele(&f, gt.end1 + 1, gt.end2) : a1;
        if (a1 < 0 || a2 < 0) {
            int_list_add(&missing, s);
            is_phased = false;
        } else {
            int_list_add(&lists[a1], s << 1);
            int_list_add(&lists[a2], (s << 1) | 1);
        }
    }

    /* HapListRep.majorAllele: the most frequent allele, the lowest index on ties. */
    int major = 0;
    for (int a = 1; a < n_alleles; ++a) {
        if (lists[a].n > lists[major].n) major = a;
    }
    int non_major = 0;
    for (int a = 0; a < n_alleles; ++a) {
        if (a != major) non_major += lists[a].n;
    }

    rec->is_phased = is_phased;
    rec->n_haps = smp->n << 1;
    rec->major_allele = major;
    rec->n_missing = missing.n;
    rec->missing = missing.v;
    if (non_major <= (smp->n >> 7)) {
        rec->kind = n_alleles == 2 ? GT_LOW_MAF_DIALLELIC : GT_LOW_MAF;
        rec->hap_list_len = util_malloc((size_t)n_alleles * sizeof *rec->hap_list_len);
        rec->hap_lists = util_malloc((size_t)n_alleles * sizeof *rec->hap_lists);
        for (int a = 0; a < n_alleles; ++a) {
            bool keep = a != major;
            rec->hap_list_len[a] = keep ? lists[a].n : 0;
            rec->hap_lists[a] = keep ? lists[a].v : NULL;
            if (!keep) free(lists[a].v);
        }
    } else {
        rec->kind = GT_BIT_ARRAY;
        rec->bits_per_allele = marker_bits_per_allele(&rec->marker);
        size_t n_bits = (size_t)rec->n_haps * (size_t)rec->bits_per_allele;
        rec->allele_bits = util_malloc(((n_bits + 63) / 64) * sizeof(uint64_t));
        memset(rec->allele_bits, 0, ((n_bits + 63) / 64) * sizeof(uint64_t));
        rec->missing_bits = util_malloc(((size_t)smp->n + 63) / 64 * sizeof(uint64_t));
        memset(rec->missing_bits, 0, ((size_t)smp->n + 63) / 64 * sizeof(uint64_t));
        for (int a = 0; a < n_alleles; ++a) {
            for (int j = 0; j < lists[a].n; ++j) {
                size_t index = (size_t)lists[a].v[j] * (size_t)rec->bits_per_allele;
                for (int k = 0; k < rec->bits_per_allele; ++k, ++index) {
                    if ((a >> k) & 1) rec->allele_bits[index >> 6] |= 1ULL << (index & 63);
                }
            }
            free(lists[a].v);
        }
        for (int j = 0; j < missing.n; ++j) {
            rec->missing_bits[missing.v[j] >> 6] |= 1ULL << (missing.v[j] & 63);
        }
    }
    free(lists);
}

void gt_rec_free(gt_rec *rec) {
    marker_free(&rec->marker);
    free(rec->missing);
    if (rec->hap_lists != NULL) {
        for (int a = 0; a < marker_n_alleles(&rec->marker); ++a) free(rec->hap_lists[a]);
    }
    free(rec->hap_lists);
    free(rec->hap_list_len);
    free(rec->allele_bits);
    free(rec->missing_bits);
    memset(rec, 0, sizeof *rec);
}

gt_rec *gt_rec_retain(gt_rec *rec) {
    ++rec->refs;
    return rec;
}

void gt_rec_release(gt_rec *rec) {
    if (rec != NULL && --rec->refs == 0) {
        gt_rec_free(rec);
        free(rec);
    }
}

int gt_rec_get(const gt_rec *rec, int hap) {
    if (hap < 0 || hap >= rec->n_haps) util_exit("java.lang.IndexOutOfBoundsException: %d", hap);
    if (rec->kind == GT_BIT_ARRAY) {
        int s = hap >> 1;
        if ((rec->missing_bits[s >> 6] >> (s & 63)) & 1) return -1;
        size_t index = (size_t)hap * (size_t)rec->bits_per_allele;
        int allele = 0;
        for (int k = 0; k < rec->bits_per_allele; ++k, ++index) {
            if ((rec->allele_bits[index >> 6] >> (index & 63)) & 1) allele |= 1 << k;
        }
        return allele;
    }
    for (int a = 0; a < marker_n_alleles(&rec->marker); ++a) {
        if (a != rec->major_allele && jarrays_search_int(rec->hap_lists[a], 0, rec->hap_list_len[a], hap) >= 0) return a;
    }
    return jarrays_search_int(rec->missing, 0, rec->n_missing, hap >> 1) >= 0 ? -1 : rec->major_allele;
}

const char *gt_rec_class_name(const gt_rec *rec) {
    switch (rec->kind) {
        case GT_LOW_MAF_DIALLELIC: return "LowMafDiallelicGTRec";
        case GT_LOW_MAF: return "LowMafGTRec";
        default: return "BitArrayGTRec";
    }
}
