/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/IntervalVcfIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/interval_it.h"

#include <stdio.h>
#include <stdlib.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/utilities.h"

typedef struct {
    sample_file_it_ops ops;  /* ours, with the record accessors of src */
    sample_file_it src;
    const chrom_interval *interval;
    void *next;              /* the lookahead */
} interval_it;

static _Noreturn void empty_interval(const chrom_interval *interval) {
    /* ChromInterval.toString */
    char start[16] = "", end[16] = "";
    if (interval->start != INT32_MIN) snprintf(start, sizeof start, "%d", interval->start);
    if (interval->end != INT32_MAX) snprintf(end, sizeof end, "%d", interval->end);
    bool bounded = interval->start != INT32_MIN || interval->end != INT32_MAX;
    util_exit("No VCF records found in the specified interval.\n"
            "Check chromosome identifier and interval: %s%s%s%s%s",
            chrom_ids_id(interval->chrom_index), bounded ? ":" : "", start, bounded ? "-" : "", end);
}

static bool contains(const interval_it *it, const void *rec) {
    const marker *m = it->src.ops->marker(rec);
    return chrom_interval_contains(it->interval, m->chrom_index, m->pos);
}

static const samples *interval_samples(const void *self) {
    const interval_it *it = self;
    return sample_file_it_samples(it->src);
}

/* IntervalVcfIt.next and readNextRecord */
static void *interval_next(void *self) {
    interval_it *it = self;
    void *cur = it->next;
    if (cur == NULL) return NULL;
    it->next = sample_file_it_next(it->src);
    if (it->next != NULL && !contains(it, it->next)) {
        it->src.ops->release(it->next);
        it->next = NULL;
    }
    return cur;
}

static void interval_close(void *self) {
    interval_it *it = self;
    if (it->next != NULL) it->src.ops->release(it->next);
    sample_file_it_close(it->src);
    free(it);
}

sample_file_it interval_it_open(sample_file_it src, const chrom_interval *interval) {
    if (interval == NULL) return src;
    interval_it *it = util_malloc(sizeof *it);
    *it = (interval_it){.ops = {.samples = interval_samples, .next = interval_next, .close = interval_close,
            .marker = src.ops->marker, .release = src.ops->release}, .src = src, .interval = interval};
    /* IntervalVcfIt.readFirstRecord */
    while (it->next == NULL) {
        void *rec = sample_file_it_next(src);
        if (rec == NULL) break;
        if (contains(it, rec)) it->next = rec;
        else src.ops->release(rec);
    }
    if (it->next == NULL) empty_interval(interval);
    return (sample_file_it){&it->ops, it};
}
