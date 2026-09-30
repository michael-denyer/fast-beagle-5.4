/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/Markers.java and
 * vcf/BasicMarker.java (equals, hashCode, toString); modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/markers.h"

#include <stdlib.h>
#include <string.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/utilities.h"

static bool span_equals(span a, span b) {
    return a.n == b.n && memcmp(a.s, b.s, (size_t)a.n) == 0;
}

bool marker_equals(const marker *a, const marker *b) {
    return a->chrom_index == b->chrom_index && a->pos == b->pos
            && span_equals(marker_alleles(a), marker_alleles(b)) && a->end == b->end;
}

static _Noreturn void marker_order_error(const marker *const *m) {
    span id[3], al[3];
    for (int j = 0; j < 3; ++j) {
        id[j] = marker_id(m[j]);
        al[j] = marker_alleles(m[j]);
    }
    util_exit("java.lang.IllegalArgumentException: markers not in chromosomal order: "
            "\n%s\t%d\t%.*s\t%.*s"
            "\n%s\t%d\t%.*s\t%.*s"
            "\n%s\t%d\t%.*s\t%.*s",
            marker_chrom(m[0]), m[0]->pos, id[0].n, id[0].s, al[0].n, al[0].s,
            marker_chrom(m[1]), m[1]->pos, id[1].n, id[1].s, al[1].n, al[1].s,
            marker_chrom(m[2]), m[2]->pos, id[2].n, id[2].s, al[2].n, al[2].s);
}

/* A hash consistent with marker_equals, for the duplicate check. */
static uint32_t marker_hash(const marker *m) {
    uint32_t h = 2166136261u;
    span al = marker_alleles(m);
    h = (h ^ (uint32_t)m->chrom_index) * 16777619u;
    h = (h ^ (uint32_t)m->pos) * 16777619u;
    for (int j = 0; j < al.n; ++j) h = (h ^ (unsigned char)al.s[j]) * 16777619u;
    h = (h ^ (uint32_t)m->end) * 16777619u;
    return h;
}

void markers_check(const marker *const *markers, int n) {
    /* Markers.checkMarkerPosOrder */
    if (n >= 2) {
        int *seen = util_malloc((size_t)n * sizeof *seen);
        int n_seen = 0;
        seen[n_seen++] = markers[0]->chrom_index;
        if (markers[1]->chrom_index != markers[0]->chrom_index) seen[n_seen++] = markers[1]->chrom_index;
        for (int j = 2; j < n; ++j) {
            int chr0 = markers[j - 2]->chrom_index, chr1 = markers[j - 1]->chrom_index, chr2 = markers[j]->chrom_index;
            if (chr0 == chr1 && chr1 == chr2) {
                int32_t pos0 = markers[j - 2]->pos, pos1 = markers[j - 1]->pos, pos2 = markers[j]->pos;
                if ((pos1 < pos0 && pos1 < pos2) || (pos1 > pos0 && pos1 > pos2)) {
                    free(seen);
                    marker_order_error(markers + j - 2);
                }
            } else if (chr1 != chr2) {
                for (int k = 0; k < n_seen; ++k) {
                    if (seen[k] == chr2) util_exit("markers on chromosome are not contiguous: %s", chrom_ids_id(chr2));
                }
                seen[n_seen++] = chr2;
            }
        }
        free(seen);
    }
    /* Markers.markerSet: open addressing over marker_equals. */
    size_t cap = 16;
    while (cap < 2 * (size_t)n) cap <<= 1;
    const marker **table = util_malloc(cap * sizeof *table);
    memset(table, 0, cap * sizeof *table);
    for (int j = 0; j < n; ++j) {
        size_t slot = marker_hash(markers[j]) & (cap - 1);
        while (table[slot] != NULL) {
            if (marker_equals(table[slot], markers[j])) {
                const marker *m = markers[j];
                span id = marker_id(m), al = marker_alleles(m);
                free(table);
                util_exit("java.lang.IllegalArgumentException: Duplicate marker: %s\t%d\t%.*s\t%.*s",
                        marker_chrom(m), m->pos, id.n, id.s, al.n, al.s);
            }
            slot = (slot + 1) & (cap - 1);
        }
        table[slot] = markers[j];
    }
    free(table);
}
