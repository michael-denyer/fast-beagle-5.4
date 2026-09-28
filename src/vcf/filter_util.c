/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/FilterUtil.java and
 * blbutil/Utilities.java (idSet); modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/filter_util.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/line_reader.h"
#include "blbutil/utilities.h"

str_set *filter_id_set(const char *path) {
    if (path == NULL) return NULL;
    str_set *set = str_set_new();
    line_reader *r = line_reader_open(path, 1);
    kstring_t line = {0, 0, NULL};
    while (line_reader_next(r, &line)) {
        /* String.trim, then StringUtil.countFields on white space (chars <= ' '). */
        const char *s = line.s;
        size_t start = 0, end = line.l;
        while (start < end && (unsigned char)s[start] <= ' ') ++start;
        while (end > start && (unsigned char)s[end - 1] <= ' ') --end;
        if (start == end) continue;
        for (size_t j = start; j < end; ++j) {
            if ((unsigned char)s[j] <= ' ') {
                util_exit("java.lang.IllegalArgumentException: line has >1 white-space delimited fields: %.*s",
                        (int)(end - start), s + start);
            }
        }
        str_set_index(set, s + start, end - start);
    }
    free(line.s);
    line_reader_close(r);
    return set;
}

bool filter_accept_marker(const str_set *exclude, const marker *m) {
    if (exclude == NULL || str_set_size(exclude) == 0) return true;
    if (marker_has_id(m)) {
        /* FilterUtil.markerIsInSet: each entry of the stored ID list. */
        span id = marker_id(m);
        const char *f = id.s, *end = id.s + id.n;
        for (;;) {
            const char *semi = memchr(f, ';', (size_t)(end - f));
            const char *f_end = semi == NULL ? end : semi;
            if (str_set_find(exclude, f, (size_t)(f_end - f)) >= 0) return false;
            if (semi == NULL) break;
            f = semi + 1;
        }
    }
    kstring_t pos_id = {0, 0, NULL};
    ksprintf(&pos_id, "%s:%d", marker_chrom(m), m->pos);
    bool accept = str_set_find(exclude, pos_id.s, pos_id.l) < 0;
    free(pos_id.s);
    return accept;
}
