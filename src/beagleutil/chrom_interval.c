/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/ChromInterval.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "beagleutil/chrom_interval.h"

#include <string.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/utilities.h"

/* Java's check for a leading zero compares the char with 0, not '0', so it
 * never rejects anything; Character.isDigit also accepts non-ASCII digits,
 * which this port does not. */
static bool is_valid_pos(const char *s, size_t start, size_t end) {
    for (size_t j = start; j < end; ++j) {
        if (s[j] < '0' || s[j] > '9') return false;
    }
    return true;
}

static int32_t parse_pos(const char *s, size_t start, size_t end) {
    int64_t v = 0;
    for (size_t j = start; j < end; ++j) {
        v = 10 * v + (s[j] - '0');
        if (v > INT32_MAX) util_exit("NumberFormatException: For input string: \"%.*s\"", (int)(end - start), s + start);
    }
    return (int32_t)v;
}

static const char *last_index_of(const char *s, size_t len, char c) {
    for (size_t j = len; j > 0; --j) {
        if (s[j - 1] == c) return s + j - 1;
    }
    return NULL;
}

bool chrom_interval_parse(const char *str, chrom_interval *out) {
    while (*str != '\0' && (unsigned char)*str <= ' ') ++str;
    size_t length = strlen(str);
    while (length > 0 && (unsigned char)str[length - 1] <= ' ') --length;
    if (length == 0) return false;
    int32_t start = INT32_MIN, end = INT32_MAX;
    const char *colon = last_index_of(str, length, ':');
    const char *hyphen = last_index_of(str, length, '-');
    long chr_delim = colon == NULL ? -1 : colon - str;
    long pos_delim = hyphen == NULL ? -1 : hyphen - str;
    size_t chrom_len;
    if (chr_delim == -1) {
        chrom_len = length;
    } else if (chr_delim == (long)length - 1) {
        chrom_len = length - 1;
    } else {
        if (pos_delim == -1 || pos_delim <= chr_delim || chr_delim == (long)length - 2
                || !is_valid_pos(str, (size_t)chr_delim + 1, (size_t)pos_delim)
                || !is_valid_pos(str, (size_t)pos_delim + 1, length)) {
            return false;
        }
        if (pos_delim > chr_delim + 1) start = parse_pos(str, (size_t)chr_delim + 1, (size_t)pos_delim);
        if ((long)length > pos_delim + 1) end = parse_pos(str, (size_t)pos_delim + 1, length);
        if (start > end) return false;
        chrom_len = (size_t)chr_delim;
    }
    out->chrom_index = chrom_ids_index(str, chrom_len);
    out->start = start;
    out->end = end;
    return true;
}

bool chrom_interval_contains(const chrom_interval *ci, int chrom_index, int32_t pos) {
    return chrom_index == ci->chrom_index && ci->start <= pos && pos <= ci->end;
}
