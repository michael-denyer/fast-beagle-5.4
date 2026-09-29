/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/GeneticMap.java,
 * vcf/PositionMap.java and vcf/PlinkGenMap.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/genetic_map.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/line_reader.h"
#include "blbutil/utilities.h"
#include "jcompat/jarrays.h"
#include "jcompat/jnum.h"
#include "vcf/marker.h"

#define MIN_END_CM_DIST 5.0

typedef struct {
    int n;
    int32_t *base_pos;
    double *gen_pos;
    char **lines;   /* kept until parsed, for error messages */
    int cap;
} chrom_map;

struct genetic_map {
    bool plink;
    double scale_factor, inv_scale_factor;  /* PositionMap */
    int n_chrom;
    chrom_map *chroms;                      /* PlinkGenMap, indexed by chromosome */
};

/* Arrays.binarySearch(double[]): equal values are ordered by their bits. */
static int search_double(const double *a, int n, double key) {
    int low = 0, high = n - 1;
    while (low <= high) {
        int mid = (int)((unsigned)(low + high) >> 1);
        double v = a[mid];
        if (v < key) {
            low = mid + 1;
        } else if (v > key) {
            high = mid - 1;
        } else {
            /* Double.doubleToLongBits, with NaN canonicalised. */
            int64_t vb, kb;
            double vv = isnan(v) ? NAN : v, kk = isnan(key) ? NAN : key;
            memcpy(&vb, &vv, sizeof vb);
            memcpy(&kb, &kk, sizeof kb);
            if (vb == kb) return mid;
            if (vb < kb) low = mid + 1;
            else high = mid - 1;
        }
    }
    return -(low + 1);
}

/* StringUtil.getFields(s): the fields of the trimmed line separated by runs of
 * characters <= ' '. Returns the count and fills up to max spans. */
static int whitespace_fields(const char *s, size_t len, span *out, int max) {
    size_t j = 0;
    int n = 0;
    while (j < len) {
        while (j < len && (unsigned char)s[j] <= ' ') ++j;
        if (j == len) break;
        size_t start = j;
        while (j < len && (unsigned char)s[j] > ' ') ++j;
        if (n < max) out[n] = (span){s + start, (int)(j - start)};
        ++n;
    }
    return n;
}

static int32_t parse_int(span f, const char *line) {
    int32_t v;
    if (!jnum_parse_int(f.s, (size_t)f.n, &v)) {
        util_exit("java.lang.NumberFormatException: For input string: \"%.*s\" (%s)", f.n, f.s, line);
    }
    return v;
}

/* PlinkGenMap.fillMapPositions */
static void fill_map_positions(chrom_map *cm) {
    cm->base_pos = util_malloc((size_t)cm->n * sizeof *cm->base_pos);
    cm->gen_pos = util_malloc((size_t)cm->n * sizeof *cm->gen_pos);
    for (int j = 0; j < cm->n; ++j) {
        const char *line = cm->lines[j];
        span f[4];
        if (whitespace_fields(line, strlen(line), f, 4) != 4) util_exit("java.lang.IllegalArgumentException: Map file format error: %s", line);
        cm->base_pos[j] = parse_int(f[3], line);
        if (!jnum_parse_double(f[2].s, (size_t)f[2].n, &cm->gen_pos[j])) {
            util_exit("java.lang.NumberFormatException: For input string: \"%.*s\"", f[2].n, f[2].s);
        }
        if (!isfinite(cm->gen_pos[j])) {
            char s[JNUM_DOUBLE_STRING_SIZE];
            jnum_double_to_string(s, cm->gen_pos[j]);
            util_exit("java.lang.IllegalArgumentException: invalid map position: %s", s);
        }
        if (j > 0) {
            if (cm->base_pos[j] == cm->base_pos[j - 1]) util_exit("java.lang.IllegalArgumentException: duplication position: %s", line);
            if (cm->base_pos[j] < cm->base_pos[j - 1] || cm->gen_pos[j] < cm->gen_pos[j - 1]) {
                util_exit("java.lang.IllegalArgumentException: map positions not in ascending order: %s", line);
            }
        }
    }
    if (cm->n > 0 && cm->gen_pos[0] == cm->gen_pos[cm->n - 1]) {
        char s[JNUM_DOUBLE_STRING_SIZE];
        jnum_double_to_string(s, cm->gen_pos[0]);
        util_exit("java.lang.IllegalArgumentException: All loci in genetic map have the same genetic position [%s]: %s",
                s, cm->lines[0]);
    }
    for (int j = 0; j < cm->n; ++j) free(cm->lines[j]);
    free(cm->lines);
    cm->lines = NULL;
}

genetic_map *genetic_map_open(const char *map_path, const chrom_interval *interval) {
    genetic_map *gm = util_malloc(sizeof *gm);
    *gm = (genetic_map){0};
    if (map_path == NULL) {
        gm->scale_factor = 1e-6;
        gm->inv_scale_factor = 1.0 / gm->scale_factor;
        return gm;
    }
    gm->plink = true;
    /* PlinkGenMap.divideByChrom, keeping only the interval's chromosome if given. */
    const char *keep = interval == NULL ? NULL : chrom_ids_id(interval->chrom_index);
    line_reader *r = line_reader_open(map_path, 1);
    kstring_t line = {0, 0, NULL};
    while (line_reader_next(r, &line)) {
        span f[4];
        int n = whitespace_fields(line.s, line.l, f, 4);
        if (n == 0) continue;
        if (n < 4) util_exit("java.lang.IllegalArgumentException: Map file format error: %s", line.s);
        if (keep != NULL && ((size_t)f[0].n != strlen(keep) || memcmp(f[0].s, keep, (size_t)f[0].n) != 0)) continue;
        int chrom = chrom_ids_index(f[0].s, (size_t)f[0].n);
        if (chrom >= gm->n_chrom) {
            gm->chroms = util_realloc(gm->chroms, (size_t)(chrom + 1) * sizeof *gm->chroms);
            memset(gm->chroms + gm->n_chrom, 0, (size_t)(chrom + 1 - gm->n_chrom) * sizeof *gm->chroms);
            gm->n_chrom = chrom + 1;
        }
        chrom_map *cm = &gm->chroms[chrom];
        if (cm->n == cm->cap) {
            cm->cap = cm->cap == 0 ? 200 : 2 * cm->cap;
            cm->lines = util_realloc(cm->lines, (size_t)cm->cap * sizeof *cm->lines);
        }
        cm->lines[cm->n++] = util_strndup(line.s, line.l);
    }
    free(line.s);
    line_reader_close(r);
    for (int c = 0; c < gm->n_chrom; ++c) fill_map_positions(&gm->chroms[c]);
    return gm;
}

static const chrom_map *check_chrom(const genetic_map *gm, int chrom) {
    if (chrom >= gm->n_chrom || gm->chroms[chrom].n == 0) {
        util_exit("java.lang.IllegalArgumentException: missing genetic map for chromosome %s", chrom_ids_id(chrom));
    }
    return &gm->chroms[chrom];
}

/* The map points used to interpolate at insertion point ins_pt, extended to 5
 * cM inside the map beyond either end (PlinkGenMap.genPos and basePos). */
static void end_points(const chrom_map *cm, int ins_pt, int *a_index, int *b_index) {
    int size_m1 = cm->n - 1;
    *a_index = ins_pt - 1;
    *b_index = ins_pt;
    if (*a_index == size_m1) {
        int p = search_double(cm->gen_pos, cm->n, cm->gen_pos[size_m1] - MIN_END_CM_DIST);
        if (p < 0) p = -p - 2;
        *a_index = p > 0 ? p : 0;
        *b_index = size_m1;
    } else if (*b_index == 0) {
        int p = search_double(cm->gen_pos, cm->n, cm->gen_pos[0] + MIN_END_CM_DIST);
        if (p < 0) p = -p - 1;
        *a_index = 0;
        *b_index = p < size_m1 ? p : size_m1;
    }
}

double genetic_map_gen_pos(const genetic_map *gm, int chrom, int32_t base_pos) {
    if (!gm->plink) return gm->scale_factor * base_pos;
    const chrom_map *cm = check_chrom(gm, chrom);
    int index = jarrays_search_int(cm->base_pos, 0, cm->n, base_pos);
    if (index >= 0) return cm->gen_pos[index];
    int a_index, b_index;
    end_points(cm, -index - 1, &a_index, &b_index);
    int32_t a = cm->base_pos[a_index], b = cm->base_pos[b_index];
    double fa = cm->gen_pos[a_index], fb = cm->gen_pos[b_index];
    /* Java's int subtraction wraps; -fwrapv gives the same. */
    return fa + (((double)(int32_t)(base_pos - a) / (double)(int32_t)(b - a)) * (fb - fa));
}

int32_t genetic_map_base_pos(const genetic_map *gm, int chrom, double gen_pos) {
    if (!gm->plink) {
        int64_t pos = jnum_round_d(gen_pos * gm->inv_scale_factor);
        if (pos > INT32_MAX) {
            util_exit("An estimated base position exceeds the maximum integer value\nIs the window parameter in cM units?");
        }
        return (int32_t)(uint32_t)(uint64_t)pos;  /* (int) of a long keeps the low 32 bits */
    }
    const chrom_map *cm = check_chrom(gm, chrom);
    int index = search_double(cm->gen_pos, cm->n, gen_pos);
    if (index >= 0) return cm->base_pos[index];
    int a_index, b_index;
    end_points(cm, -index - 1, &a_index, &b_index);
    double a = cm->gen_pos[a_index], b = cm->gen_pos[b_index];
    int32_t fa = cm->base_pos[a_index], fb = cm->base_pos[b_index];
    double interp = fa + ((gen_pos - a) / (b - a)) * (int32_t)(fb - fa);
    if (interp >= INT32_MAX) {
        util_exit("An estimated base position exceeds the maximum integer value\n"
                "Are the window parameter and the genetic map in cM units?");
    }
    return (int32_t)(uint32_t)(uint64_t)jnum_round_d(interp);  /* (int) of a long keeps the low 32 bits */
}

void genetic_map_free(genetic_map *gm) {
    for (int c = 0; c < gm->n_chrom; ++c) {
        free(gm->chroms[c].base_pos);
        free(gm->chroms[c].gen_pos);
    }
    free(gm->chroms);
    free(gm);
}
