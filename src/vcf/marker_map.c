/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/MarkerMap.java,
 * vcf/GeneticMap.java and vcf/Steps.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/marker_map.h"

#include <math.h>
#include <stdlib.h>

#include "blbutil/int_list.h"
#include "blbutil/utilities.h"

static double mean_single_base_gen_dist(const genetic_map *gm, const marker *const *markers, int n) {
    const marker *a = markers[0];
    const marker *b = markers[n - 1];
    if (a->chrom_index != b->chrom_index) util_exit("java.lang.IllegalArgumentException: inconsistent data");
    if (a->pos == b->pos) {
        util_exit("java.lang.IllegalArgumentException: Window has only one position: CHROM=%s POS=%d",
                marker_chrom(a), a->pos);
    }
    double dist = fabs(genetic_map_gen_pos(gm, b->chrom_index, b->pos) - genetic_map_gen_pos(gm, a->chrom_index, a->pos))
            / abs(b->pos - a->pos);
    return dist > 1e-8 ? dist : 1e-8;
}

static void set_gen_pos(marker_map *mm, double *gen_pos, int n) {
    mm->n = n;
    mm->gen_pos = gen_pos;
    mm->gen_dist = util_malloc((size_t)n * sizeof *mm->gen_dist);
    mm->gen_dist[0] = 0.0f;
    for (int j = 1; j < n; ++j) {
        float d = (float)(gen_pos[j] - gen_pos[j - 1]);
        mm->gen_dist[j] = d < 1e-7f ? 1e-7f : d;
    }
}

void marker_map_init(marker_map *mm, const genetic_map *gm, const marker *const *markers, int n) {
    double min_gen_dist = mean_single_base_gen_dist(gm, markers, n);
    double *gen_pos = util_malloc((size_t)n * sizeof *gen_pos);
    int chrom = markers[0]->chrom_index;
    gen_pos[0] = genetic_map_gen_pos(gm, chrom, markers[0]->pos);
    double last_map_pos = gen_pos[0];
    for (int j = 1; j < n; ++j) {
        double map_pos = genetic_map_gen_pos(gm, chrom, markers[j]->pos);
        double dist = map_pos - last_map_pos;
        gen_pos[j] = gen_pos[j - 1] + (dist > min_gen_dist ? dist : min_gen_dist);
        last_map_pos = map_pos;
    }
    set_gen_pos(mm, gen_pos, n);
}

void marker_map_restrict(marker_map *dst, const marker_map *src, const int *indices, int n) {
    double *gen_pos = util_malloc((size_t)n * sizeof *gen_pos);
    for (int j = 0; j < n; ++j) {
        if (j > 0 && indices[j] <= indices[j - 1]) util_exit("java.lang.IllegalArgumentException: %d", indices[j]);
        gen_pos[j] = src->gen_pos[indices[j]];
    }
    set_gen_pos(dst, gen_pos, n);
}

void marker_map_free(marker_map *mm) {
    free(mm->gen_pos);
    free(mm->gen_dist);
}

void steps_init(steps *s, const marker_map *map, float min_step) {
    int_list ends = {0};
    int end = 0;
    while (end < map->n) {
        double min_gen_pos = map->gen_pos[end] + (double)min_step;
        ++end;
        while (end < map->n && map->gen_pos[end] < min_gen_pos) ++end;
        int_list_add(&ends, end);
    }
    s->n = ends.n;
    s->ends = ends.v;
}

void steps_free(steps *s) {
    free(s->ends);
}
