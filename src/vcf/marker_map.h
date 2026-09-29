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
#ifndef VCF_MARKER_MAP_H
#define VCF_MARKER_MAP_H

#include "vcf/genetic_map.h"
#include "vcf/marker.h"

/* The cM position of each marker in a window, strictly increasing, and the
 * distance from the previous marker (0 for the first). */
typedef struct {
    int n;
    double *gen_pos;
    float *gen_dist;
} marker_map;

/* MarkerMap.create(genMap, meanSingleBaseGenDist(genMap, markers), markers):
 * successive markers are at least the window's mean cM per base apart. Exits
 * as Java does when the markers span two chromosomes or one position. */
void marker_map_init(marker_map *mm, const genetic_map *gm, const marker *const *markers, int n);
/* MarkerMap.restrict(indices): indices must be increasing. */
void marker_map_restrict(marker_map *dst, const marker_map *src, const int *indices, int n);
void marker_map_free(marker_map *mm);

/* Steps: consecutive markers grouped so that each step spans at least min_step
 * cM from its first marker, except possibly the last. */
typedef struct {
    int n;
    int *ends;   /* exclusive end marker of each step */
} steps;

void steps_init(steps *s, const marker_map *map, float min_step);
void steps_free(steps *s);

/* Steps.start(step): the first marker of a step */
static inline int steps_start(const steps *s, int step) {
    return step == 0 ? 0 : s->ends[step - 1];
}

#endif
