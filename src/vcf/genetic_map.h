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
#ifndef VCF_GENETIC_MAP_H
#define VCF_GENETIC_MAP_H

#include <stdint.h>

#include "beagleutil/chrom_interval.h"

/* Converts between base positions and cM. Without a map file, 1 cM is 1 Mb
 * (PositionMap). With a PLINK map, positions are interpolated linearly, and
 * beyond either end of a chromosome's map the line through the end point and
 * the map point 5 cM inside it is used (PlinkGenMap). */
typedef struct genetic_map genetic_map;

/* GeneticMap.geneticMap(file, chromInt): map_path may be NULL. With an
 * interval, only its chromosome is read from the map file. */
genetic_map *genetic_map_open(const char *map_path, const chrom_interval *interval);
double genetic_map_gen_pos(const genetic_map *gm, int chrom_index, int32_t base_pos);
int32_t genetic_map_base_pos(const genetic_map *gm, int chrom_index, double gen_pos);
void genetic_map_free(genetic_map *gm);

#endif
