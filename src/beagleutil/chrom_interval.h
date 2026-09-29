/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/ChromInterval.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BEAGLEUTIL_CHROM_INTERVAL_H
#define BEAGLEUTIL_CHROM_INTERVAL_H

#include <stdbool.h>
#include <stdint.h>

typedef struct {
    int chrom_index;
    int32_t start;  /* INT32_MIN when absent */
    int32_t end;    /* INT32_MAX when absent */
} chrom_interval;

/* ChromInterval.parse: false where Java returns null. Registers the chromosome. */
bool chrom_interval_parse(const char *str, chrom_interval *out);
bool chrom_interval_contains(const chrom_interval *ci, int chrom_index, int32_t pos);

#endif
