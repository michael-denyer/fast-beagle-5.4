/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/ChromIds.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BEAGLEUTIL_CHROM_IDS_H
#define BEAGLEUTIL_CHROM_IDS_H

#include <stddef.h>

/* The process-wide chromosome index, assigned in first-seen order. Java assigns
 * indices from parallel parsing threads, so its order can differ when several
 * new chromosomes share a buffer; indices are only compared for equality, so
 * the order is not observable. */
int chrom_ids_index(const char *id, size_t len);
const char *chrom_ids_id(int index);

#endif
