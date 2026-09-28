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
#ifndef VCF_FILTER_UTIL_H
#define VCF_FILTER_UTIL_H

#include <stdbool.h>

#include "blbutil/str_set.h"
#include "vcf/marker.h"

/* Utilities.idSet: the non-blank trimmed lines of a file, one identifier each,
 * or NULL when path is NULL. */
str_set *filter_id_set(const char *path);

/* FilterUtil.markerFilter: false when any ';'-separated ID of the marker, or
 * its CHROM:POS, is in the exclusion set. A NULL set accepts every marker. */
bool filter_accept_marker(const str_set *exclude, const marker *m);

#endif
