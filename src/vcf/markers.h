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
#ifndef VCF_MARKERS_H
#define VCF_MARKERS_H

#include <stdbool.h>

#include "vcf/marker.h"

/* Marker.equals: same chromosome, position, alleles and INFO/END value; the ID
 * is not compared. */
bool marker_equals(const marker *a, const marker *b);

/* The checks in the Markers constructor: a position that is a local maximum or
 * minimum among three consecutive markers on one chromosome, a chromosome that
 * reappears after another, or a duplicate marker, exits with Java's message. */
void markers_check(const marker *const *markers, int n);

#endif
