/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/Samples.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_SAMPLES_H
#define VCF_SAMPLES_H

#include <stdbool.h>

typedef struct {
    int n;
    char **ids;
    bool *is_diploid;
} samples;

/* new Samples(ids, isDiploid): takes ownership of both arrays, exits on an
 * empty identifier and warns on standard error about duplicates. */
void samples_init(samples *s, int n, char **ids, bool *is_diploid);
void samples_free(samples *s);

#endif
