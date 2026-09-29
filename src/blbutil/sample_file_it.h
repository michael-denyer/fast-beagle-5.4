/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/SampleFileIt.java and
 * blbutil/FileIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BLBUTIL_SAMPLE_FILE_IT_H
#define BLBUTIL_SAMPLE_FILE_IT_H

#include "vcf/marker.h"
#include "vcf/samples.h"

/* The records of a file of samples. Records are untyped because C has no
 * generics: each reader's records have one type (gt_rec or ref_gt_rec), and
 * marker and release take a record of that type, standing in for Java's
 * `E extends GTRec` and for ownership. */
typedef struct {
    const samples *(*samples)(const void *it);
    /* The next record, or NULL at the end. The caller owns it (release). */
    void *(*next)(void *it);
    void (*close)(void *it);
    const marker *(*marker)(const void *rec);
    void (*release)(void *rec);
} sample_file_it_ops;

typedef struct {
    const sample_file_it_ops *ops;
    void *it;
} sample_file_it;

static inline const samples *sample_file_it_samples(sample_file_it f) {
    return f.ops->samples(f.it);
}

static inline void *sample_file_it_next(sample_file_it f) {
    return f.ops->next(f.it);
}

static inline void sample_file_it_close(sample_file_it f) {
    f.ops->close(f.it);
}

#endif
