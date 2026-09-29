/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/ChromIds.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "beagleutil/chrom_ids.h"

#include <pthread.h>

#include "blbutil/str_set.h"
#include "blbutil/utilities.h"

/* Shared by the reference reader thread and the main thread. */
static str_set *ids;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

int chrom_ids_index(const char *id, size_t len) {
    if (len == 0) util_exit("id.isEmpty()");
    pthread_mutex_lock(&lock);
    if (ids == NULL) ids = str_set_new();
    int index = str_set_index(ids, id, len);
    pthread_mutex_unlock(&lock);
    return index;
}

const char *chrom_ids_id(int index) {
    pthread_mutex_lock(&lock);
    const char *id = str_set_get(ids, index);
    pthread_mutex_unlock(&lock);
    return id;
}
