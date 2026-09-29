/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) ints/IntIntMap.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "blbutil/int_int_map.h"

#include <stdlib.h>

#include <htslib/khash.h>

#include "blbutil/utilities.h"

KHASH_MAP_INIT_INT(int_int, int)

struct int_int_map {
    khash_t(int_int) *h;
};

int_int_map *int_int_map_new(void) {
    int_int_map *m = util_malloc(sizeof *m);
    m->h = kh_init(int_int);
    return m;
}

void int_int_map_free(int_int_map *m) {
    kh_destroy(int_int, m->h);
    free(m);
}

void int_int_map_clear(int_int_map *m) {
    kh_clear(int_int, m->h);
}

int int_int_map_get(const int_int_map *m, int key, int sentinel) {
    khiter_t k = kh_get(int_int, m->h, key);
    return k == kh_end(m->h) ? sentinel : kh_value(m->h, k);
}

void int_int_map_put(int_int_map *m, int key, int value) {
    int absent;
    khiter_t k = kh_put(int_int, m->h, key, &absent);
    if (absent < 0) util_exit(PROGRAM ": out of memory");
    kh_value(m->h, k) = value;
}

void int_int_map_remove(int_int_map *m, int key) {
    khiter_t k = kh_get(int_int, m->h, key);
    if (k != kh_end(m->h)) kh_del(int_int, m->h, k);
}
