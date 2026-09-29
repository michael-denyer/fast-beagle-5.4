/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) ints/IntIntMap.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BLBUTIL_INT_INT_MAP_H
#define BLBUTIL_INT_INT_MAP_H

/* A map from int keys to int values. It has no iteration, so its hash order
 * cannot reach the output. */
typedef struct int_int_map int_int_map;

int_int_map *int_int_map_new(void);
void int_int_map_free(int_int_map *m);
void int_int_map_clear(int_int_map *m);
/* IntIntMap.get(key, sentinel) */
int int_int_map_get(const int_int_map *m, int key, int sentinel);
void int_int_map_put(int_int_map *m, int key, int value);
void int_int_map_remove(int_int_map *m, int key);

#endif
