/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "blbutil/str_set.h"

#include <stdlib.h>
#include <string.h>

#include <htslib/khash.h>

#include "blbutil/utilities.h"

KHASH_MAP_INIT_STR(str_index, int)

struct str_set {
    khash_t(str_index) *map;
    char **items;
    int size, cap;
};

str_set *str_set_new(void) {
    str_set *set = util_malloc(sizeof *set);
    set->map = kh_init(str_index);
    if (set->map == NULL) util_oom();
    set->items = NULL;
    set->size = set->cap = 0;
    return set;
}

int str_set_find(const str_set *set, const char *s, size_t len) {
    char stack[256];
    char *key = len < sizeof stack ? stack : util_malloc(len + 1);
    memcpy(key, s, len);
    key[len] = '\0';
    khiter_t k = kh_get(str_index, set->map, key);
    if (key != stack) free(key);
    return k == kh_end(set->map) ? -1 : kh_val(set->map, k);
}

int str_set_index(str_set *set, const char *s, size_t len) {
    int index = str_set_find(set, s, len);
    if (index >= 0) return index;
    if (set->size == set->cap) {
        set->cap = set->cap == 0 ? 8 : 2 * set->cap;
        set->items = util_realloc(set->items, (size_t)set->cap * sizeof *set->items);
    }
    char *key = util_strndup(s, len);
    int absent;
    khiter_t k = kh_put(str_index, set->map, key, &absent);
    if (absent < 0) util_oom();
    kh_val(set->map, k) = set->size;
    set->items[set->size] = key;
    return set->size++;
}

int str_set_size(const str_set *set) {
    return set->size;
}

const char *str_set_get(const str_set *set, int index) {
    return set->items[index];
}

void str_set_free(str_set *set) {
    if (set == NULL) return;
    for (int j = 0; j < set->size; ++j) free(set->items[j]);
    free(set->items);
    kh_destroy(str_index, set->map);
    free(set);
}
