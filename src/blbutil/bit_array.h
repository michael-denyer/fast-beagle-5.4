/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/BitArray.java and the
 * storePhasing and allele methods of phase/RevPbwtPhaser.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BLBUTIL_BIT_ARRAY_H
#define BLBUTIL_BIT_ARRAY_H

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "blbutil/utilities.h"

/* A zeroed array of n_bits bits. */
static inline uint64_t *bit_array_new(size_t n_bits) {
    size_t size = ((n_bits + 63) >> 6) * sizeof(uint64_t);
    uint64_t *bits = util_malloc(size);
    memset(bits, 0, size);
    return bits;
}

static inline int bit_array_get(const uint64_t *bits, size_t index) {
    return (int)((bits[index >> 6] >> (index & 63)) & 1);
}

static inline void bit_array_set(uint64_t *bits, size_t index) {
    bits[index >> 6] |= (uint64_t)1 << (index & 63);
}

/* BitArray.copyFrom(src, from, to): bits [from, to) of src into dst, a word
 * at a time. */
static inline void bit_array_copy_range(uint64_t *dst, const uint64_t *src, int from, int to) {
    if (from >= to) return;
    int start_word = from >> 6;
    int end_word = (to - 1) >> 6;
    uint64_t start_mask = ~(uint64_t)0 << (from & 63);
    uint64_t end_mask = ~(uint64_t)0 >> ((unsigned)-to & 63);
    if (start_word == end_word) {
        uint64_t mask = start_mask & end_mask;
        dst[start_word] ^= (dst[start_word] ^ src[start_word]) & mask;
        return;
    }
    dst[start_word] ^= (dst[start_word] ^ src[start_word]) & start_mask;
    for (int j = start_word + 1; j < end_word; ++j) dst[j] = src[j];
    dst[end_word] ^= (dst[end_word] ^ src[end_word]) & end_mask;
}

/* BitArray.equal(other, from, to) */
static inline bool bit_array_equal_range(const uint64_t *a, const uint64_t *b, int from, int to) {
    if (from >= to) return true;
    int start_word = from >> 6;
    int end_word = (to - 1) >> 6;
    uint64_t start_mask = ~(uint64_t)0 << (from & 63);
    uint64_t end_mask = ~(uint64_t)0 >> ((unsigned)-to & 63);
    if (start_word == end_word) return ((a[start_word] ^ b[start_word]) & start_mask & end_mask) == 0;
    if (((a[start_word] ^ b[start_word]) & start_mask) != 0) return false;
    for (int j = start_word + 1; j < end_word; ++j) {
        if (a[j] != b[j]) return false;
    }
    return ((a[end_word] ^ b[end_word]) & end_mask) == 0;
}

/* Markers.allele(hapBits, marker): marker m's allele in a haplotype whose
 * markers start at bits hap_bits[m]. */
static inline int bit_array_allele(const uint64_t *bits, const int *hap_bits, int m) {
    int allele = 0;
    for (int b = hap_bits[m], k = 0; b < hap_bits[m + 1]; ++b, ++k) allele |= bit_array_get(bits, (size_t)b) << k;
    return allele;
}

/* Markers.setAllele(marker, allele, bitList): writes marker m's allele into
 * bits [hap_bits[m], hap_bits[m + 1]), setting and clearing. */
static inline void bit_array_set_allele(uint64_t *bits, const int *hap_bits, int m, int allele) {
    for (int b = hap_bits[m], k = 0; b < hap_bits[m + 1]; ++b, ++k) {
        uint64_t mask = (uint64_t)1 << (b & 63);
        if ((allele >> k) & 1) bits[b >> 6] |= mask;
        else bits[b >> 6] &= ~mask;
    }
}

/* BitArray.hash(from, to): the words overlapping bits [from, to), masked to
 * the range but not shifted, folded with XOR, then Long.hashCode. The value
 * depends on where the range sits within its words. */
static inline int32_t bit_array_hash(const uint64_t *words, int from, int to) {
    if (from == to) return 0;
    int start_word = from >> 6;
    int end_word = (to - 1) >> 6;
    uint64_t start_mask = ~(uint64_t)0 << (from & 63);
    uint64_t end_mask = ~(uint64_t)0 >> ((unsigned)-to & 63);
    uint64_t v;
    if (start_word == end_word) {
        v = words[start_word] & start_mask & end_mask;
    } else {
        v = words[start_word] & start_mask;
        for (int j = start_word + 1; j < end_word; ++j) v ^= words[j];
        v ^= words[end_word] & end_mask;
    }
    return (int32_t)(uint32_t)(v ^ (v >> 32));
}

/* Alleles of n_haps haplotypes packed bits_per_allele bits each, low bit first. */
static inline uint64_t *allele_bits_store(const int *alleles, int n_haps, int bits_per_allele) {
    uint64_t *bits = bit_array_new((size_t)n_haps * (size_t)bits_per_allele);
    size_t bit = 0;
    for (int h = 0; h < n_haps; ++h) {
        for (int j = 0; j < bits_per_allele; ++j, ++bit) {
            if ((alleles[h] >> j) & 1) bit_array_set(bits, bit);
        }
    }
    return bits;
}

static inline int allele_bits_get(const uint64_t *bits, int hap, int bits_per_allele) {
    size_t bit = (size_t)hap * (size_t)bits_per_allele;
    int allele = 0;
    for (int j = 0; j < bits_per_allele; ++j, ++bit) allele |= bit_array_get(bits, bit) << j;
    return allele;
}

#endif
