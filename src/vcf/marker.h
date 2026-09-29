/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) vcf/BasicMarker.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_MARKER_H
#define VCF_MARKER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* A substring: s[0..n) with no terminator. */
typedef struct {
    const char *s;
    int n;
} span;

/* Java's BasicMarker. `fields` holds the ID list joined by ';' when there is
 * one, then REF '\t' ALT unless those match an entry of SNV_PERMS, which
 * field_info indexes instead. REF '\t' ALT is the record's text, which is how
 * BasicMarker prints its parsed alleles. */
typedef struct {
    int chrom_index;
    int32_t pos;
    int32_t end;        /* the INFO END value, -1 when absent */
    int n_alleles;
    uint16_t field_info;
    char *fields;       /* NULL when empty */
    size_t fields_len;
} marker;

/* BasicMarker(String vcfRecord). Exits with Java's message on malformed input. */
void marker_parse(marker *m, const char *rec, size_t len);
void marker_free(marker *m);

const char *marker_chrom(const marker *m);
int marker_n_alleles(const marker *m);
/* Marker.bitsPerAllele: the bits a non-missing allele needs, 0 for one allele. */
int marker_bits_per_allele(const marker *m);
bool marker_has_id(const marker *m);
span marker_id(const marker *m);         /* the IDs joined by ';', "." when there are none */
span marker_alleles(const marker *m);    /* REF '\t' ALT */

#endif
