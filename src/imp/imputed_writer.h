/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) imp/ImputedVcfWriter.java,
 * imp/RefHapHash.java and WindowWriter.printImputed; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef IMP_IMPUTED_WRITER_H
#define IMP_IMPUTED_WRITER_H

#include "imp/imp_ls.h"
#include "main/window_writer.h"

/* WindowWriter.printImputed: the imputed records for reference markers
 * [start, end), one cluster at a time. Each marker's allele probabilities
 * come from the state probabilities of its cluster, interpolated towards the
 * next cluster between clusters; genotyped markers keep the phased alleles.
 * Writes trace seam T5d. */
void imputed_writer_print(window_writer *ww, const imp_data *id, const state_probs *sp, const par *p, int start, int end);

#endif
