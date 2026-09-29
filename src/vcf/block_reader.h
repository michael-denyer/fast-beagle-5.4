/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/BlockLineReader.java and
 * the parseLines step of vcf/RefIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef VCF_BLOCK_READER_H
#define VCF_BLOCK_READER_H

#include <htslib/kstring.h>

#include "blbutil/line_reader.h"
#include "vcf/ref_gt_rec.h"
#include "vcf/vcf_header.h"

/* The data lines of a reference VCF, read a batch at a time on a reader
 * thread and parsed on n_threads threads from a parser thread, so that
 * reading, parsing and the caller's work overlap. Records come back in file
 * order. Blank lines are not skipped. */
typedef struct block_reader block_reader;

/* Starts the reader and parser threads. reader is positioned after the header, which
 * vcf_header_read left in line; the block_reader takes over line's buffer and
 * reads reader until close, and header must outlive it. */
block_reader *block_reader_open(line_reader *reader, kstring_t line, const vcf_header *header, int n_threads);
/* The next parsed record, or NULL at the end of the file. The caller owns it
 * (ref_gt_rec_release). */
ref_gt_rec *block_reader_next(block_reader *r);
/* Stops the reader and parser threads and releases the records not yet returned. */
void block_reader_close(block_reader *r);

#endif
