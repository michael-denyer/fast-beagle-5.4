/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) blbutil/InputIt.java and
 * blbutil/BGZipIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BLBUTIL_LINE_READER_H
#define BLBUTIL_LINE_READER_H

#include <stdbool.h>

#include <htslib/kstring.h>

/* InputIt.fromBGZipFile: the lines of a text file. A name ending in ".gz" or
 * ".bgz" is decompressed. A bgzip file splits lines at '\n' and drops one
 * trailing '\r' (BGZipIt); any other file splits at "\n", "\r" or "\r\n"
 * (BufferedReader.readLine), and a bgzip file drops text after its last '\n'.
 * Malformed UTF-8 becomes U+FFFD, as Java's decoding does. */
typedef struct line_reader line_reader;

line_reader *line_reader_open(const char *path, int n_threads);  /* threads decompress BGZF input */
bool line_reader_next(line_reader *r, kstring_t *line);  /* false at end of file */
/* line_reader_next without the UTF-8 sanitising; the caller runs
 * jutf8_sanitize on the line before anything reads its text. */
bool line_reader_next_raw(line_reader *r, kstring_t *line);
const char *line_reader_name(const line_reader *r);      /* File.getName() */
void line_reader_close(line_reader *r);

#endif
