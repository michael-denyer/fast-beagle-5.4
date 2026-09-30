/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef BGEN_BGEN_FILES_H
#define BGEN_BGEN_FILES_H

#include <stdio.h>

/* Removal of partial BGEN output at process exit. A run creates at most the
 * three members (.bgen, .info, .sample); exit removes the ones it created
 * unless bgen_files_complete ran first. Any thread's exit() can run the
 * removal, so creation and removal are serialized, and removal is terminal:
 * later opens fail. */

/* Once, before starting any thread that can exit. */
void bgen_files_register_cleanup(void);
/* fopen(path, "wb"), recording path for removal; path is borrowed until
 * exit or completion. Returns NULL on an open error or after exit cleanup or
 * completion. */
FILE *bgen_files_open(const char *path);
/* After every member is written and closed: exit keeps them. */
void bgen_files_complete(void);

#endif
