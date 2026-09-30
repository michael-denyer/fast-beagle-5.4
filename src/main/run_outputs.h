/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef MAIN_RUN_OUTPUTS_H
#define MAIN_RUN_OUTPUTS_H

#include "main/par.h"

/* fast-beagle's output destinations. Paths are immutable and borrowed by
 * their format writers. */
typedef struct run_outputs run_outputs;
typedef enum {
    RUN_OUTPUT_VCF, RUN_OUTPUT_LOG, RUN_OUTPUT_BGEN, RUN_OUTPUT_INFO,
    RUN_OUTPUT_SAMPLE, RUN_OUTPUT_TBI, RUN_OUTPUT_COUNT
} run_output_file;

/* Check enabled destinations against inputs without opening files. Called at
 * the existing parameter validation point. */
void run_outputs_check(const par *p);
run_outputs *run_outputs_new(const par *p);
/* Optional destinations are NULL when disabled. out must outlive the path. */
const char *run_outputs_path(const run_outputs *out, run_output_file file);
void run_outputs_free(run_outputs *out);

#endif
