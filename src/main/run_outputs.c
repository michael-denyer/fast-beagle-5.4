/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "main/run_outputs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "blbutil/utilities.h"

/* Also the preflight order, which determines the first collision reported. */
static const char *const suffixes[RUN_OUTPUT_COUNT] = {
    ".vcf.gz", ".log", ".bgen", ".info", ".sample", ".vcf.gz.tbi"
};

struct run_outputs {
    char *paths[RUN_OUTPUT_COUNT];
};

static bool enabled(const par *p, run_output_file file) {
    if (file == RUN_OUTPUT_TBI) return p->tbi;
    if (file == RUN_OUTPUT_BGEN || file == RUN_OUTPUT_INFO || file == RUN_OUTPUT_SAMPLE)
        return p->bgen != BGEN_NONE;
    return true;
}

static char *destination(const par *p, run_output_file file) {
    size_t size = strlen(p->out) + strlen(suffixes[file]) + 1;
    char *path = util_malloc(size);
    snprintf(path, size, "%s%s", p->out, suffixes[file]);
    return path;
}

/* java.io.File normalization applies only to comparison copies. Writers use
 * the original spelling, including repeated and trailing slashes. */
static char *file_path(char *s) {
    size_t k = 0;
    for (size_t j = 0; s[j] != '\0'; ++j) {
        if (s[j] != '/' || k == 0 || s[k - 1] != '/') s[k++] = s[j];
    }
    if (k > 1 && s[k - 1] == '/') --k;
    s[k] = '\0';
    return s;
}

static void check_file(const par *p, run_output_file file) {
    char *output = file_path(destination(p, file));
    struct stat output_stat;
    bool exists = stat(output, &output_stat) == 0;
    const char *inputs[] = {p->ref, p->gt, p->map, p->excludesamples, p->excludemarkers, p->ped, p->truth};
    for (size_t j = 0; j < sizeof inputs / sizeof *inputs; ++j) {
        if (inputs[j] == NULL) continue;
        char *input = file_path(util_strndup(inputs[j], strlen(inputs[j])));
        bool same_path = strcmp(output, input) == 0;
        if (file == RUN_OUTPUT_VCF && j < 2 && same_path)
            util_exit("ERROR: VCF output file equals input file: %s", input);
        struct stat input_stat;
        if (same_path || (exists && stat(input, &input_stat) == 0
                && output_stat.st_dev == input_stat.st_dev && output_stat.st_ino == input_stat.st_ino)) {
            util_exit(PROGRAM ": output file %s equals input file %s", output, input);
        }
        free(input);
    }
    free(output);
}

void run_outputs_check(const par *p) {
    struct stat st;
    if (stat(p->out, &st) == 0 && S_ISDIR(st.st_mode))
        util_exit("ERROR: \"out\" parameter cannot be a directory: \"%s\"", p->out);
    for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file)
        if (enabled(p, file)) check_file(p, file);
}

run_outputs *run_outputs_new(const par *p) {
    run_outputs *out = util_malloc(sizeof *out);
    *out = (run_outputs){0};
    for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file)
        if (enabled(p, file)) out->paths[file] = destination(p, file);
    return out;
}

const char *run_outputs_path(const run_outputs *out, run_output_file file) {
    return out->paths[file];
}

void run_outputs_free(run_outputs *out) {
    for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file) free(out->paths[file]);
    free(out);
}
