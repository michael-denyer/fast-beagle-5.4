/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) main/Par.java and
 * blbutil/Validate.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "main/par.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bgen/plink2_num.h"
#include "blbutil/utilities.h"
#include "jcompat/jmath.h"
#include "jcompat/jnum.h"

/* Validate.argsToMap, kept in argument order; `used` marks removed keys.
 * Values point into argv, so the parameters taken from them outlive the map. */
typedef struct {
    int n;
    char **key;
    const char **value;
    bool *used;
} args_map;

static void args_to_map(args_map *m, int argc, char **argv) {
    m->n = argc;
    m->key = util_malloc((size_t)argc * sizeof *m->key);
    m->value = util_malloc((size_t)argc * sizeof *m->value);
    m->used = util_malloc((size_t)argc * sizeof *m->used);
    for (int j = 0; j < argc; ++j) {
        const char *arg = argv[j];
        const char *eq = strchr(arg, '=');
        if (eq == NULL) util_exit("missing delimiter character (=): %s", arg);
        if (eq == arg) util_exit("missing key in key-value pair: %s", arg);
        if (eq[1] == '\0') util_exit("missing value in key-value pair: %s", arg);
        m->key[j] = util_strndup(arg, (size_t)(eq - arg));
        m->value[j] = eq + 1;
        m->used[j] = false;
        for (int k = 0; k < j; ++k) {
            if (strcmp(m->key[k], m->key[j]) == 0) util_exit("duplicate arguments: %s", m->key[j]);
        }
    }
}

/* Map.remove: the value for key, or NULL. */
static const char *take(args_map *m, const char *key) {
    for (int j = 0; j < m->n; ++j) {
        if (!m->used[j] && strcmp(m->key[j], key) == 0) {
            m->used[j] = true;
            return m->value[j];
        }
    }
    return NULL;
}

static const char *string_arg(args_map *m, const char *key, bool required) {
    const char *v = take(m, key);
    if (v == NULL && required) util_exit("missing %s argument", key);
    return v;
}

/* Validate.getFile */
static const char *file_arg(args_map *m, const char *key, bool required) {
    const char *name = string_arg(m, key, required);
    if (name == NULL) return NULL;
    struct stat st;
    const char *err = NULL;
    if (stat(name, &st) != 0) err = "File does not exist";
    else if (S_ISDIR(st.st_mode)) err = "File cannot be a directory";
    if (err != NULL) util_exit("\nError     :  %s\nFilename  :  %s", err, name);
    return name;
}

static void check_range(const char *key, int64_t v, int64_t min, int64_t max) {
    if (v < min) util_exit("Error in \"%s\" argument: value=%lld < %lld", key, (long long)v, (long long)min);
    if (v > max) util_exit("Error in \"%s\" argument: value=%lld > %lld", key, (long long)v, (long long)max);
}

static int int_arg(args_map *m, const char *key, int def, int min, int max) {
    const char *v = take(m, key);
    if (v == NULL) return def;
    int32_t x;
    if (!jnum_parse_int(v, strlen(v), &x)) util_exit("%s is not a number", v);
    check_range(key, x, min, max);
    return x;
}

static int64_t long_arg(args_map *m, const char *key, int64_t def, int64_t min, int64_t max) {
    const char *v = take(m, key);
    if (v == NULL) return def;
    int64_t x;
    if (!jnum_parse_long(v, strlen(v), &x)) util_exit("%s is not a number", v);
    check_range(key, x, min, max);
    return x;
}

static float float_arg(args_map *m, const char *key, float def, float min, float max) {
    const char *v = take(m, key);
    if (v == NULL) return def;
    float f;
    if (!jnum_parse_float(v, strlen(v), &f)) util_exit("%s is not a number", v);
    if (isnan(f)) util_exit("Error in \"%s\" argument: value=NaN", key);
    if (f < min || f > max) {
        char fs[JNUM_DOUBLE_STRING_SIZE], bound[JNUM_DOUBLE_STRING_SIZE];
        jnum_float_to_string(fs, f);
        jnum_float_to_string(bound, f < min ? min : max);
        util_exit("Error in \"%s\" argument: value=%s %c %s", key, fs, f < min ? '<' : '>', bound);
    }
    return f;
}

/* A number as plink2 parses its arguments (ScanadvDouble), or NAN when absent. */
static double plink2_number_arg(args_map *m, const char *key) {
    const char *v = take(m, key);
    if (v == NULL) return NAN;
    double x;
    const char *end = plink2_scanadv_double(v, &x);
    if (end == NULL || *end != '\0') util_exit("%s is not a number", v);
    return x;
}

static bool boolean_arg(args_map *m, const char *key, bool def) {
    const char *v = take(m, key);
    if (v == NULL) return def;
    if (strcasecmp(v, "true") == 0 || strcasecmp(v, "t") == 0) return true;
    if (strcasecmp(v, "false") == 0 || strcasecmp(v, "f") == 0) return false;
    util_exit("%s is not \"true\" or \"false\"", v);
}

/* java.io.File's path, normalized in place: repeated slashes become one and a
 * trailing slash is dropped. File.equals compares these paths. */
static char *file_path(char *s) {
    size_t k = 0;
    for (size_t j = 0; s[j] != '\0'; ++j) {
        if (s[j] != '/' || k == 0 || s[k - 1] != '/') s[k++] = s[j];
    }
    if (k > 1 && s[k - 1] == '/') --k;
    s[k] = '\0';
    return s;
}

/* Check log, VCF, BGEN and tabix destinations before their writers open files.
 * Main.checkOutputPrefix compares only the VCF path with ref=, then gt=; every
 * other collision is fast-beagle's. stat also detects aliases through
 * relative paths, symlinks and hard links. */
static void check_output_file(const par *p, const char *suffix) {
    size_t size = strlen(p->out) + strlen(suffix) + 1;
    char *output = util_malloc(size);
    snprintf(output, size, "%s%s", p->out, suffix);
    file_path(output);
    struct stat output_stat;
    bool exists = stat(output, &output_stat) == 0;
    bool vcf = strcmp(suffix, ".vcf.gz") == 0;
    const char *inputs[] = {p->ref, p->gt, p->map, p->excludesamples, p->excludemarkers, p->ped, p->truth};
    for (size_t j = 0; j < sizeof inputs / sizeof *inputs; ++j) {
        if (inputs[j] == NULL) continue;
        char *input = file_path(util_strndup(inputs[j], strlen(inputs[j])));
        bool same_path = strcmp(output, input) == 0;
        if (vcf && j < 2 && same_path) util_exit("ERROR: VCF output file equals input file: %s", input);
        struct stat input_stat;
        if (same_path || (exists && stat(input, &input_stat) == 0
                && output_stat.st_dev == input_stat.st_dev && output_stat.st_ino == input_stat.st_ino)) {
            util_exit(PROGRAM ": output file %s equals input file %s", output, input);
        }
        free(input);
    }
    free(output);
}

/* Main.checkOutputPrefix and Main.parameters, after Par has read every argument. */
static void check_parameters(const par *p) {
    struct stat st;
    if (stat(p->out, &st) == 0 && S_ISDIR(st.st_mode)) {
        util_exit("ERROR: \"out\" parameter cannot be a directory: \"%s\"", p->out);
    }
    check_output_file(p, ".vcf.gz");
    check_output_file(p, ".log");
    if (p->bgen != BGEN_NONE) {
        check_output_file(p, ".bgen");
        check_output_file(p, ".info");
        check_output_file(p, ".sample");
    }
    if (p->tbi) check_output_file(p, ".vcf.gz.tbi");
    if (p->window < 1.1 * p->overlap) {
        util_exit("ERROR: The \"window\" parameter must be at least 1.1 times the \"overlap\" parameter");
    }
}

float par_err(const par *p, int n_haps) {
    return p->err >= 0 ? p->err : par_li_stephens_p_mismatch(n_haps);
}

float par_li_stephens_p_mismatch(int n_haps) {
    double theta = 1 / (jmath_log(n_haps) + 0.5);
    return (float)(theta / (2 * (theta + n_haps)));
}

void par_parse(par *p, int argc, char **argv) {
    const int IMAX = INT32_MAX;
    const float FMIN = FLT_TRUE_MIN, FMAX = FLT_MAX;
    memset(p, 0, sizeof *p);
    p->argc = argc;
    p->argv = argv;
    args_map m;
    args_to_map(&m, argc, argv);

    p->gt = file_arg(&m, "gt", true);
    p->ref = file_arg(&m, "ref", false);
    p->out = string_arg(&m, "out", true);
    p->ped = file_arg(&m, "ped", false);
    p->map = file_arg(&m, "map", false);
    const char *chrom = string_arg(&m, "chrom", false);
    if (chrom != NULL) {
        p->has_chrom_int = chrom_interval_parse(chrom, &p->chrom_int);
        if (!p->has_chrom_int) util_exit("Invalid chrom parameter: %s", chrom);
    }
    p->excludesamples = file_arg(&m, "excludesamples", false);
    p->excludemarkers = file_arg(&m, "excludemarkers", false);

    p->burnin = int_arg(&m, "burnin", 3, 1, IMAX);
    p->iterations = int_arg(&m, "iterations", 12, 1, IMAX);
    p->phase_states = int_arg(&m, "phase-states", 280, 1, IMAX);
    p->step_scale = float_arg(&m, "step-scale", 3.0f, FMIN, FMAX);
    p->rare = float_arg(&m, "rare", 0.002f, FMIN, 0.5f);

    p->impute = boolean_arg(&m, "impute", true);
    p->imp_states = int_arg(&m, "imp-states", 1600, 1, IMAX);
    p->imp_segment = float_arg(&m, "imp-segment", 6.0f, FMIN, FMAX);
    p->imp_step = float_arg(&m, "imp-step", 0.1f, FMIN, FMAX);
    p->imp_nsteps = int_arg(&m, "imp-nsteps", 7, 1, IMAX);
    p->cluster = float_arg(&m, "cluster", 0.005f, 0.0f, FMAX);
    p->ap = boolean_arg(&m, "ap", false);
    p->gp = boolean_arg(&m, "gp", false);

    p->em = boolean_arg(&m, "em", true);
    p->ne = float_arg(&m, "ne", 100000.0f, FMIN, FMAX);
    p->err = float_arg(&m, "err", -FMIN, -FMIN, FMAX);
    p->window = float_arg(&m, "window", 40.0f, FMIN, FMAX);
    p->overlap = float_arg(&m, "overlap", 2.0f, FMIN, FMAX);
    p->buffer = float_arg(&m, "buffer", 1.0f, FMIN, FMAX);
    p->seed = long_arg(&m, "seed", -99999, INT64_MIN, INT64_MAX);
    int raw_nthreads = int_arg(&m, "nthreads", IMAX, 1, IMAX);
    long n_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    p->no_nthreads = raw_nthreads == IMAX;
    p->nthreads = p->no_nthreads ? (n_cpus > 0 ? (int)n_cpus : 1) : raw_nthreads;

    p->truth = file_arg(&m, "truth", false);
    p->trace = string_arg(&m, "trace", false);
    const char *bgen = string_arg(&m, "bgen", false);
    if (bgen == NULL) p->bgen = BGEN_NONE;
    else if (strcmp(bgen, "plink2") == 0) p->bgen = BGEN_PLINK2;
    else if (strcmp(bgen, "phased") == 0) p->bgen = BGEN_PHASED;
    else util_exit("Error in \"bgen\" argument: %s (expected plink2 or phased)", bgen);
    int bgen_bits = int_arg(&m, "bgen-bits", 0, 1, 16);
    if (bgen_bits != 0 && p->bgen == BGEN_NONE) util_exit("bgen-bits needs bgen=plink2 or bgen=phased");
    p->bgen_bits = bgen_bits != 0 ? bgen_bits : 8;
    int bgen_chr_set = int_arg(&m, "bgen-chr-set", 0, 1, 95);
    if (bgen_chr_set != 0 && p->bgen != BGEN_PLINK2) util_exit("bgen-chr-set needs bgen=plink2");
    p->bgen_chr_set = bgen_chr_set != 0 ? bgen_chr_set : 22;
    double min_dr2 = plink2_number_arg(&m, "bgen-min-dr2");
    double min_maf = plink2_number_arg(&m, "bgen-min-maf");
    if (p->bgen != BGEN_PLINK2 && !(isnan(min_dr2) && isnan(min_maf))) {
        util_exit("bgen-min-dr2 and bgen-min-maf need bgen=plink2");
    }
    p->has_bgen_min_dr2 = !isnan(min_dr2);
    p->bgen_min_dr2 = p->has_bgen_min_dr2 ? min_dr2 : 0.0;
    p->bgen_min_maf = isnan(min_maf) ? 0.0 : min_maf;
    if (p->bgen_min_maf < 0.0 || p->bgen_min_maf > 1.0) {
        util_exit("Error in \"bgen-min-maf\" argument: value=%g is not in [0, 1]", p->bgen_min_maf);
    }
    p->tbi = boolean_arg(&m, "tbi", false);

    /* Validate.confirmEmptyMap. Java lists the keys in HashMap order. */
    int unused = 0;
    for (int j = 0; j < m.n; ++j) unused += !m.used[j];
    if (unused > 0) {
        fprintf(stderr, "Error: unrecognized parameter%s", unused == 1 ? ":" : "s:");
        for (int j = 0; j < m.n; ++j) {
            if (!m.used[j]) fprintf(stderr, " %s=%s", m.key[j], m.value[j]);
        }
        fputc('\n', stderr);
        exit(1);
    }
    for (int j = 0; j < m.n; ++j) free(m.key[j]);
    free(m.key);
    free(m.value);
    free(m.used);
    check_parameters(p);
}
