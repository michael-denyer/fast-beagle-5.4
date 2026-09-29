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
#ifndef MAIN_PAR_H
#define MAIN_PAR_H

#include <stdbool.h>
#include <stdint.h>

#include "beagleutil/chrom_interval.h"

/* fast-beagle only: the bgen= output mode */
typedef enum { BGEN_NONE, BGEN_PLINK2, BGEN_PHASED } bgen_mode;

/* Beagle's command-line parameters, with Java's defaults and ranges. */
typedef struct {
    /* data parameters; NULL when absent */
    const char *gt, *ref, *out, *ped, *map, *excludesamples, *excludemarkers, *truth;
    bool has_chrom_int;
    chrom_interval chrom_int;

    /* phasing parameters */
    int burnin, iterations, phase_states;
    float step_scale, rare;

    /* imputation parameters */
    bool impute, ap, gp;
    int imp_states, imp_nsteps;
    float imp_segment, imp_step, cluster;

    /* general parameters */
    bool em;
    float ne, err, window, overlap, buffer;
    int nthreads;
    bool no_nthreads;      /* nthreads= was absent */
    int64_t seed;

    /* the arguments as given, for the log's command line */
    int argc;
    char **argv;

    /* fast-beagle only: directory for trace seams, or NULL */
    const char *trace;

    /* fast-beagle only: BGEN output, its bits per probability, its plink2
     * --chr-set autosome count, and its plink2 --extract-if-info "DR2 >= x"
     * and --maf filters; bgen_min_maf 0 means no filter */
    bgen_mode bgen;
    int bgen_bits;
    int bgen_chr_set;
    bool has_bgen_min_dr2;
    double bgen_min_dr2, bgen_min_maf;

    /* fast-beagle only: also write <out>.vcf.gz.tbi */
    bool tbi;
} par;

/* Parses key=value arguments as Par(String[]) does, exiting on invalid input. */
void par_parse(par *p, int argc, char **argv);

/* Par.err(nHaps): the err parameter, or the Li-Stephens value if absent */
float par_err(const par *p, int n_haps);
/* Par.liStephensPMismatch(nHaps) */
float par_li_stephens_p_mismatch(int n_haps);

#endif
