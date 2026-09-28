/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) phase/PbwtPhaser.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/pbwt_phaser.h"

#include <stdbool.h>
#include <stdlib.h>

#include "blbutil/parallel.h"
#include "blbutil/utilities.h"
#include "phase/fwd_pbwt_phaser.h"

/* The positions are strictly increasing, so Java's binary search insertion
 * point is the first position at or above pos. */
static int ins_pt(const marker_map *map, double pos) {
    int lo = 0, hi = map->n;
    while (lo < hi) {
        int mid = (int)((unsigned)(lo + hi) >> 1);
        if (map->gen_pos[mid] < pos) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

int_list pbwt_phaser_windows(const fixed_phase_data *fpd, int nthreads) {
    const marker_map *map = &fpd->stage1_map;
    double total_cm = map->gen_pos[map->n - 1] - map->gen_pos[0];
    double overlap_cm = 1.5;
    double per_thread = (total_cm - overlap_cm) / nthreads;
    double advance_cm = 2 * overlap_cm > per_thread ? 2 * overlap_cm : per_thread;
    int_list windows = {0};
    int start = 0;
    int end = ins_pt(map, map->gen_pos[start] + overlap_cm + advance_cm) + 1;
    while (end < map->n) {
        int_list_add(&windows, start);
        int_list_add(&windows, end);
        start = ins_pt(map, map->gen_pos[end] - overlap_cm) - 1;
        end = ins_pt(map, map->gen_pos[end] + advance_cm) + 1;
    }
    int_list_add(&windows, start);
    int_list_add(&windows, map->n);
    return windows;
}

/* PbwtPhaser.indices: each sample's missing genotypes, and its heterozygotes
 * at or after the stage-1 overlap other than its first heterozygote. */
static void indices(const fixed_phase_data *fpd, int_list *miss, int_list *hets) {
    int n_samples = fpd_n_targ_haps(fpd) >> 1;
    bool *seen_het = util_malloc((size_t)n_samples * sizeof *seen_het);
    for (int s = 0; s < n_samples; ++s) {
        miss[s] = (int_list){0};
        hets[s] = (int_list){0};
        seen_het[s] = false;
    }
    for (int m = 0; m < fpd->n_stage1; ++m) {
        for (int s = 0; s < n_samples; ++s) {
            int a1 = fpd_targ_allele(fpd, m, s << 1);
            int a2 = fpd_targ_allele(fpd, m, (s << 1) | 1);
            if (a1 < 0 || a2 < 0) {
                int_list_add(&miss[s], m);
            } else if (a1 != a2) {
                if (m >= fpd->stage1_overlap && seen_het[s]) int_list_add(&hets[s], m);
                else seen_het[s] = true;
            }
        }
    }
    free(seen_het);
}

/* The heterozygote in [start, overlap_end) nearest copy_start at which the
 * previous window's labels are aligned, or -1. */
static int alignment_het(const int_list *hets, int start, int copy_start, int overlap_end) {
    if (hets->n == 0) return -1;
    int index = 0;
    while (index < hets->n && hets->v[index] < copy_start) ++index;
    if (index == hets->n || (hets->v[index] >= overlap_end && index > 0)) index -= 1;
    int het = hets->v[index];
    return start <= het && het < overlap_end ? het : -1;
}

/* PbwtPhaser.copyHaps: Java swaps the two rows it writes into when the labels
 * switch, which leaves the rows' earlier markers in place. */
static void copy_haps(const fwd_pbwt_phaser *fp, int **hap1, int **hap2, const int_list *hets, int s_start, int s_end,
        int overlap_end) {
    int copy_start = (int)((unsigned)(fp->start + overlap_end) >> 1);
    for (int s = s_start; s < s_end; ++s) {
        int h1 = s << 1;
        int h2 = h1 | 1;
        bool swap = false;
        if (fp->start > 0) {
            int het = alignment_het(&hets[s], fp->start, copy_start, overlap_end);
            swap = het >= 0 && hap1[s][het] == fwd_pbwt_phaser_allele(fp, het, h2)
                    && hap2[s][het] == fwd_pbwt_phaser_allele(fp, het, h1);
        }
        int *row1 = swap ? hap2[s] : hap1[s];
        int *row2 = swap ? hap1[s] : hap2[s];
        for (int m = copy_start; m < fp->end; ++m) {
            row1[m] = fwd_pbwt_phaser_allele(fp, m, h1);
            row2[m] = fwd_pbwt_phaser_allele(fp, m, h2);
        }
    }
}

/* PbwtPhaser.pbwtPhasers and initPhase: the windows' phasers are built in
 * parallel, then each block of samples copies its haplotypes from every
 * window in order. */
typedef struct {
    const fixed_phase_data *fpd;
    const int_list *windows;
    int64_t seed;
    fwd_pbwt_phaser *fps;
    int **hap1, **hap2;
    const int_list *hets, *miss;
    sample_phase *phase;
    int n_samples, block;
} init_ctx;

static void phaser_task(void *worker, int j) {
    const init_ctx *c = worker;
    fwd_pbwt_phaser_init(&c->fps[j], c->fpd, c->windows->v[2 * j], c->windows->v[2 * j + 1],
            (int64_t)((uint64_t)c->seed + (uint64_t)j));
}

static void samples_task(void *worker, int b) {
    const init_ctx *c = worker;
    int s_start = b * c->block;
    int s_end = s_start + c->block < c->n_samples ? s_start + c->block : c->n_samples;
    for (int s = s_start; s < s_end; ++s) {
        c->hap1[s] = util_malloc((size_t)c->fpd->n_stage1 * sizeof **c->hap1);
        c->hap2[s] = util_malloc((size_t)c->fpd->n_stage1 * sizeof **c->hap2);
    }
    int overlap_end = 0;
    for (int j = 0; j < c->windows->n >> 1; ++j) {
        copy_haps(&c->fps[j], c->hap1, c->hap2, c->hets, s_start, s_end, overlap_end);
        overlap_end = c->fps[j].end;
    }
    for (int s = s_start; s < s_end; ++s) {
        sample_phase_init(&c->phase[s], s, c->fpd, c->hap1[s], c->hap2[s], c->hets[s].v, c->hets[s].n, c->miss[s].v,
                c->miss[s].n);
        free(c->hap1[s]);
        free(c->hap2[s]);
    }
}

sample_phase *pbwt_phaser_init_phase(const fixed_phase_data *fpd, int nthreads, int64_t seed) {
    int n_samples = fpd_n_targ_haps(fpd) >> 1;
    int_list *miss = util_malloc((size_t)n_samples * sizeof *miss);
    int_list *hets = util_malloc((size_t)n_samples * sizeof *hets);
    indices(fpd, miss, hets);
    int_list windows = pbwt_phaser_windows(fpd, nthreads);
    int n_windows = windows.n >> 1;
    init_ctx c = {fpd, &windows, seed, util_malloc((size_t)n_windows * sizeof *c.fps),
            util_malloc((size_t)n_samples * sizeof *c.hap1), util_malloc((size_t)n_samples * sizeof *c.hap2), hets, miss,
            util_malloc((size_t)n_samples * sizeof *c.phase), n_samples, 128};
    int n_blocks = (n_samples + c.block - 1) / c.block;
    parallel_for(parallel_threads(nthreads, n_windows), n_windows, &c, 0, phaser_task);
    parallel_for(parallel_threads(nthreads, n_blocks), n_blocks, &c, 0, samples_task);
    for (int j = 0; j < n_windows; ++j) fwd_pbwt_phaser_free(&c.fps[j]);
    for (int s = 0; s < n_samples; ++s) {
        free(hets[s].v);
        free(miss[s].v);
    }
    sample_phase *phase = c.phase;
    free(c.fps);
    free(c.hap1);
    free(c.hap2);
    free(windows.v);
    free(hets);
    free(miss);
    return phase;
}
