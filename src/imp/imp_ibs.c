/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) imp/CodedSteps.java and
 * imp/ImpIbs.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "imp/imp_ibs.h"

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/int_list.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"

/* CodedSteps.stepStarts: the first step is half of imp-step. Cluster
 * positions increase strictly, so the insertion point is a lower bound. */
static int_list step_starts(const imp_data *id, const par *p) {
    int_list starts = {0};
    int_list_add(&starts, 0);
    double step = p->imp_step;
    double next_pos = id->pos[0] + step / 2;
    int index = 0;
    for (;;) {
        while (index < id->n_clusters && id->pos[index] < next_pos) ++index;
        if (index >= id->n_clusters) break;
        int_list_add(&starts, index);
        next_pos = id->pos[index] + step;
    }
    return starts;
}

/* CodedSteps.codeStep: sequences over the step's clusters, numbered from 1
 * in order of first appearance among the target haplotypes; 0 for a
 * reference haplotype no target haplotype matches. */
static int code_step(const imp_data *id, int start, int end, int *seq) {
    int n_ref = id->n_ref_haps;
    int n_haps = id->n_haps;
    for (int h = 0; h < n_haps; ++h) seq[h] = 1;
    int seq_cnt = 2;
    for (int c = start; c < end; ++c) {
        int n_alleles = id->hap_to_seq[c].n_seq;
        size_t size = (size_t)seq_cnt * (size_t)n_alleles;
        int *seq_map = util_malloc(size * sizeof *seq_map);
        memset(seq_map, 0, size * sizeof *seq_map);
        seq_cnt = 1;
        for (int h = n_ref; h < n_haps; ++h) {
            int index = n_alleles * seq[h] + imp_data_seq(id, c, h);
            if (seq_map[index] == 0) seq_map[index] = seq_cnt++;
            seq[h] = seq_map[index];
        }
        for (int h = 0; h < n_ref; ++h) {
            if (seq[h] != 0) seq[h] = seq_map[seq[h] * n_alleles + imp_data_seq(id, c, h)];
        }
        free(seq_map);
    }
    return seq_cnt;
}

typedef struct {
    int_list *v;
    int n, cap;
} list_list;

static void list_list_add(list_list *l, int_list x) {
    if (l->n == l->cap) {
        l->cap = l->cap == 0 ? 8 : 2 * l->cap;
        l->v = util_realloc(l->v, (size_t)l->cap * sizeof *l->v);
    }
    l->v[l->n++] = x;
}

typedef struct {
    const imp_data *id;
    int64_t seed;
    int n_haps_per_step;
    imp_ibs *ib;
    int step;
} ibs_ctx;

/* The number of reference haplotypes in an increasing list. */
static int n_ref_in(const int_list *l, int n_ref_haps) {
    int lo = 0, hi = l->n;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (l->v[mid] < n_ref_haps) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static const ibs_set *new_set(ibs_ctx *c, const int *haps, int n) {
    imp_ibs *ib = c->ib;
    ibs_set *set = &ib->owned[c->step][ib->n_owned[c->step]++];
    set->n = n;
    set->haps = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *set->haps);
    if (n > 0) memcpy(set->haps, haps, (size_t)n * sizeof *haps);
    return set;
}

static void set_result(ibs_ctx *c, const int_list *child, int first_targ, const ibs_set *set) {
    for (int j = first_targ; j < child->n; ++j) c->ib->ibs[c->step][child->v[j] - c->id->n_ref_haps] = set;
}

/* ImpIbs.partition: the parent's haplotypes grouped by sequence, a group per
 * target sequence in order of first appearance, each in increasing order. */
static list_list partition(const int_list *parent, int first_targ, const int *seq, int n_seq) {
    int *slot = util_malloc((size_t)n_seq * sizeof *slot);
    for (int s = 0; s < n_seq; ++s) slot[s] = -1;
    list_list children = {0};
    for (int k = first_targ; k < parent->n; ++k) {
        int s = seq[parent->v[k]];
        if (slot[s] < 0) {
            slot[s] = children.n;
            list_list_add(&children, (int_list){0});
        }
    }
    for (int k = 0; k < parent->n; ++k) {
        int s = seq[parent->v[k]];
        if (slot[s] >= 0) int_list_add(&children.v[slot[s]], parent->v[k]);
    }
    free(slot);
    return children;
}

/* ImpIbs.ibsHaps: the child's reference haplotypes and a random subset of the
 * parent's other reference haplotypes, up to n_haps_per_step, sorted. */
static const ibs_set *ibs_haps(ibs_ctx *c, const int_list *parent, const int_list *child, int n_child_ref) {
    int n_ref_haps = c->id->n_ref_haps;
    int_list combined = {0};
    for (int j = 0; j < n_child_ref; ++j) int_list_add(&combined, child->v[j]);
    int size = c->n_haps_per_step - n_child_ref;
    jrandom r;
    jrandom_init(&r, (int64_t)((uint64_t)c->seed + (uint64_t)(int64_t)parent->v[0]));
    /* ImpIbs.uniqToParent */
    int_list uniq = {0};
    int n_child_ref_m1 = n_child_ref - 1;
    int n_parent_ref = n_ref_in(parent, n_ref_haps);
    int ci = 0;
    int c_val = child->v[ci];
    for (int p = 0; p < n_parent_ref; ++p) {
        int p_val = parent->v[p];
        while (c_val < p_val && ci < n_child_ref_m1) c_val = child->v[++ci];
        if (p_val != c_val) int_list_add(&uniq, p_val);
    }
    /* ImpIbs.randomSubset */
    if (uniq.n < size) size = uniq.n;
    util_shuffle(uniq.v, uniq.n, size, &r);
    for (int j = 0; j < size; ++j) int_list_add(&combined, uniq.v[j]);
    if (combined.n > 1) qsort(combined.v, (size_t)combined.n, sizeof *combined.v, util_compare_ints);
    const ibs_set *set = new_set(c, combined.v, combined.n);
    free(combined.v);
    free(uniq.v);
    return set;
}

static void free_lists(list_list *l) {
    for (int j = 0; j < l->n; ++j) free(l->v[j].v);
    free(l->v);
}

/* ImpIbs.getIbsHaps for one step. */
static void step_ibs_haps(ibs_ctx *c, int n_steps_to_merge) {
    imp_ibs *ib = c->ib;
    const imp_data *id = c->id;
    int n_ref = id->n_ref_haps;
    int index = c->step;
    /* initPartition: every haplotype, grouped like a parent holding them all */
    int_list all = {0};
    for (int h = 0; h < id->n_haps; ++h) int_list_add(&all, h);
    list_list children = partition(&all, n_ref, ib->step_seq[index], ib->step_n_seq[index]);
    free(all.v);
    list_list next_parents = {0};
    for (int j = 0; j < children.n; ++j) {
        int n_ref_in_child = n_ref_in(&children.v[j], n_ref);
        if (n_ref_in_child <= c->n_haps_per_step) {
            set_result(c, &children.v[j], n_ref_in_child, new_set(c, children.v[j].v, n_ref_in_child));
            free(children.v[j].v);
        } else {
            list_list_add(&next_parents, children.v[j]);
        }
    }
    free(children.v);
    for (int i = 1; i < n_steps_to_merge; ++i) {
        list_list parents = next_parents;
        next_parents = (list_list){0};
        for (int j = 0; j < parents.n; ++j) {
            const int_list *parent = &parents.v[j];
            list_list kids = partition(parent, n_ref_in(parent, n_ref), ib->step_seq[index + i], ib->step_n_seq[index + i]);
            for (int k = 0; k < kids.n; ++k) {
                int n_child_ref = n_ref_in(&kids.v[k], n_ref);
                if (n_child_ref <= c->n_haps_per_step) {
                    set_result(c, &kids.v[k], n_child_ref, ibs_haps(c, parent, &kids.v[k], n_child_ref));
                    free(kids.v[k].v);
                } else {
                    list_list_add(&next_parents, kids.v[k]);
                }
            }
            free(kids.v);
        }
        free_lists(&parents);
    }
    /* finalUpdateResults */
    for (int j = 0; j < next_parents.n; ++j) {
        const int_list *child = &next_parents.v[j];
        int n_ref_in_child = n_ref_in(child, n_ref);
        int *list = util_malloc((size_t)(n_ref_in_child > 0 ? n_ref_in_child : 1) * sizeof *list);
        memcpy(list, child->v, (size_t)n_ref_in_child * sizeof *list);
        int n = n_ref_in_child;
        if (c->n_haps_per_step < n) {
            jrandom r;
            jrandom_init(&r, (int64_t)((uint64_t)c->seed + (uint64_t)(int64_t)child->v[0]));
            util_shuffle(list, n, n, &r);
            n = c->n_haps_per_step;
            qsort(list, (size_t)n, sizeof *list, util_compare_ints);
        }
        set_result(c, child, n_ref_in_child, new_set(c, list, n));
        free(list);
    }
    free_lists(&next_parents);
}

static void trace(const imp_ibs *ib, const imp_data *id) {
    kstring_t s = {0, 0, NULL};
    kputs("steps", &s);
    for (int j = 0; j < ib->n_steps; ++j) {
        uint64_t h = TRACE_FNV_BASIS;
        for (int hap = 0; hap < id->n_haps; ++hap) h = trace_fold(h, (uint32_t)ib->step_seq[j][hap]);
        ksprintf(&s, "\t%d:%d:%" PRIx64, ib->step_starts[j], ib->step_n_seq[j], h);
    }
    trace_line("T5b", "%s", s.s);
    for (int j = 0; j < ib->n_steps; ++j) {
        for (int t = 0; t < ib->n_targ_haps; ++t) {
            const ibs_set *set = ib->ibs[j][t];
            s.l = 0;
            ksprintf(&s, "ibs\t%d\t%d\t", j, t);
            for (int k = 0; k < set->n; ++k) {
                if (k > 0) kputc(',', &s);
                kputw(set->haps[k], &s);
            }
            trace_line("T5b", "%s", s.s);
        }
    }
    free(s.s);
}

void imp_ibs_init(imp_ibs *ib, const imp_data *id, const par *p) {
    int_list starts = step_starts(id, p);
    ib->n_steps = starts.n;
    ib->step_starts = starts.v;
    ib->n_targ_haps = id->n_targ_haps;
    ib->step_seq = util_malloc((size_t)ib->n_steps * sizeof *ib->step_seq);
    ib->step_n_seq = util_malloc((size_t)ib->n_steps * sizeof *ib->step_n_seq);
    for (int j = 0; j < ib->n_steps; ++j) {
        int end = j + 1 < ib->n_steps ? ib->step_starts[j + 1] : id->n_clusters;
        ib->step_seq[j] = util_malloc((size_t)id->n_haps * sizeof **ib->step_seq);
        ib->step_n_seq[j] = code_step(id, ib->step_starts[j], end, ib->step_seq[j]);
    }
    int n_steps_per_segment = (int)jnum_round_f(p->imp_segment / p->imp_step);
    if (n_steps_per_segment == 0) util_exit("java.lang.ArithmeticException: / by zero");
    ibs_ctx c = {id, p->seed, p->imp_states / n_steps_per_segment, ib, 0};
    ib->ibs = util_malloc((size_t)ib->n_steps * sizeof *ib->ibs);
    ib->owned = util_malloc((size_t)ib->n_steps * sizeof *ib->owned);
    ib->n_owned = util_malloc((size_t)ib->n_steps * sizeof *ib->n_owned);
    for (int j = 0; j < ib->n_steps; ++j) {
        ib->ibs[j] = util_malloc((size_t)id->n_targ_haps * sizeof **ib->ibs);
        /* each target haplotype lies in exactly one set */
        ib->owned[j] = util_malloc((size_t)id->n_targ_haps * sizeof **ib->owned);
        ib->n_owned[j] = 0;
        c.step = j;
        int n_to_merge = p->imp_nsteps < ib->n_steps - j ? p->imp_nsteps : ib->n_steps - j;
        step_ibs_haps(&c, n_to_merge);
    }
    if (trace_on()) trace(ib, id);
}

void imp_ibs_free(imp_ibs *ib) {
    for (int j = 0; j < ib->n_steps; ++j) {
        for (int k = 0; k < ib->n_owned[j]; ++k) free(ib->owned[j][k].haps);
        free(ib->owned[j]);
        free((void *)ib->ibs[j]);
        free(ib->step_seq[j]);
    }
    free(ib->owned);
    free(ib->n_owned);
    free((void *)ib->ibs);
    free(ib->step_seq);
    free(ib->step_n_seq);
    free(ib->step_starts);
}
