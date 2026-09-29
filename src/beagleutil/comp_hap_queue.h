/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) beagleutil/CompHapSegment.java, with
 * the priority queue operations Beagle uses on it, and the queue
 * state and updateHeadOfQ that phase/BasicPhaseStates.java,
 * phase/LowFreqPhaseStates.java and imp/ImpStates.java each repeat;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef BEAGLEUTIL_COMP_HAP_QUEUE_H
#define BEAGLEUTIL_COMP_HAP_QUEUE_H

#include <stdbool.h>
#include <stddef.h>

#include "blbutil/int_int_map.h"

/* A reference haplotype segment copied into one composite haplotype. */
typedef struct {
    int hap;
    int start_step;
    int last_ibs_step;
    int comp_hap_index;
} comp_hap_segment;

/* A binary min-heap on last_ibs_step. Segments with equal keys must leave in
 * the same order as in Java. It holds pointers to segments it does not own. */
typedef struct {
    comp_hap_segment **es;
    int size;
    int cap;
} comp_hap_queue;

void comp_hap_queue_init(comp_hap_queue *q, int capacity);
void comp_hap_queue_offer(comp_hap_queue *q, comp_hap_segment *s);
/* NULL when empty. */
comp_hap_segment *comp_hap_queue_poll(comp_hap_queue *q);
static inline comp_hap_segment *comp_hap_queue_peek(const comp_hap_queue *q) {
    return q->size > 0 ? q->es[0] : NULL;
}
static inline void comp_hap_queue_clear(comp_hap_queue *q) {
    q->size = 0;
}
void comp_hap_queue_free(comp_hap_queue *q);

/* The queue of at most max_states segments, the segments it points into, and
 * the step at which each queued haplotype was last an IBS neighbour.
 * Only comp_hap_queue.c reads these fields. */
typedef struct {
    comp_hap_queue q;
    comp_hap_segment *segs;
    int_int_map *last_ibs_step;
    int max_states;
} comp_hap_tracker;

void comp_hap_tracker_init(comp_hap_tracker *t, int max_states);
void comp_hap_tracker_free(comp_hap_tracker *t);
void comp_hap_tracker_clear(comp_hap_tracker *t);
/* One completed observation. index == -1 means an existing haplotype was
 * seen again. Otherwise index now holds the observed haplotype from end_step.
 * old_hap == -1 means a new composite haplotype; otherwise its retired segment
 * held old_hap over [start_step, end_step). Callers map steps to their markers. */
typedef struct {
    int index;
    int old_hap;
    int start_step, end_step;
} comp_hap_change;

/* Observes hap at step, including lazy queue updates, any retirement and the
 * map update. Retires only when full or stale by min_steps; INT_MAX disables
 * staleness. Repeated observations keep the same composite haplotype. */
comp_hap_change comp_hap_tracker_observe(comp_hap_tracker *t, int hap, int step, int min_steps);
/* Random fallback only, after an empty IBS scan: adds a whole haplotype,
 * including repeated draws. Do not resume observations until _clear. */
int comp_hap_tracker_seed(comp_hap_tracker *t, int hap);
int comp_hap_tracker_size(const comp_hap_tracker *t);
/* The unfinished segment of composite haplotype index, until the next
 * observation or clear. It runs from start_step to the end of the window. */
const comp_hap_segment *comp_hap_tracker_segment(const comp_hap_tracker *t, int index);

#endif
