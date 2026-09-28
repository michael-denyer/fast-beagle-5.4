/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.4 (29Oct24) main/RunStats.java and the printing
 * methods of blbutil/Utilities.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
/* clock_gettime and localtime_r are POSIX, which glibc hides under -std=c11.
 * macOS shows them anyway, and this macro would hide its ru_maxrss. */
#ifndef __APPLE__
#define _POSIX_C_SOURCE 200809L
#endif

#include "main/run_stats.h"

#include <htslib/kstring.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

#include "blbutil/utilities.h"
#include "jcompat/jnum.h"

/* The first two lines of Main.SHORT_HELP, naming this program. Beagle's third
 * line points to a usage text that this program does not print. */
static const char BANNER[] = "fast-beagle: a C port of beagle.29Oct24.c8e.jar (version 5.4)\n"
                             "Copyright (C) 2014-2022 Brian L. Browning\n";

int64_t run_stats_nanos(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000 + ts.tv_nsec;
}

/* Utilities.duoPrint */
__attribute__((format(printf, 2, 3))) static void duo_print(run_stats *rs, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    va_list ap2;
    va_copy(ap2, ap);
    vfprintf(stdout, fmt, ap);
    vfprintf(rs->log, fmt, ap2);
    va_end(ap2);
    va_end(ap);
}

/* stdout is line buffered (run_stats_open), so only the log needs flushing. */
static void flush(run_stats *rs) {
    fflush(rs->log);
}

/* String.format("%,d", n) in Locale.US */
static const char *grouped(char *buf, int64_t n) {
    char digits[24];
    int len = snprintf(digits, sizeof digits, "%lld", (long long)llabs(n));
    char *p = buf;
    if (n < 0) *p++ = '-';
    for (int j = 0; j < len; ++j) {
        if (j > 0 && (len - j) % 3 == 0) *p++ = ',';
        *p++ = digits[j];
    }
    *p = '\0';
    return buf;
}

/* Utilities.timeStamp: "hh:mm a z 'on' dd MMM yyyy" in Locale.US */
static const char *time_stamp(char *buf, size_t size) {
    time_t now = time(NULL);
    struct tm tm;
    localtime_r(&now, &tm);
    strftime(buf, size, "%I:%M %p %Z on %d %b %Y", &tm);
    return buf;
}

/* Utilities.elapsedNanos: "H hours M minutes S seconds" */
static const char *elapsed(char *buf, size_t size, int64_t nanos) {
    int64_t seconds = jnum_round_d(nanos / 1000000000.0);
    int n = 0;
    if (seconds >= 3600) {
        int64_t hours = seconds / 3600;
        n += snprintf(buf + n, size - n, "%lld %s", (long long)hours, hours == 1 ? "hour " : "hours ");
        seconds %= 3600;
    }
    if (seconds >= 60) {
        int64_t minutes = seconds / 60;
        n += snprintf(buf + n, size - n, "%lld %s", (long long)minutes, minutes == 1 ? "minute " : "minutes ");
        seconds %= 60;
    }
    snprintf(buf + n, size - n, "%lld %s", (long long)seconds, seconds == 1 ? "second" : "seconds");
    return buf;
}

/* RunStats.duoPrintNanos, with the message padded as String.format("%1$-31s") */
static void print_nanos(run_stats *rs, const char *message, int64_t nanos) {
    char buf[96];
    duo_print(rs, "%-31s%s\n", message, elapsed(buf, sizeof buf, nanos));
}

/* String.format("%1$7.1e", (double)x) for x >= 0. Java's Formatter rounds the
 * shortest decimal digits of the double half up, where printf rounds the
 * exact binary value half to even. */
static const char *sci1(char *buf, size_t size, float x) {
    char s[JNUM_DOUBLE_STRING_SIZE];
    jnum_double_to_string(s, x);
    char *e = strchr(s, 'E');
    int exp10 = 0;
    if (e != NULL) {
        exp10 = atoi(e + 1);
        *e = '\0';
    }
    char dig[JNUM_DOUBLE_STRING_SIZE + 2] = {0};
    int n_dig = 0, n_seen = 0, point = -1, first = -1;
    for (const char *c = s; *c != '\0'; ++c) {
        if (*c == '.') {
            point = n_seen;
            continue;
        }
        if (first < 0 && *c != '0') first = n_seen;
        if (first >= 0) dig[n_dig++] = *c;
        ++n_seen;
    }
    if (point < 0) point = n_seen;
    if (first < 0) {
        snprintf(buf, size, "%7s", "0.0e+00");
        return buf;
    }
    exp10 += point - first - 1;
    int d0 = dig[0] - '0';
    int d1 = n_dig > 1 ? dig[1] - '0' : 0;
    if (n_dig > 2 && dig[2] >= '5' && ++d1 == 10) {
        d1 = 0;
        if (++d0 == 10) {
            d0 = 1;
            ++exp10;
        }
    }
    char t[16];
    snprintf(t, sizeof t, "%d.%de%c%02d", d0, d1, exp10 < 0 ? '-' : '+', abs(exp10));
    snprintf(buf, size, "%7s", t);
    return buf;
}

void run_stats_open(run_stats *rs, const par *p, const char *program) {
    *rs = (run_stats){.par = p, .start_nanos = run_stats_nanos()};
    kstring_t path = KS_INITIALIZE;
    ksprintf(&path, "%s.log", p->out);
    rs->log = fopen(path.s, "w");
    if (rs->log == NULL) util_exit("Error opening %s", path.s);
    ks_free(&path);
    util_exit_log(rs->log);
    /* Java's System.out flushes at each line; a pipe would otherwise hold the
     * progress lines until the run ends. */
    setvbuf(stdout, NULL, _IOLBF, 0);

    char ts[64];
    duo_print(rs, "%s", BANNER);
    duo_print(rs, "Start time: %s\n", time_stamp(ts, sizeof ts));
    duo_print(rs, "\nCommand line: %s\n", program);
    for (int j = 0; j < p->argc; ++j) duo_print(rs, "  %s\n", p->argv[j]);
    if (p->no_nthreads) duo_print(rs, "  nthreads=%d\n", p->nthreads);
    if (p->map == NULL) duo_print(rs, "\nNo genetic map is specified: using 1 cM = 1 Mb\n");
    flush(rs);
}

void run_stats_sample_summary(run_stats *rs, int n_ref_samples, int n_targ_samples) {
    char buf[32];
    duo_print(rs, "\n");
    duo_print(rs, "Reference samples: %20s\n", grouped(buf, n_ref_samples));
    duo_print(rs, "Study     samples: %20s\n", grouped(buf, n_targ_samples));
    flush(rs);
}

static const marker *window_marker(const window *w, bool ref, int m) {
    return ref ? &w->ref[m]->marker : &w->targ[m]->marker;
}

void run_stats_window_update(run_stats *rs, const window *w) {
    bool ref = rs->par->ref != NULL;
    int n_markers = ref ? w->n_ref : w->n_targ;
    const marker *first = window_marker(w, ref, 0);
    const marker *last = window_marker(w, ref, n_markers - 1);
    const char *chr = marker_chrom(first);
    duo_print(rs, "\nWindow %d [", w->index);
    if (strcmp(chr, ".") != 0) duo_print(rs, "%s:", chr);
    duo_print(rs, "%d-", (int)first->pos);
    if (strcmp(chr, marker_chrom(last)) != 0) duo_print(rs, "%s:", marker_chrom(last));
    duo_print(rs, "%d]\n", (int)last->pos);
    char buf[32];
    if (ref) duo_print(rs, "Reference markers: %20s\n", grouped(buf, n_markers));
    duo_print(rs, "Study     markers: %20s\n", grouped(buf, w->n_targ));
    flush(rs);
}

void run_stats_stage1(run_stats *rs, const phase_data *pd, int64_t nanos) {
    const par *p = rs->par;
    if (pd->it == p->burnin && p->em) {
        /* PhaseData.ne() */
        int64_t ne = (int64_t)ceil(25 * pd->recomb_intensity * pd->fpd->n_haps);
        char err[16];
        duo_print(rs, "\n");
        duo_print(rs, "%-31s%lld\n", "Estimated ne:", (long long)ne);
        duo_print(rs, "%-31s%s\n", "Estimated err:", sci1(err, sizeof err, pd->p_mismatch));
    }
    rs->total_phase_nanos += nanos;
    bool burnin = pd->it < p->burnin;
    int it = burnin ? pd->it : pd->it - p->burnin;   /* counted from 0 within its stage */
    if (it == 0) duo_print(rs, "\n");
    char msg[48];
    snprintf(msg, sizeof msg, "%s iteration %d:", burnin ? "Burnin " : "Phasing", it + 1);
    print_nanos(rs, msg, nanos);
    flush(rs);
}

void run_stats_stage2(run_stats *rs, int64_t nanos) {
    rs->total_phase_nanos += nanos;
    print_nanos(rs, "Low frequency phasing:", nanos);
    flush(rs);
}

void run_stats_imputation(run_stats *rs, int64_t nanos) {
    rs->total_impute_nanos += nanos;
    duo_print(rs, "\n");
    print_nanos(rs, "Imputation time:", nanos);
    flush(rs);
}

/* getrusage's max resident set size in bytes: Linux reports kilobytes and
 * macOS bytes. */
static int64_t max_rss_bytes(const struct rusage *ru) {
#ifdef __APPLE__
    return ru->ru_maxrss;
#else
    return (int64_t)ru->ru_maxrss * 1024;
#endif
}

static int64_t timeval_nanos(struct timeval tv) {
    return (int64_t)tv.tv_sec * 1000000000 + (int64_t)tv.tv_usec * 1000;
}

void run_stats_close(run_stats *rs, int64_t n_targ_markers, int64_t n_markers) {
    int64_t total = run_stats_nanos() - rs->start_nanos;
    char buf[32];
    duo_print(rs, "\n");
    duo_print(rs, "Cumulative Statistics:\n\n");
    if (n_targ_markers != n_markers) duo_print(rs, "Reference markers: %20s\n", grouped(buf, n_markers));
    duo_print(rs, "Study     markers: %20s\n\n", grouped(buf, n_targ_markers));
    if (rs->total_phase_nanos > 1000) print_nanos(rs, "Haplotype phasing time:", rs->total_phase_nanos);
    if (rs->total_impute_nanos > 0) print_nanos(rs, "Imputation time:", rs->total_impute_nanos);
    print_nanos(rs, "Total time:", total);
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    print_nanos(rs, "CPU time:", timeval_nanos(ru.ru_utime) + timeval_nanos(ru.ru_stime));
    duo_print(rs, "%-31s%s MB\n", "Max memory:", grouped(buf, max_rss_bytes(&ru) / (1024 * 1024)));
    char ts[64];
    duo_print(rs, "\nEnd time: %s\n", time_stamp(ts, sizeof ts));
    duo_print(rs, "fast-beagle finished\n");
    util_exit_log(NULL);
    if (fclose(rs->log) != 0) util_exit("Error writing %s.log", rs->par->out);
}
