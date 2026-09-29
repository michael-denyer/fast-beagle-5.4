/* Reads the output of JcompatFixtures.math() and prints it again, recomputing
 * every result with jmath. */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "jcompat/jmath.h"

static double d(const char *hex) {
    uint64_t b = strtoull(hex, NULL, 16);
    double x;
    memcpy(&x, &b, sizeof x);
    return x;
}

static void print_d(double x) {
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    printf("%" PRIx64 "\n", b);
}

static void print_f(float x) {
    uint32_t b;
    memcpy(&b, &x, sizeof b);
    printf("%" PRIx32 "\n", b);
}

int main(void) {
    char op[32], a[32], b[32], result[32];
    char line[256];
    while (fgets(line, sizeof line, stdin)) {
        if (sscanf(line, "%31s %31s %31s %31s", op, a, b, result) < 3) {
            fprintf(stderr, "bad line: %s", line);
            return 1;
        }
        if (strcmp(op, "strict-log") == 0) {
            printf("%s %s ", op, a);
            print_d(jmath_log(d(a)));
        } else if (strcmp(op, "strict-log10") == 0) {
            printf("%s %s ", op, a);
            print_d(jmath_log10(d(a)));
        } else if (strcmp(op, "strict-expm1") == 0) {
            printf("%s %s ", op, a);
            print_d(jmath_expm1(d(a)));
        } else if (strcmp(op, "strict-pow") == 0) {
            printf("%s %s %s ", op, a, b);
            print_d(jmath_pow(d(a), d(b)));
        } else if (strcmp(op, "site-pmismatch") == 0) {
            int n = atoi(a);
            double theta = 1 / ((jmath_log(n) + 0.5));
            printf("%s %s ", op, a);
            print_f((float)(theta / (2 * (theta + n))));
        } else if (strcmp(op, "site-maxnseq") == 0) {
            double exponent = 2 * jmath_log10(atoi(a)) + 1;
            printf("%s %s %" PRId64 "\n", op, a, (int64_t)floor(jmath_pow(2.0, exponent)));
        } else if (strcmp(op, "site-unphased") == 0) {
            printf("%s %s %s ", op, a, b);
            print_f((float)jmath_pow(atoi(a), -1.0 / atoi(b)));
        } else if (strcmp(op, "site-precomb") == 0) {
            printf("%s %s %s ", op, a, b);
            print_f((float)-jmath_expm1(d(a) * d(b)));
        } else {
            fprintf(stderr, "unknown op: %s\n", op);
            return 1;
        }
    }
    return 0;
}
