/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "jnum.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "jutf8.h"

int32_t jnum_d2i(double x) {
    if (isnan(x)) return 0;
    if (x >= 2147483647.0) return INT32_MAX;
    if (x <= -2147483648.0) return INT32_MIN;
    return (int32_t)x;
}

int64_t jnum_d2l(double x) {
    if (isnan(x)) return 0;
    if (x >= 9223372036854775808.0) return INT64_MAX;
    if (x <= -9223372036854775808.0) return INT64_MIN;
    return (int64_t)x;
}

/* Math.round Javadoc: the closest integer, ties towards positive infinity.
 * x - floor(x) is exact, so the tie test needs no x + 0.5, which can round.
 * Doubles >= 2^52 are integers with remainder 0; the casts saturate. */
static double round_half_up(double x) {
    double f = floor(x);
    return x - f >= 0.5 ? f + 1.0 : f;
}

int32_t jnum_round_f(float x) {
    return jnum_d2i(round_half_up((double)x));
}

int64_t jnum_round_d(double x) {
    return jnum_d2l(round_half_up(x));
}

/* DecimalFormat has its own nonfinite strings. For finite values, printf
 * rounds the exact binary value, ties to even, as DecimalFormat does. */
void jnum_format_fixed(char *buf, size_t size, double x, int digits) {
    if (!isfinite(x)) {
        const char *s = isnan(x) ? "NaN" : x < 0 ? "-\xe2\x88\x9e" : "\xe2\x88\x9e";
        snprintf(buf, size, "%s", s);
    } else {
        snprintf(buf, size, "%.*f", digits, x);
    }
}

void jnum_format_hash2(char *buf, size_t size, double x) {
    jnum_format_fixed(buf, size, x, 2);
    char *end = buf + strlen(buf);
    while (end[-1] == '0') *--end = '\0';
    if (end[-1] == '.') end[-1] = '\0';
}

void jnum_format_sci1(char *buf, size_t size, double x) {
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
        return;
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
}

/* Checks Java's FloatingDecimal grammar and copies the number, without its
 * suffix, into a NUL-terminated buffer for strtod or strtof, which round the
 * same way. Returns 1 for a number, 2 for NaN, 3 for +Infinity, 4 for
 * -Infinity, and 0 for invalid text. */
static int java_number_text(const char *s, size_t len, char **text) {
    size_t start = 0, end = len;
    while (start < end && (unsigned char)s[start] <= ' ') ++start;
    while (end > start && (unsigned char)s[end - 1] <= ' ') --end;
    const char *p = s + start, *e = s + end;
    bool neg = p < e && *p == '-';
    if (p < e && (*p == '+' || *p == '-')) ++p;
    if ((size_t)(e - p) == 3 && memcmp(p, "NaN", 3) == 0) return 2;
    if ((size_t)(e - p) == 8 && memcmp(p, "Infinity", 8) == 0) return neg ? 4 : 3;
    const char *q = p;
    if (e - q >= 2 && q[0] == '0' && (q[1] == 'x' || q[1] == 'X')) {
        q += 2;
        int digits = 0;
        while (q < e && isxdigit((unsigned char)*q)) ++q, ++digits;
        if (q < e && *q == '.') {
            ++q;
            while (q < e && isxdigit((unsigned char)*q)) ++q, ++digits;
        }
        if (digits == 0 || q == e || (*q != 'p' && *q != 'P')) return 0;
        ++q;
        if (q < e && (*q == '+' || *q == '-')) ++q;
        if (q == e || !isdigit((unsigned char)*q)) return 0;
        while (q < e && isdigit((unsigned char)*q)) ++q;
    } else {
        int digits = 0;
        while (q < e && isdigit((unsigned char)*q)) ++q, ++digits;
        if (q < e && *q == '.') {
            ++q;
            while (q < e && isdigit((unsigned char)*q)) ++q, ++digits;
        }
        if (digits == 0) return 0;
        if (q < e && (*q == 'e' || *q == 'E')) {
            ++q;
            if (q < e && (*q == '+' || *q == '-')) ++q;
            if (q == e || !isdigit((unsigned char)*q)) return 0;
            while (q < e && isdigit((unsigned char)*q)) ++q;
        }
    }
    const char *number_end = q;
    if (q < e && strchr("fFdD", *q) != NULL) ++q;
    if (q != e) return 0;
    size_t n = (size_t)(number_end - (s + start));
    *text = malloc(n + 1);
    if (*text == NULL) return 0;
    memcpy(*text, s + start, n);
    (*text)[n] = '\0';
    return 1;
}

bool jnum_parse_double(const char *s, size_t len, double *out) {
    char *text = NULL;
    switch (java_number_text(s, len, &text)) {
        case 1: *out = strtod(text, NULL); free(text); return true;
        case 2: *out = NAN; return true;
        case 3: *out = INFINITY; return true;
        case 4: *out = -INFINITY; return true;
        default: return false;
    }
}

bool jnum_parse_float(const char *s, size_t len, float *out) {
    char *text = NULL;
    switch (java_number_text(s, len, &text)) {
        case 1: *out = strtof(text, NULL); free(text); return true;
        case 2: *out = NAN; return true;
        case 3: *out = INFINITY; return true;
        case 4: *out = -INFINITY; return true;
        default: return false;
    }
}

/* Every BMP code point c with Character.digit(c, 10) == 0 in JDK 21; each is
 * the start of a run of ten digits. */
static const uint16_t DECIMAL_ZEROS[] = {
    0x0030, 0x0660, 0x06f0, 0x07c0, 0x0966, 0x09e6, 0x0a66, 0x0ae6, 0x0b66, 0x0be6,
    0x0c66, 0x0ce6, 0x0d66, 0x0de6, 0x0e50, 0x0ed0, 0x0f20, 0x1040, 0x1090, 0x17e0,
    0x1810, 0x1946, 0x19d0, 0x1a80, 0x1a90, 0x1b50, 0x1bb0, 0x1c40, 0x1c50, 0xa620,
    0xa8d0, 0xa900, 0xa9d0, 0xa9f0, 0xaa50, 0xabf0, 0xff10,
};

static int decimal_digit(int32_t c) {
    size_t n = sizeof DECIMAL_ZEROS / sizeof DECIMAL_ZEROS[0];
    for (size_t i = 0; i < n; ++i) {
        if (c < DECIMAL_ZEROS[i]) break;
        if (c - DECIMAL_ZEROS[i] < 10) return (int)(c - DECIMAL_ZEROS[i]);
    }
    return -1;
}

bool jnum_parse_long(const char *s, size_t len, int64_t *out) {
    if (len == 0) return false;
    size_t i = 0;
    int32_t first = jutf8_next_bmp(s, len, &i);
    bool negative = false;
    int64_t limit = -INT64_MAX;
    if (first < '0') {
        if (first == '-') {
            negative = true;
            limit = INT64_MIN;
        } else if (first != '+') {
            return false;
        }
        if (i == len) return false;
    } else {
        i = 0;
    }
    int64_t multmin = limit / 10;
    int64_t result = 0;
    while (i < len) {
        int digit = decimal_digit(jutf8_next_bmp(s, len, &i));
        if (digit < 0 || result < multmin) return false;
        result *= 10;
        if (result < limit + digit) return false;
        result -= digit;
    }
    *out = negative ? result : -result;
    return true;
}

bool jnum_parse_int(const char *s, size_t len, int32_t *out) {
    int64_t v;
    if (!jnum_parse_long(s, len, &v) || v < INT32_MIN || v > INT32_MAX) return false;
    *out = (int32_t)v;
    return true;
}

/* Whether the decimal in s reads back as x, in the precision to_string prints. */
static bool round_trips(const char *s, double x, bool as_float) {
    return as_float ? (double)strtof(s, NULL) == x : strtod(s, NULL) == x;
}

/* The p-digit decimal nearest x, as digits[0..p) and the exponent of its first
 * digit, if it rounds back to x. Otherwise the p-digit decimal next to it on
 * the other side of x, which is the nearest on that side, if that rounds back
 * to x (the rounding interval is wider above a power of two). */
static bool shortest_candidate(double x, bool as_float, int p, char *digits, int *exp10) {
    char s[40];
    snprintf(s, sizeof s, "%.*e", p - 1, x);
    int e = atoi(strchr(s, 'e') + 1);
    int n = 0;
    for (const char *c = s; *c != 'e'; ++c) {
        if (isdigit((unsigned char)*c)) digits[n++] = *c;
    }
    if (!round_trips(s, x, as_float)) {
        double nearest = strtod(s, NULL);
        int step = nearest < x ? 1 : -1;
        int j = p - 1;
        if (step > 0) {
            while (j >= 0 && digits[j] == '9') digits[j--] = '0';
            if (j < 0) {
                digits[0] = '1';
                ++e;
            } else {
                ++digits[j];
            }
        } else {
            while (j >= 0 && digits[j] == '0') digits[j--] = '9';
            --digits[j];
            if (digits[0] == '0') {
                memmove(digits, digits + 1, (size_t)p - 1);
                digits[p - 1] = '9';
                --e;
            }
        }
        snprintf(s, sizeof s, "%.*se%d", p, digits, e - p + 1);
        if (!round_trips(s, x, as_float)) return false;
    }
    *exp10 = e;
    return true;
}

/* Double.toString and Float.toString share one algorithm (JDK 19 and later):
 * the shortest decimal that reads back as x, at least two digits, written
 * plainly for 1e-3 <= |x| < 1e7 and as d.dddE<n> otherwise. */
static void to_string(char *buf, double x, bool as_float) {
    if (isnan(x)) {
        strcpy(buf, "NaN");
        return;
    }
    if (signbit(x)) *buf++ = '-';
    x = fabs(x);
    if (isinf(x)) {
        strcpy(buf, "Infinity");
        return;
    }
    if (x == 0.0) {
        strcpy(buf, "0.0");
        return;
    }
    char digits[20];
    int p = 1, e;
    while (!shortest_candidate(x, as_float, p, digits, &e)) ++p;
    if (p == 1) shortest_candidate(x, as_float, p = 2, digits, &e);
    while (p > 1 && digits[p - 1] == '0') --p;
    if (x >= 1e-3 && x < 1e7) {
        int n = 0;
        if (e < 0) {
            buf[n++] = '0';
            buf[n++] = '.';
            for (int j = -1; j > e; --j) buf[n++] = '0';
            memcpy(buf + n, digits, (size_t)p);
            n += p;
        } else {
            for (int j = 0; j <= e; ++j) buf[n++] = j < p ? digits[j] : '0';
            buf[n++] = '.';
            if (p > e + 1) {
                memcpy(buf + n, digits + e + 1, (size_t)(p - e - 1));
                n += p - e - 1;
            } else {
                buf[n++] = '0';
            }
        }
        buf[n] = '\0';
    } else {
        snprintf(buf, JNUM_DOUBLE_STRING_SIZE - 1, "%c.%.*sE%d", digits[0], p > 1 ? p - 1 : 1, p > 1 ? digits + 1 : "0", e);
    }
}

void jnum_double_to_string(char *buf, double x) {
    to_string(buf, x, false);
}

void jnum_float_to_string(char *buf, float x) {
    to_string(buf, x, true);
}
