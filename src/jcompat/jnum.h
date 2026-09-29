/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef JCOMPAT_JNUM_H
#define JCOMPAT_JNUM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Java's narrowing casts: NaN becomes 0 and out-of-range values saturate,
 * where the C casts are undefined. */
int32_t jnum_d2i(double x);
int64_t jnum_d2l(double x);

/* Float.floatToRawIntBits */
static inline uint32_t jnum_float_bits(float f) {
    uint32_t bits;
    memcpy(&bits, &f, sizeof bits);
    return bits;
}

/* Math.round: the closest integer, with ties rounded towards positive infinity. */
int32_t jnum_round_f(float x);
int64_t jnum_round_d(double x);

/* java.text.DecimalFormat patterns "0.00" and "0.0000" (jnum_format_fixed) and
 * "#.##" (jnum_format_hash2), with its default HALF_EVEN rounding. Exact for
 * Beagle's inputs: floats in [0, 1] widened to double, and j/100.0. There
 * DecimalFormat equals half-even rounding of the exact binary value, which is
 * what printf does (checked for all 1,065,353,217 such floats). For general
 * doubles DecimalFormat rounds a shortest-digits string instead and can differ
 * in the last place, so these functions are not a general DecimalFormat.
 * Both write at most `size` bytes including the terminator. */
void jnum_format_fixed(char *buf, size_t size, double x, int digits);
void jnum_format_hash2(char *buf, size_t size, double x);

/* String.format("%1$7.1e", x) in Locale.US for finite x >= 0. Java's Formatter
 * rounds the shortest decimal digits of x half up, where printf rounds the
 * exact binary value half to even. */
void jnum_format_sci1(char *buf, size_t size, double x);

/* Double.parseDouble and Float.parseFloat: surrounding characters <= ' ' are
 * ignored; the text is an optional sign and then "NaN", "Infinity", a decimal
 * number with optional exponent, or a hex number with a required binary
 * exponent, and a decimal or hex number may end in f, F, d or D. Returns false
 * where Java throws NumberFormatException. */
bool jnum_parse_double(const char *s, size_t len, double *out);
bool jnum_parse_float(const char *s, size_t len, float *out);

/* Integer.parseInt on the UTF-8 decoding of s[0..len): an optional sign, then
 * the decimal digits of any script in the BMP. Returns false where Java throws
 * NumberFormatException, including on malformed UTF-8, which Java decodes to
 * U+FFFD. */
bool jnum_parse_int(const char *s, size_t len, int32_t *out);

/* Long.parseLong, with the same rules as jnum_parse_int. */
bool jnum_parse_long(const char *s, size_t len, int64_t *out);

/* Double.toString as JDK 19 and later print it: the shortest decimal that
 * rounds to x (or, if that has one digit, the closest one with at most two),
 * written plainly for 1e-3 <= |x| < 1e7 and as d.dddE<n> otherwise. buf needs
 * JNUM_DOUBLE_STRING_SIZE bytes. */
#define JNUM_DOUBLE_STRING_SIZE 32
void jnum_double_to_string(char *buf, double x);

/* Float.toString by the same rules, deciding shortness in float precision. */
void jnum_float_to_string(char *buf, float x);

#endif
