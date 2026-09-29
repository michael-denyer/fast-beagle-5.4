/*
 * bgen_quantise: each probability times max rounded down, then one more unit
 * to the largest remainders, lower alleles first on ties; the last allele's
 * value is left for the reader to infer.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "bgen/bgen_writer.h"

static int failures;

static void check(const char *name, const float *p, int n, uint32_t max, const uint32_t *want) {
    uint32_t got[8];
    bgen_quantise(p, n, max, got);
    for (int a = 0; a + 1 < n; ++a) {
        if (got[a] != want[a]) {
            printf("FAIL %s: allele %d is %u, want %u\n", name, a, got[a], want[a]);
            ++failures;
            return;
        }
    }
    printf("pass %s\n", name);
}

int main(void) {
    check("certain first allele", (const float[]){1.0f, 0.0f}, 2, 255, (const uint32_t[]){255});
    check("certain last allele", (const float[]){0.0f, 1.0f}, 2, 255, (const uint32_t[]){0});
    /* 127.5 and 127.5: one unit left, the tie goes to allele 0 */
    check("even split", (const float[]){0.5f, 0.5f}, 2, 255, (const uint32_t[]){128});
    /* As floats, 25.50000038, 76.50000304 and 153.00000608: floors sum to 254
     * and the unit goes to allele 1, whose remainder is largest */
    check("three alleles", (const float[]){0.1f, 0.3f, 0.6f}, 3, 255, (const uint32_t[]){25, 77});
    /* 85.0000025 each: floors already sum to 255 */
    check("thirds", (const float[]){1.0f / 3, 1.0f / 3, 1.0f / 3}, 3, 255, (const uint32_t[]){85, 85});
    /* 1.02 and 253.98: floors sum to 254, the unit goes to allele 1 */
    check("larger remainder wins", (const float[]){0.004f, 0.996f}, 2, 255, (const uint32_t[]){1});
    /* 4 bits: 10.4999998 and 4.5000002, the unit goes to allele 1 */
    check("max 15", (const float[]){0.7f, 0.3f}, 2, 15, (const uint32_t[]){10});
    /* 16 bits: 45874.4992 and 19660.5008, the unit goes to allele 1 */
    check("max 65535", (const float[]){0.7f, 0.3f}, 2, 65535, (const uint32_t[]){45874});
    /* 1 bit: 0.6 and 0.4 floor to 0, the unit goes to the larger */
    check("max 1", (const float[]){0.6f, 0.4f}, 2, 1, (const uint32_t[]){1});
    check("max 1 minor allele", (const float[]){0.4f, 0.6f}, 2, 1, (const uint32_t[]){0});
    /* NaN counts as 0 and takes no unit, so every unit goes to allele 0 */
    check("nan", (const float[]){NAN, NAN}, 2, 255, (const uint32_t[]){255});
    return failures != 0;
}
