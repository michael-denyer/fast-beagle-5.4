/* Reads the output of JcompatFixtures.numbers() and prints it again, recomputing
 * every result with jnum. */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "jcompat/jnum.h"

int main(void) {
    char op[32], in[32], line[256], buf[64];
    while (fgets(line, sizeof line, stdin)) {
        if (sscanf(line, "%31s %31s", op, in) != 2) {
            fprintf(stderr, "bad line: %s", line);
            return 1;
        }
        uint64_t bits = strtoull(in, NULL, 16);
        double x;
        memcpy(&x, &bits, sizeof x);
        uint32_t fbits = (uint32_t)bits;
        float fx;
        memcpy(&fx, &fbits, sizeof fx);
        printf("%s %s ", op, in);
        if (strcmp(op, "d2i") == 0) printf("%" PRId32 "\n", jnum_d2i(x));
        else if (strcmp(op, "d2l") == 0) printf("%" PRId64 "\n", jnum_d2l(x));
        else if (strcmp(op, "roundf") == 0) printf("%" PRId32 "\n", jnum_round_f(fx));
        else if (strcmp(op, "roundd") == 0) printf("%" PRId64 "\n", jnum_round_d(x));
        else {
            if (strcmp(op, "hash2") == 0) jnum_format_hash2(buf, sizeof buf, x);
            else if (strcmp(op, "fixed2") == 0) jnum_format_fixed(buf, sizeof buf, x, 2);
            else if (strcmp(op, "fixed4") == 0) jnum_format_fixed(buf, sizeof buf, x, 4);
            else if (strcmp(op, "sci1") == 0) jnum_format_sci1(buf, sizeof buf, x);
            else if (strcmp(op, "tostring") == 0) jnum_double_to_string(buf, x);
            else if (strcmp(op, "ftostring") == 0) jnum_float_to_string(buf, fx);
            else {
                fprintf(stderr, "unknown op: %s\n", op);
                return 1;
            }
            printf("%s\n", buf);
        }
    }
    return 0;
}
