/* Marker validation errors must be catchable by the read-ahead reader. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"
#include "vcf/markers.h"

typedef struct {
    const marker *const *markers;
    int n;
} check_ctx;

static void check_markers(void *arg) {
    const check_ctx *c = arg;
    markers_check(c->markers, c->n);
}

static void check(const char *name, const char *const *lines, int n, const char *expected) {
    marker records[3];
    const marker *markers[3];
    for (int j = 0; j < n; ++j) {
        marker_parse(&records[j], lines[j], strlen(lines[j]));
        markers[j] = &records[j];
    }
    check_ctx c = {markers, n};
    char *error = util_try(check_markers, &c);
    if (expected == NULL) assert(error == NULL);
    else assert(error != NULL && strcmp(error, expected) == 0);
    free(error);
    for (int j = 0; j < n; ++j) marker_free(&records[j]);
    printf("PASS %s\n", name);
}

int main(void) {
    const char *ascending[] = {
        "1\t100\tfirst\tA\tC\t.\t.\t.\tGT",
        "1\t200\tmiddle\tA\tC\t.\t.\t.\tGT",
        "1\t300\tlast\tA\tC\t.\t.\t.\tGT"
    };
    const char *descending[] = {ascending[2], ascending[1], ascending[0]};
    const char *duplicate[] = {ascending[0], "1\t100\tdup%id\tA\tC\t.\t.\t.\tGT"};
    const char *disordered[] = {ascending[0], ascending[2], ascending[1]};
    check("empty markers", NULL, 0, NULL);
    check("ascending markers", ascending, 3, NULL);
    check("descending markers", descending, 3, NULL);
    check("duplicate marker is caught with its original message", duplicate, 2,
            "java.lang.IllegalArgumentException: Duplicate marker: 1\t100\tdup%id\tA\tC");
    check("marker order error is caught with its original message", disordered, 3,
            "java.lang.IllegalArgumentException: markers not in chromosomal order: \n1\t100\tfirst\tA\tC"
            "\n1\t300\tlast\tA\tC\n1\t200\tmiddle\tA\tC");
    return 0;
}
