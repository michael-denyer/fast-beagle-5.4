/* The SNV allele table matches MarkerUtils.snvPerms(). The test includes the
 * implementation to read the static table. */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "vcf/marker.c"

/* FNV-1a 64 of the table's rows, each followed by '\n'. Regenerate with
 *   python3 -c 'from itertools import permutations as P
 *   r=sorted([p[0]+"\t"+",".join(p[1:]) for p in P("*ACGT") if p[0]!="*"]+[b+"\t." for b in "ACGT"])
 *   h=0xcbf29ce484222325
 *   for b in "".join(s+"\n" for s in r).encode(): h=((h^b)*0x100000001b3)%2**64
 *   print(hex(h))' */
#define SNV_PERMS_FNV 0x2c5556d88a25060u

int main(void) {
    init_snv_perms();
    assert(n_snv_perms == N_SNV_PERMS);
    uint64_t h = 0xcbf29ce484222325u;
    for (int i = 0; i < N_SNV_PERMS; ++i) {
        for (const char *p = snv_perms[i]; *p != '\0'; ++p) h = (h ^ (unsigned char)*p) * 0x100000001b3u;
        h = (h ^ '\n') * 0x100000001b3u;
    }
    assert(strcmp(snv_perms[0], "A\t*,C,G,T") == 0);
    assert(strcmp(snv_perms[N_SNV_PERMS - 1], "T\tG,C,A,*") == 0);
    assert(h == SNV_PERMS_FNV);
    printf("snv_perms: %d rows match\n", N_SNV_PERMS);
    return 0;
}
