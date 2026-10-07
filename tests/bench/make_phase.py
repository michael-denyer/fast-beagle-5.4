"""Turn a 1000 Genomes high-coverage panel into an unphased phasing input.

Usage: make_phase.py <panel.vcf.gz> <out.vcf.gz> [<first>-<last>]

Keeps every sample and every record except symbolic structural variants and
repeats of an earlier (POS, REF, ALT), which Beagle rejects, and removes the
phase from each genotype. With <first>-<last>, keeps only the records whose
POS is in that range. chrX needs 2781480-155701382, the region outside the
pseudoautosomal regions, where each sample has one ploidy.
"""

import gzip
import sys


def main():
    src, out = sys.argv[1], sys.argv[2]
    first, last = (int(x) for x in sys.argv[3].split("-")) if len(sys.argv) > 3 else (0, float("inf"))
    seen = set()
    with gzip.open(src, "rt") as f, gzip.open(out, "wt", compresslevel=1) as dst:
        for line in f:
            if line.startswith("##"):
                if line.startswith(("##fileformat", "##contig=", "##FORMAT=<ID=GT")):
                    dst.write(line)
                continue
            if line.startswith("#"):
                dst.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            key = (cols[1], cols[3], cols[4])
            if cols[4].startswith("<") or key in seen or not first <= int(cols[1]) <= last:
                continue
            seen.add(key)
            dst.write("\t".join(cols[:7] + [".", "GT"] + [g.replace("|", "/") for g in cols[9:]]) + "\n")


if __name__ == "__main__":
    main()
