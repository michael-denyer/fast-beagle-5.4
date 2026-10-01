"""Turn a 1000 Genomes high-coverage panel into an unphased phasing input.

Usage: make_phase.py <panel.vcf.gz> <out.vcf.gz>

Keeps every sample and every record except symbolic structural variants and
repeats of an earlier (POS, REF, ALT), which Beagle rejects, and removes the
phase from each genotype.
"""

import gzip
import sys


def main():
    src, out = sys.argv[1], sys.argv[2]
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
            if cols[4].startswith("<") or key in seen:
                continue
            seen.add(key)
            dst.write("\t".join(cols[:7] + [".", "GT"] + [g.replace("|", "/") for g in cols[9:]]) + "\n")


if __name__ == "__main__":
    main()
