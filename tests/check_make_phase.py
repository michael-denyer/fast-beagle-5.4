#!/usr/bin/env python3
"""Check the records and genotypes that tests/bench/make_phase.py keeps."""

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "bench/make_phase.py"
HEADER = "##fileformat=VCFv4.2\n##source=panel\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfemale\tmale\n"
# Males are diploid in the pseudoautosomal regions (100 and 900 here) and
# haploid between them.
RECORDS = (
    (100, "A", "C", "0|1", "1|0"),
    (200, "A", "C", "0|1", "1"),
    (300, "A", "<DEL>", "0|1", "0"),
    (400, "A", "C", "1|1", "0"),
    (400, "A", "C", "1|1", "0"),
    (400, "A", "G", "1|0", "1"),
    (500, "A", "C", "0|0", "1"),
    (900, "A", "C", "1|0", "0|1"),
)


def make_phase(*args):
    with tempfile.TemporaryDirectory() as tmp:
        panel, out = Path(tmp) / "panel.vcf.gz", Path(tmp) / "out.vcf.gz"
        with gzip.open(panel, "wt") as f:
            f.write(HEADER)
            for pos, ref, alt, female, male in RECORDS:
                f.write(f"chrX\t{pos}\tid\t{ref}\t{alt}\t.\tPASS\tAC=1\tGT\t{female}\t{male}\n")
        subprocess.run([sys.executable, str(SCRIPT), str(panel), str(out), *args], check=True, timeout=30)
        with gzip.open(out, "rt") as f:
            return [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]


class MakePhase(unittest.TestCase):
    def test_whole_panel(self):
        rows = make_phase()
        self.assertEqual(
            [(r[1], r[4]) for r in rows],
            [("100", "C"), ("200", "C"), ("400", "C"), ("400", "G"), ("500", "C"), ("900", "C")],
        )
        self.assertEqual(rows[0][7:], [".", "GT", "0/1", "1/0"])

    def test_range_is_inclusive(self):
        rows = make_phase("200-500")
        self.assertEqual([(r[1], r[4]) for r in rows], [("200", "C"), ("400", "C"), ("400", "G"), ("500", "C")])
        self.assertEqual({r[10] for r in rows}, {"0", "1"}, "a kept male genotype is diploid")

    def test_range_inside_the_panel(self):
        self.assertEqual([r[1] for r in make_phase("201-499")], ["400", "400"])


if __name__ == "__main__":
    unittest.main()
