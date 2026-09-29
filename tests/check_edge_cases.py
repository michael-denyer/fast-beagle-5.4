#!/usr/bin/env python3
"""Check output collisions, empty windows, and non-finite imputation output."""

import gzip
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BEAGLE = Path(os.environ.get("BEAGLE", ROOT / "build/beagle")).resolve()
SUFFIXES = (".vcf.gz", ".log", ".bgen", ".sample", ".info", ".vcf.gz.tbi")


def vcf(samples, rows):
    header = (
        '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    )
    return header + "".join(f"20\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n" for pos, gts in rows)


TARGET = vcf(["S1", "S2"], [(1000, ["0|0", "1|1"]), (2000, ["0|0", "1|1"]), (3000, ["0|0", "1|1"])])


class EdgeCases(unittest.TestCase):
    def run_beagle(self, folder, *args):
        result = subprocess.run(
            [BEAGLE, *args, "nthreads=1"], cwd=folder, capture_output=True, text=True, timeout=20, check=False
        )
        self.assertNotIn("Sanitizer", result.stderr)
        self.assertNotIn("runtime error:", result.stderr)
        return result

    def refused_unchanged(self, folder, *args):
        before = {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}
        result = self.run_beagle(folder, *args)
        after = {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}
        self.assertEqual(after, before, "a rejected output collision changed files")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("output file equals input file", result.stderr)

    def test_path_aliases(self):
        for alias in ("dot", "parent", "absolute", "symlink", "hardlink", "directory-symlink"):
            with self.subTest(alias=alias), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                (folder / "input.vcf.gz").write_bytes(gzip.compress(TARGET.encode()))
                if alias == "dot":
                    out = "./input"
                elif alias == "parent":
                    (folder / "sub").mkdir()
                    out = "sub/../input"
                elif alias == "absolute":
                    out = str(folder / "input")
                elif alias == "directory-symlink":
                    (folder / "link").symlink_to(folder, target_is_directory=True)
                    out = "link/input"
                else:
                    output = folder / "result.vcf.gz"
                    if alias == "symlink":
                        output.symlink_to("input.vcf.gz")
                    else:
                        os.link(folder / "input.vcf.gz", output)
                    out = "result"
                self.refused_unchanged(folder, "gt=input.vcf.gz", f"out={out}")

    def test_output_suffixes_and_inputs(self):
        for suffix in SUFFIXES:
            for key in ("gt", "ref", "map", "excludesamples", "excludemarkers", "ped", "truth"):
                with self.subTest(suffix=suffix, key=key), tempfile.TemporaryDirectory() as tmp:
                    folder = Path(tmp)
                    (folder / "target.vcf").write_text(TARGET)
                    (folder / ("input" + suffix)).write_text(TARGET)
                    args = [f"{key}=input{suffix}", "out=./input", "bgen=phased", "tbi=true"]
                    if key != "gt":
                        args.append("gt=target.vcf")
                    self.refused_unchanged(folder, *args)

    def test_unused_outputs_are_allowed(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "input.bgen").write_text(TARGET)
            result = self.run_beagle(folder, "gt=input.bgen", "out=input")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((folder / "input.bgen").read_text(), TARGET)

    def test_empty_window(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "target.vcf").write_text(TARGET)
            result = self.run_beagle(folder, "gt=target.vcf", "out=result", "window=0.0000001", "overlap=0.00000001")
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertIn("java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0", result.stderr)

    def test_nan_imputation(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "target.vcf").write_text(vcf(["S1"], [(1000, ["1|1"]), (3000, ["1|1"])]))
            (folder / "ref.vcf").write_text(
                vcf(["R1", "R2"], [(1000, ["0|0", "0|0"]), (2000, ["0|1", "1|0"]), (3000, ["0|0", "0|0"])])
            )
            result = self.run_beagle(folder, "gt=target.vcf", "ref=ref.vcf", "out=result", "err=0")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(folder / "result.vcf.gz", "rt") as stream:
                record = next(line.rstrip() for line in stream if line.startswith("20\t2000\t"))
            self.assertEqual(record, "20\t2000\t.\tA\tC\t.\tPASS\tDR2=0.00;AF=NaN;IMP\tGT:DS\t0|0:0")


if __name__ == "__main__":
    unittest.main()
