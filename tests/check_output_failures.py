#!/usr/bin/env python3
"""Check output/input collisions and numeric output boundaries."""

import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BEAGLE = Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else ROOT / "build/beagle"
PROGRAM = "fast-beagle"
OPTIONAL_SUFFIXES = (".bgen", ".info", ".sample", ".vcf.gz.tbi")
SUFFIXES = (".vcf.gz", ".log", *OPTIONAL_SUFFIXES)
VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
    "1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0|1\t1|0\t0|0\n"
    "1\t200\tv2\tA\tC\t.\t.\t.\tGT\t0|1\t1|0\t0|0\n"
)
INPUTS = {
    "gt": VCF,
    "ref": VCF,
    "map": "1 v1 0 100\n1 v2 1 200\n",
    "excludesamples": "S3\n",
    "excludemarkers": "absent-marker\n",
    "truth": VCF,
    "ped": "",
}


def run(args):
    return subprocess.run(
        [str(BEAGLE), "nthreads=2", "seed=-99999", *[f"{k}={v}" for k, v in args.items()]],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


class OutputFailures(unittest.TestCase):
    def collision(self, field, suffix, mode, alias="direct"):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF)
            prefix = root / "result"
            dest = Path(str(prefix) + suffix)
            contents = INPUTS[field].encode()
            if suffix == ".vcf.gz":
                contents = gzip.compress(contents)
            source = dest
            if alias in ("symlink", "hardlink"):
                source = root / "source"
            source.write_bytes(contents)
            if alias == "symlink":
                dest.symlink_to(source)
            elif alias == "hardlink":
                os.link(source, dest)
            elif alias == "relative":
                (root / "sub").mkdir()
                prefix = root / "sub" / ".." / "result"
            args = {"gt": target, "out": prefix, field: source, "tbi": "true"}
            if mode is not None:
                args["bgen"] = mode
            proc = run(args)
            self.assertEqual(proc.returncode, 1, proc.stderr)
            # Beagle refuses only a VCF path that equals the gt= or ref= path.
            if suffix == ".vcf.gz" and field in ("gt", "ref") and alias == "direct":
                self.assertIn(f"ERROR: VCF output file equals input file: {source}", proc.stderr)
            else:
                self.assertIn(f"{PROGRAM}: output file {prefix}{suffix} equals input file ", proc.stderr)
            self.assertEqual(source.read_bytes(), contents, "input was overwritten")
            if alias == "symlink":
                self.assertTrue(dest.is_symlink(), "input link was removed")
            for ext in SUFFIXES:
                if ext != suffix:
                    self.assertFalse(Path(str(prefix) + ext).exists(), f"unexpected {ext}")

    def test_input_collisions(self):
        for mode in ("phased", "plink2"):
            for field in INPUTS:
                for suffix in SUFFIXES:
                    with self.subTest(mode=mode, field=field, suffix=suffix):
                        self.collision(field, suffix, mode)

    def test_alias_collisions(self):
        for alias in ("relative", "symlink", "hardlink"):
            for suffix in SUFFIXES:
                with self.subTest(alias=alias, suffix=suffix):
                    self.collision("gt", suffix, "phased", alias)

    def test_log_collision_without_bgen(self):
        self.collision("gt", ".log", None)

    def test_disabled_outputs_do_not_collide(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF)
            for suffix in OPTIONAL_SUFFIXES:
                with self.subTest(suffix=suffix):
                    prefix = root / suffix[1:]
                    exclude = Path(str(prefix) + suffix)
                    exclude.write_text("S3\n")
                    proc = run({"gt": target, "out": prefix, "excludesamples": exclude})
                    self.assertEqual(proc.returncode, 0, proc.stderr)
                    self.assertEqual(exclude.read_text(), "S3\n")

    def test_large_population_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                "1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\n"
                "1\t20\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\n"
            )
            proc = run({"gt": target, "out": root / "result", "ne": "1e20", "burnin": 1, "iterations": 1})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("9223372036854775807", (root / "result.log").read_text())
            self.assertNotIn("runtime error:", proc.stderr)
            self.assertNotIn("Sanitizer", proc.stderr)

    def test_nonfinite_phased_probabilities(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "result"
            proc = run(
                {
                    "gt": ROOT / "data/target.thin.vcf.gz",
                    "ref": ROOT / "data/ref.vcf.gz",
                    "out": prefix,
                    "err": 0,
                    "bgen": "phased",
                }
            )
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn(f"{PROGRAM}: bgen=phased: cannot encode a non-finite allele probability", proc.stderr)
            self.assertNotIn("runtime error:", proc.stderr)
            self.assertNotIn("Sanitizer", proc.stderr)
            for suffix in (".bgen", ".info", ".sample"):
                self.assertFalse(Path(str(prefix) + suffix).exists(), f"partial {suffix}")


if __name__ == "__main__":
    unittest.main()
