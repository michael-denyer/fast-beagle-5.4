#!/usr/bin/env python3
"""Check completed records and captured probabilities through the writer."""

import gzip
import subprocess
import tempfile
from pathlib import Path

from check_bgen_phased import decode, read_bgen

ROOT = Path(__file__).resolve().parent.parent


def main():
    with tempfile.TemporaryDirectory() as tmp:
        outputs = []
        for mode in range(3):
            out = Path(tmp) / str(mode)
            subprocess.run([ROOT / "build/output/record_fixture", out, str(mode)], check=True)
            with gzip.open(f"{out}.vcf.gz", "rt") as f:
                outputs.append([s.rstrip("\n").split("\t") for s in f if s[0] != "#"])
        assert outputs[0] == outputs[1] == outputs[2]
        phased, imputed, observed, mono, five = outputs[0]
        assert phased[6:] == ["PASS", ".", "GT", "0|1", "1", "0|1"]
        assert imputed[6] == "PASS" and imputed[7].endswith(";IMP")
        assert imputed[8] == "GT:DS:AP1:AP2:GP"
        assert imputed[9].split(":")[:4] == ["0|1", "1.21", "0.33", "0.88"]
        assert imputed[10:] == ["0:0.25:0.25", "0|0:0:0:0:1,0,0"]
        assert observed[7].endswith(";END=301")
        assert observed[9:] == [
            "2|0:0,1:0,1:0,0:0,0,0,1,0,0",
            "1:1,0:1,0",
            "0|0:0,0:0,0:0,0:1,0,0,0,0,0",
        ]
        assert mono[7] == "IMP" and mono[9:] == ["0|0:1", "0", "0|0:1"]
        assert five[9].startswith("0|0:0,0,0,0:0,0,0,0:0,0,0,0:1,")
        print("PASS phased, imputed, observed, haploid and multiallelic records")

        ids, variants = read_bgen(Path(tmp) / "2.bgen")
        assert ids == ["D1", "H", "D2"] and len(variants) == 5
        top, probs = decode(variants[1][-1], 3, 2)
        assert top == 65535
        # At 16 bits, rounding to the VCF's two decimals before capture fails.
        for got, want in zip(probs[0], [0.3338, 0.8766], strict=True):
            assert abs(got[1] / top - want) < 1 / top
        assert abs(probs[1][0][1] / top - 0.25) < 1 / top
        assert probs[2] == [[top, 0], [top, 0]]
        _, probs = decode(variants[2][-1], 3, 3)
        assert probs == [
            [[0, 0, top], [top, 0, 0]],
            [[0, top, 0]],
            [[top, 0, 0], [top, 0, 0]],
        ]
        print("PASS BGEN captures probabilities before buffer reuse and VCF rounding")


if __name__ == "__main__":
    main()
