#!/usr/bin/env python3
"""Failed oracle runs must not replace the recorded logs, even after a good run.

A successful recording replaces tests/logs/ entirely, dropping removed cases.
"""

import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VCF = "##fileformat=VCFv4.2\n"
FAKE_JAVA = """#!/usr/bin/env python3
import gzip
import os
import sys
from pathlib import Path

out = Path(next(arg[4:] for arg in sys.argv[1:] if arg.startswith("out=")))
Path(str(out) + ".log").write_text("new " + out.name + "\\n")
vcf = "wrong\\n" if out.name == os.environ.get("WRONG_CASE") else "##fileformat=VCFv4.2\\n"
Path(str(out) + ".vcf.gz").write_bytes(gzip.compress(vcf.encode()))
sys.exit(73 if out.name == os.environ.get("FAIL_CASE") else 0)
"""


class LogRecording(unittest.TestCase):
    def check_recording(self, fail_case="", wrong_case=""):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tests = root / "tests"
            logs = tests / "logs"
            logs.mkdir(parents=True)
            for name in ("check-log.sh", "cases.sh"):
                shutil.copy2(ROOT / "tests" / name, tests / name)
            # These cases have no input fixtures; keep fixture downloads out of this harness.
            fetch = tests / "fetch-fixtures.sh"
            fetch.write_text("#!/bin/bash\nexit 0\n")
            fetch.chmod(0o755)
            expected_hash = hashlib.sha256(VCF.encode()).hexdigest()[:16]
            (tests / "oracle-cases.txt").write_text(
                f"first {expected_hash} - ignored=true\nsecond {expected_hash} - ignored=true\n"
            )
            for name in ("first", "second", "removed"):
                (logs / f"{name}.log").write_text(f"old {name}\n")
            java = root / "java"
            java.write_text(FAKE_JAVA)
            java.chmod(0o755)
            env = {
                **os.environ,
                "PATH": f"{root}:{os.environ['PATH']}",
                "RECORD": "1",
                "CASES": "",
                "FAIL_CASE": fail_case,
                "WRONG_CASE": wrong_case,
            }
            proc = subprocess.run(
                ["bash", str(tests / "check-log.sh"), "java"],
                env=env,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            failed = bool(fail_case or wrong_case)
            self.assertEqual(proc.returncode, int(failed), proc.stdout + proc.stderr)
            for name in ("first", "second"):
                want = f"{'old' if failed else 'new'} {name}\n"
                self.assertEqual((logs / f"{name}.log").read_text(), want)
            # A case removed from the table loses its log only when recording succeeds.
            self.assertEqual((logs / "removed.log").exists(), failed)

    def test_failed_run_preserves_recordings(self):
        for name in ("first", "second"):
            with self.subTest(case=name):
                self.check_recording(fail_case=name)

    def test_wrong_vcf_preserves_recordings(self):
        self.check_recording(wrong_case="second")

    def test_success_publishes_recordings(self):
        self.check_recording()


if __name__ == "__main__":
    unittest.main()
