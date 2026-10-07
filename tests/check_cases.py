#!/usr/bin/env python3
"""Exercise preparation and expected exits through the actual case runners."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run(root, command, *, env=None, ok=True):
    result = subprocess.run(command, cwd=root, env=env, text=True, capture_output=True, check=False)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "tests").mkdir()
        for name in (
            "cases.sh",
            "fetch-fixtures.sh",
            "oracle-cases.txt",
            "trace-cases.txt",
            "bgen-cases.txt",
            "run-trace.sh",
            "check-trace.sh",
            "check-oracle.sh",
        ):
            shutil.copy2(ROOT / "tests" / name, root / "tests" / name)
        shutil.copytree(ROOT / "data", root / "data")
        manifest = root / "data/.fixtures.sha256"
        manifest.unlink(missing_ok=True)
        prepare = ["bash", "-c", "ROOT=$PWD; source tests/cases.sh"]

        # A cache made by an older recipe has the old sentinel but no new inputs.
        missing = root / "data/test.common.vcf.gz"
        missing.unlink()
        run(root, prepare)
        assert missing.is_file() and manifest.is_file()
        print("PASS old fixture cache regenerated")

        stamp = manifest.stat().st_mtime_ns
        run(root, prepare)
        assert manifest.stat().st_mtime_ns == stamp
        print("PASS complete fixture cache reused")

        missing.unlink()
        run(root, prepare)
        assert missing.is_file()
        print("PASS deleted fixture restored despite existing manifest")

        recipe = root / "tests/fetch-fixtures.sh"
        with recipe.open("a") as f:
            f.write("\n# Changed recipe for the cache regression.\n")
        run(root, prepare)
        assert manifest.stat().st_mtime_ns != stamp
        print("PASS changed recipe regenerated")

        fake_bin = root / "bin"
        fake_bin.mkdir()
        java = fake_bin / "java"
        java.write_text("#!/bin/bash\nexit 73\n")
        java.chmod(0o755)
        env = dict(os.environ, PATH=f"{fake_bin}:{os.environ['PATH']}")
        missing.unlink()
        failed = run(root, prepare, env=env, ok=False).stderr
        assert not manifest.exists()
        assert "FAIL java cannot run" in failed and "FAIL fixtures" in failed, failed
        print("PASS failed preparation stops runner and invalidates manifest")
        tla = str(ROOT / "tests/check-tla.sh")
        assert (
            run(root, [tla], env=env, ok=False).stderr
            == "FAIL java cannot run, and TLC needs it: put a JDK first on PATH\n"
        )
        assert run(root, [tla, "BlockReader Nope"], env=env, ok=False).stdout == "FAIL unknown spec BlockReader Nope\n"
        print("PASS model check stops once without java or for an unknown spec")
        run(root, prepare)
        print("PASS interrupted fixture generation recovers")

        # Both implementations produce the same trace and the same wrong exit.
        # The comparison must still fail; an intentional refusal must pass.
        fake = """#!/bin/bash
for arg in "$@"; do
  case "$arg" in *trace=*) dir=${arg#*=} ;; esac
done
printf 'same trace\\n' > "$dir/T1a.txt"
exit 1
"""
        java.write_text(fake)
        (root / "build").mkdir()
        binary = root / "build/beagle"
        binary.write_text(fake)
        binary.chmod(0o755)
        env["CASES"] = "gt"
        run(root, ["tests/check-trace.sh", "T1a"], env=env, ok=False)
        print("PASS equal unexpected failures rejected")
        env["CASES"] = "edge-markers"
        run(root, ["tests/check-trace.sh", "T1a"], env=env)
        print("PASS intentional refusal accepted")

        # A Beagle that writes the same VCF on every run and exits FAKE_EXIT.
        writes_vcf = root / "writes-vcf"
        writes_vcf.write_text("""#!/bin/bash
for arg in "$@"; do
  case "$arg" in out=*) out=${arg#*=} ;; esac
done
printf '##fileformat=VCFv4.2\\n#CHROM\\n' | gzip > "$out.vcf.gz"
exit "${FAKE_EXIT:-0}"
""")
        writes_vcf.chmod(0o755)
        verdicts = run(
            root,
            [
                "bash",
                "-c",
                f"""ROOT=$PWD; source tests/cases.sh
judge() {{  # label expect nthreads
  case_verdict "$2" gt=@test.vcf.gz "$PWD/build/$1" "$3" {writes_vcf}
  echo "$1 $? $VERDICT"
}}
{writes_vcf} out="$PWD/build/h"
h=$(vcf_hash build/h.vcf.gz)
judge threads-1 "1:$h,2:$h,18:0123456789abcdef" 1
judge threads-18 "1:$h,2:$h,18:0123456789abcdef" 18
judge exit-0 "$h" 2
FAKE_EXIT=3 judge exit-3 "$h" 2
FAKE_EXIT=3 judge refusal "exit=3" 2
""",
            ],
        ).stdout
        rc = {line.split()[0]: line.split()[1] for line in verdicts.splitlines()}
        assert rc == {"threads-1": "0", "threads-18": "1", "exit-0": "0", "exit-3": "1", "refusal": "0"}, verdicts
        assert "want exit=0 0123456789abcdef" in verdicts, verdicts
        print("PASS verdict rejects a wrong hash at 18 threads and an unexpected exit")

        selection = ["bash", "-c", "ROOT=$PWD; source tests/cases.sh; check_selection tests/oracle-cases.txt"]
        run(root, selection)
        run(root, selection, env=dict(os.environ, CASES="gt"))
        typo = run(root, selection, env=dict(os.environ, CASES="gt no-such-case"), ok=False).stdout
        assert "FAIL no-such-case is not a case" in typo and "FAIL gt " not in typo, typo
        run(root, ["tests/check-oracle.sh", "true"], env=dict(os.environ, CASES="no-such-case"), ok=False)
        print("PASS a CASES name missing from the tables fails the runner")

        # A Beagle that prints the refusal, leaves the files named by LEAVE and exits FAKE_EXIT.
        refuses = root / "refuses"
        refuses.write_text("""#!/bin/bash
echo "ERROR: refused" >&2
for f in ${LEAVE:-}; do : > "$f"; done
exit "${FAKE_EXIT:-1}"
""")
        refuses.chmod(0o755)
        verdicts = run(
            root,
            [
                "bash",
                "-c",
                f"""ROOT=$PWD; source tests/cases.sh
judge() {{  # label
  rm -f build/r.bgen
  {refuses} 2> build/r.log
  refused build/r.log $? "ERROR: refused" build/r.bgen
  echo "$1 $? $VERDICT"
}}
FAKE_EXIT=2 judge crash-after-message
LEAVE=build/r.bgen judge leftover-bgen
judge clean
""",
            ],
        ).stdout
        rc = {line.split()[0]: line.split()[1] for line in verdicts.splitlines()}
        assert rc == {"crash-after-message": "1", "leftover-bgen": "1", "clean": "0"}, verdicts
        assert "exit=2, want 1" in verdicts and "left r.bgen" in verdicts, verdicts
        print("PASS refusal rejects a wrong exit and a leftover file")

        # check-sanitizers.sh with its builds replaced by the fake Beagle.
        fake_make = root / "make-bin/make"
        fake_make.parent.mkdir()
        fake_make.write_text(f'#!/bin/bash\nmkdir -p "$2/build"\ncp {writes_vcf} "$2/build/beagle"\n')
        fake_make.chmod(0o755)
        shutil.copy2(ROOT / "tests/check-sanitizers.sh", root / "tests")
        env = dict(os.environ, PATH=f"{fake_make.parent}:{os.environ['PATH']}", CASES="gt gt-x", NTHREADS="1")
        runs = run(root, ["tests/check-sanitizers.sh"], env=env, ok=False).stdout
        modes = {tuple(line.split()[1:4]) for line in runs.splitlines() if line.startswith(("PASS", "FAIL"))}
        assert ("gt", "nthreads=1", "plink2") in modes and ("gt-x", "nthreads=1", "phased") in modes, runs
        assert ("gt-x", "nthreads=1", "plink2") not in modes, runs
        print("PASS nonautosome case skipped under bgen=plink2 by the sanitizer runs")


if __name__ == "__main__":
    main()
