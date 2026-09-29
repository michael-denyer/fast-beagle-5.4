# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["hypothesis==6.168.1"]
# [tool.uv]
# exclude-newer = "2026-09-26T00:00:00Z"
# ///
"""Differential fuzzing: the C binary against the Java release jar.

Usage: uv run --python 3.12 --script tests/check_fuzz.py [--examples N] [--invalid-examples N] [--random]

Hypothesis generates small inputs: a target VCF, and sometimes a reference
panel, a genetic map and options (windows, iterations, states, imputation,
gp/ap, bgen= for the C run). Haplotypes are mosaics of a few founders, so
the panels have IBS segments. Both tools run with the same seed. Each input
must either succeed in both with the same VCF (without ##source and
##filedate), or fail in both: Java with an exception or ERROR line, C with
exit code 1 and a message line (as every C error exit does). A C crash
or abort never counts as failing alike. Where the two messages differ the
input still passes, and the summary counts it under "messages differ".
A failure is shrunk to a small input, whose files are copied to
build/fuzz-fail with the two commands.

The invalid-parameter examples change one parameter of such an input so that
Main.parameters or Par rejects it (out naming an input or a directory, window
below 1.1 times overlap, a value outside its bounds or not a number, an
unknown parameter). Java and C must both exit 1 with the same message (Java
may prefix it with the exception class) and leave the input files unchanged.

First, each directory in tests/fuzz-regressions (inputs from earlier
failures, with their arguments in args.txt) must give the C result named in
expect.txt: an exit code and a message, or exit code 0 and the VCF hash.

Without --random the examples are the same on every run (the gate's mode).
BEAGLE overrides the C binary (default build/beagle), JAR the jar (default
data/beagle.29Oct24.c8e.jar) and JAVA the java command.
"""

import argparse
import gzip
import hashlib
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from hypothesis import HealthCheck, given, reproduce_failure, settings
from hypothesis import strategies as st

ROOT = Path(__file__).resolve().parent.parent
BEAGLE = os.environ.get("BEAGLE", str(ROOT / "build" / "beagle"))
JAR = os.environ.get("JAR", str(ROOT / "data" / "beagle.29Oct24.c8e.jar"))
JAVA = os.environ.get("JAVA", "java")
FAIL_DIR = ROOT / "build" / "fuzz-fail"
REGRESSIONS = ROOT / "tests" / "fuzz-regressions"
# bgen= must fail when no record is left for the BGEN (as plink2 does), after
# writing the VCF.
BGEN_EMPTY = "no variants remaining after the bgen filters"
BASES = "ACGT"
OUTCOMES = Counter()


@st.composite
def scenarios(draw):
    """Sizes, rates and options drawn by Hypothesis; genotypes come from a PRNG seeded by it."""
    has_ref = draw(st.booleans())
    chrom = draw(st.sampled_from(["20", "chr7", "X"]))
    s = {
        "rng_seed": draw(st.integers(0, 2**32 - 1)),
        "chrom": chrom,
        "has_ref": has_ref,
        "n_ref": draw(st.integers(2, 40)) if has_ref else 0,
        "n_targ": draw(st.integers(1, 30)),
        "n_markers": draw(st.integers(1, 400)),
        "span_mb": draw(st.sampled_from([0.01, 1, 5, 20])),
        "n_founders": draw(st.integers(1, 8)),
        "switch": draw(st.sampled_from([0.0, 0.01, 0.1])),
        "multi": draw(st.sampled_from([0.0, 0.05, 0.3])),
        "missing": draw(st.sampled_from([0.0, 0.01, 0.2])),
        "phased": draw(st.sampled_from([0.0, 0.5, 1.0])),
        "haploid": draw(st.sampled_from([0.0, 0.3])) if chrom == "X" else 0.0,
        "targ_every": draw(st.integers(1, 10)) if has_ref else 1,
        "extra_targ": draw(st.booleans()) if has_ref else False,
        "map": draw(st.booleans()),
        "opts": {},
    }
    opts = s["opts"]
    opts["seed"] = draw(st.integers(-99999, 99999))
    opts["nthreads"] = draw(st.sampled_from([1, 2, 3, 5]))
    if draw(st.booleans()):
        opts["burnin"] = draw(st.integers(1, 3))
        opts["iterations"] = draw(st.integers(1, 4))
    if draw(st.booleans()):
        opts["phase-states"] = draw(st.integers(1, 40))
    if draw(st.booleans()):
        window = draw(st.sampled_from([0.5, 1, 2, 5]))
        opts["window"] = window
        opts["overlap"] = draw(st.sampled_from([0.1, 0.3, window / 2]))
    if has_ref:
        if draw(st.booleans()):
            opts["imp-states"] = draw(st.integers(1, 60))
        if draw(st.booleans()):
            opts["impute"] = draw(st.sampled_from(["true", "false"]))
        if draw(st.booleans()):
            opts["cluster"] = draw(st.sampled_from([0.001, 0.005, 0.05]))
    if draw(st.booleans()):
        opts["gp"] = "true"
    if draw(st.booleans()):
        opts["ap"] = "true"
    # bgen= is C only; its VCF must still match Java's. plink2 mode rejects chrX.
    modes = ["none", "phased"] + ([] if chrom == "X" else ["plink2"])
    s["bgen"] = draw(st.sampled_from(modes))
    return s


# Parameters outside Par's bounds, one per parameter Par checks with a bound.
OUT_OF_BOUNDS = [
    ("window", "0"),
    ("overlap", "0"),
    ("buffer", "0"),
    ("burnin", "0"),
    ("iterations", "0"),
    ("phase-states", "0"),
    ("imp-states", "0"),
    ("imp-segment", "0"),
    ("imp-step", "0"),
    ("imp-nsteps", "0"),
    ("initial-lr", "0.5"),  # Beagle 5.5 only: Beagle 5.4 refuses it as unknown
    ("step-scale", "0"),
    ("rare", "0.6"),
    ("cluster", "-0.5"),
    ("ne", "0"),
    ("err", "-1"),
    ("nthreads", "0"),
    ("window-markers", "99999"),  # Beagle 5.5 only, as initial-lr
    ("seed", "1e3"),
    ("window", "abc"),
    ("nthreads", "2.5"),
]

INVALID = [
    ("out-equals-gt", "/"),
    ("out-equals-gt", "//"),
    ("out-directory",),
    ("unknown",),
    *[("window-overlap", o, o * f) for o in (0.1, 0.5, 1, 2, 5) for f in (0.25, 1.0, 1.05)],
    *[("out-of-bounds", key, value) for key, value in OUT_OF_BOUNDS],
    ("out-equals-ref", "/"),
    ("out-equals-ref", "//"),
]


def invalid_scenarios(invalid):
    """Scenarios with one of the INVALID changes."""
    base = scenarios().filter(lambda s: s["has_ref"]) if invalid[0] == "out-equals-ref" else scenarios()
    return base.map(lambda s: {**s, "invalid": invalid})


def mosaic(rng, founders, n_markers, switch):
    f = rng.randrange(len(founders))
    hap = []
    for m in range(n_markers):
        if rng.random() < switch:
            f = rng.randrange(len(founders))
        hap.append(founders[f][m])
    return hap


def write_vcf(path, chrom, positions, alleles, samples, genotypes):
    with gzip.open(path, "wt") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={chrom}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for m, pos in enumerate(positions):
            ref, *alt = alleles[m]
            fields = [chrom, str(pos), f"m{pos}", ref, ",".join(alt) if alt else ".", ".", "PASS", ".", "GT"]
            out.write("\t".join(fields + [g[m] for g in genotypes]) + "\n")


def build(s, d):
    """Writes the scenario's inputs into d and returns Beagle's arguments."""
    rng = random.Random(s["rng_seed"])
    n = s["n_markers"]
    span = max(n, int(s["span_mb"] * 1e6))
    positions = sorted(rng.sample(range(1, span + 1), n))
    alleles = []
    for _ in range(n):
        k = 3 if rng.random() < s["multi"] else 2
        alleles.append(rng.sample(BASES, k))
    founders = []
    for _ in range(s["n_founders"]):
        founders.append([0 if rng.random() < 0.7 else rng.randrange(1, len(alleles[m])) for m in range(n)])

    def hap():
        return mosaic(rng, founders, n, s["switch"])

    args = [f"gt={d}/targ.vcf.gz", f"out={d}/out"]
    if s["has_ref"]:
        ref_gt = [[f"{a}|{b}" for a, b in zip(hap(), hap())] for _ in range(s["n_ref"])]
        write_vcf(d / "ref.vcf.gz", s["chrom"], positions, alleles, [f"R{i}" for i in range(s["n_ref"])], ref_gt)
        args.append(f"ref={d}/ref.vcf.gz")
    targ_markers = list(range(0, n, s["targ_every"]))
    targ_pos = [positions[m] for m in targ_markers]
    targ_alleles = [alleles[m] for m in targ_markers]
    if s["extra_targ"]:  # a marker the reference lacks
        p = positions[-1] + 17
        targ_markers.append(-1)
        targ_pos.append(p)
        targ_alleles.append(["A", "G"])
    targ_gt = []
    for _ in range(s["n_targ"]):
        haploid = rng.random() < s["haploid"]
        h1, h2 = hap(), hap()
        g = []
        for m in targ_markers:
            a, b = (h1[m], h2[m]) if m >= 0 else (rng.randrange(2), rng.randrange(2))
            if rng.random() < s["missing"]:
                g.append("." if haploid else "./.")
            elif haploid:
                g.append(str(a))
            else:
                g.append(f"{a}{'|' if rng.random() < s['phased'] else '/'}{b}")
        targ_gt.append(g)
    write_vcf(d / "targ.vcf.gz", s["chrom"], targ_pos, targ_alleles, [f"T{i}" for i in range(s["n_targ"])], targ_gt)
    if s["map"]:
        with open(d / "plink.map", "w") as out:
            cm = 0.0
            for p in sorted(set(positions + targ_pos)):
                cm += rng.random() * 0.5
                out.write(f"{s['chrom']}\t.\t{cm:.6f}\t{p}\n")
        args.append(f"map={d}/plink.map")
    args += [f"{k}={v}" for k, v in s["opts"].items()]
    return args


def vcf_hash(path):
    try:
        with gzip.open(path, "rt") as f:
            lines = [ln for ln in f if not ln.startswith(("##source", "##filedate"))]
    except (OSError, EOFError):
        return None
    return hashlib.sha256("".join(lines).encode()).hexdigest()[:16]


TIMEOUT = 120  # seconds; every generated input runs in a few


# Java Beagle can keep running after its main thread throws, while another
# thread is still alive; that counts as exit code 1.
MAIN_THREW = 'Exception in thread "main"'
# Java reports a failure on an exception line or an "ERROR" line.
JAVA_ERROR = re.compile(r"^.*(?:\bERROR\b|Exception).*$", re.MULTILINE)
# Par also refuses an unknown parameter with an "Error:" line.
JAVA_REFUSAL = re.compile(r"^.*(?:\bERROR\b|^Error: |Exception).*$", re.MULTILINE)


def run(cmd, out):
    """Exit code (None on timeout), VCF hash and the output."""
    log = Path(f"{out}.run.log")
    with open(log, "w") as f:
        proc = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        start = time.monotonic()
        while proc.poll() is None:
            time.sleep(0.1)
            text = log.read_text(errors="replace")
            # The prefix and exception message are separate writes. Killing
            # after the prefix alone can lose the message under a slow JVM.
            main_exception = text.find(MAIN_THREW)
            if main_exception >= 0 and "\n" in text[main_exception:]:
                proc.kill()
                proc.wait()
                return 1, None, text[-2000:]
            if time.monotonic() - start > TIMEOUT:
                proc.kill()
                proc.wait()
                return None, None, f"timed out after {TIMEOUT} s\n{text[-2000:]}"
    return proc.returncode, vcf_hash(f"{out}.vcf.gz"), log.read_text(errors="replace")


def error_messages(log_j, log_c, java_error=JAVA_ERROR):
    """Java's error line and C's, or None where either is missing.

    Every C error exit (util_exit and four direct exit(1) calls) prints one
    line and exits 1, so C's message is its last line.
    """
    java = java_error.search(log_j)
    c = log_c.strip().splitlines()
    if java is None or not c:
        return None
    return java.group(0).strip().removeprefix(MAIN_THREW).strip(), c[-1].strip()


def check(s):
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        args = build(s, d)
        java = [JAVA, "-jar", JAR, *args]
        rc_j, h_j, log_j = run(java, d / "out")
        c_args = [a.replace(f"out={d}/out", f"out={d}/c") for a in args]
        if s["bgen"] != "none":
            c_args.append(f"bgen={s['bgen']}")
        c = [BEAGLE, *c_args]
        rc_c, h_c, log_c = run(c, d / "c")
        bgen_empty = rc_j == 0 and rc_c == 1 and BGEN_EMPTY in log_c and h_c == h_j is not None
        # Both fail alike when Java prints an error line and C exits 1 with a
        # message. A C signal (negative exit) or abort is never agreement.
        messages = error_messages(log_j, log_c) if rc_j not in (None, 0) and rc_c == 1 else None
        both_succeed = rc_j == 0 and rc_c == 0 and h_j == h_c and h_j is not None
        ok = bgen_empty or messages is not None or both_succeed
        if bgen_empty:
            OUTCOMES["same VCF; no record left for the BGEN, so bgen= fails"] += 1
        elif both_succeed:
            OUTCOMES["both succeed, same VCF"] += 1
        elif messages is not None and messages[0] == messages[1]:
            OUTCOMES[f"both fail: {messages[0][:90]}"] += 1
        elif messages is not None:
            # Counted, not failed: both fail with a message, worded differently.
            OUTCOMES[f"both fail, messages differ: java {messages[0][:60]!r}, C {messages[1][:60]!r}"] += 1
        if not ok:
            if FAIL_DIR.exists():
                shutil.rmtree(FAIL_DIR)
            shutil.copytree(d, FAIL_DIR)
            (FAIL_DIR / "commands.txt").write_text(
                " ".join(java).replace(str(d), str(FAIL_DIR)) + "\n" + " ".join(c).replace(str(d), str(FAIL_DIR)) + "\n"
            )
            raise AssertionError(
                f"java exit={rc_j} hash={h_j}; C exit={rc_c} hash={h_c}\n"
                f"inputs in {FAIL_DIR}\n--- java ---\n{log_j[-2000:]}\n--- C ---\n{log_c[-2000:]}"
            )


def invalid_args(s, d):
    """The scenario's arguments with its invalid change applied."""
    args = build(s, d)
    kind, *v = s["invalid"]
    if kind in ("out-equals-gt", "out-equals-ref"):
        stem = "targ" if kind == "out-equals-gt" else "ref"
        return [a for a in args if not a.startswith("out=")] + [f"out={d}{v[0]}{stem}"]
    if kind == "out-directory":
        (d / "outdir").mkdir()
        return [a for a in args if not a.startswith("out=")] + [f"out={d}/outdir"]
    if kind == "window-overlap":
        overlap, window = v
        return [a for a in args if not a.startswith(("window=", "overlap="))] + [
            f"window={window}",
            f"overlap={overlap}",
        ]
    if kind == "out-of-bounds":
        key, value = v
        return [a for a in args if not a.startswith(f"{key}=")] + [f"{key}={value}"]
    return args + ["foo=1"]


def check_invalid(s):
    """Java and C must both refuse the run with a message and leave the inputs unchanged."""
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        args = invalid_args(s, d)
        inputs = {p: p.read_bytes() for p in d.glob("*.vcf.gz")}
        java = [JAVA, "-jar", JAR, *args]
        rc_j, _, log_j = run(java, d / "run-java")
        c = [BEAGLE, *args]
        rc_c, _, log_c = run(c, d / "run-c")
        messages = error_messages(log_j, log_c, JAVA_REFUSAL) if rc_j == 1 and rc_c == 1 else None
        # Java prefixes an exception's message with its class.
        same = messages is not None and messages[0].endswith(messages[1])
        kind = s["invalid"][0]
        changed = [p.name for p, data in inputs.items() if p.read_bytes() != data]
        if messages is None or not same or changed:
            raise AssertionError(
                f"{s['invalid']}: java exit={rc_j}, C exit={rc_c}, same message: {same}, inputs changed: {changed}\n"
                f"{' '.join(c)}\n--- java ---\n{log_j[-2000:]}\n--- C ---\n{log_c[-2000:]}"
            )
        OUTCOMES[f"invalid {kind}: both refuse, same message"] += 1


def check_runner():
    """Java prints its exception prefix and message in separate writes."""
    message = "java.lang.IllegalArgumentException: complete test message"
    child = f"""
import sys, time
sys.stderr.write({(MAIN_THREW + " ")!r})
sys.stderr.flush()
time.sleep(0.4)
sys.stderr.write('java.lang.IllegalArgumentException: ')
sys.stderr.flush()
time.sleep(0.4)
sys.stderr.write('complete test message\\n')
sys.stderr.flush()
# Beagle's non-daemon workers can keep the process alive after main throws.
time.sleep(10)
"""
    with tempfile.TemporaryDirectory() as tmp:
        rc, _, log = run([sys.executable, "-c", child], Path(tmp) / "partial-exception")
        assert rc == 1 and message in log, f"runner truncated the exception: {log!r}"
    print("PASS runner waits for the complete Java exception line")
    # Beagle 5.4 can print the same window bounds several times in a row.
    child = """
import time
for w in (3, 4, 5):
    print(f"Window {w} [20:1000-2000]", flush=True)
    time.sleep(0.3)
"""
    with tempfile.TemporaryDirectory() as tmp:
        rc, _, log = run([sys.executable, "-c", child], Path(tmp) / "repeated-window")
        assert rc == 0, f"runner stopped a run that repeats window bounds: {log!r}"
    print("PASS runner lets a run repeat window bounds")


def check_regressions():
    for case in sorted(d for d in REGRESSIONS.iterdir() if d.is_dir()):
        args = (case / "args.txt").read_text().split()
        want_rc, want = (case / "expect.txt").read_text().rstrip("\n").split(" ", 1)
        with tempfile.TemporaryDirectory() as tmp:
            cmd = [BEAGLE, f"gt={case}/targ.vcf", f"out={tmp}/c", *args]
            if (case / "ref.vcf").exists():
                cmd.append(f"ref={case}/ref.vcf")
            rc, vcf, log = run(cmd, Path(tmp) / "c")
        # Exit code 0 names the VCF hash, any other a message.
        ok = rc == int(want_rc) and (vcf == want if rc == 0 else want in log)
        if not ok:
            sys.exit(f"FAIL regression {case.name}: exit {rc} hash {vcf}, want {want_rc} and {want!r}\n{log[-2000:]}")
        print(f"PASS regression {case.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--examples", type=int, default=100)
    ap.add_argument("--invalid-examples", type=int, default=2, help="examples per invalid parameter change")
    ap.add_argument("--random", action="store_true", help="new examples on every run")
    ap.add_argument("--reproduce", metavar="BLOB", help="rerun the failure Hypothesis printed as this blob")
    a = ap.parse_args()
    for path in (BEAGLE, JAR):
        if not Path(path).exists():
            sys.exit(f"FAIL {path} is missing")
    check_runner()
    check_regressions()
    fuzz = settings(
        deadline=None,
        derandomize=not a.random,
        database=None,
        suppress_health_check=[HealthCheck.too_slow],
        print_blob=True,
    )
    if a.examples:
        test = settings(fuzz, max_examples=a.examples)(given(scenarios())(check))
        if a.reproduce:
            test = reproduce_failure("6.168.1", a.reproduce.encode())(test)
        test()
    if a.invalid_examples and not a.reproduce:
        for invalid in INVALID:
            settings(fuzz, max_examples=a.invalid_examples)(given(invalid_scenarios(invalid))(check_invalid))()
    for outcome, n in OUTCOMES.most_common():
        print(f"{n:5d}  {outcome}")
    mode = "random" if a.random else "fixed"
    print(
        f"PASS {a.examples} fuzz examples and {a.invalid_examples} for each of {len(INVALID)}"
        f" invalid parameter changes ({mode})"
    )


if __name__ == "__main__":
    main()
