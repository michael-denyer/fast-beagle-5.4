#!/bin/bash
# The gate's checks, in order, run in the checkout <root>. tests/check-local.sh
# runs them natively and on Linux x86_64 in docker, and CI runs them on each
# runner. Prints one pass or FAIL line per check; each check's output goes to
# build/check-<name>.log, or build/check-c-<name>.log in the C tier.
#
# Usage: [GATE_TIER=c] [GATE_FUZZ=random] tests/gate-steps.sh <root>
# Needs Java 21, htslib and uv; the full tier also needs PLINK2 naming the
# pinned plink2 binary (see tests/check-bgen.sh). GATE_TIER=c runs only the checks that compare the C
# binary against recorded results; Java then only builds the bref3 fixtures.
# It skips every check that runs Java or the jar alongside it, the
# fixture-cache check, the sanitizers and the TLA+ model. It checks the
# BGEN output against the hashes in tests/bgen-hashes.txt instead of
# plink2, and runs the saved fuzz regressions but no new fuzz
# examples. The default is the full gate. GATE_FUZZ=random fuzzes new examples
# instead of the fixed 200.
# shellcheck disable=SC2329  # java_build, oracle_trace and trace_threads run through step
set -uo pipefail
cd "$1" || exit 1
SEAMS="T1a T1b T1c T1d T2 T2b T3a T3b0 T3b1 T3b T3c T3d T3e T4a T4b T4c T4d T5a T5b T5c T5d"
# Seams whose content depends on nthreads, rechecked at 1 and 18 threads on the
# cases whose windows are long enough for the thread count to split them: the
# cases with per-thread hashes.
THREAD_SEAMS="T3b0 T3b1 T3b T3c T3d T3e T4a T4b T4c T4d"

mkdir -p build
fail=0
tier=${GATE_TIER:-full}
case $tier in full|c) ;; *) echo "GATE_TIER must be full or c, not $tier"; exit 2 ;; esac
bgen_oracle=live logs=build/check-
[ "$tier" = c ] && bgen_oracle=recorded logs=build/check-c-
fuzz_args=(--examples 200)
case ${GATE_FUZZ:-fixed} in
  fixed) ;;
  random) fuzz_args+=(--random) ;;
  *) echo "GATE_FUZZ must be fixed or random, not $GATE_FUZZ"; exit 2 ;;
esac
step() {  # name command...
  local name=$1; shift
  if "$@" > "$logs$name.log" 2>&1; then
    echo "  pass  $name"
  else
    echo "  FAIL  $name ($logs$name.log)"; fail=1
  fi
}
full_step() {  # name command...: runs only in the full tier
  if [ "$tier" = full ]; then step "$@"; else echo "  skip  $1 (full tier only)"; fi
}

java_build() {
  mkdir -p build/classes \
    && find java/src -name '*.java' > build/java-sources.txt \
    && javac -nowarn -d build/classes @build/java-sources.txt
}

oracle_trace() {
  mkdir -p build/trace-oracle \
    && tests/check-oracle.sh java -Dbeagle.trace=build/trace-oracle -cp build/java-trace/classes main.Main
}

trace_threads() {
  local t thread_cases
  # shellcheck source=cases.sh
  thread_cases=$(ROOT=$PWD; source tests/cases.sh; cases | while read -r name expect _; do
    if thread_dependent "$expect"; then printf '%s ' "$name"; fi; done) || return 1
  for t in 1 18; do
    # shellcheck disable=SC2086  # the seam list splits into arguments
    NTHREADS=$t CASES="$thread_cases" tests/check-trace.sh $THREAD_SEAMS || return 1
  done
}

step fixtures tests/fetch-fixtures.sh
full_step cases python3 tests/check_cases.py
full_step jcompat make check-jcompat
step tracker make check-tracker
step interval make check-interval
full_step oracle-jar tests/check-oracle.sh java -ea -jar data/beagle.29Oct24.c8e.jar
full_step failures-jar tests/check-failures.sh java -ea -jar data/beagle.29Oct24.c8e.jar
full_step java-build java_build
full_step oracle-source tests/check-oracle.sh java -ea -cp build/classes main.Main
full_step java-trace make java-trace
full_step oracle-trace oracle_trace
step c-build make build/beagle
step oracle-c tests/check-oracle.sh build/beagle
step failures-c tests/check-failures.sh build/beagle
full_step log tests/check-log.sh
step piece-size make check-piece-size
step bgen-unit make check-bgen-unit
step records make check-records
step vcf-index make check-vcf-index
step tbi make check-tbi
step bgen env BGEN_ORACLE="$bgen_oracle" tests/check-bgen.sh
# shellcheck disable=SC2086  # the seam list splits into arguments
full_step trace tests/check-trace.sh $SEAMS
full_step sanitizers tests/check-sanitizers.sh
full_step tla tests/check-tla.sh
full_step fuzz uv run --python 3.12 --script tests/check_fuzz.py "${fuzz_args[@]}"
if [ "$tier" = c ]; then
  step fuzz-regressions uv run --python 3.12 --script tests/check_fuzz.py --examples 0 --invalid-examples 0
fi
full_step trace-threads trace_threads
exit $fail
