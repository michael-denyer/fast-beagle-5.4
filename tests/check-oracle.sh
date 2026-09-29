#!/bin/bash
# Run a Beagle implementation on every case in tests/oracle-cases.txt and
# compare its output to the recorded hash.
#
# Usage: tests/check-oracle.sh <command...>
#   tests/check-oracle.sh java -ea -jar data/beagle.29Oct24.c8e.jar
#   tests/check-oracle.sh java -ea -cp build/classes main.Main
#   tests/check-oracle.sh build/beagle            (the C binary, once it exists)
# NTHREADS overrides the thread counts tried (default "1 2 18").
# CASES restricts the run to the named cases (default all).
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

check_selection "$ROOT/tests/oracle-cases.txt" || exit 1
fail=0
while read -r name expect _ args; do
  selected "$name" || continue
  for t in ${NTHREADS:-1 2 18}; do
    out="$OUT/$name.t$t"
    if case_verdict "$expect" "$args" "$out" "$t" "$@"; then
      echo "PASS $name nthreads=$t $VERDICT"
    else
      echo "FAIL $name nthreads=$t $VERDICT"; tail -5 "$out.run.log"; fail=1
    fi
  done
done < <(cases)
exit $fail
