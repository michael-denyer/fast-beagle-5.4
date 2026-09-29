#!/bin/bash
# Runs a Beagle implementation on each oracle case at nthreads=2 and compares
# its <out>.log with Beagle's recorded log in tests/logs/<case>.log. Both are
# masked: timings, timestamps, the data directory, out= and the lines that name
# the program, and fast-beagle's added CPU time and max memory lines are
# dropped; every other line must match. For fast-beagle, standard output must
# also equal the log, and a run that fails on a malformed reference must print
# one error message and end its log with it, even when every parse worker fails
# at once.
#
# Usage: tests/check-log.sh <command...>
#   tests/check-log.sh java -ea -jar data/beagle.29Oct24.c8e.jar
#   tests/check-log.sh build/beagle
# RECORD=1 rewrites tests/logs/ from the command's logs; record with the jar.
# CASES restricts the run to the named cases (default all).
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
LOGS="$ROOT/tests/logs"
BEAGLE=("$@")

mask() {  # log
  sed -E \
    -e '/^Enter "java -jar /d' \
    -e '/^(CPU time|Max memory): /d' \
    -e 's/^(beagle\.29Oct24\.c8e\.jar|fast-beagle: a C port of beagle\.29Oct24\.c8e\.jar) \(version 5\.4\)$/<banner>/' \
    -e 's/^(Start|End) time: .*/\1 time: <time>/' \
    -e 's/^Command line: .*/Command line: <program>/' \
    -e 's/^  out=.*/  out=<out>/' \
    -e "s#=$DATA/#=@/#" \
    -e 's/^(beagle\.29Oct24\.c8e\.jar|fast-beagle) finished$/<program> finished/' \
    -e 's/^(.{31})([0-9]+ hours? )?([0-9]+ minutes? )?[0-9]+ seconds?$/\1<elapsed>/' \
    "$1"
}

check_selection "$ROOT/tests/oracle-cases.txt" || exit 1
if [ -n "${RECORD:-}" ]; then
  [ -z "${CASES:-}" ] || { echo "FAIL RECORD=1 needs every case"; exit 1; }
  mkdir "$OUT/logs" || exit 1
fi
fail=0
while read -r name expect _ args; do
  selected "$name" || continue
  out="$OUT/$name"
  if ! case_verdict "$expect" "$args" "$out" 2 "${BEAGLE[@]}"; then
    echo "FAIL $name: $VERDICT"; tail -3 "$out.run.log"; fail=1
  elif [ ! -s "$out.log" ]; then
    echo "FAIL $name: no log"; tail -3 "$out.run.log"; fail=1
  elif [ -n "${RECORD:-}" ]; then
    if mask "$out.log" > "$OUT/logs/$name.log"; then
      echo "PASS $name recorded $(wc -l < "$out.log" | tr -d ' ') lines"
    else
      echo "FAIL $name: could not mask the log"; fail=1
    fi
  elif [ "$1" != java ] && ! cmp -s "$out.log" "$out.run.log"; then
    echo "FAIL $name: standard output differs from the log"
    diff "$out.log" "$out.run.log" | head -10; fail=1
  elif ! diff "$LOGS/$name.log" <(mask "$out.log") > "$out.diff" 2>&1; then
    echo "FAIL $name: the log differs from tests/logs/$name.log"; head -20 "$out.diff"; fail=1
  else
    echo "PASS $name $(wc -l < "$out.log" | tr -d ' ') lines"
  fi
done < <(cases)

# bad-alleles makes every parse worker fail at once.
check_bad_ref() {  # label ref
  local label=$1 out="$OUT/$1" rc
  "${BEAGLE[@]}" gt="$DATA/target.vcf.gz" ref="$2" out="$out" nthreads=18 > /dev/null 2> "$out.err"
  rc=$?
  if [ "$rc" -ne 1 ] || [ "$(wc -l < "$out.err")" -ne 1 ] || [ "$(tail -1 "$out.log")" != "$(cat "$out.err")" ]; then
    echo "FAIL $label: exit=$rc, log ends '$(tail -1 "$out.log")', stderr:"; head -c 600 "$out.err"; echo; fail=1
  else
    echo "PASS $label: the log ends with '$(cat "$out.err")'"
  fi
}
if [ "$1" != java ] && [ -z "${CASES:-}" ]; then
  echo junk | gzip > "$OUT/junk.vcf.gz"
  check_bad_ref bad-ref "$OUT/junk.vcf.gz"
  gzip -dc "$DATA/ref.vcf.gz" | awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$NF = "7|7"; print}' | gzip > "$OUT/bad-alleles.ref.vcf.gz"
  check_bad_ref bad-alleles "$OUT/bad-alleles.ref.vcf.gz"
fi
if [ -n "${RECORD:-}" ] && [ "$fail" -eq 0 ]; then
  rm -rf "$LOGS" && cp -R "$OUT/logs" "$LOGS" || fail=1
fi
exit $fail
