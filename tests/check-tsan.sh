#!/bin/bash
# Builds build/beagle under ThreadSanitizer in build/tsan, so the normal build
# is untouched, then runs every oracle case at 18 threads, and the
# thread-dependent cases again with trace= (parallel_for workers write trace
# seams). Every run must pass the oracle verdict with no ThreadSanitizer
# report. It runs on macOS only: ThreadSanitizer cannot start in the gate's
# docker leg, because it needs an address layout that the x86_64 emulation
# gives only with ASLR off, and docker forbids turning ASLR off.
#
# Usage: tests/check-tsan.sh
# CASES restricts the run to the named cases (default all).
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"

if [ "$(uname -s)" != Darwin ]; then
  echo "ThreadSanitizer runs on the macOS host only (see the header)"
  exit 1
fi
TSAN="$ROOT/build/tsan"
mkdir -p "$TSAN"
rsync -a --delete --exclude build "$ROOT/Makefile" "$ROOT/java" "$ROOT/src" "$ROOT/third_party" "$ROOT/tests" "$TSAN/"
# make does not rebuild objects when only CFLAGS change.
make -C "$TSAN" clean > /dev/null
make -C "$TSAN" CFLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=thread" build/beagle check-block-reader \
  > "$TSAN/make.log" 2>&1 || { echo "FAIL ThreadSanitizer build ($TSAN/make.log)"; exit 1; }
echo "PASS ThreadSanitizer build"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

check_selection "$ROOT/tests/oracle-cases.txt" || exit 1
fail=0
while read -r name expect _ args; do
  selected "$name" || continue
  for traced in no yes; do
    [ $traced = yes ] && ! thread_dependent "$expect" && continue
    out="$OUT/$name.tsan.$traced"
    trace=()
    [ $traced = yes ] && { mkdir -p "$out.trace"; trace=(trace="$out.trace"); }
    case_verdict "$expect" "$args" "$out" 18 "$TSAN/build/beagle" ${trace[@]+"${trace[@]}"}
    ok=$?
    if grep -q ThreadSanitizer "$out.run.log"; then
      echo "FAIL $name nthreads=18 traced=$traced: ThreadSanitizer report"; grep -A20 ThreadSanitizer "$out.run.log" | head -40; fail=1
    elif [ $ok -eq 0 ]; then
      echo "PASS $name nthreads=18 traced=$traced tsan"
    else
      echo "FAIL $name nthreads=18 traced=$traced tsan $VERDICT"; tail -5 "$out.run.log"; fail=1
    fi
  done
done < <(cases)
exit $fail
