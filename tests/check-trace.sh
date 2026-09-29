#!/bin/bash
# Run the Java trace build and the C binary on every case and compare the named
# trace seams byte for byte, then each case's exit code and the seam files C
# writes. Needs `make build/beagle java-trace` first.
#
# Usage: tests/check-trace.sh <seam...>      e.g. tests/check-trace.sh T1a
# CASES restricts the run to the named cases (default all); SKIP leaves cases out.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
fail=0
"$ROOT/tests/run-trace.sh" "$OUT/java" java '-Dbeagle.trace={trace}' -cp "$ROOT/build/java-trace/classes" main.Main > "$OUT/java.txt" || fail=1
"$ROOT/tests/run-trace.sh" "$OUT/c" "$ROOT/build/beagle" 'trace={trace}' > "$OUT/c.txt" || fail=1
# Each case must exit as Java's does (Java exits 1 on edge-markers and
# imp-chroms-noimpute), and C must write no seam file that Java does not.
while read -r name exit_java _; do
  exit_c=$(awk -v n="$name" '$1 == n {print $2}' "$OUT/c.txt")
  if [ "$exit_java" = "$exit_c" ]; then
    echo "PASS $name $exit_java"
  else
    echo "FAIL $name java $exit_java, C ${exit_c:-no run}"
    tail -5 "$OUT/java/$name/out.run.log" "$OUT/c/$name/out.run.log"
    fail=1
  fi
done < "$OUT/java.txt"
for c in "$OUT"/c/*/*.txt; do
  [ -e "$c" ] || continue
  file=${c#"$OUT/c/"}
  if [ ! -e "$OUT/java/$file" ]; then
    echo "FAIL $file written by C only"
    fail=1
  fi
done

compared=" "
for dir in "$OUT"/java/*/; do
  name=$(basename "$dir")
  for seam in "$@"; do
    for java in "$dir/$seam.txt" "$dir/$seam"-*.txt; do
      [ -e "$java" ] || continue
      file=$(basename "$java")
      compared="$compared$seam "
      if cmp -s "$java" "$OUT/c/$name/$file"; then
        echo "PASS $name $file $(wc -l < "$java" | tr -d ' ') lines"
      else
        echo "FAIL $name $file"
        diff "$java" "$OUT/c/$name/$file" 2>&1 | head -5
        fail=1
      fi
    done
  done
done
# A seam with nothing to compare fails: the Java trace build may be missing,
# or the seam renamed.
for seam in "$@"; do
  if [[ "$compared" != *" $seam "* ]]; then
    echo "FAIL no $seam files from the Java trace build"
    fail=1
  fi
done
exit $fail
