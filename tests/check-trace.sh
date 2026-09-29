#!/bin/bash
# Run the Java trace build and the C binary on every case and compare the named
# trace seams byte for byte, then each case's exit code and the seam files C
# writes. Needs `make build/beagle java-trace` first.
#
# Each seam of a case is judged under the parity ratchet (tests/c54-ratchet.sh)
# with the key "trace <case> <seam> nthreads=<n>": it differs when one of its
# files differs from Java's or C writes one that Java does not.
#
# Usage: tests/check-trace.sh <seam...>      e.g. tests/check-trace.sh T1a
# CASES restricts the run to the named cases (default all); SKIP leaves cases out.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=c54-ratchet.sh
source "$ROOT/tests/c54-ratchet.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
T=${NTHREADS:-2}
fail=0
"$ROOT/tests/run-trace.sh" "$OUT/java" java '-Dbeagle.trace={trace}' -cp "$ROOT/build/java-trace/classes" main.Main > "$OUT/java.txt" || fail=1
"$ROOT/tests/run-trace.sh" "$OUT/c" "$ROOT/build/beagle" 'trace={trace}' > "$OUT/c.txt" || fail=1
# Each case must exit as Java's does (Java exits 1 on edge-markers and
# imp-chroms-noimpute).
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

# The seam a trace file belongs to: T3d-0-bwd.txt is T3d's.
seam_of() {  # file
  local seam=${1%.txt}
  echo "${seam%%-*}"
}

# "<case> <seam>" per line, for every seam compared and every seam that differs.
: > "$OUT/judged.txt"
: > "$OUT/differs.txt"
for c in "$OUT"/c/*/*.txt; do
  [ -e "$c" ] || continue
  file=${c#"$OUT/c/"}
  if [ ! -e "$OUT/java/$file" ]; then
    name=${file%%/*}
    pair="$name $(seam_of "${file#*/}")"
    echo "$pair" >> "$OUT/judged.txt"
    echo "$pair" >> "$OUT/differs.txt"
    echo "$(miss_label "trace $pair nthreads=$T") $file written by C only"
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
      echo "$name $seam" >> "$OUT/judged.txt"
      if cmp -s "$java" "$OUT/c/$name/$file"; then
        echo "PASS $name $file $(wc -l < "$java" | tr -d ' ') lines"
      else
        echo "$(miss_label "trace $name $seam nthreads=$T") $name $file"
        diff "$java" "$OUT/c/$name/$file" 2>&1 | head -5
        echo "$name $seam" >> "$OUT/differs.txt"
      fi
    done
  done
done
# A listed seam of a case that ran, with no file from either build, matches.
while read -r name seam; do
  if grep -q "^$name " "$OUT/java.txt" && [[ " $* " == *" $seam "* ]]; then
    echo "$name $seam" >> "$OUT/judged.txt"
  fi
done < <(awk -v t="nthreads=$T" '$1 == "trace" && $4 == t {print $2, $3}' "$C54_PENDING")
while read -r pair; do
  status=0
  grep -qxF "$pair" "$OUT/differs.txt" && status=1
  ratchet "trace $pair nthreads=$T" $status || fail=1
done < <(sort -u "$OUT/judged.txt")
# A seam with nothing to compare fails: the Java trace build may be missing,
# or the seam renamed.
for seam in "$@"; do
  if [[ "$compared" != *" $seam "* ]]; then
    echo "FAIL no $seam files from the Java trace build"
    fail=1
  fi
done
exit $fail
