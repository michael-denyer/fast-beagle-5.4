#!/bin/bash
# Run Java Beagle and build/beagle alternately on the chr20 benchmark and
# print one line per run: tool, round, wall s, CPU s (user + sys), peak RSS
# GB, and the output hash as tests/check-oracle.sh computes it. Alternating
# runs share the machine's thermal and load state, so compare them in pairs.
#
# Usage: tests/bench/bench.sh <dir from fetch-chr20.sh> <rounds> [nthreads]
# JAVA overrides the java command (default java -Xmx32g).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
DIR=$1 ROUNDS=$2 THREADS=${3:-18}
JAR="$ROOT/data/beagle.29Oct24.c8e.jar"
read -r -a JAVA <<< "${JAVA:-java -Xmx32g}"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
case "$(uname)" in
  Darwin) TIME=(/usr/bin/time -l) ;;
  *) TIME=(/usr/bin/time -v) ;;
esac

run() {  # name round command...
  local name=$1 round=$2; shift 2
  "${TIME[@]}" "$@" ref="$DIR/ref.vcf.gz" gt="$DIR/target.vcf.gz" map="$DIR/chr20.map" \
    out="$OUT/$name" nthreads="$THREADS" > "$OUT/$name.stdout" 2> "$OUT/$name.time"
  local hash
  hash=$(gzip -dc "$OUT/$name.vcf.gz" | grep -v -e '^##source' -e '^##filedate' | shasum -a 256 | cut -c1-16)
  awk -v n="$name" -v r="$round" -v h="$hash" '
    / real / {wall = $1; cpu = $3 + $5}
    /maximum resident set size/ {rss = $1 / 1e9}
    /Elapsed \(wall clock\)/ {k = split($NF, t, ":"); wall = t[k] + 60 * t[k - 1] + (k == 3 ? 3600 * t[1] : 0)}
    /User time|System time/ {cpu += $NF}
    /Maximum resident set size/ {rss = $NF / 1e6}
    END {printf "%-5s %s %6.2f %7.1f %5.2f %s\n", n, r, wall, cpu, rss, h}' "$OUT/$name.time"
}

echo "tool  round wall_s   cpu_s rss_gb hash"
for i in $(seq 1 "$ROUNDS"); do
  run java "r${i}a" "${JAVA[@]}" -jar "$JAR"
  run c "r${i}a" "$ROOT/build/beagle"
  run c "r${i}b" "$ROOT/build/beagle"
  run java "r${i}b" "${JAVA[@]}" -jar "$JAR"
done
