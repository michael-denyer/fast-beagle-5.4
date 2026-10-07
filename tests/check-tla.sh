#!/bin/bash
# Model-checks the TLA+ specs in tla/ with TLC over a matrix of constants.
# docs/testing.md says what each spec models. Each run checks the spec's
# invariants, no deadlock, and its liveness properties. Needs Java; downloads
# tla2tools.jar TLA_VERSION into data/, checked against TLA_SHA256.
#
# Usage: tests/check-tla.sh [Spec ...]   (default: every spec in SPECS)
# shellcheck disable=SC2329  # the function named after each spec is called through SPECS
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SPECS=(ParallelOrdered BlockReader SlidingWindow FatalExit BgenCleanup)
[ $# -gt 0 ] || set -- "${SPECS[@]}"
for spec in "$@"; do
  known=0
  for s in "${SPECS[@]}"; do [ "$s" != "$spec" ] || known=1; done
  [ $known = 1 ] || { echo "FAIL unknown spec $spec"; exit 1; }
done
TLA_VERSION=v1.7.4
TLA_SHA256=936a262061c914694dfd669a543be24573c45d5aa0ff20a8b96b23d01e050e88
JAR="$ROOT/data/tla2tools-$TLA_VERSION.jar"
SHA=$(command -v sha256sum || echo "shasum -a 256")
# As in tests/fetch-fixtures.sh: run java rather than look for it.
java -version > /dev/null 2>&1 || {
  echo "FAIL java cannot run, and TLC needs it: put a JDK first on PATH" >&2
  exit 1
}

mkdir -p "$ROOT/data"
if [ ! -f "$JAR" ]; then
  curl -LSf -o "$JAR.part" "https://github.com/tlaplus/tlaplus/releases/download/$TLA_VERSION/tla2tools.jar" \
    || { echo "FAIL download tla2tools.jar $TLA_VERSION"; exit 1; }
  mv "$JAR.part" "$JAR"
fi
got=$($SHA "$JAR" | cut -d' ' -f1)
[ "$got" = "$TLA_SHA256" ] || { echo "FAIL $JAR has sha256 $got, want $TLA_SHA256"; exit 1; }

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
fail=0
run=0

# check <Spec> <label> <invariants> <properties> <<< "<constant assignments>":
# one TLC run of tla/<Spec>.tla.
check() {
  local spec=$1 label=$2
  run=$((run + 1))
  { echo "CONSTANTS"; cat; echo "SPECIFICATION Spec"; echo "INVARIANTS $3"; echo "PROPERTIES $4"; } > "$OUT/MC.cfg"
  cp "$ROOT/tla/$spec.tla" "$OUT/"
  if java -XX:+UseParallelGC -cp "$JAR" tlc2.TLC -workers 4 -metadir "$OUT/states-$run" \
      -config "$OUT/MC.cfg" "$OUT/$spec.tla" > "$OUT/tlc.log" 2>&1 \
      && grep -q "No error has been found" "$OUT/tlc.log"; then
    echo "PASS $spec $label $(grep -o '[0-9,]* distinct states found' "$OUT/tlc.log" | head -1)"
  else
    echo "FAIL $spec $label: $(grep -m1 "^Error:" "$OUT/tlc.log")"; grep -v "^\s*$" "$OUT/tlc.log" | tail -60; fail=1
  fi
}

# worker_set <n>: w1,...,w<n>
worker_set() {
  seq -f "w%g" 1 "$1" | paste -sd, -
}

ParallelOrdered() {
  # workers items window
  for model in "1 3 1" "2 4 1" "2 5 2" "3 5 2" "3 6 3" "2 3 5" "3 8 3" "4 7 2"; do
    read -r nw ni win <<< "$model"
    check ParallelOrdered "workers=$nw items=$ni window=$win" \
      "TypeOK WindowBound NoOverwrite ConsumesOwnItem InOrder" "Terminates AllConsumed" <<EOF
  Workers = {$(worker_set "$nw")}
  NItems = $ni
  Window = $win
EOF
  done
}

BlockReader() {
  local invariants="TypeOK SlotsPartition ReaderSlotExclusive InOrder"
  # The code has BLOCK_READER_SLOTS = 4; Slots = 1 is the boundary where a
  # missing broadcast stalls the pipeline.
  for slots in 1 2 3 4; do
    for nb in 0 1 3 5; do
      check BlockReader "slots=$slots batches=$nb close=any" \
        "$invariants" "CloseTerminates WaitEnds Progress" <<EOF
  Slots = $slots
  NBatches = $nb
  MayClose = TRUE
EOF
      check BlockReader "slots=$slots batches=$nb close=never" \
        "$invariants" "WaitEnds SeesSentinel" <<EOF
  Slots = $slots
  NBatches = $nb
  MayClose = FALSE
EOF
    done
  done
}

SlidingWindow() {
  # windows failat
  for model in "1 {}" "1 {1}" "2 {}" "2 {1,2}" "3 {}" "3 {1,2,3}" "4 {}" "4 {1,2,3,4}" "4 {3}"; do
    read -r nw failat <<< "$model"
    check SlidingWindow "windows=$nw failat=$failat" \
      "TypeOK InOrder AheadNotGiven NoLeak ErrorAfterPrefix ErrorNoAhead" \
      "CloseTerminates NextReturns ExitOnlyOnError" <<EOF
  NWindows = $nw
  FailAt = $failat
EOF
  done
}

FatalExit() {
  for nw in 1 2; do
    check FatalExit "workers=$nw" \
      "TypeOK" "Terminates LockReleased MainNotBlocked DeferredRaised" <<EOF
  Workers = {$(worker_set "$nw")}
EOF
  done
}

BgenCleanup() {
  check BgenCleanup "members=3" \
    "TypeOK NoPartialAfterCleanup CompletedSurvive" "Terminates" <<EOF
  Members = {bgen, info, sample}
EOF
}

for spec in "$@"; do "$spec"; done
exit $fail
