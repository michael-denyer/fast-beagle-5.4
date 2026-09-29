#!/bin/bash
# Checks tests/gate-tier.sh on changed-path lists with known tiers, including a
# documentation-only list, which must run no tier, a new test script, which
# must run the full tier without an entry of its own, and an empty list, which
# must fail rather than pick a tier.
#
# Usage: tests/check-gate-tier.sh
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
fail=0

expect() {  # tier paths...
  local want=$1 got; shift
  got=$(printf '%s\n' "$@" | "$ROOT/tests/gate-tier.sh")
  if [ "$got" = "$want" ]; then
    echo "PASS $want: $*"
  else
    echo "FAIL $*: got '$got', want $want"; fail=1
  fi
}

expect c src/phase/phase_ls.c
expect c src/main/run_stats.c src/vcf/sliding_window.h
expect none docs/usage.md README.md
expect none docs/logo.jpg
expect none docs/how-it-works.webp images/diagram.svg LICENSE
expect c src/phase/phase_ls.c docs/usage.md README.md
expect c third_party/libdeflate/lib/utils.c
expect full src/jcompat/jnum.c
expect full tests/check-log.sh
expect full tests/check-something-new.sh
expect full tests/logs/gt.log
expect full src/main/main.c tests/bgen-hashes.txt
expect full Makefile
expect full java/trace.patch
expect full tla/ParallelOrdered.tla
expect full .github/workflows/gate.yml
expect full docs/testing.md tests/gate-steps.sh

if "$ROOT/tests/gate-tier.sh" < /dev/null > /dev/null 2>&1; then
  echo "FAIL an empty path list picked a tier"; fail=1
else
  echo "PASS an empty path list fails"
fi
exit $fail
