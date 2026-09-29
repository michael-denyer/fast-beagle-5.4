#!/bin/bash
# Runs build/beagle with tbi=true on oracle cases at 1 and 18 threads. The VCF
# must have the oracle hash and match a run without tbi=true, and
# <out>.vcf.gz.tbi must be byte-identical to the index htslib's
# tbx_index_build3 (tabix -p vcf) builds from that VCF.
#
# Usage: tests/check-tbi.sh    (after make build/beagle build/output/vcf_index_test)
# CASES overrides the cases run.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
BEAGLE=$ROOT/build/beagle
TABIX="$ROOT/build/output/vcf_index_test tabix"
# Several chromosomes, INFO END=, multiallelic records, and phasing only with
# 17 BGZF blocks.
CASES=${CASES:-imp-chroms imp-bref3-end imp-bref3-multi gt-map}

check_selection "$ROOT/tests/oracle-cases.txt" || exit 1
fail=0
while read -r name expect _ args; do
  selected "$name" || continue
  for t in 1 18; do
    out="$OUT/$name.t$t"
    if ! case_verdict "$expect" "$args tbi=true" "$out" "$t" "$BEAGLE"; then
      echo "FAIL $name nthreads=$t $VERDICT"; tail -5 "$out.run.log"; fail=1; continue
    fi
    case_run "$args" "$out.plain" "$t" "$BEAGLE"
    if [ "$(vcf_hash "$out.plain.vcf.gz")" != "$(vcf_hash "$out.vcf.gz")" ]; then
      echo "FAIL $name nthreads=$t: the VCF differs from a run without tbi=true"; fail=1
    elif [ -e "$out.plain.vcf.gz.tbi" ]; then
      echo "FAIL $name nthreads=$t: a run without tbi=true wrote an index"; fail=1
    elif ! $TABIX "$out.vcf.gz" "$out.tabix.tbi"; then
      echo "FAIL $name nthreads=$t: tabix cannot index the VCF"; fail=1
    elif ! cmp "$out.vcf.gz.tbi" "$out.tabix.tbi"; then
      echo "FAIL $name nthreads=$t: the index differs from tabix's"; fail=1
    else
      echo "PASS $name nthreads=$t $VERDICT, $(($(wc -c < "$out.vcf.gz"))) bytes, index matches tabix"
    fi
  done
done < <(cases)
exit $fail
