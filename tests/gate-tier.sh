#!/bin/bash
# Prints the gate tier a change needs, given its changed paths one per line on
# standard input. none: every path is documentation only (docs/, Markdown,
# images and LICENSE), which no check reads. c: every other path is under src/
# (except src/jcompat/) or third_party/, which the C tier checks. full: any
# other path, including every test script and table, so a new full-tier check
# needs no entry here.
#
# Usage: git diff --name-only origin/main | tests/gate-tier.sh
set -uo pipefail
tier=none n=0
while IFS= read -r path; do
  n=$((n + 1))
  case $path in
    docs/* | *.md | *.png | *.svg | *.jpg | *.jpeg | *.gif | *.webp | LICENSE) ;;
    src/jcompat/*) tier=full ;;
    src/* | third_party/*) [ "$tier" = none ] && tier=c ;;
    *) tier=full ;;
  esac
done
[ "$n" -gt 0 ] || { echo "gate-tier.sh: no changed paths on standard input" >&2; exit 1; }
echo "$tier"
