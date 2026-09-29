#!/bin/bash
# The parity ratchet of the Beagle 5.4 port. tests/c54-pending.txt lists the
# checks of build/beagle whose result still differs from Beagle 5.4's, one key
# per line (kind, then case). The check scripts judge each such result here:
# - an unlisted key that differs fails, as it did before the list existed;
# - a listed key that differs prints "pending <key>" and passes;
# - a listed key that now matches fails, so each fix removes its line and the
#   list only shrinks.
# Sourced by tests/cases.sh. Run as a script (tests/c54-ratchet.sh <key>
# <status>), it judges one key for tests/check_fuzz.py.
# shellcheck shell=bash
C54_PENDING="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/c54-pending.txt"

# True when tests/c54-pending.txt lists the key.
pending() {  # key
  awk -v key="$1" '$0 == key {found = 1} END {exit !found}' "$C54_PENDING"
}

# The label of a result that differs: "differs" while the key is listed, FAIL
# otherwise.
miss_label() {  # key
  if pending "$1"; then echo differs; else echo FAIL; fi
}

# Judges a result of build/beagle, status 0 when it matched Beagle 5.4. Returns
# non-zero when the check must fail.
ratchet() {  # key status
  if ! pending "$1"; then
    return "$2"
  elif [ "$2" -eq 0 ]; then
    echo "FAIL $1 matches Beagle 5.4: remove $1 from tests/c54-pending.txt"
    return 1
  fi
  echo "pending $1"
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  ratchet "$@"
fi
