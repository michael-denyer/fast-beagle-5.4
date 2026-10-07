# Sourced by the scripts that run Beagle on the case tables tests/*-cases.txt:
# how a case row becomes a run, and how a run is judged. The caller sets ROOT.
# shellcheck shell=bash
DATA="$ROOT/data"
SEED=-99999
SHA=$(command -v sha256sum || echo "shasum -a 256")
# Preparation errors must stop even callers that deliberately run without -e.
"$ROOT/tests/fetch-fixtures.sh" --ensure >&2 || { echo "FAIL fixtures: tests/fetch-fixtures.sh --ensure" >&2; exit 1; }

# The case rows (name expect tags args) of the named tables, default the
# oracle's. expect is the recorded hash or exit=<status>.
cases() {  # [table...]
  awk '!/^#/ && NF' "${@:-$ROOT/tests/oracle-cases.txt}"
}

# True when CASES is unset or names the case.
selected() {  # name
  [ -z "${CASES:-}" ] || [[ " $CASES " == *" $1 "* ]]
}

# Fails when CASES names a case the tables lack, so that a mistyped filter
# cannot select nothing and pass.
check_selection() {  # table...
  local known name bad=0
  known=" $(cases "$@" | awk '{printf "%s ", $1}')"
  for name in ${CASES:-}; do
    [[ "$known" == *" $name "* ]] || { echo "FAIL $name is not a case in $*"; bad=1; }
  done
  return $bad
}

has_tag() {  # tags tag
  [[ ",$1," == *",$2,"* ]]
}

# A case's hash field is one hash, or thread:hash pairs separated by commas
# when Beagle's output depends on the thread count.
thread_dependent() {  # expect
  [[ "$1" == *:* ]]
}

want_hash() {  # hashes nthreads
  if thread_dependent "$1"; then
    tr ',' '\n' <<< "$1" | awk -F: -v t="$2" '$1 == t {print $2}'
  else
    echo "$1"
  fi
}

vcf_hash() {  # vcf.gz
  gzip -dc "$1" | grep -v '^##source\|^##filedate' | $SHA | cut -c1-16
}

# Runs command with the case's arguments (@ as the data directory), then out,
# seed and nthreads, writing its standard output and error to <out>.run.log,
# since Beagle writes <out>.log itself.
case_run() {  # args out nthreads command...
  local args=$1 out=$2 t=$3; shift 3
  # shellcheck disable=SC2086  # args is the case's argument list; splitting it is the point
  "$@" ${args//@/$DATA/} out="$out" seed=$SEED nthreads="$t" > "$out.run.log" 2>&1
}

# Runs a case and judges it: the run must exit as expect says, and a case with
# a hash must exit 0 and write the VCF with the hash for nthreads. Sets VERDICT
# to the exit status and hash, followed on failure by what was wanted.
case_verdict() {  # expect args out nthreads command...
  local expect=$1 out=$3 t=$4 rc want; shift
  case_run "$@"
  rc=$?
  if [[ "$expect" == exit=* ]]; then
    VERDICT="exit=$rc"
    want=$expect
  else
    VERDICT="exit=$rc $(vcf_hash "$out.vcf.gz")"
    want="exit=0 $(want_hash "$expect" "$t")"
  fi
  [ "$VERDICT" = "$want" ] || { VERDICT="$VERDICT, want $want"; return 1; }
}

# Judges a run that Beagle must refuse: exit 1, the message in its log, and no
# regular file at any of the paths (a caller may put a directory there to make
# the run fail). A crash after the message exits otherwise. Sets VERDICT to
# what failed.
refused() {  # log rc message [path-that-must-be-absent...]
  local log=$1 rc=$2 message=$3 path; shift 3
  VERDICT=
  [ "$rc" -eq 1 ] || { VERDICT="exit=$rc, want 1: $(tail -1 "$log")"; return 1; }
  grep -qF -- "$message" "$log" || { VERDICT="output lacks '$message': $(tail -1 "$log")"; return 1; }
  for path; do
    [ ! -f "$path" ] || { VERDICT="left $(basename "$path")"; return 1; }
  done
}
