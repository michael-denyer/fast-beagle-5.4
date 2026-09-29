#!/bin/bash
# Runs a Beagle implementation on arguments that Beagle refuses. Each run must
# exit 1 and print Java's message, and a refused run must leave its input
# files unchanged and write no output.
#
# Usage: tests/check-failures.sh <command...>
#   tests/check-failures.sh java -ea -jar data/beagle.29Oct24.c8e.jar
#   tests/check-failures.sh build/beagle
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

fail=0
failed() { echo "FAIL $*"; fail=1; }

# ABSENT names an output the refused run must not write.
check() {  # name message args...
  local name=$1 message=$2; shift 2
  "${BEAGLE[@]}" "$@" seed=$SEED nthreads=2 > "$OUT/$name.log" 2>&1
  refused "$OUT/$name.log" $? "$message" ${ABSENT:+"$ABSENT"} || { failed "$name $VERDICT"; return; }
  echo "PASS $name: $message"
}

unchanged() {  # name file hash
  [ "$($SHA "$2" | cut -c1-16)" = "$3" ] || failed "$1 changed $(basename "$2")"
}

BEAGLE=("$@")
cp "$DATA/target.vcf.gz" "$OUT/in.vcf.gz"
cp "$DATA/ref.vcf.gz" "$OUT/ref.vcf.gz"
in_hash=$($SHA "$OUT/in.vcf.gz" | cut -c1-16)
ref_hash=$($SHA "$OUT/ref.vcf.gz" | cut -c1-16)
mkdir "$OUT/dir"

# Main.parameters and Main.checkOutputPrefix. File.equals compares normalized
# paths, so a doubled slash still names the input, and the message prints the
# input's normalized path.
norm=$(tr -s / <<< "$OUT")
check out-equals-gt "ERROR: VCF output file equals input file: $norm/in.vcf.gz" gt="$OUT/in.vcf.gz" out="$OUT/in"
unchanged out-equals-gt "$OUT/in.vcf.gz" "$in_hash"
check out-equals-gt-slashes "ERROR: VCF output file equals input file: $norm/in.vcf.gz" gt="$OUT/in.vcf.gz" out="$OUT//in"
unchanged out-equals-gt-slashes "$OUT/in.vcf.gz" "$in_hash"
check gt-slashes "ERROR: VCF output file equals input file: $norm/in.vcf.gz" gt="$OUT//in.vcf.gz" out="$OUT/in"
unchanged gt-slashes "$OUT/in.vcf.gz" "$in_hash"
check out-equals-ref "ERROR: VCF output file equals input file: $norm/ref.vcf.gz" \
  ref="$OUT/ref.vcf.gz" gt="$DATA/target.thin.vcf.gz" out="$OUT/ref"
unchanged out-equals-ref "$OUT/ref.vcf.gz" "$ref_hash"
ABSENT="$OUT/dir.vcf.gz" check out-directory "ERROR: \"out\" parameter cannot be a directory" \
  gt="$DATA/target.vcf.gz" out="$OUT/dir"
ABSENT="$OUT/window.vcf.gz" check window-overlap \
  "ERROR: The \"window\" parameter must be at least 1.1 times the \"overlap\" parameter" \
  gt="$DATA/target.vcf.gz" out="$OUT/window" window=1 overlap=1
# A window shorter than the marker spacing holds no marker, and BasicGT
# indexes the empty marker array.
check empty-window "java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0" \
  gt="$DATA/target.thin.vcf.gz" window=0.0000001 overlap=0.00000001 out="$OUT/empty-window"

# VcfRecGTParser reads a one-character allele as c - '0' and longer ones with
# Integer.parseInt, so a lone Arabic-Indic one (U+0661) is allele 1585.
gzip -dc "$DATA/target.vcf.gz" | LC_ALL=C awk 'BEGIN {OFS = "\t"} /^#/ {print; next}
  !done {$10 = sprintf("%c%c", 217, 161) substr($10, 2); done = 1} {print}' | gzip > "$OUT/one-char.vcf.gz"
check one-char-allele "ERROR: Invalid allele [$(printf '\331\241')]" gt="$OUT/one-char.vcf.gz" out="$OUT/one-char-out"

# BasicMarker exits on a malformed marker. Each file edits the 100th record of
# target.vcf.gz (22:20008093, REF G, ALT A), so only one record fails:
# parallel parsing makes the message of the first of several failures vary.
marker_check() {  # name message field value
  local name=$1 message=$2
  gzip -dc "$DATA/target.vcf.gz" | awk -v f="$3" -v v="$4" 'BEGIN {OFS = "\t"} /^#/ {print; next} ++n == 100 {$f = v} {print}' \
    | gzip > "$OUT/$name.vcf.gz"
  check "$name" "$message" gt="$OUT/$name.vcf.gz" out="$OUT/$name-out"
}
marker_check end-before-pos "ERROR: invalid INFO:END field at 22:20008093 [20008092]" 8 END=20008092
marker_check ref-not-acgtn "ERROR: REF field is not a sequence of A, C, T, G, or N characters at 22:20008093 [R]" 4 R
marker_check duplicate-alleles "ERROR: duplicate allele at 22:20008093 [G, G]" 5 G

# Beagle 5.5's parameters, which Beagle 5.4's Par does not know.
check initial-lr "Error: unrecognized parameter: initial-lr=100" \
  gt="$DATA/target.vcf.gz" out="$OUT/initial-lr" initial-lr=100
check window-markers "Error: unrecognized parameter: window-markers=100000" \
  gt="$DATA/target.vcf.gz" out="$OUT/window-markers" window-markers=100000

# A 2.5 cM gap in the target markers, longer than the window. Imputing, Beagle
# 5.4 fails on the window 22:20049939-20068772, which has no target markers.
# Phasing only, it fails on a window with one position.
check imp-gap-map "22:20049939-20068772" \
  ref="$DATA/ref.vcf.gz" gt="$DATA/target.gap.vcf.gz" map="$DATA/map.map" window=1.5 overlap=0.5 out="$OUT/gap"
check imp-gap-map-noimp "java.lang.IllegalArgumentException: Window has only one position: CHROM=22 POS=20029411" \
  ref="$DATA/ref.vcf.gz" gt="$DATA/target.gap.vcf.gz" map="$DATA/map.map" window=1.5 overlap=0.5 impute=false out="$OUT/gap-noimp"

# PlinkGenMap throws IllegalArgumentException, and prints a genetic position
# with Double.toString.
map_check() {  # name message map-lines...
  local name=$1 message=$2; shift 2
  printf '%s\n' "$@" > "$OUT/$name.map"
  check "$name" "java.lang.IllegalArgumentException: $message" gt="$DATA/target.vcf.gz" map="$OUT/$name.map" out="$OUT/$name"
}
map_check map-same-cm "All loci in genetic map have the same genetic position [0.480887]: 22 . 0.480887 20000000" \
  "22 . 0.480887 20000000" "22 . 0.480887 20100000"
map_check map-same-cm-sci "All loci in genetic map have the same genetic position [1.0E-4]: 22 . 0.0001 20000000" \
  "22 . 0.0001 20000000" "22 . 0.00010 20100000"
map_check map-infinite "invalid map position: Infinity" "22 . 0.1 20000000" "22 . +Infinity 20100000"
map_check map-duplicate "duplication position: 22 . 0.2 20000000" "22 . 0.1 20000000" "22 . 0.2 20000000"
map_check map-descending "map positions not in ascending order: 22 . 0.05 20100000" "22 . 0.1 20000000" "22 . 0.05 20100000"
map_check map-short-line "Map file format error: 22 . 0.1" "22 . 0.1"
map_check map-long-line "Map file format error: 22 . 0.1 20000000 x" "22 . 0.1 20000000 x" "22 . 0.2 20100000"
map_check map-missing-chrom "missing genetic map for chromosome 22" "21 . 0.1 20000000" "21 . 0.2 20100000"

# Bref3Reader.readMarker indexes SNV_PERMS with alleleCode >> 2, so allele code
# byte 0x80 on the first marker (20000086, rs138720731) is permutation -32.
LC_ALL=C perl -0777 -pe 's/(\x01\x31\x2d\x56\x01\x00\x0brs138720731)./$1\x80/s or die "no first marker\n"' \
  "$DATA/ref.bref3" > "$OUT/snv-code.bref3" || failed "bref3-snv-code cannot edit $DATA/ref.bref3"
check bref3-snv-code "java.lang.ArrayIndexOutOfBoundsException" \
  ref="$OUT/snv-code.bref3" gt="$DATA/target.thin.vcf.gz" out="$OUT/snv-code"

# A bref3 header with 2^30 samples, truncated after four empty IDs; the sample
# count doubled overflows int. Java first allocates String[2^30], so it reads
# to EOF with a heap of 4 GB or more and runs out of memory with less. Only the
# port's result, EOF, is fixed.
if [ "$1" != java ]; then
  perl -e 'print pack("N", 2055763188), pack("n", 1), "x", pack("N", 0x40000000), "\0\0" x 4' > "$OUT/samples.bref3"
  check bref3-sample-count "java.io.EOFException" \
    ref="$OUT/samples.bref3" gt="$DATA/target.thin.vcf.gz" out="$OUT/samples"
fi

# ImpIbs divides imp-states by round(imp-segment / imp-step), which is 0 here.
check imp-segment "java.lang.ArithmeticException: / by zero" \
  ref="$DATA/ref.vcf.gz" gt="$DATA/target.thin.vcf.gz" out="$OUT/segment" imp-segment=0.01

exit $fail
