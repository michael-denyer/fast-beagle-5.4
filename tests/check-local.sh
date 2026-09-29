#!/bin/bash
# The pre-merge gate, run on this machine in place of CI: the lint hooks once,
# then the checks in tests/gate-steps.sh (the list CI runs) natively (arm64 on
# an Apple Silicon Mac) and on Linux x86_64 in docker, each in the full tier and
# then in the C tier that CI runs on pull requests, without plink2.
#
# Usage: tests/check-local.sh
# Needs Java 21, htslib, bgzip and prek on the host, and a running docker (colima)
# that can run linux/amd64 images and mount this checkout. The bgen check needs
# plink2 v2.0.0-a.7.8 (see docs/testing.md): PLINK2_ARM64 and PLINK2_LINUX override the
# default paths of the macOS arm64 and Linux x86_64 binaries. It also needs uv on
# the host; the Linux leg installs uv UV_VERSION, checked against UV_LINUX_SHA256.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PLINK2_ARM64=${PLINK2_ARM64:-$HOME/Claude/beagle-tools/plink2-20260919/mac/plink2}
PLINK2_LINUX=${PLINK2_LINUX:-$HOME/Claude/beagle-tools/plink2-20260919/linux/plink2}
UV_VERSION=0.11.32
UV_LINUX_SHA256=aab924fd522efd06f1c5f3b93a243864fc453132c94b2dc49f1371b528a4b967

fail=0
mkdir -p "$ROOT/build"
if (cd "$ROOT" && prek run --all-files) > "$ROOT/build/check-lint.log" 2>&1; then
  echo "  pass  lint"
else
  echo "  FAIL  lint (build/check-lint.log)"; fail=1
fi
echo "host $(uname -m)"
PLINK2="$PLINK2_ARM64" "$ROOT/tests/gate-steps.sh" "$ROOT" || fail=1
echo "host $(uname -m) c tier"
env -u PLINK2 GATE_TIER=c "$ROOT/tests/gate-steps.sh" "$ROOT" || fail=1

# Linux x86_64: a copy of the checkout without build outputs, in a directory
# the docker VM can mount (colima mounts the home directory).
LINUX_DIR="$HOME/.cache/fast-beagle-check/$(basename "$ROOT")"
mkdir -p "$LINUX_DIR"
rsync -a --delete --exclude build --exclude data --exclude .git "$ROOT/" "$LINUX_DIR/"
echo "linux x86_64 (docker)"
docker run --rm --platform linux/amd64 -v "$LINUX_DIR":/w -w /w -v "$PLINK2_LINUX":/opt/plink2:ro -e PLINK2=/opt/plink2 \
    -e UV_VERSION="$UV_VERSION" -e UV_LINUX_SHA256="$UV_LINUX_SHA256" eclipse-temurin:21-jdk bash -c '
  apt-get update >/dev/null 2>&1 && apt-get install -y gcc make libhts-dev tabix patch curl perl rsync python3 >/dev/null 2>&1 \
    || { echo "  FAIL  apt-get"; exit 1; }
  curl -LSf -o /tmp/uv.tar.gz "https://github.com/astral-sh/uv/releases/download/$UV_VERSION/uv-x86_64-unknown-linux-gnu.tar.gz" \
    && echo "$UV_LINUX_SHA256  /tmp/uv.tar.gz" | sha256sum -c >/dev/null && tar -xzf /tmp/uv.tar.gz -C /usr/local/bin --strip-components=1 \
    || { echo "  FAIL  uv install"; exit 1; }
  tests/gate-steps.sh /w; full=$?
  echo "linux x86_64 c tier"
  env -u PLINK2 GATE_TIER=c tests/gate-steps.sh /w || exit 1
  exit $full' || fail=1

[ $fail -eq 0 ] && echo "ALL PASS" || echo "FAILED"
exit $fail
