#!/bin/bash
# Builds fast-beagle against the htslib in the host environment. The Makefile
# appends its fixed floating-point flags after the compiler's CFLAGS.
set -euo pipefail
make -j"${CPU_COUNT}" install PREFIX="${PREFIX}" HTSLIB_PREFIX="${PREFIX}"
