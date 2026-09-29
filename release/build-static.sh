#!/bin/bash
# Builds a portable fast-beagle binary with htslib linked statically, checks
# that it loads no shared htslib, and packages it with the licenses as
# dist/fast-beagle-<version>-<os>-<arch>.tar.gz. Run from the root of a source
# tree. The release workflow runs it in a manylinux_2_28 container on Linux,
# so the Linux binary needs glibc 2.28 or later.
#
# Usage: release/build-static.sh <version>
set -euo pipefail
version=$1
# htslib is built without libcurl, GCS and S3 support, plugins, CRAM's bzip2
# and lzma codecs and libdeflate, so it needs only zlib, which every system has.
HTSLIB_VERSION=1.24
HTSLIB_SHA256=28a8de191381c7a97a35675ceac76fa1ea95e7b678d6a2e9d600a7874e4077de

case $(uname -s) in
  Linux) os=linux jobs=$(nproc) ;;
  Darwin) os=macos jobs=$(sysctl -n hw.ncpu); export MACOSX_DEPLOYMENT_TARGET=11.0 ;;
  *) echo "build-static.sh: unsupported system $(uname -s)" >&2; exit 1 ;;
esac
case $(uname -m) in
  x86_64) arch=x86_64 ;;
  aarch64 | arm64) arch=arm64 ;;
  *) echo "build-static.sh: unsupported machine $(uname -m)" >&2; exit 1 ;;
esac

make clean
work=$PWD/build/static
htslib_src=$work/htslib-$HTSLIB_VERSION
htslib_prefix=$work/htslib
mkdir -p "$work"
curl -fSL -o "$work/htslib.tar.bz2" \
  "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2"
if command -v sha256sum > /dev/null; then sum=sha256sum; else sum="shasum -a 256"; fi
echo "$HTSLIB_SHA256  $work/htslib.tar.bz2" | $sum -c -
tar -xjf "$work/htslib.tar.bz2" -C "$work"
(
  cd "$htslib_src"
  ./configure CFLAGS="-g -O2 -fPIC" --disable-libcurl --disable-gcs --disable-s3 --disable-plugins \
    --disable-bz2 --disable-lzma --without-libdeflate
  make -j"$jobs" lib-static
)
# Only the static library goes in the prefix, so -lhts cannot pick a shared one.
mkdir -p "$htslib_prefix/include/htslib" "$htslib_prefix/lib"
cp "$htslib_src"/htslib/*.h "$htslib_prefix/include/htslib/"
cp "$htslib_src/libhts.a" "$htslib_prefix/lib/"

# HTSLIB_PREFIX is empty so that the Makefile adds no rpath to the build
# directory. The static htslib needs zlib after it on the link line.
LDLIBS=-lz make -j"$jobs" HTSLIB_PREFIX= CFLAGS="-O2 -g -I$htslib_prefix/include" LDFLAGS="-L$htslib_prefix/lib" \
  build/beagle

if [ "$os" = linux ]; then deps=$(ldd build/beagle); else deps=$(otool -L build/beagle); fi
echo "$deps"
if grep libhts <<< "$deps"; then
  echo "build-static.sh: build/beagle loads a shared htslib" >&2
  exit 1
fi

name=fast-beagle-$version-$os-$arch
stage=$work/$name
mkdir -p "$stage" dist
cp build/beagle "$stage/fast-beagle"
cp LICENSE NOTICE.md "$stage/"
cp "$htslib_src/LICENSE" "$stage/LICENSE.htslib"
tar -czf "dist/$name.tar.gz" -C "$work" "$name"
echo "dist/$name.tar.gz"
