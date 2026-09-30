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
# htslib is built without libcurl, GCS and S3 support, plugins and CRAM's bzip2
# and lzma codecs, so it needs only libdeflate, built here, and zlib, which
# every system has. libdeflate compresses BGZF about 1.5 times as fast as zlib.
# libdeflate 1.25 is the version that third_party/libdeflate holds.
HTSLIB_VERSION=1.24
HTSLIB_SHA256=28a8de191381c7a97a35675ceac76fa1ea95e7b678d6a2e9d600a7874e4077de
LIBDEFLATE_VERSION=1.25
LIBDEFLATE_SHA256=fed5cd22f00f30cc4c2e5329f94e2b8a901df9fa45ee255cb70e2b0b42344477

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
libdeflate_src=$work/libdeflate-$LIBDEFLATE_VERSION
mkdir -p "$work"
curl -fSL -o "$work/htslib.tar.bz2" \
  "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2"
curl -fSL -o "$work/libdeflate.tar.gz" \
  "https://github.com/ebiggers/libdeflate/releases/download/v$LIBDEFLATE_VERSION/libdeflate-$LIBDEFLATE_VERSION.tar.gz"
if command -v sha256sum > /dev/null; then sum=sha256sum; else sum="shasum -a 256"; fi
echo "$HTSLIB_SHA256  $work/htslib.tar.bz2" | $sum -c -
echo "$LIBDEFLATE_SHA256  $work/libdeflate.tar.gz" | $sum -c -
tar -xjf "$work/htslib.tar.bz2" -C "$work"
tar -xzf "$work/libdeflate.tar.gz" -C "$work"
# Only static libraries go in the prefix, so -lhts and -ldeflate cannot pick
# shared ones.
mkdir -p "$htslib_prefix/include/htslib" "$htslib_prefix/lib"
(
  cd "$libdeflate_src"
  # libdeflate builds from its lib/ sources alone. Both lib/arm and lib/x86
  # have a cpu_features.c, so each object is named after its whole path.
  for src in lib/*.c lib/*/*.c; do
    ${CC:-cc} -O2 -fPIC -c -o "${src//\//_}.o" "$src"
  done
  ar rcs "$htslib_prefix/lib/libdeflate.a" ./*.o
  cp libdeflate.h "$htslib_prefix/include/"
)
(
  cd "$htslib_src"
  ./configure CFLAGS="-g -O2 -fPIC" CPPFLAGS="-I$htslib_prefix/include" LDFLAGS="-L$htslib_prefix/lib" \
    --disable-libcurl --disable-gcs --disable-s3 --disable-plugins --disable-bz2 --disable-lzma --with-libdeflate
  make -j"$jobs" lib-static
)
cp "$htslib_src"/htslib/*.h "$htslib_prefix/include/htslib/"
cp "$htslib_src/libhts.a" "$htslib_prefix/lib/"

# HTSLIB_PREFIX is empty so that the Makefile adds no rpath to the build
# directory. The static htslib needs libdeflate and zlib after it on the link
# line. build/beagle already holds third_party/libdeflate's compressor, so the
# link takes from libdeflate.a only what that copy lacks: the decompressor,
# crc32 and the gzip wrappers.
LDLIBS="-ldeflate -lz" make -j"$jobs" HTSLIB_PREFIX= CFLAGS="-O2 -g -I$htslib_prefix/include" \
  LDFLAGS="-L$htslib_prefix/lib" build/beagle

if [ "$os" = linux ]; then deps=$(ldd build/beagle); else deps=$(otool -L build/beagle); fi
echo "$deps"
if grep -E 'libhts|libdeflate' <<< "$deps"; then
  echo "build-static.sh: build/beagle loads a shared htslib or libdeflate" >&2
  exit 1
fi

name=fast-beagle-$version-$os-$arch
stage=$work/$name
mkdir -p "$stage" dist
cp build/beagle "$stage/fast-beagle"
cp LICENSE NOTICE.md "$stage/"
cp "$htslib_src/LICENSE" "$stage/LICENSE.htslib"
cp "$libdeflate_src/COPYING" "$stage/LICENSE.libdeflate"
tar -czf "dist/$name.tar.gz" -C "$work" "$name"
echo "dist/$name.tar.gz"
