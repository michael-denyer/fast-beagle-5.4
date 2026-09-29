# Releases

A release of this repository is a Beagle 5.4 (29Oct24) edition of fast-beagle, versioned 5.4.x. [fast-beagle-5.5](https://github.com/michael-denyer/fast-beagle-5.5) releases the 5.5.x versions of the same `fast-beagle` package. Each release publishes a GitHub Release with:

- `fast-beagle-<version>.tar.gz`, the source tarball
- `fast-beagle-<version>-linux-x86_64.tar.gz`, `fast-beagle-<version>-linux-arm64.tar.gz` and `fast-beagle-<version>-macos-arm64.tar.gz`, each holding a `fast-beagle` binary with htslib linked statically, and the licenses
- `SHA256SUMS`, the SHA-256 of every file above

The bioconda package `fast-beagle` builds from the source tarball with the recipe in `release/conda/`.

## Cut a release

1. Merge every change for the release into `main`, with the gate passing.
2. Tag the `main` commit `v5.4.<n>` and push the tag. The release workflow fails on any other tag name.

   ```bash
   git tag v5.4.0 origin/main
   git push origin v5.4.0
   ```

3. Wait for the release workflow run on the tag. It publishes the GitHub Release once every build and check has passed. If a job fails, nothing is published. Fix the cause on `main`, then delete the tag and push it again on the fixed commit.
4. Submit the recipe to [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) as `recipes/fast-beagle/`. The run's `conda-recipe` artifact holds `recipe.yaml` with the release's version and the tarball's sha256 filled in, with `build.sh` and `small.vcf`. After the first submission, bioconda's autobump bot opens the pull request for each new release.

## What the release workflow proves

`.github/workflows/release.yml` runs on a pushed `v*` tag, on pull requests that change the workflow, `release/` or the `Makefile`, and by hand. Pull requests and manual runs do everything but publish. They use the version in `release/conda/recipe.yaml`.

- The `source` job makes the tarball with `git archive`, with every path under `fast-beagle-<version>/`. It fills in the recipe's version and sha256 from the tarball. Every other job builds from this tarball.
- The `conda` jobs build the recipe with rattler-build on `linux-64`, `linux-aarch64` and `osx-arm64`, from conda-forge and bioconda, as bioconda does. The build uses conda's compilers and its `CFLAGS` (`-O3` on `linux-aarch64`), and the Makefile appends its fixed floating-point flags after them. The Makefile's `-Werror` applies, so a warning from conda's compilers fails the build. rattler-build runs the recipe's test, a phasing run on `small.vcf`. The job then installs the package into a new environment and runs `tests/check-oracle.sh` on the installed `fast-beagle`, so every oracle case must match its hash at 1, 2 and 18 threads.
- The `binaries` jobs run `release/build-static.sh`. It builds htslib 1.24, the version that bioconda packages, from its checksum-verified release tarball, without libcurl, GCS and S3 support, plugins, bzip2, lzma and libdeflate, so that it needs only zlib. It links htslib statically and fails if the binary loads a shared htslib. The Linux binaries build in `manylinux_2_28` containers pinned by digest, so they need glibc 2.28 or later, and run on the runner itself for the check. Each job runs `tests/check-oracle.sh` on the binary unpacked from its package.
- The `publish` job runs only on a pushed tag, after every other job has passed. It writes `SHA256SUMS` and creates the GitHub Release with the default `GITHUB_TOKEN`, which has `contents: write` in this job only.

rattler-build 0.76.1, micromamba 2.9.0-0 and htslib 1.24 are pinned by version and SHA-256 in the workflow and in `release/build-static.sh`.
