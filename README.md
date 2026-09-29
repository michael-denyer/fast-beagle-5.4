# Fast Beagle v5.4

![C11](https://img.shields.io/badge/C-C11-00599C?logo=c&logoColor=white)
![Python test scripts](https://img.shields.io/badge/Python-test_scripts-3776AB?logo=python&logoColor=white)
![Beagle 5.4 byte-identical](https://img.shields.io/badge/Beagle_5.4-byte--identical-2E7D32)
![macOS](https://img.shields.io/badge/macOS-supported-D32F2F?logo=apple&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-supported-FCC624?logo=linux&logoColor=black)
![License GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue)

<img src="docs/logo.jpg" width="144" align="right" alt="Beagle with a DNA helix">

**1.7× to 2.8× faster, with over 2× less CPU time and up to 4.4× less memory than Java Beagle 5.4**, and byte-identical output. The gain grows with scale, from a 321-sample public benchmark to a production imputation of about 100,000 samples. [Performance](docs/perf-baseline.md) has both runs.

fast-beagle is a standalone C port of [Beagle 5.4](https://faculty.washington.edu/browning/beagle/beagle.html) (29Oct24), the genotype phasing and imputation tool. Its output is byte-identical to the Java release `beagle.29Oct24.c8e.jar` run with the same `nthreads=`. It adds BGEN v1.2 output and a tabix index written in the same pass as the VCF.

This repository is the Beagle 5.4 edition of [fast-beagle](https://github.com/michael-denyer/fast-beagle), which ports Beagle 5.5 (27Feb25). The two editions share their design, their checks and the outputs they add. Each one writes the output of its own Beagle release, and Beagle 5.4 and 5.5 phase differently. Use the edition that matches the Beagle version whose output you need.

## How it works

Beagle phases each target sample's genotypes into two haplotypes, then imputes the markers the target lacks from a reference panel. It works along the chromosome in overlapping windows. fast-beagle runs the same steps and writes the same output, and can also write BGEN.

![Four steps. 1, input genotypes: observed allele pairs from a target VCF, with some markers untyped. 2, phase: estimate which alleles belong on each chromosome, giving two haplotypes. 3, impute: use a reference panel to estimate the untyped variants, with probabilities. 4, output: a phased VCF, dosages and quality scores, and optional BGEN v1.2. Reference haplotypes are optional for phasing and required for imputation. An optional genetic map provides recombination distances.](docs/how-it-works.webp)

## Install

Build fast-beagle from source. It needs a C11 compiler, `make` and htslib.

```bash
git clone https://github.com/michael-denyer/fast-beagle-5.4.git
brew install htslib              # macOS
sudo apt-get install libhts-dev  # Debian and Ubuntu
make -C fast-beagle-5.4
```

`make` builds the binary `build/beagle`.

## Getting started

Run `build/beagle` with Beagle's `key=value` arguments:

```bash
build/beagle gt=target.vcf.gz ref=ref.vcf.gz map=plink.map out=result
```

The run writes `result.vcf.gz`. Without `ref=` it phases the target only. Add `bgen=plink2` or `bgen=phased` to write `result.bgen`, `result.sample` and `result.info` as well. Add `tbi=true` to write `result.vcf.gz.tbi`.

## Details

- [Parameters and output files](docs/usage.md)
- [How fast-beagle differs from Beagle 5.4](docs/beagle-divergences.md)
- [Architecture and source layout](docs/architecture.md)
- [Checks and the pre-merge gate](docs/testing.md)
- [Byte identity with Beagle 5.4](docs/byte-identity.md)
- [Performance](docs/perf-baseline.md)
- [Domain terms](CONTEXT.md)

## Contributing

Bug reports and pull requests are welcome. Before you open a pull request, run the [lint hooks](docs/testing.md#run-the-lint-hooks) and the [pre-merge gate](docs/testing.md#run-the-pre-merge-gate). A change to the output must keep every oracle hash.

## License

fast-beagle is licensed under [GPL-3.0-or-later](LICENSE), the same license as Beagle, because it is a derivative work of Beagle 5.4 and 5.5 (Copyright (C) 2014-2024 Brian L. Browning). [NOTICE.md](NOTICE.md) lists the third-party code in this repository and its licenses.
