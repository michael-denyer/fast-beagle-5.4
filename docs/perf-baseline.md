# Performance

This page compares `build/beagle` with Java Beagle 5.4 (`beagle.29Oct24.c8e.jar`). Both tools write the same VCF, so the comparison is of time and memory only. The gain grows with the size of the imputation.

| Run | Target samples | Threads | Faster | Less CPU time | Less max memory |
|---|---|---|---|---|---|
| [1000 Genomes chr20](#1000-genomes-chr20-benchmark) | 321 | 18 | 1.7× | 2.0× | 3.2× |
| [Production scale](#production-scale-run) | about 100,000 | 60 | 2.8× | 2.1× | 4.4× |

The chr20 figures are the ratio of the Java median to the C median over 6 runs per tool on one day. The chr20 benchmark uses public data, so anyone can rerun it.

## Production-scale run

This run imputed one chromosome for about 100,000 target samples from a bref3 reference. It ran on a Databricks node of type `Standard_E64ds_v6`, which has 64 vCPUs (Intel Xeon Platinum 8573C) and 512 GiB of memory, with Databricks Runtime 16.4 LTS on Ubuntu 24.04. Both tools ran with `nthreads=60`, and both wrote the same VCF. The C build was commit `c7e7aa8`, built in an `ubuntu:24.04` container. Java ran with `-Xmx500g`. The figures cover the Beagle step alone, one run per tool on 2026-09-29, taken from each process's resource usage. The input data is not public, so this run cannot be reproduced from this repository.

| | Java | C | Java / C |
|---|---|---|---|
| Wall time | 773 s | 277 s | 2.8× faster |
| CPU time | 27,251 s | 13,219 s | 2.1× less |
| Max memory | 367.5 GB | 83.9 GB | 4.4× less |

Java kept about 60 cores busy for its first 480 s, then spent about 285 s on 1 to 3 cores.

The C run did not write BGEN. With BGEN output added, the same step took 336 s and 17,097 s of CPU time at 138.8 GB max memory, wrote the same VCF, and wrote an 8.2 GB `.bgen` file.

## 1000 Genomes chr20 benchmark

### Input

`tests/bench/fetch-chr20.sh <dir>` downloads the 1000 Genomes high-coverage chr20 phased panel and Beagle's GRCh38 genetic map. It checks their SHA-256 and derives three files with `tests/bench/make_chr20.py`:

- `target.vcf.gz`: 321 samples (every 10th panel sample), unphased, at 15,879 markers (every 10th biallelic SNV with 0.05 <= AF <= 0.95).
- `ref.vcf.gz`: the other 2,881 samples at 1,642,181 markers (every panel record except symbolic structural variants and repeats of an earlier position, REF and ALT).
- `chr20.map`: the chr20 PLINK map with chromosome `chr20`.

```bash
tests/bench/fetch-chr20.sh ~/beagle-bench
```

### Method

`tests/bench/bench.sh <dir> <rounds> [nthreads]` runs Java and C in the order Java, C, C, Java in each round, so both tools see the same thermal and load state. It records wall time, CPU time (user plus system) and max memory (maximum resident set size) with `/usr/bin/time`. It also records the output hash as `tests/check-oracle.sh` computes it: the first 16 hex digits of the SHA-256 of the decompressed VCF without its `##source` and `##filedate` lines. The script runs the jar at `data/beagle.29Oct24.c8e.jar`, which `tests/fetch-fixtures.sh` downloads, and `build/beagle` as `make` builds it (the Makefile's default `CFLAGS`, `-O2 -g`, plus its fixed flags).

```bash
make
JAVA="/opt/homebrew/opt/openjdk@21/bin/java -Xmx32g" tests/bench/bench.sh ~/beagle-bench 3
```

### Result

These runs used commit `80f140d` on 2026-09-29, with 18 threads. Every run of both tools wrote the VCF with hash `41b92f8850e50b23`.

| Round | Java wall | Java CPU | Java max memory | C wall | C CPU | C max memory |
|---|---|---|---|---|---|---|
| 1 | 31.6 s, 32.2 s | 470 s, 482 s | 8.53 GB, 8.39 GB | 18.1 s, 18.2 s | 227 s, 228 s | 2.65 GB, 2.65 GB |
| 2 | 31.3 s, 34.3 s | 466 s, 485 s | 8.33 GB, 8.27 GB | 19.8 s, 19.7 s | 238 s, 237 s | 2.65 GB, 2.64 GB |
| 3 | 34.2 s, 32.5 s | 511 s, 488 s | 8.24 GB, 9.15 GB | 19.7 s, 18.8 s | 247 s, 238 s | 2.64 GB, 2.65 GB |

| | Java | C | C / Java |
|---|---|---|---|
| Wall time, median (range) | 32.3 s (31.3 to 34.3) | 19.2 s (18.1 to 19.8) | 0.60 |
| CPU time, median (range) | 484 s (466 to 511) | 237 s (227 to 247) | 0.49 |
| Max memory, median (range) | 8.4 GB (8.2 to 9.2) | 2.65 GB (2.64 to 2.65) | 0.32 |

Each cell lists the tool's two runs in the order they ran. Each round ran as `tests/bench/bench.sh ~/beagle-bench 1`, so that the machine could be checked before each round: on AC power, and no other job using more than 2 cores. Two more rounds wrote the same hash but are left out, because another job started during them.

Java's max memory depends on `-Xmx` and the garbage collector. These runs used `-Xmx32g` and the default collector.

The machine was an Apple M5 Pro (6 super and 12 performance cores) with 64 GB and macOS 27.0. Java was OpenJDK 21.0.12.1 (Homebrew). The C build used Apple clang 21.0.0 and htslib 1.24.
