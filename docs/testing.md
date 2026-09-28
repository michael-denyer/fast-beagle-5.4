# Checks and the pre-merge gate

Most checks compare fast-beagle with Beagle 5.4 or with a tool whose output it must match. The unmodified Beagle 5.4 Java source in `java/src/` and the release jar `beagle.29Oct24.c8e.jar` are the oracle. The C engine is still the Beagle 5.5 port, so the checks that compare it with Beagle 5.4 go through the [parity ratchet](#parity-ratchet). [Byte identity with Beagle 5.5](byte-identity.md) summarises what these checks proved for the Beagle 5.5 port.

## Check an implementation

`tests/check-oracle.sh` takes the command of any implementation. To check the release jar:

```bash
tests/check-oracle.sh java -ea -jar data/beagle.29Oct24.c8e.jar
```

To check fast-beagle, run `tests/check-oracle.sh build/beagle`.

## Fixtures

`tests/fetch-fixtures.sh` downloads Beagle's public test data and release jars, verifies checksums, and derives the fixtures that the case tables run. The fixtures include:

- the reference/target split
- chrX with haploid samples
- multiallelic markers
- missing genotypes
- a genetic map on which the region spans 5 cM (5 windows with `window=1.5 overlap=0.5`)
- a target with a 50 kb gap in its markers
- exclusion lists
- bref3 references, one of them with END= on SNV and indel records
- a bref3 edge case (IDs with a 4-byte character, an invalid byte and two entries, a marker with no ALT allele, a symbolic allele)
- Beagle 5.4's marker and record text: FILTER values and END= in INFO, END= with leading zeros or given twice, and ID lists with a `.` entry
- a flat genetic map whose plateau gives genetic distances below 1e-7 cM
- a 20-sample reference whose sequence coding reaches its limit
- an exclusion list of reference samples only, which Beagle 5.4 does not apply to a bref3 reference
- for `tests/check-bgen.sh`, the chrX split moved to chromosome 22 and the reference/target split moved to chromosome 38

## Oracle hashes

`tests/check-oracle.sh` runs an implementation on the fixtures at 1, 2 and 18 threads and compares its output with the hashes in `tests/oracle-cases.txt`. `NTHREADS` overrides the thread counts.

- `CASES="gt imp"` restricts the run. A name that is not in the case tables fails the run. The same holds for the sanitizer, trace and BGEN runners.
- Beagle 5.4's output depends on the thread count when a window is longer than 4.5 cM. Such cases record one hash per thread count, as `thread:hash` pairs.

## Case tables

Every case table row holds a name, the expected outcome, tags and Beagle's arguments. The expected outcome is a hash, or `exit=N` where the row records no hash. The tags (`nonautosome`, `multiallelic`) tell `tests/check-bgen.sh` how plink2 reads the output.

- `tests/cases.sh` reads the tables. It holds the one verdict that `tests/check-oracle.sh`, `tests/check-sanitizers.sh` and `tests/check-bgen.sh` apply to a run. The verdict requires the expected exit and, for a case with a hash, the hash for the run's thread count.
- `tests/cases.sh` also holds the one refusal rule that `tests/check-failures.sh` and `tests/check-bgen.sh` apply. The rule requires exit 1, the expected message, and no file at the paths the caller lists.
- `tests/check_cases.py` tests the verdict and the refusal rule with a fake Beagle.

## Log and console output

`tests/check-log.sh` runs Java Beagle and `build/beagle` on each oracle case and compares their `<out>.log` files. It masks the timings, the start and end times, `out=` and the lines that name the program, and drops fast-beagle's `CPU time:` and `Max memory:` lines. Every other line must match. fast-beagle's standard output must also equal its log, and a run that fails on a malformed reference must end its log with the error message. `CASES` restricts the cases and `NTHREADS` sets the thread count, default 2. The full gate tier runs it.

## Refused inputs

`tests/check-failures.sh` runs an implementation on arguments and inputs that Beagle refuses. Each run must exit 1 with Java's message. The cases are:

- an `out=` that names the `gt=` or `ref=` file or a directory. A refused `out=` must leave the input unchanged and write no VCF.
- `window` less than 1.1 times `overlap`
- `imp-segment` below half of `imp-step`
- genetic map files that Beagle rejects:
  - genetic positions that are all equal
  - an infinite position
  - a repeated base position
  - descending positions
  - lines with too few or too many fields
  - no map for the target's chromosome
- a one-character non-ASCII GT allele
- a marker that Beagle 5.4 rejects: END= before POS, a REF allele with a character other than A, C, G, T or N, and a repeated allele
- the Beagle 5.5 parameters `initial-lr=` and `window-markers=`, which Beagle 5.4 does not know
- a bref3 SNV allele code whose permutation index is negative
- a bref3 header whose sample count overflows when doubled. This case runs for the C build only, because the jar's result depends on its heap size.

## Compare trace seams

`java/trace.patch` holds trace hooks for the Beagle 5.5 Java source and does not apply to the Beagle 5.4 source in `java/src/`. Until it is rebased, the gate prints `skip  <step> (5.4 trace patch pending)` for `java-trace`, `oracle-trace`, `trace` and `trace-threads`.

`java/trace.patch` holds trace hooks for the Java source. `make java-trace` applies the patch to a copy in `build/java-trace/`. When that build runs with `-Dbeagle.trace=<dir>`, it writes each trace seam to `<dir>/<seam>.txt`. To change the hooks, edit a patched copy and regenerate the patch with `diff -ruN` against `java/src`.

- `tests/run-trace.sh` runs an implementation on every case with tracing on, so two runs can be compared with `diff -r`.
- `tests/check-trace.sh <seam...>` runs the Java trace build and `build/beagle` on every case, including the input edge cases in `tests/trace-cases.txt`, and compares the named trace seams.

### Seam formats

Fields are tab-separated. Floating-point values print as raw bits in hex (`float` 8 digits at most, `double` 16). A digest is FNV-1a over 32-bit words. Lists are comma-separated, and `-` marks an empty list. The hook for each seam is in `java/trace.patch`, and the `Trace` class there documents each record.

| Seam | File | Records |
|---|---|---|
| T1a | `T1a-target`, `T1a-ref` | Per input record: CHROM, POS, the stored IDs joined by `;` (`.` for none), REF, the ALT alleles joined by `,` (`.` for none), the allele count, and the END value as an int (empty for none). |
| T1b | `T1b-target`, `T1b-ref` | Per sample: ID and `diploid` or `haploid`. A bref3 reference lists every sample, because Beagle 5.4 does not apply `excludesamples=` to bref3. |
| T1c | `T1c-target` | Per target record: record class, phase flag, allele per haplotype (`.` for missing). |
| T1d | `T1d-ref` | Per reference record: record class, major allele, sequence group and sequence-to-allele map for a `SeqCodedRefGTRec` (`-` otherwise), allele per haplotype. A `group` line gives each group's haplotype-to-sequence map once. |
| T2 | `T2` | Per window: index, chromosome, last-window flag, `endCm` (double bits), end position, marker counts, overlap and splice indices; then the positions and the target-to-marker map. |
| T2b | `T2b` | Per window: fixed phasing data, per marker the position, distance, previous stage-1 marker and weight, and the carriers of each allele; then the stage-1 markers, their distances and the step ends. |
| T3a | `T3a` | Per window: stage-1 minor allele frequencies, IBS2 markers, step starts, and per sample its IBS2 segments. |
| T3b0 | `T3b0-<start>` | Per initial PBWT sub-window: its marker range and seed, then the alleles at each marker. |
| T3b1 | `T3b1` | Per window: `pd` seed, recombination intensity, mismatch probability; `pRecomb` per stage-1 marker; `leaveUnph`, the per-sample proportion of unphased heterozygotes left unphased each iteration; then one sample line per target sample (below). |
| T3b | `T3b` | After each stage-1 iteration: `it` iteration count and swap rate (double bits), then one sample line per target sample. |
| T3c | `T3c` | Per parameter update: the mismatch and switch data in summation order, then the estimated and stored mismatch probability and recombination intensity. |
| T3d | `T3d-<it>-fwd`, `T3d-<it>-bwd` | Per IBS search: candidate count, steps per batch, overlap steps, each step's haplotype-to-sequence map, and each step's IBS haplotypes. |
| T3e | `T3e` | Per stage-1 iteration: `it` iteration index, then one line per target sample (below). |
| T4a | `T4a` | Per window: the stage-2 IBS haplotypes, forward then backward, per step and target haplotype. |
| T4b | `T4b` | Per target haplotype: stage-2 state count and a digest of the states and state probabilities. |
| T4c | `T4c` | Per stage-2 marker: the haplotypes that carry each allele, `M` for the implicit major allele. |
| T4d | `T4d` | Per target sample: a digest of the stage-2 phase and imputation probabilities. |
| T5a | `T5a` | Per window: imputation cluster count, haplotype counts, mismatch rate; per cluster its ranges, position, mismatch and switch probabilities and a digest of its sequence map; a digest of the interpolation weights. |
| T5b | `T5b` | Per window: the imputation steps, then each step's IBS reference haplotypes per target haplotype. |
| T5c | `T5c` | Per target haplotype: state count, a digest of the states and matches, and a digest of the kept state probabilities. |
| T5d | `T5d` | Per imputed record: cluster, reference marker, a digest of the allele probabilities, and per ALT allele the dose sums and DR2. |

A stage-1 sample line in T3b1 and T3b is `sample`, the sample index, the alleles of each haplotype at each stage-1 marker, the sizes of the sample's marker clusters, the unphased heterozygote markers, and the missing-genotype markers. A heterozygote is unphased when its phase relative to the previous heterozygote is still open. Beagle 5.4 has no cluster types.

A T3e sample line records one sample's `PhaseBaum1` update. A sample with no unphased heterozygote and no missing genotype skips the HMM and prints `sample <s> 0 0`. Other samples print these fields:

1. `sample` and the sample index.
2. The counts of unphased heterozygotes and missing genotypes before the update.
3. The cluster count and the HMM state count.
4. A digest of the state mismatches (rows 0 to 2, then each cluster, then each state), followed by the state alleles at each missing-genotype cluster.
5. A digest of the cluster recombination probabilities.
6. A digest of the emission pair `(1 - size*pMismatch, size*pMismatch)` at each HMM step. The backward pass comes first, then the forward pass.
7. The likelihood ratio of each unphased heterozygote, in marker order.
8. The count of phase switches.
9. The likelihood-ratio threshold, or `-` when none applies. None applies in burn-in or when no heterozygote is unphased.
10. The unphased heterozygote markers after the update.

The 5.5 build of these seams differed in a few fields. Its T1a printed the raw ID, REF, ALT and first `END=` text. Its T3b1 `pd` line and T3b `it` line carried the per-iteration LR threshold. Its sample lines carried `type:size` clusters, with no unphased or missing lists. It had no `leaveUnph` line and no T3e. Record class names in T1c and T1d are the 5.4 classes, for example `SeqCodedRefGTRec` for 5.5's `HapRefGTRec`.

## Piece size check

`make check-piece-size` (`tests/check-piece-size.sh`) checks the imputation writer's split of a long cluster. The writer splits a long cluster into work items of at most `PIECE_RECORDS` reference markers (500). `PIECE_RECORDS` is a tuning value that must not change the output.

The check builds `build/beagle-piece1` with one marker per work item. It runs `build/beagle-piece1` and `build/beagle` on every oracle case with `ref=` at 2 and 18 threads. Both must write the same VCF and the same trace seam T5d. The VCF alone would miss a split cluster whose pieces do not use the whole cluster's hash, because rounding to printed precision hides the difference. T5d shows it.

## Unit checks

- `make check-jcompat` compares each Java library reproduction in `src/jcompat/` against output printed by real Java (`tests/jcompat/JcompatFixtures.java`).
- `make check-interval` tests `src/vcf/interval_it.c` over an in-memory record source (`tests/vcf/interval_it_test.c`).
- `make check-records`: `tests/output/record_fixture.c` writes phased, imputed, genotyped, haploid and multiallelic records through the window writer with no BGEN and in both `bgen=` modes. `tests/check_records.py` requires the same VCF from all three runs and the expected VCF fields. It also requires phased BGEN probabilities captured before the VCF rounds them.
- `make check-tracker` tests the composite haplotype tracker in `src/beagleutil/comp_hap_queue.c` through the interface every caller uses (`tests/beagleutil/tracker_test.c`).
- `make check-bgen-unit` tests:
  - the scaling rule, exactly (`tests/bgen/quantise_test.c`)
  - the packing of values at 1 to 16 bits (`tests/bgen/pack_test.c`)
  - the reading of INFO entries behind the `.info` and `bgen-min-dr2`, including keys that prefix longer keys (`tests/bgen/info_test.c`)

## Tabix index checks

`make check-vcf-index` writes synthetic VCFs with one and four BGZF threads and requires `src/main/vcf_index.c` to write the same `.tbi` bytes as htslib's `tbx_index_build3` (`tests/output/vcf_index_test.c`). The synthetic VCFs cover:

- three chromosomes
- lines longer than a block
- lines and a header ending exactly on a block boundary
- a multi-base REF
- INFO END= variants
- a header-only file

`make check-tbi` (`tests/check-tbi.sh`) runs `build/beagle` with `tbi=true` on four oracle cases at 1 and 18 threads. The VCF must keep its oracle hash and match a run without `tbi=true`. The `.tbi` must match `tbx_index_build3`'s byte for byte.

## BGEN checks

`tests/check-bgen.sh` runs `build/beagle` with `bgen=plink2` on every case in `tests/oracle-cases.txt` and `tests/bgen-cases.txt`, with and without the filters. The cases in `tests/bgen-cases.txt` have no oracle hash and run on fixtures moved to chromosomes 22 and 38.

Some cases run again at 1, 3, 12 and 16 bits and with `bgen-chr-set=38`. The `bgen-chr-set=38` runs use the `imp` fixture moved to chromosome 38, and `imp-chroms`, whose chromosome 23 is then an autosome. The script runs plink2 with the same `bits=` on the resulting VCF and compares the `.bgen` and `.sample` with `cmp`. It checks that the cases tagged `nonautosome` fail in both tools.

It also runs every case with `bgen=phased`, and some cases with `bgen=phased` at other bit depths. `tests/check_bgen_phased.py` then checks the BGEN against the VCF for:

- the same records
- certain GT alleles where the VCF has no DS
- each haplotype's GT allele as its most probable
- dosages within one stored unit per haplotype of DS, AP1 and AP2, at any bit depth

plink2 must load the BGEN of every autosomal case.

In both modes, `tests/check_bgen_reader.py` checks that [bgen-reader](https://pypi.org/project/bgen-reader/) 4.0.9 decodes every `.bgen` to the same variants, ploidy, phased flag and probabilities as our own decoder. The live checks need [uv](https://docs.astral.sh/uv/) on `PATH` and fail without it. `tests/check_bgen_info.py` checks that the `.info` has one row per BGEN variant and that each row copies the variant's VCF record.

Failed runs must leave no `.bgen`, `.sample` or `.info`. The rule covers runs that fail when the writer opens, when the VCF cannot be opened, and after records are written. `bgen-bits` values outside 1 to 16 or without `bgen=` must fail with a message. Every such refusal must exit 1.

### Recorded BGEN hashes

`tests/bgen-hashes.txt` records the SHA-256 prefixes of the `.bgen`, `.sample` and `.info` of every run above. Every run must match its row. `BGEN_ORACLE=recorded tests/check-bgen.sh` checks these hashes instead of running plink2, bgen-reader and the Python checks, so it needs neither plink2 nor uv. It still checks each run's exit status, its VCF oracle hash and the refusals, and the C tier of the gate runs it that way. After a change to the BGEN output, rewrite the table with the live checks:

```bash
RECORD=1 PLINK2=/path/to/plink2 tests/check-bgen.sh
```

The script writes the table only when plink2, bgen-reader and our own decoder accept every file.

### Pinned plink2 build

The live checks need `PLINK2` naming the pinned plink2 binary, and fail without it. The pinned build is PLINK v2.0.0-a.7.8 (19 Sep 2026), source tag `v2.0.0-a.7.8`, commit `d293523ae31bfd10f3e8461084ee53e398705ee5`:

- macOS arm64: `plink2_mac_arm64_20260919.zip`, SHA-256 `04be4b865ad4b79a61b61ba2eb7468b4dbf62674dc5f7c94484fe0e5471539e4`.
- Linux x86_64: `plink2_linux_x86_64_20260919.zip`, SHA-256 `7ded2a083cf6863e997099d1adf5e4b86ad6b5112d2934e2b020cfd77f9bdaeb`.

## Sanitizers

`tests/check-sanitizers.sh` builds `build/beagle` and the unit tests of `make check-bgen-unit`, `make check-records` and `make check-tracker` with AddressSanitizer and UndefinedBehaviorSanitizer in `build/san`. It runs those tests. It then runs every oracle case at 1 and 2 threads as VCF only, with `bgen=plink2` and with `bgen=phased`. The `bgen=plink2` runs skip the cases tagged `nonautosome`. Every run must exit 0 with the oracle hash and no sanitizer report.

- On Linux, LeakSanitizer also runs, so any memory still allocated at a normal exit fails the check. LeakSanitizer does not support macOS arm64.
- On macOS, the script then builds `build/beagle` with ThreadSanitizer in `build/tsan`. It runs every oracle case at 18 threads, and runs the cases with per-thread hashes again with `trace=`. ThreadSanitizer runs on macOS only, because it cannot start under the x86_64 emulation of the gate's docker leg.

## Differential fuzzing

`tests/check_fuzz.py` runs differential fuzzing with [Hypothesis](https://hypothesis.readthedocs.io/).

- It generates small target VCFs. Some examples also get a reference panel, a genetic map and options (windows, iterations, states, imputation, `gp`/`ap`, and `bgen=` for the C run). It runs the release jar and `build/beagle` with the same seed. Both must fail, or both must write the same VCF. Where Java repeats without end a window that cannot advance, `build/beagle` must instead exit with its "does not advance" error.
- Invalid-parameter examples then change one parameter of a generated input so that Beagle refuses it. Both tools must exit 1 with the same message and leave the inputs unchanged. Java may prefix the message with the exception class. The change is one of:
  - `out` naming an input or a directory
  - `window` below 1.1 times `overlap`
  - a value out of bounds or not a number
  - an unknown parameter
- The script shrinks a failure to a small input and saves it in `build/fuzz-fail` with both commands. Inputs from past failures go in `tests/fuzz-regressions/`, which runs first. Each one's `expect.txt` holds the C result recorded from Beagle 5.4: an exit code and a message, or exit code 0 and the VCF hash. Beagle 5.4 finishes `window-stall`, which Beagle 5.5 repeats without end. `seq-coder-full` is a reference record whose allele count reaches the sequence limit of a 3-sample panel: Beagle 5.4 throws at 1 thread and hangs after throwing at 2.
- The full tier runs a fixed set of 200 examples (about 30 s per 100 on an M5) and 2 for each invalid-parameter change. The nightly CI run fuzzes 200 new examples. `uv run --python 3.12 --script tests/check_fuzz.py --examples 1000 --random` tries new ones.

## Model check the pipelined writer

`tests/check-tla.sh` model-checks [tla/ParallelOrdered.tla](../tla/ParallelOrdered.tla), the protocol of `parallel_ordered` (the pipelined imputed writer). It runs TLC from tla2tools.jar v1.7.4 for several worker counts, item counts and windows. It checks that:

- at most `window` items are claimed but not consumed
- no slot is overwritten before it is consumed
- items are consumed in order
- there is no deadlock
- every run finishes with every item consumed

The model represents condition-variable waits with wait sets and spurious wakeups, so a lost wakeup fails the check.

## Benchmark

`tests/bench/` holds the 1000 Genomes chr20 benchmark. `fetch-chr20.sh <dir>` downloads and derives the inputs, and `bench.sh <dir> <rounds>` times Java and C alternately. [perf-baseline.md](perf-baseline.md) has the method and the current result.

## Parity ratchet

`tests/c54-pending.txt` lists the checks of `build/beagle` whose result still differs from Beagle 5.4's, one key per line: `oracle <case>`, `failure <case>`, `log <case>`, `fuzz-regression <directory>`, `fuzz` for the generated examples, and `fuzz-invalid <change>`. `tests/c54-ratchet.sh` holds the one rule that `tests/check-oracle.sh`, `tests/check-failures.sh`, `tests/check-log.sh` and `tests/check_fuzz.py` apply to such a result:

- A key that is not listed must match, as before.
- A listed key that differs prints `pending <key>` and passes.
- A listed key that matches fails with `remove <key> from tests/c54-pending.txt`, so the list only shrinks.

`tests/check-oracle.sh` and `tests/check-failures.sh` apply the ratchet to every command except `java`. The sanitizer, tabix and BGEN checks run oracle cases too. For a listed `oracle <case>` they require exit 0 and any hash, and keep their own checks.

## Run the pre-merge gate

`tests/gate-steps.sh` holds the gate's list of checks. It needs Java 21, htslib and uv. The full tier also needs `PLINK2` naming the [pinned plink2 build](#pinned-plink2-build).

The script has two tiers. The full tier, the default, runs every check. `GATE_TIER=c` runs the C tier, which compares `build/beagle` against recorded results only. Java then only builds the bref3 fixtures. The C tier skips every check that runs Java or the jar next to the C binary (`jcompat`, `oracle-jar`, `failures-jar`, `java-build`, `oracle-source`, `java-trace`, `oracle-trace`, `log`, `trace`, `fuzz` and `trace-threads`), the fixture-cache check `cases`, the sanitizers and the TLA+ model check. Its `bgen` step checks the BGEN output against `tests/bgen-hashes.txt` instead of plink2 ([recorded hashes](#recorded-bgen-hashes)). It runs the saved fuzz regressions in `tests/fuzz-regressions/` as `fuzz-regressions`. It prints a `skip` line for each check it leaves out. `GATE_FUZZ=random` makes the full tier fuzz 200 new examples instead of the fixed 200.

`tests/check-local.sh` is the pre-merge gate. It runs the lint hooks once, then every check in `tests/gate-steps.sh` natively and on Linux x86_64 in docker. It prints one pass or fail line per check.

```bash
tests/check-local.sh
```

Two GitHub Actions workflows define the same checks:

- `.github/workflows/gate.yml` runs `tests/gate-steps.sh` on GitHub-hosted `macos-latest` (arm64) and `ubuntu-latest` (x86_64) runners. Pushes to `main` and pull requests run the C tier. A pull request runs the full tier if it changes `java/`, `src/jcompat/`, `tla/`, the workflow, or the scripts and tables the Java-side checks use. Manual runs use the full tier with random fuzz examples. The nightly run does the same, and skips itself when `main` has not moved since the last nightly run that passed.
- `.github/workflows/lint.yml` runs every hook in `.pre-commit-config.yaml` on the macOS runner.

## Run the lint hooks

`.pre-commit-config.yaml` runs actionlint, zizmor, shellcheck, ruff, typos, markdownlint and the standard file checks. It skips `third_party/`, `java/src/` and `src/jcompat/fdlibm/`, which stay as upstream wrote them. Install [prek](https://github.com/j178/prek) (`brew install prek`), then:

```bash
prek install              # run the hooks on every commit
prek run --all-files      # run them on the whole tree
```
