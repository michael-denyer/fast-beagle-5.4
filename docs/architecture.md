# Architecture and source layout

fast-beagle is a file-by-file C port of Beagle 5.4 (29Oct24). Each ported source file names the Beagle 5.4 Java file whose behaviour it implements. The [checks](testing.md) hold every file's behaviour to Beagle 5.4's. [CONTEXT.md](../CONTEXT.md) defines the domain terms.

## Each window runs through four stages

```mermaid
graph LR
    IN["gt=, ref=, map="] --> READ["vcf/, bref/<br/>read records, cut windows"]
    READ --> PHASE["phase/<br/>phase the window"]
    PHASE -->|"reference markers<br/>the target lacks"| IMP["imp/<br/>impute"]
    PHASE -->|"nothing to impute"| WRITE["main/window_writer.c"]
    IMP --> WRITE
    WRITE --> VCF["out.vcf.gz"]
    WRITE --> BGEN["bgen/<br/>out.bgen, .sample, .info"]

    classDef input fill:#1565C0,stroke:#0D47A1,stroke-width:3px,color:#FFFFFF,font-weight:bold
    classDef read fill:#00838F,stroke:#004D40,stroke-width:3px,color:#FFFFFF,font-weight:bold
    classDef phase fill:#6A1B9A,stroke:#38006B,stroke-width:3px,color:#FFFFFF,font-weight:bold
    classDef impute fill:#E65100,stroke:#8C2F00,stroke-width:3px,color:#FFFFFF,font-weight:bold
    classDef write fill:#C62828,stroke:#7F0000,stroke-width:3px,color:#FFFFFF,font-weight:bold
    classDef output fill:#2E7D32,stroke:#1B5E20,stroke-width:3px,color:#FFFFFF,font-weight:bold
    class IN input
    class READ read
    class PHASE phase
    class IMP impute
    class WRITE write
    class VCF,BGEN output
    linkStyle default stroke:#9E9E9E,stroke-width:2px
```

`src/main/main.c` reads the parameters (`src/main/par.c`) and runs each window through these stages in turn. `src/blbutil/`, `src/beagleutil/` and `src/jcompat/` hold the utilities that every stage uses.

The port reads both input files and splits them into overlapping windows along the genetic map. It phases every window in four parts: fixed phasing inputs, the initial PBWT phase, the stage-1 iterations and stage 2. It then imputes the reference markers the target lacks. It writes Beagle's VCF with DS, DR2, AF and optional AP1/AP2 and GP fields.

## Threaded loops

The port threads the loops that dominate run time:

- reference record parsing per block
- the initial PBWT phase per sub-window and per block of samples
- step coding per step
- stage-1 phasing and parameter estimation per sample
- stage 2 per sample
- the IBS neighbour searches per batch of steps
- imputation per haplotype
- the imputed writer, which builds records on worker threads and prints them in order

`nthreads=` sets both the thread count and the partitions that Beagle's output depends on, so the output matches Beagle run with the same thread count. [tla/ParallelOrdered.tla](../tla/ParallelOrdered.tla) models the protocol of `parallel_ordered`, the pipelined imputed writer ([model check](testing.md#model-check-the-pipelined-writer)).

The imputation writer splits a long cluster into work items of at most `PIECE_RECORDS` reference markers (500). `PIECE_RECORDS` is a tuning value that must not change the output ([piece size check](testing.md#piece-size-check)).

## Source layout

- `src/`: the C port. Its directories mirror the Java packages (`blbutil`, `beagleutil`, `vcf`, `main`, …), and each file names the Java file it was ported from.
- `src/jcompat/`: C reproductions of the Java library behaviour Beagle depends on. `make check-jcompat` compares each against output printed by real Java (`tests/jcompat/JcompatFixtures.java`).
- `src/vcf/interval_it.c`: the `chrom=` interval that both the target and the reference records pass through. It skips records before the interval, and the first record after the interval ends the stream. It reads the source at most one record past the interval. With no interval, every record passes through.
- `src/beagleutil/comp_hap_queue.c`: the composite haplotype tracker.
- `src/bgen/`: the BGEN writer for `bgen=plink2` and `bgen=phased` ([BGEN output](usage.md#bgen-output)).
- `src/main/vcf_index.c`: the tabix index writer for `tbi=true`.
- `third_party/libdeflate/`: the files of libdeflate 1.25 (MIT, `COPYING`) that zlib compression needs. They come from plink2's source tree at tag v2.0.0-a.7.8. plink2 compresses BGEN with this library at level 6. zlib's own deflate gives different bytes.
- `java/src/`: the unmodified Beagle 5.4 Java source, used as the oracle.
- `java/trace.patch`: trace hooks for the Java source ([trace seams](testing.md#compare-trace-seams)).
- `tests/`: the checks ([checks and the pre-merge gate](testing.md)).
- `tla/`: the TLA+ model of the pipelined imputed writer.

## Tabix index in the same pass

`src/main/vcf_index.c` records each line's chromosome, interval and uncompressed end as the VCF is written. It then reads the BGZF block headers of the closed file to place those ends, so it works with the multithreaded writer.
