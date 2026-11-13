# fgkl and Lorikeet against the callers' kernels

Issue #1391.

**Verdicts:**

- **Lorikeet: leave.** Its licence is AGPL-3.0, and the project refuses GPL-family source (htsjdk-rs decision 0023). It has had no commit since 2024-01, and it does not claim byte identity.
- **fgkl: not adopted for now, blocked by a run nobody has done.** It meets the licence and provenance guard. Its own README says double precision matched GATK's Java in 0.1.0, and that 0.2.0's suffix sharing changes some last digits. The bit-for-bit comparison is ready to run. It has not run because auto mode refused to build and run third-party code on this machine. Until it runs, the port keeps its own kernels.

## The table

| | fgkl | Lorikeet |
| --- | --- | --- |
| What | PairHMM, partially determined PairHMM and Smith-Waterman kernels for GATK's native bindings: scalar, NEON, AVX2 and AVX-512 | HaplotypeCaller re-implemented in Rust |
| Source read | `fulcrumgenomics/fgkl` at b31732c, v0.2.0 | `rhysnewell/Lorikeet` at 3fbbdec (2024-01-15), after v0.6.0rc1 |
| Licence | MIT | `LICENSE` says AGPL-3.0, while `Cargo.toml` says GPL-3.0 |
| Against the guard | Acceptable. The kernel crates `fgkl-pairhmm` and `fgkl-smithwaterman` have no dependencies, and their scalar `reference` modules port GATK's Apache-2.0 Java | **Refused** (GPL family) |
| Maintenance | Active: 0.1.0, then 0.2.0, from Fulcrum Genomics | Nothing since 2024-01, still a release candidate |
| What it claims | Double precision byte-identical to GATK's Java `LoglessPairHMM` in 0.1.0. In 0.2.0, sharing changes "rare last digits" (one record in 680,000 on a centromere shard). Smith-Waterman reproduces `SmithWatermanJavaAligner` "exactly" | Speed on par with the Java. No byte claim |
| Measured against the oracle | **Not run** (see below) | Not measured, and not needed for this verdict |
| Speed against the port | Not run | — |
| Verdict | Blocked by the bit-for-bit run, then adopt or leave by its result | Leave. Ideas only, and nothing copied |

## What was measured

**The port's own kernel is bit-identical to GATK 4.7.0.0 on 20,220 pairs.** That is the target any
replacement would have to match. The Java side is `LoglessPairHMM.computeReadLikelihoodGivenHaplotypeLog10`,
called directly from a harness in GATK's own package. Its raw `doubleToRawLongBits` are compared,
because `--pair-hmm-results-file` prints `%e` and hides the last digits.

| Corpus | Pairs | Port identical to Java | Cells |
| --- | --- | --- | --- |
| Real HaplotypeCaller pairs: 3 haplotype sets × 99 reads on `svbp_reads.bam`, from the reference's own HC run | 220 | 220 | 4.2 M |
| Synthetic, seeded (`make_pairs.py`, seed 1391). Haplotypes 1–620 bp and reads 1–255 bp. 3% mismatches, 0.5% N, an indel in 1 read in 5. Base qualities 0–6, 60 or 93, or 2–41 | 20,000 | 20,000 | 305.8 M |

The port's speed, measured on an arm64 Mac with `/usr/bin/time -l` and three runs of 10 passes over
the synthetic set (3.06 G cells): 120.8–124.4 G cycles, about 40 cycles per cell, scalar. fgkl's
SIMD lanes would have to beat that figure, and only on rows where the bits are equal.

Read without running it, fgkl's scalar reference (`crates/pairhmm/src/reference.rs`, `model.rs`)
performs GATK's operations in GATK's order:

- the `2^1020 / haplotype length` initial condition;
- priors of `1 - e` and `e / 3`, and `N` matching anything;
- the match, insertion and deletion recursion;
- the last row's `m + x` summed;
- `log10` of the sum minus the initial condition's `log10`;
- the transition table built through `approximateLog10SumLog10` and `log1p`.

That fits the 0.1.0 claim, but it is not a measurement.

## What is left, and how to run it

1. **Allow building and running fgkl's source on this machine, or run it in CI.** Auto mode refused
   `cargo build` on a probe crate that depends on the clone, because that runs third-party code.
   Nothing here goes around the refusal.
2. Run `probe check pairs.txt bits.txt synth.txt synth_bits.txt`. For each backend (scalar, NEON
   here; AVX2 and AVX-512 on x86-64 CI), in double and float, one pair per call, it prints how many
   results are identical to the oracle's bits and the worst ULP distance. It does the same for 16×16
   batches with `share_suffixes` off and on.
3. Run `probe sw smith_waterman.txt.gz`: every row of the Smith-Waterman golden through fgkl, which
   must produce the same offset and CIGAR.
4. **The rule from the issue: a kernel that differs by one ULP is not adopted.** Two configurations
   are already excluded:
   - `share_suffixes = true`, which by fgkl's own account moves a last digit;
   - single precision, which differs from GKL's output on every shard fgkl compared.

   That leaves double precision without sharing. Only if it is identical on all 20,220 pairs and
   every SW row does the speed bench (`probe bench`, with the same alternating `/usr/bin/time` runs)
   decide.

The probe crate, the Java harness (`PairHmmBits.java`), the real and synthetic corpora with their
Java bits, and the runner scripts are in the session-notes folder `fgkl-1391/`. Outside the
repository, `synth.txt` is regenerated byte for byte by `make_pairs.py synth.txt`.

## Why Lorikeet is left

- The project's licence guard refuses GPL-family source, so no Lorikeet code can be copied or linked.
- It reimplements the HaplotypeCaller algorithm and does not port its bytes. Its output would be a
  third opinion, not a reference, and the oracle already gives the first one.
- Its PairHMM is the `gkl` crate 0.1.1 (philipc/gkl-rs, Apache-2.0), a binding to Intel GKL. GKL is
  the vectorised kernel #251 already compares with `LoglessPairHMM`, so it adds nothing new to
  measure.
- As a source of test cases its value is small. The port's corpus already holds the reference's own
  HaplotypeCaller rows.
