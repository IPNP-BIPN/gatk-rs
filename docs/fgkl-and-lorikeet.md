# fgkl and Lorikeet against the callers' kernels

Issue #1391.

**Verdicts:**

- **fgkl Smith-Waterman: adopt.** Its aligner gives the golden's offset and CIGAR on all 704
  alignments, on both scalar and NEON. It agrees with the port on 32,000 synthetic alignments, and
  it is **14× faster**.
- **fgkl PairHMM: don't adopt. Speed up the port's own kernel instead.** On x86-64 Linux, where CI
  and users run, no configuration is bit-identical to the oracle. Its best one, scalar double with
  one haplotype per call, differs on 7 of 20,000 synthetic pairs (worst 256 ULP). Every one of the
  7 comes from its final `f64::log10` (glibc), where the port uses Java's (`jmath`). On macOS, a table
  fix made it exact on all 20,220 pairs, but only because Apple's `log10` happens to agree with
  Java's there; on Linux that fix changes nothing. Upstream closed the report
  (fulcrumgenomics/fgkl#25) and the fix (#26) as won't-fix: bit-exactness with a JVM is not its
  goal, and its tolerance is 10⁻⁹ in log10. Where its bits were equal, it ran **5.3×** to **16×**
  faster than the port without SIMD, so that speed is within reach of a kernel that keeps the
  port's operations.
- **Lorikeet: leave.** It is AGPL-3.0, has had no commit since 2024-01, and makes no byte claim.

The rule from the issue: a kernel that differs by one ULP is not adopted. Everything recommended
above differs by 0 ULP on every input measured.

## The table

| | fgkl | Lorikeet |
| --- | --- | --- |
| What | PairHMM, partially determined PairHMM and Smith-Waterman kernels for GATK's native bindings: scalar, NEON, AVX2, AVX-512 | HaplotypeCaller re-implemented in Rust |
| Source read | `fulcrumgenomics/fgkl` at b31732c, v0.2.0 | `rhysnewell/Lorikeet` at 3fbbdec (2024-01-15), after v0.6.0rc1 |
| Licence | MIT | `LICENSE` says AGPL-3.0, while `Cargo.toml` says GPL-3.0 |
| Against the guard | Acceptable. `fgkl-pairhmm` and `fgkl-smithwaterman` have no dependencies, and their scalar `reference` modules port GATK's Apache-2.0 Java | **Refused** (GPL family) |
| Maintenance | Active: 0.1.0, then 0.2.0, from Fulcrum Genomics | Nothing since 2024-01, still a release candidate |
| Its claim | Double precision byte-identical to `LoglessPairHMM` in 0.1.0; sharing changes rare last digits. Smith-Waterman reproduces `SmithWatermanJavaAligner` exactly | Speed on par with the Java. No byte claim |
| PairHMM against the oracle | See below: never 0 ULP on x86-64 Linux (its final `log10`); on macOS 0 ULP only with scalar double, a table fix, and one haplotype per call | Not measured |
| Smith-Waterman against the oracle | 704/704 alignments of the golden, plus 32,000/32,000 synthetic ones against the port | Uses the `gkl` crate (see below) |
| Speed where equal | PairHMM 5.3× and 16× (macOS only); Smith-Waterman 14× | — |
| Verdict | Smith-Waterman: adopt. PairHMM: don't; speed up the port's kernel | Leave |

## How it was measured

The oracle's answers come from the pinned GATK 4.7.0.0. A harness in GATK's own package calls
`LoglessPairHMM.computeReadLikelihoodGivenHaplotypeLog10` and prints its
`doubleToRawLongBits`, because `--pair-hmm-results-file` prints `%e` and hides the last digits.
There are two corpora:

- **Real:** 220 pairs from 3 haplotype sets × 99 reads that the reference's own HaplotypeCaller
  formed on `svbp_reads.bam`.
- **Synthetic:** 20,000 pairs from `make_pairs.py` with seed 1391. Haplotypes are 1–620 bp and
  reads 1–255 bp. There are 3% mismatches, 0.5% `N`, and an indel in 1 read in 5. Base qualities are
  0–6, 60 or 93, or 2–41.

The port is bit-identical to the oracle on all 20,220 pairs, so batched results, which have no
oracle answer of their own, are compared with the port. A read-haplotype pair's answer in GATK does
not depend on the batch it came in.

Smith-Waterman is compared with the port's `smith_waterman.txt.gz` golden:

- 46 pairs × 4 parameter sets × 4 overhang strategies;
- the 2,000 first synthetic pairs, aligned as haplotype against read, with every parameter set and
  strategy, against the port.

The probe crate lives outside the repository (session notes, `fgkl-1391/`). Every timing is three
alternating runs of `/usr/bin/time -l` on an Apple M4 Max, as cycles elapsed. Wall time was noisy
because other sessions shared the machine; cycles were not.

## PairHMM: the measurements

### On x86-64 Linux: the verdict

Measured for #1433, in a `linux/amd64` container against the same oracle bits (x86-64 Java):

| Configuration | Real (220) | Synthetic (20,000) |
| --- | --- | --- |
| The port | 220, 0 ULP | 20,000, 0 ULP |
| fgkl 0.2.0, scalar double, one haplotype per call | 220 | 19,993 (7 differ, worst 256 ULP) |
| The same with the table fix below | 220 | 19,993 (the same 7) |

On Linux, `10f64.powf(x)` already compiles to `pow`, so the table fix changes nothing. The 7 are
the last step: the port's own recursion with only its final `log10` swapped for `f64::log10`
reproduces exactly those 7. The port uses `jmath::math::log10`, a port of Java's.

Upstream's answer on fulcrumgenomics/fgkl#25: every difference is within its documented tolerance
(10⁻⁹ in log10), HotSpot's `Math.pow` and `Math.log10` give different bits on x86-64 and on
aarch64, and matching one JVM's bits is not a goal. The oracle, and every golden, is x86-64 Java,
so that is the JVM the port has to match.

### On macOS (aarch64): what led to the first verdict

The rest of this section was measured on an Apple M4 Max. It explains the three causes, but its
0 ULP result does not hold on Linux.

As published (0.2.0), against the oracle's bits:

| Configuration | Real (220) | Synthetic (20,000) |
| --- | --- | --- |
| The port | 220, 0 ULP | 20,000, 0 ULP |
| fgkl `reference` (its scalar port of the Java) | 220 | 19,962 (38 differ, worst 128 ULP) |
| fgkl scalar double, one pair per call | 220 | 19,962 (38 differ, worst 128 ULP) |
| fgkl NEON double, one pair per call | 219 (worst 8 ULP) | 19,918 (82 differ, worst 128 ULP) |
| fgkl scalar or NEON float | 30 | about 7,650 (worst 5.5·10¹⁰ ULP) |

### Three causes, each isolated

1. **The tables: `10f64.powf(x)`.** `model.rs` builds the error-probability and match-to-match
   tables with a literal base of 10, and the compiler does not compute that as `pow(10, x)`. The
   port builds its tables with `std::hint::black_box(10f64).powf(x)` for this reason. Swapping the
   final `log10` for Java's (`jmath`) changes nothing: still 38 differ. Adding `black_box` to the
   four `powf` calls in a scratch copy makes fgkl's reference 20,220/20,220 identical on macOS. On
   Linux the compiler already calls `pow`, and the final `log10` differs instead (see above).

2. **Fused multiply-add in every SIMD backend.** Vector `mul_add` is `vfmaq_f64` on NEON and
   `_mm*_fmadd_pd` on AVX2 and AVX-512. It rounds once where Java's `a * b + c` rounds twice. The
   scalar backend does not fuse. With the tables fixed:

   | With the table fix | Real (220) | Synthetic (20,000) |
   | --- | --- | --- |
   | fgkl `reference` | 220 | 20,000 |
   | fgkl scalar double, one pair per call | 220 | 20,000 |
   | fgkl NEON double, one pair per call | 219 (worst 8 ULP) | 19,933 (67 differ, worst 16 ULP) |

   NEON, AVX2 and AVX-512 double are therefore out under the rule, and so is float.

3. **The prefix shared between haplotypes.** A call with several haplotypes always shares forward
   columns over their common prefix. `Config.share_suffixes` controls only the suffix half.
   `kernel.rs:641` rescales the shared state by `src.hap_len / hap_len`. GATK seeds each haplotype
   with `2^1020 / its own length`, so the rescale rounds. With the tables fixed and sharing off:
   - every read of a chunk against one haplotype differs on 0 of 6,400 pairs;
   - one read against 16 haplotypes differs on 4 of 6,400, each by 1 ULP, and each is exact when
     computed alone;
   - 16 × 16 batches differ on 180 of 320,000 (worst 512 ULP), and on 211 with sharing on.

   **The exact shape on macOS is one haplotype per call, with all of a region's reads in it**:
   **320,000/320,000** synthetic and **3,472/3,472** real batched pairs identical to the port,
   0 ULP.

### Speed, where the bits are equal (macOS)

| Shape | The port | fgkl scalar double (tables fixed) | Ratio |
| --- | --- | --- | --- |
| One pair per call, synthetic set × 3 (917 M cells) | 36.3–37.4 G cycles | 6.91–7.01 G | **5.3×** |
| Reads of each chunk of 16 against one haplotype per call, 2,000 pairs (32,000 likelihoods) | 18.8–19.1 G | 1.15–1.19 G | **16×** |

The port runs at about 40 cycles per cell. fgkl's scalar kernel does the same operations in the
same order, apart from its last `log10`, so 5–16× is within reach of the port's own kernel without
SIMD and without losing a bit: no matrices allocated per pair, and the tables built once per read.

## Smith-Waterman: the measurements

| Input | Result |
| --- | --- |
| The golden: 46 pairs × 4 parameter sets × 4 strategies, scalar and NEON | 704/704 alignments identical (offset and CIGAR). The 32 empty-reference and empty-alternate rows are refused by both sides; GATK throws `IllegalArgumentException` and fgkl returns `sequences must be non-empty` |
| 2,000 synthetic pairs × 16 settings (32,000 alignments), against the port | 32,000/32,000 identical. Both sides wrote the same 558,874 characters |
| Speed on those 32,000 alignments | The port 9.35–9.43 G cycles; fgkl (NEON) 0.67 G: **14×** |

Smith-Waterman is integer arithmetic, so FMA and the tables do not touch it.

## What adopting would take

Smith-Waterman only. It is integer arithmetic, so the platform's libm does not touch it. #1433
adds `fgkl-smithwaterman` at a pinned commit, the way bwa-mem4 is pinned. `Aligner::new` picks the
widest backend the CPU has, so x86-64 CI runs AVX2 or AVX-512 against the goldens and covering
arrays of HaplotypeCaller and Mutect2.

PairHMM stays the port's. Making the port's kernel faster is separate work.

## Why Lorikeet is left

- The project's licence guard refuses GPL-family source.
- It reimplements the HaplotypeCaller algorithm and does not port its bytes.
- Its PairHMM and Smith-Waterman are the `gkl` crate 0.1.1 (philipc/gkl-rs, Apache-2.0), a binding
  to Intel GKL's AVX C kernels. Those are the vectorised kernels #251 compares with
  `LoglessPairHMM`, so there is nothing new to measure.
- The port's corpus already holds the reference's own HaplotypeCaller rows, so Lorikeet adds little
  as a source of test cases.
