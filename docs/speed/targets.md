# Targets, chosen from the baseline

Milestone S, S.5 (#112). A path is on this list because the baseline ([baseline.md](baseline.md))
says it is slow, not because it looks slow. Each entry carries its measured ratio, where the time
goes according to a profile, the optimisation and the argument for why it cannot move a byte
(written before the change, not after), and the suite that proves it did not.

The baseline names four tools whose steady ratio reaches 1. Two are targets; two are not yet.

## 1. `jmath`'s correctly rounded logarithm (`ModelSegments`, `VariantRecalibrator`)

**Ratio.** `ModelSegments` 18.4 steady, 4.83 cold (the only tool slower than the reference even
with the JVM's start-up charged to it); `VariantRecalibrator` 4.63 steady.

**Where the time goes.** A native profile of `ModelSegments`' baseline row (macOS `sample`): 94% of
samples in `jmath::log::ln_dd`, the 22-term double-double series behind `jmath::math::log`;
`VariantRecalibrator`'s row is dominated by the same function. Both call `Math.log` where the
reference has a hardware intrinsic, a dozen times per het per sampler evaluation in the first
case, per Gaussian per variant in the second.

**The optimisation.** In htsjdk-rs (IPNP-BIPN/htsjdk-rs#243, decision 0045 there): a fast phase,
a table, an exact reduction and a short polynomial with an error bound, answering only when
twice the bound rounds to one `f64`, and falling back to the series otherwise.

**Why it cannot move a byte.** `Math.log` is correctly rounded (decision 0006), so it has one
answer per input; the fast phase returns an answer only where it is provably that one, and the
series, unchanged, decides everything else. Checked on 200 million inputs (99.94% decided fast,
none differing) and on two million more on every CI run.

**The suites that prove it.** The `jmath` suite and its hard-to-round corpus in htsjdk-rs; here,
every conformance suite and covering array on the bump's PR, `model-segments` and
`VariantRecalibrator`'s array among them, all green; locally, the four tools' baseline rows give
identical outputs before and after.

**Result, natively (M-series, median of three, identical outputs).**

| Tool | before | after | |
|---|---:|---:|---:|
| `ModelSegments` | 5.19 s | 0.46 s | 11.3x |
| `VariantRecalibrator` | 0.358 s | 0.114 s | 3.1x |
| `HaplotypeBasedVariantRecaller` | 68 ms | 59 ms | 1.15x |
| `FlowPairHMMAlignReadsToHaplotypes` | 87 ms | 87 ms | 1.0x |

**Result on real x86-64** (`Speed` run 37923376066 on this change, `tools/speed/current.json`,
now STATUS.md's cost column; the baseline stays in `baseline.json`):

| Tool | steady, before | steady, after | cold, before | cold, after | port, before | port, after |
|---|---:|---:|---:|---:|---:|---:|
| `ModelSegments` | 18.39 | **1.50** | 4.83 | **0.38** | 9.67 s | 1.26 s |
| `VariantRecalibrator` | 4.63 | **0.75** | 0.37 | 0.06 | 1.00 s | 0.15 s |
| `HaplotypeBasedVariantRecaller` | 1.30 | 1.05 | 0.07 | 0.06 | 158 ms | 125 ms |
| `FlowPairHMMAlignReadsToHaplotypes` | 1.05 | 0.98 | 0.09 | 0.09 | 120 ms | 183 ms |

`ModelSegments` is now cheaper than the reference cold and within a factor of 1.5 of it warm;
`VariantRecalibrator` is cheaper than a warm JVM. The last two rows moved inside the noise, in
both directions, which is what an unchanged path does. The other 131 tools moved a median 7% cold
and 10% steady, at most 45%, the same spread two runs of one commit showed in the baseline.

## 2. The shipped binary's profile, `dist` (W12, P0)

**Profile.** Not one tool: the whole binary. The default `release` profile compiles each crate in
sixteen units and stops inlining at crate boundaries, so a hot call from `gatk-tools` into
`gatk-engine`, `htsjdk-*` or `jmath` is a real call. `gatk-tools` alone is 488,000 lines in one crate.

**The optimisation.** `[profile.dist]`: `inherits = "release"`, `lto = "fat"`, `codegen-units = 1`.
Every job that runs the port against the reference builds with it; `release` stays the test
profile.

**Why it cannot move a byte.** `byte-neutrality.md`, "What is admitted under classification": no
fast-math flag exists to be propagated, and Rust never contracts, so wider inlining changes which
calls are made, not which operations. The one way it could move a byte, a host transcendental
folded at compile time once its argument becomes a constant across crates, is exactly what the
arrays check, because they run the `dist` binary.

**The suites that prove it.** Every covering array, `thread_independence.py` and the end-to-end
`Speed` rows, all on `target/dist/gatk`; the in-process suites unchanged on `release`.

**Result.** Section "P0 and P3, measured" below.

## 3. A global allocator (W12, P0)

**Profile.** Every streaming tool allocates per record (a decoded BAM record, a VCF line's
fields); glibc's malloc serialises its arenas under threads and keeps freed pages per arena.

**The optimisation.** `mimalloc` as `#[global_allocator]` in `crates/gatk-cli/src/main.rs` only.

**Why it cannot move a byte.** An allocator chooses addresses, and an address reaches no output
the suites check (ASLR already moves every address on every run; the one known case is masked by
layout, `pointers-that-reach-the-output.md`). The relative order of allocations is what is new,
and the `determinism` job and every end-to-end test run the binary, so its allocator, under a
second environment.

**The suites that prove it.** Every end-to-end test and covering array, and `determinism`. Kept only
if peak RSS does not rise against glibc on any measured tool, since the gate forbids trading RSS
for CPU.

**Result: admitted, not landed by W12.** On the macOS measurement below it saves 3 to 13% of CPU on
top of `dist`, but PrintReads' peak RSS rises 8% (2,478 MB against 2,294 MB), which is the condition
above failing on the one host measured. That host's baseline is Apple's allocator, not glibc, so it
does not settle Linux. #1474 carries the same `#[global_allocator]` line; it is decided there, on
`Speed at scale`'s Linux rows, against this criterion. The guard admits the line in `main.rs`
either way.

## 4. PairHMM across pairs, dispatched at run time (W12, P3)

**Profile.** HaplotypeCaller and Mutect2 spend most of their compute in
`pair_hmm::reads_likelihoods_given_haplotype_log10`, which already computes four reads side by
side against one haplotype. Without `target-cpu` (refused), x86-64 compiles those four lanes to
SSE2, two per register; GKL, which GATK calls, uses AVX-512 or AVX2.

**The optimisation.** The same body, generic over the lane count, compiled three times: portable
(four lanes, SSE2 or NEON), `#[target_feature(enable = "avx2")]` (four lanes, one 256-bit register),
`#[target_feature(enable = "avx512f")]` (eight lanes, one 512-bit register), chosen once by
`is_x86_feature_detected!`.

**Why it cannot move a byte.** One read per lane, each lane doing the single pair's multiplies and
adds in `LoglessPairHMM`'s order; a lane-wise IEEE 754 multiply or add is the scalar one, bit for
bit, and nothing crosses lanes (the final sum stays scalar, per read, in the reference's order).
`fma` is not enabled and the guard refuses it in this file. Identity holds whether or not LLVM
vectorises the body, so it does not rest on the autovectoriser.

**The suites that prove it.** `reads_side_by_side_give_each_reads_own_bits` forces every kernel
the host has against `read_likelihood_given_haplotype_log10`, which `pair_hmm_port` holds to the
reference's transcription; then every HaplotypeCaller and Mutect2 suite and array.

## P0 and P3, measured (W12)

Native Apple M-series (aarch64, 16 cores, shared with other jobs), toolchain 1.97.1. Input: the
chr20 smoke slice (NA12878 30x, chr20:30-32 Mb, `tools/speed/scale.json`). Three rounds,
alternating binaries; CPU is user + sys, median; RSS is the peak over the rounds. Every output was
the same across all binaries and rounds, once the `##GATKCommandLine` line (it carries the run's
date) was left out.

| Tool | `release` CPU | `dist` | `dist` + mimalloc | RSS `dist` | RSS + mimalloc |
|---|---:|---:|---:|---:|---:|
| HaplotypeCaller, 300 kb | 11.08 s | 0.98x | 0.88x | 1.00x | 0.96x |
| Mutect2, 100 kb, `--ignore-itr-artifacts` | 2.85 s | 0.96x | 0.84x | 1.00x | 0.97x |
| PrintReads, whole slice | 5.97 s | 1.01x | 0.97x | 1.00x | **1.08x** |
| SelectVariants, gnomAD chr20, SNPs | 29.19 s | 0.88x | 0.82x | 1.00x | 1.01x |

`dist` costs build time: 704 s against 164 s for `release` on this host, with a 3.9 GB peak (it
was OOM-killed in a 7.75 GB Docker VM that other jobs shared). It makes the binary smaller,
22.1 MB against 27.9 MB. Wall time on this shared host moved by up to 2x between rounds of the same
binary, so it is not reported.

**P3 on aarch64.** The portable kernel is the one that was already there, so as expected it is
unchanged: HaplotypeCaller 10.96 s against 11.07 s (median of three, identical output). A sampled
HaplotypeCaller run puts 78% of its busy samples in
`reads_likelihoods_given_haplotype_log10_with`, which bounds what the x86-64 kernels can win.
**x86-64 is pending.** Rosetta exposes neither AVX2 nor AVX-512, so the AVX2 and AVX-512 kernels
have been checked only by their assembly (8 `vmulpd` and 4 `vaddpd` per cell, `ymm` and `zmm`
respectively) until CI's test log reports which kernels each runner held to the single pair and
`Speed`/`Speed at scale` time them.

**Found while measuring, owned elsewhere.** Without `--ignore-itr-artifacts`, Mutect2 on 100 kb
ran for more than ten minutes with nearly all its CPU in `memmove`. The runner's `contig_bases`
clones the whole contig (64 MB for chr20) for every read the palindrome clip sees. #1474 shares the
bases through an `Rc`.

## Not targets yet

`HaplotypeBasedVariantRecaller` (1.30 steady) and `FlowPairHMMAlignReadsToHaplotypes` (1.05) are
at or near one, and the baseline's own noise says a single tool can move by up to half between two
runs on different runners. Neither is "slow" by a margin that survives that. Their rows take 60 to
90 ms natively, too short for a sampling profiler to attribute, so there is also no profile to
choose an optimisation from. They re-enter this list when a larger input, still covered by a
golden, gives them a run long enough to profile and a ratio that clears the noise.

## Out of scope, as S.2 says

Anything that would reorder a reduction, contract an FMA or swap a `jmath` call for the host's,
however tempting a profile makes it look. Target 1 is admissible precisely because it changes how
the correctly rounded value is found and never which value it is.
