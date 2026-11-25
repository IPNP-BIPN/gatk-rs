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

## 6. The four fixture losers, again (W12 §3b)

**Ratio.** After entry 1, on x86-64 (`tools/speed/current.json`): `FlowPairHMMAlignReadsToHaplotypes`
2.83 steady, `ModelSegments` 1.51, `VariantRecalibrator` 1.01, `HaplotypeBasedVariantRecaller` 0.97;
peak RSS 0.04x to 0.09x the reference's on all four, so only CPU is at stake.

**Where the time goes, and the optimisation.** Profiles by macOS `sample`, on longer runs of the
same rows where the fixture is too short to sample (more MCMC samples, or the input repeated):

| Tool | profile | change |
|---|---|---|
| `FlowPairHMMAlignReadsToHaplotypes` | `jmath::math::log10` per cell, on a few hundred distinct inputs | a per-thread direct-mapped cache of `log10` by the input's bits |
| `ModelSegments` | `hetLogLikelihood`'s logarithms and log-gammas, a dozen per het per sampler evaluation | three per-thread caches keyed by their inputs' bits: the whole answer (a segment's hets repeat counts), the terms the outlier probability does not enter (its sampler moves `pi` alone), the outlier term's count-only parts |
| `VariantRecalibrator` | every resource scanned whole per variant; then `Math.log(10)` per datum per Gaussian, an array per evaluation, a temporary matrix per scatter term | the resources indexed by contig and start, in file order; `log(10)` once; the cross terms and scatter terms summed in place |
| `HaplotypeBasedVariantRecaller` | the one-read PairHMM | the engine's batched recursion, several reads side by side |

**Why it cannot move a byte.** Each cache holds the value of a pure function of exactly the bits it
is keyed by, so a hit returns the bits the computation returns, and a collision only recomputes.
The `ModelSegments` terms are cached as the operands they were, so every sum keeps its order.
The resource index returns the scan's records in the scan's order, filtered by the same
predicate. The EM's in-place sums are the same multiplies and adds in the same order (a term
added as formed is what `plusEquals` of the term matrix does element by element). The batched
PairHMM does each lane's own unfused operations.

**The suites that prove it.** `pair_hmm_port` (batched against one-read, bit for bit); new unit
tests: the flow `log10` cache against `jmath`, the `ModelSegments` likelihood against the
uncached body kept under `cfg(test)` over more keys than slots (a planted reordering fails it),
the resource index against the scan. The tools' suites, and locally every output file of the four
fixture rows, and of `ModelSegments` at 2000 samples, identical to main's. `ModelSegments`' sampler
is otherwise covered only by the x86 conformance job: its outputs differ from an **arm64** JVM's
on main already (its `Math.exp` is not x86's), so a local oracle run cannot vouch for it.

**Result, natively (M-series, one binary each side, minimum CPU of five, alternated twice,
identical outputs).** x86-64 pending: the `Speed` run on this branch.

| Tool | main | this change | |
|---|---:|---:|---:|
| `FlowPairHMMAlignReadsToHaplotypes` | 169 ms | 21 ms | 8.0x |
| `ModelSegments` | 490 ms | 216 ms | 2.3x |
| `ModelSegments`, 2000 samples | 10.04 s | 4.16 s | 2.4x |
| `VariantRecalibrator` | 110 ms | 56 ms | 2.0x |
| `HaplotypeBasedVariantRecaller` | 20 ms (340M instructions) | 16.5 ms (242M) | 1.2x |
## 7. The allocator, measured against glibc (W12, P0)

**Profile.** Entry 3 (#1541) measured mimalloc against Apple's allocator only: 3 to 13% less CPU,
and PrintReads' peak RSS 8% higher. Its own criterion is against glibc, which is what the shipped
Linux binary links.

**The optimisation.** The same `#[global_allocator]` line in `crates/gatk-cli/src/main.rs`.

**Why it cannot move a byte.** Entry 3's argument: an allocator chooses addresses, and no address
reaches an output the suites check.

**The suites that prove it.** Entry 3's. The decision rule, written before measuring: kept only if,
on Linux against glibc, no measured tool's peak RSS rises and the CPU of each falls or holds.

## 8. Profile-guided optimisation of the shipped binary (W12, P0)

**Profile.** The binary's hot paths (the PairHMM, BGZF, the VCF and BAM codecs, the read
walkers) are a small part of a 28 MB binary built without knowing which branches are taken.

**The optimisation.** `-Cprofile-generate`, a training run over the speed rows and the chr20 smoke
slice, `llvm-profdata merge`, `-Cprofile-use`. Only the build changes.

**Why it cannot move a byte.** A profile steers inlining, block layout and branch weights, never the
meaning of an operation. Rust has no fast-math, so LLVM may not reassociate, contract or otherwise
change a floating-point result whatever it inlines; integer and memory semantics are the language's.
What could move bytes is a miscompile, and the suites run the PGO binary to catch one, as they would
any compiler change.

**The suites that prove it.** Every end-to-end test and covering array, run on the PGO binary, and
the smoke slice's outputs compared with the non-PGO binary's. Kept only if CPU falls with no RSS rise.

## Not targets yet

Superseded by entry 6 for these two, kept for the record. `HaplotypeBasedVariantRecaller` (1.30
steady) and `FlowPairHMMAlignReadsToHaplotypes` (1.05) are
at or near one, and the baseline's own noise says a single tool can move by up to half between two
runs on different runners. Neither is "slow" by a margin that survives that. Their rows take 60 to
90 ms natively, too short for a sampling profiler to attribute, so there is also no profile to
choose an optimisation from. They re-enter this list when a larger input, still covered by a
golden, gives them a run long enough to profile and a ratio that clears the noise.

## Out of scope, as S.2 says

Anything that would reorder a reduction, contract an FMA or swap a `jmath` call for the host's,
however tempting a profile makes it look. Target 1 is admissible precisely because it changes how
the correctly rounded value is found and never which value it is.
