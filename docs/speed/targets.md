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
