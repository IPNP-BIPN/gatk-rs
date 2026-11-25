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

## Scale J. The Spark read tools held every read (`PrintReadsSpark`, `ExtractOriginalAlignmentRecordsByNameSpark`, `FlagStatSpark`, the four Spark metrics tools)

**Ratio.** At scale (chr20 30x BAM, linux/arm64, 6 GiB container), each one killed at the limit:

| Tool | GATK CPU | GATK RSS |
|---|---:|---:|
| `PrintReadsSpark` | 71 s | 2.22 GiB |
| `ExtractOriginalAlignmentRecordsByNameSpark` | 27 s | 2.19 GiB |
| `FlagStatSpark` | 28 s | 1.59 GiB |
| `QualityScoreDistributionSpark` | 34 s | 1.58 GiB |
| `CollectQualityYieldMetricsSpark` | 28 s | 1.49 GiB |
| `MeanQualityByCycleSpark` | 24 s | 1.45 GiB |
| `CollectBaseDistributionByCycleSpark` | 31 s | 1.29 GiB |

Natively on macOS, with no limit, main holds 54 GiB (`PrintReadsSpark`) and 35 to 36 GiB (the
others).

**Where the memory goes.** `read_walker::traverse` collects the whole traversal before the tool
looks at a read. #1473 streams the writer, but its callers still hand it that list.

**The optimisation.** `traverse_iter`, the same traversal a read at a time.
- `FlagStatSpark` and the metrics tools (`spark_metric_text`) fold the reads as they come. The
  `spark_metrics` folds now take any `Borrow<BamRecord>`.
- `PrintReadsSpark` and `ExtractOriginalAlignmentRecordsByNameSpark` hand #1473's writer one
  partition as the reads come (`write_spark_streamed`), which is what `parts_of` makes of the list.
  A sharded output writes its parts as it goes, and a CRAM refuses its version only after the
  traversal, so both keep the list.

**Why it cannot move a byte.** The same reads reach the same fold or the same writer in the same
order. A traversal refusal is still the first one, and leaves no output, because the writer
renames its scratch file only at the end. A refusal of the writer's own first reads the rest of
the traversal, so a later traversal refusal still wins, as when the list was read whole. The
metrics' `--chart` refusal still comes after the traversal.

**The suites that prove it.** `spark-walkers` and the tools' covering arrays. Natively, old against
new: on chr20 the outputs are identical, and on a BAM with a malformed record in the middle every
tool gives the same exit and the same message and leaves no output. The same holds for an output
in a missing directory, with a good input and with a bad one.

**Result.** macOS native, chr20, this branch (with #1473 and #1482 under it) against main:

| Tool | main RSS | RSS | main CPU | CPU | identical |
|---|---:|---:|---:|---:|---|
| `PrintReadsSpark` (`.bai`, `.sbi`) | 53.8 GiB | 40 MB | 96 s | 43 s | `.bam` `8dafdf0c`, `.bai`, `.sbi` |
| `PrintReadsSpark` to `.sam`, 1 Mb | | 31 MB | | 6 s | `f079c14f` |
| `ExtractOriginalAlignmentRecordsByNameSpark` | 36.4 GiB | 42 MB | 40 s | 19 s | `.bam` `17819ff0`, `.bai`, `.sbi` |
| `FlagStatSpark` | 36.4 GiB | 30 MB | 42 s | 19 s | `FlagStat`'s, `5434284f` |
| `MeanQualityByCycleSpark` | 35.5 GiB | 29 MB | 43 s | 25 s | `6b5a1c63` |
| `QualityScoreDistributionSpark` | 36.3 GiB | 30 MB | 48 s | 29 s | `fab3cd23` |
| `CollectBaseDistributionByCycleSpark` | 36.3 GiB | 30 MB | 47 s | 24 s | `df9fe6d5` |
| `CollectQualityYieldMetricsSpark` | 36.3 GiB | 33 MB | 42 s | 25 s | `9b2d70ed` |

The metrics md5s are of the file without its command-line and start-time lines. Against GATK the
RSS is 0.02x. The CPU is 0.6x to 0.7x for the first three and about 1x for the metrics, where the
rest is decoding (the engine's, W12).

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
