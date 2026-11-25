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

## Scale W. `VariantRecalibrator` held its driving VCFs whole

**Ratio.** W12's 40x row (104k calls, 86k training sites, 54k known sites, allele-specific QD, MQ,
FS and SOR, no `-L`): GATK arm64 takes 9.8 to 12.5 s CPU and 1.48 GB at its cgroup peak; the port
on this base, macOS native, 110.6 s CPU and 487 MB (0.33x RSS, about 10x CPU).

**Where the memory goes.** `multi_variant::open` read every `--variant` file into a `Vec` of
records before the traversal looked at the first one, then merged references into that `Vec`.

**The optimisation.** With no `-L`/`-XL` and VCF inputs, the driving files are opened for their
headers only and merged as they are read (`Merging`, the same queue as `merge_into`; one input
passes through), so a record is held only while its start is the queue's. Anything that stops
the streamed traversal (a reader's refusal, a datum's, an unwritable tranches file) starts the
whole run again with the inputs read whole, so it refuses in its own words and order. Two
things are put back first: the process's random stream, which the first attempt has drawn from,
and the tranches file, which it has emptied.

**Why it cannot move a byte.** The streamed merge pops the same records in the same order as the
whole merge (`merging_is_merging_order`, 2000 random cases), the datum loop is the same, and every
other path is the whole run, from the same random state and the same files.

**The suites that prove it.** `variant-recalibrator`, `apply-vqsr-two-modes`. Natively, old against
new: the 26 rows of the covering array, and 66 runs over a malformed record mid-file, an unsorted
file, a truncated file, no records, BGZF, a bad `QD`, two inputs, the same input twice, a missing
input and a good input with a bad one, each with and without `-L chr1` and with the tranches file
absent, already there, or in a directory that is not. The first version failed 6 of the 26 rows
(the restart reused the advanced random stream) and 6 of the edge runs (a tranches file left
behind); both are fixed and covered by those runs.

**Result.** macOS native, the 40x row: 487 MB to 329 MB, 110.6 s to 51.0 s CPU (the whole read
was the larger half of the CPU): 0.22x RSS of GATK, CPU still about 4.5x on this base, where
#1554 takes it to 1.85 s. Same `vr.recal` (73fa4dc355a0642ead8e76ac98a09011, GATK's too) and
`vr.tranches` (755319ca15df339fcd4cd8db4ce0990f) as #1554.

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
