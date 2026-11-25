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

## Scale K. The locus walkers held every read (`Pileup`, `PileupSpark`, `CallableLoci`, `CollectF1R2Counts`, `ASEReadCounter`, `GetPileupSummaries`)

**Ratio.** At scale (chr20 30x BAM, the arguments of `tools/speed/scale.json`, linux/arm64, 6 GiB
container):

| Tool | GATK CPU | GATK RSS | port |
|---|---:|---:|---|
| `CollectF1R2Counts` | 271 s | 1.77 GiB | killed at 4.76 GiB |
| `CallableLoci` | 135 s | 1.04 GiB | killed at 4.01 GiB |
| `ASEReadCounter` | 128 s | 0.58 GiB | killed at 3.62 GiB |
| `Pileup` (10 Mb) | 122 s | 0.46 GiB | killed at 5.15 GiB |
| `GetPileupSummaries` | 30 s | 3.80 GiB | refused `-V` and `-L` `.vcf.gz` |
| `PileupSpark` (10 Mb) | 83 s | 2.90 GiB | killed at 5.98 GiB |

**Where the memory goes.** `read_walker::traverse` collects every read, and
`locus_walker::traverse` every alignment context over them, before the tool looks at one;
`Pileup` and `CallableLoci` also held each contig's whole reference. Natively, on 5 Mb of chr20,
main holds 3 to 27 GB. `ASEReadCounter` and `GetPileupSummaries` also scanned every variant at
every locus, which never finished on the whole contig, and `PileupSpark` ran the locus traversal
over every read once per 10 kb shard.

**The optimisation.** `locus_walker::for_each_streamed`: the reads come from `traverse_iter` and
are cut into segments of 4096 at a strictly increasing start. A read that reaches a segment's cut
is carried into the next one, the contexts at or past the cut wait for it, and those below the
previous cut were already given. `OverlapCursor` and `EmptyLociCursor` are `-L`'s overlap and the
empty loci over the same stream. It holds only without downsampling (`--max-depth-per-sample 0`,
every tool's default); with it, the whole-list path runs as before. The reference is read a
megabase at a time (`BaseWindow`), the sites are indexed by start, with each contig's longest span,
and `CallableLoci`'s writer takes a locus at a time. `PileupSpark` streams over its shards as one
interval list, and still refuses no read or a read longer than a shard through the whole list.
`GetPileupSummaries` inflates a `.vcf.gz` given to `-V`, as `ASEReadCounter` already did (`-L`'s
is #1564), and `CollectF1R2Counts` resolves each read group's sample once, where it compared the
header's twelve 40-character ids for every pileup element (40% of its samples).

**Why it cannot move a byte.** A segment's contexts are the whole traversal's at the same loci:
every read that covers a locus below the cut started before it, so it is in this segment or carried
into it, and the pileup at that locus is built from the same reads in the same order. The window's
bases are the contig query's. Any refusal of the streamed path (a read, a locus, an unsorted input,
a tool's own) drops its output and runs the whole-list path, which refuses exactly as before.

**The suites that prove it.** The locus walker's suites and the six tools' covering arrays, and
two tests that compare the segmented traversal with the whole one: the fixture configurations
for every segment size from 1 to 1000, and a synthetic BAM whose cigars (deletions, insertions,
clips, skips) reach across the cuts, with and without empty loci. Natively, old against new, on
5 Mb of chr20 and on edge inputs (a malformed record, an unsorted BAM, a reference without the
contig, no reference, `-L`, an output in a missing directory): the same exit, message and output.

**Result.** macOS native, chr20, this branch against main:

| Tool | input | main RSS | RSS | main CPU | CPU |
|---|---|---:|---:|---:|---:|
| `Pileup` | 5 Mb | 11.6 GB | 146 MB | 9.5 s | 8.8 s |
| `CallableLoci` | 5 Mb | 26.7 GB | 145 MB | 21.9 s | 12.4 s |
| `CollectF1R2Counts` | 5 Mb | 11.0 GB | 85 MB | 23.6 s | 15.8 s |
| `ASEReadCounter` (a het VCF) | 100 kb | 365 MB | 148 MB | 25.0 s | 0.2 s |
| `GetPileupSummaries` (`-L` hapmap) | 5 Mb | 3.24 GB | 84 MB | 3.7 s | 3.2 s |
| `PileupSpark` | 200 kb | 2.50 GB | 131 MB | 4.5 s | 1.7 s |
| `PileupSpark` | 100 kb of reads, whole contig | 870 MB | 103 MB | 416 s | 0.2 s |

Every output is main's byte for byte (`CollectF1R2Counts`' archive by its inflated files, which
are also identical on the whole contig before and after the sample lookup). On the whole contig,
against GATK at scale (the GATK column of the first table):

| Tool | RSS | CPU | RSS ratio | CPU ratio |
|---|---:|---:|---:|---:|
| `Pileup` (10 Mb) | 192 MB | 24.2 s | 0.41 | 0.20 |
| `CallableLoci` | 157 MB | 85.0 s | 0.15 | 0.63 |
| `CollectF1R2Counts` | 117 MB | 214.7 s | 0.06 | 0.79 |
| `ASEReadCounter` (hapmap) | 212 MB | 88.0 s | 0.36 | 0.69 |
| `GetPileupSummaries` | 304 MB | 37.4 s | 0.08 | 1.25 |
| `PileupSpark` (10 Mb) | 144 MB | 33.3 s | 0.05 | 0.40 |

The RSS left is the reader's (about 100 MB on this base, a third of that with #1482) and, for
`ASEReadCounter` and `GetPileupSummaries`, the whole variant file decoded. The CPU left is decoding and the pileup itself, the
engine's (W12).

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
