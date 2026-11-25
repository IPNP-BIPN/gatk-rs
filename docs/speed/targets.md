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

## 9. HaplotypeCaller's allocations: the pileup iterator, the clipper, reference confidence

**Where the time goes.** A native profile of `HaplotypeCaller -ERC GVCF` on NA12878 chr20:10-11 Mb
(one thread, glibc kept, entry 7 of #1554): about a tenth of the samples in `malloc` and `free`.
The first caller is the locus iterator (`Contexts::next` with what it inlines), which allocates
at every reference base whether or not anything changed there: two `Vec`s of downsamplers and
buckets in `collect_pending_reads` (one reservoir of the target's capacity per sample when
downsampling is on) even when no read starts at that base, a fresh `Vec` of kept states in
`update_read_states`, a pileup `Vec` grown by doubling, and the contig name copied three times.
Next come whole-read copies: `ReadClipper` copied its read on construction and again in
`clip_read` before the first op copied it a third time, and every clip in finalisation goes
through it; `calculate_ref_confidence` copied the region's evidence twice (once to filter it by
sample, once to sort it) for every GVCF region.

**The optimisation.** `collect_pending_reads` returns before building anything when the next
pending read does not start at the boundary; `update_read_states` compacts the states in place;
the pileup `Vec` is sized to the number of states; the contig name is copied twice, not three
times. `ReadClipper` borrows its read, so a clip makes the one copy its first op makes; the
evidence is filtered and sorted as references, and the locus iterator takes any iterator of
references rather than a slice of records.

**Why it cannot move a byte.** The early return takes the path the loop already took when no read
matched: nothing submitted, no random draw (the reservoir draws only in `submit`), every bucket
empty, nothing admitted. The in-place compaction keeps the surviving states in the same order and
steps each machine once, as the copy did; on an error the iteration ends either way. Sorting
references with the same stable sort and the same comparator gives the same order as sorting
copies. A copy that nothing mutates is the record it was copied from. Capacities and the number of
copies of a string are not observable in any output.

**The suites that prove it.** The locus-iterator and read-state suites in `gatk-engine`; the
HaplotypeCaller and Mutect2 arrays; locally, the GVCF and VCF of the profiled region identical to
the unchanged binary's.

**Result** (one thread, NA12878 chr20, user CPU and peak RSS, three interleaved rounds each,
outputs identical to the unchanged binary's, `GATKCommandLine` aside). Linux arm64, glibc,
rustc 1.97.1:

| Row | before | after | RSS before | RSS after |
|---|---:|---:|---:|---:|
| `HaplotypeCaller -ERC GVCF`, 10-11 Mb | 13.41-13.79 s | 12.31-12.42 s (-8.6%) | 97 MiB | 96 MiB |
| `HaplotypeCaller`, 10-11 Mb | 7.59-7.66 s | 7.19-7.24 s (-5.3%) | 95 MiB | 94 MiB |
| `Mutect2`, 10-10.2 Mb | 2.47-2.54 s | 2.36-2.43 s (-4.4%) | 150 MiB | 150 MiB |

On macOS, whose allocator is dearer, the GVCF row goes from 11.9-12.9 s to 10.0-10.8 s.

## 10. Resolving a read's group: ids compared from the end

**Where the time goes.** `read_group::resolve`, under every `sample_name` / `sample_of` (once per
read when reads are admitted, partitioned or paired, and once per pileup element in some
walkers), scans the header's `@RG` lines comparing ids. NA12878's twelve ids are about forty
characters sharing a long common prefix, so each comparison is a `memcmp` that runs nearly to the
end before failing: in the GVCF profile of entry 9, `memcmp` is more than half of `sample_name`'s
samples (about 2% of the run); in `CollectF1R2Counts`, which resolves per pileup element, the whole
lookup was 40% of the main thread (reported by the scale session).

**The optimisation.** The scan tests the length first and then compares the bytes from the last
one backwards, where ids that share a prefix differ.

**Why it cannot move a byte.** Two strings are equal exactly when they have the same length and
the same bytes; the order the bytes are examined in does not change the answer. The scan still
visits the groups in header order and returns the first match, as `find` did.

**The suites that prove it.** The read-group tests in `gatk-engine`, with a test that the
reverse comparison agrees with `==` on equal ids, a prefix, a suffix and a one-byte change; every
tool's array through `sample_name`.

**Result** (macOS, one thread, three interleaved rounds, outputs identical): the HaplotypeCaller
GVCF row of entry 9 (9.51 s to 9.46-9.52 s) and the Mutect2 row (1.98-2.00 s to 2.00-2.01 s) do
not move beyond noise, as a 2% share predicts. The gain is for walkers that resolve per pileup
element, `CollectF1R2Counts` first; its number belongs to the scale session's row.

## 10. Resolving a read's group by comparing ids from the end: refused

Tried and measured, not landed. NA12878's read-group ids share a long prefix, and in the
HaplotypeCaller GVCF profile `memcmp` was about half of `sample_name`'s samples, so
`read_group::resolve` compared ids from the last byte (same equality, same first match). It moved
nothing: HaplotypeCaller GVCF and Mutect2 within noise (macOS, three rounds), and
`CollectF1R2Counts` on chr20:10-15 Mb 20.5 s against 19.8 s (scale session, three alternated runs).
What costs there is the lookup per pileup element as a whole (the `RG` tag, `Attributes::get("SM")`,
the `String`), which an id-to-sample table built once per header removes: 19.8 s to 15.4 s
(-22%), on the scale session's branch.

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
