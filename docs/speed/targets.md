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

## Scale A. A reads input read whole for its dictionary (`CheckReferenceCompatibility`, the variant walkers given `-I`, `GenomicsDBImport`, the interpreter tools)

**Ratio.** Measured at scale (the every-tool pass of `tools/speed/scale.json`, NA12878 chr20 30x,
968 MB BAM, linux/arm64 in the pinned container, 6 GiB limit, port at `c/streaming-guard`):
`CheckReferenceCompatibility` peak RSS 0.28 GiB for GATK, the port killed at the 6 GiB limit.

**Where the memory goes.** `runners.rs`'s `reads_dictionary` reads the whole BAM with `fs::read`
and inflates all of it with `decompress_all`, to answer two questions that need its first bytes:
is it block-compressed, and does it start with the BAM magic. The file is then opened again for
its header. Every tool that only wants a reads input's dictionary goes through it.

**The optimisation.** The sniff of `read_walker_startup_with`: `file_start` for the raw start,
`ReadsDataSource::check` (the same whole-file inflate, streamed, discarded block by block) then
`inflated_start` for the inflated start, and `open_checked` so the check is not run twice.

**Why it cannot move a byte.** The dictionary is the header's `@SQ` records, read by the same
`read_header`. The two questions get the same answers: the start is a prefix of the bytes read
before, and `check` refuses exactly the files `decompress_all` refused (it is `decompress_all` over
256-block batches), which both map to an empty dictionary. A file that cannot be read is refused
with the same `cannot_read` message.

**The suites that prove it.** `check-reference-compatibility`, `multi-variant-walker`,
`annotate-vcf-with-bam-depth`, `genomicsdb-import`, `interpreter-tools`, and the covering arrays
of the tools that reach it.

**Result.** `CheckReferenceCompatibility` on the chr20 BAM: peak RSS killed at 6 GiB before,
93 MB after (0.33x GATK's 0.28 GiB), and the table is GATK's byte for byte. CPU stays at 6.6 s
against GATK's 2 s: `check` still inflates the whole file, which GATK never does for a tool that
reads only the header. That is the reader's, engine-wide, and not this change's.

## Scale B. `FlagStat` collected every read before counting one

**Ratio.** At scale (chr20 30x BAM, linux/arm64, 6 GiB container): GATK 0.64 GiB peak RSS and
20 s CPU, the port killed at the limit. Natively on macOS, with no limit: 19.19 GiB.

**Where the memory goes.** `flag_stat` called `read_walker::traverse`, which collects the whole
traversal into a `Vec<BamRecord>`, then folded it into `FlagStatus`.

**The optimisation.** `traverse_iter`, the same traversal a read at a time, folded as it comes.

**Why it cannot move a byte.** `FlagStatus::add` sees the same reads in the same order; the report
is written only after the last one, so a refusal partway is the same `ReadsError`, printed the same
way, with nothing written before it.

**The suites that prove it.** `counting-walkers`, and `FlagStat`'s covering array.

**Result.** macOS native, chr20: 19.19 GiB to 0.10 GiB peak RSS, report identical (md5
`5434284f...`). CPU 35 s against GATK's 20 s: decoding, the engine's (W12), not this change's.

## Scale C. `UnmarkDuplicates` held the reads and the BAM it wrote (and `rewrite_to` for its archetype)

**Ratio.** At scale (chr20 30x BAM, linux/arm64, 6 GiB container): GATK 0.59 GiB peak RSS and
88 s CPU, the port killed at the limit. Natively on macOS, with no limit: 21.92 GiB.

**Where the memory goes.** `unmark_duplicates_with` collected the traversal, cleared 0x400 on
each read, and `write_records_with` encoded them into a `Vec<u8>` that `write_bam` then wrote,
asking `is_coordinate_sorted` of it by inflating it whole once more.

**The optimisation.** `sam_output::rewrite_to`: each read of `traverse_iter` through the tool's
`apply`, into `write_bam_file` (scratch file renamed on success, MD5 hashed while writing). The
`.bai` is kept under the rule `write_bam` applied, now read off the header (`sam_output::indexes`,
the same `SO:coordinate` test on the same `@HD` line). Every tool of the "rewrite a read, write it"
archetype (`AddOriginalAlignmentTags`, `RevertBaseQualityScores`, ...) can move onto it; this
change moves the one whose `apply` cannot refuse.

**Why it cannot move a byte.** The same reads reach the same `BamWriter` in the same order with the
same level and deflater, so the stream is the same; a traversal refusal leaves no file, as before.
A unit test checks `rewrite_to` against `write_records_with` and that a refusal partway leaves no
file.

**The suites that prove it.** `record-transform`, and `UnmarkDuplicates`' covering array.

**Result.** macOS native, chr20: 21.92 GiB to 0.10 GiB peak RSS; the BAM is the old port's byte for
byte (md5 `7e673ef0...`, both written to the same path). CPU 83 s against GATK's 88 s, inside the
noise of one run and short of 0.5x: decoding and deflating, the engine's (W12).

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
