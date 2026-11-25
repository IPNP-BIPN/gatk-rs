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

## 5. `BaseRecalibrator`'s counting on several threads (W12)

Numbered after #1541's entries 2 to 4, which this branch does not carry yet.

**Where the time goes.** A native profile of `BaseRecalibrator` on the chr20 smoke slice, on
#1515's streaming traversal (macOS `sample`): most samples are in counting, which is
`process_read`'s covariates and the table updates under it (`update_tables`, `increment_datum`,
`get_checked`), and the rest is decoding and the allocator. GATK counts on one thread, so the
counting alone is the share a second core can take, and splitting only the preparation would leave
the table updates, the largest part, on one thread.

**The optimisation.** The traversal's thread decodes, fetches each contig's bases once and hands
batches of 1,024 reads round robin to counting threads, one per spare core and at most six. Each
counting thread has an engine and tables of its own. At the end every thread's cells are added into
one engine, and `finalizeData`, the quantization and the report run on it as before.

**Why it cannot move a byte.**

- The tables are indexed by key (`NestedIntegerArray`), so the order in which cells arrive is
  nothing; the report walks them by key.
- Without BAQ, each base adds one observation and an error of exactly 0 or 1, so a cell holds whole
  numbers, its mismatch count scaled by 100,000, and these sums are exact far beyond any read
  count; adding several threads' sums gives the one thread's sum to the bit. The empirical quality
  is computed from the sums, after the merge, as before.
- With BAQ on, the errors are fractional, so their sum depends on its order. Counting then stays on
  one thread.
- GATK stops at the first read that fails. So does this: each failure carries its read's position,
  and the earliest one is returned, whether the traversal, the reference or a counting thread
  raised it.

**The suites that prove it.** `base_recalibrator`'s golden suite, now run on one thread and on three
threads handed one read at a time, the finest split there is (a planted bug that drops one thread's
share fails it); every array that runs the binary against the reference, on runners with several
cores; and the smoke slice below, identical output on every run.

**What the threads first cost, and what removed it.** Counting on seven threads first spent more
CPU and memory than one: the same instructions took 17% more cycles, and the RSS rose by 53 MB a
thread, whatever the batch size. Each engine built BAQ's 49 MB emission table, which a run without
BAQ never reads. Three fixes, each byte-neutral by construction, made the threads cheaper than
#1515's one thread:

- the BAQ table is built the first time a read is BAQ'ed (`OnceCell`);
- the per-read covariate matrix is one buffer with the same `[event][position][covariate]`
  indexing and bounds, not one allocation per event and position (450 for a 150-base read);
- a read with nothing to clip (no soft clip, its adaptor boundary outside it, no low-quality tail)
  is used as it is, where three clips each handed back a full copy of it.

The VCF known sites are also decoded a line at a time, keeping only each record's span: the whole
file and its 1.4 million records with their INFO were held together, 2.8 GB for dbSNP at chr20.
The same refusals come in the same order: eleven malformed and edge-case files (no `#CHROM`, a bad
`POS`, missing columns, a non-UTF-8 byte after a bad record, a bad genotype, CRLF, no final newline)
give the old binary's status, message and table.

**Measured (M-series, 12 performance and 4 efficiency cores, shared with other work; chr20 smoke
slice, `-L chr20:30000000-35000000`, dbSNP and Mills as known sites; release builds).** The table is
byte-identical on every run, and to GATK's.

| | GATK 4.7.0.0 | #1515 | this branch | vs GATK |
|---|---:|---:|---:|---:|
| CPU (user + sys) | 18.7, 20.5 s | 16.3, 17.3, 16.3 s | 12.4, 12.6, 12.3 s | 0.63x |
| peak memory | 662, 720 MB | 3,049 to 3,051 MB | 243, 251, 248 MB | 0.36x |
| wall | 41.3, 44.7 s | 15.8, 16.9, 15.7 s | 6.7, 6.3, 6.2 s | 0.15x |
| instructions retired | | 324.5 G | 200.6 G | |

GATK ran in the arm64 oracle image (`java -jar`, JDK 17) on the same host, with the bgzipped known
sites and their indexes; its CPU and peak memory are the container's cgroup counters, so its
memory includes the page cache and overstates its heap. On one counting thread this branch takes
10.9 to 11.3 s of CPU and 214 MB, so the six threads cost about 12% more CPU for 1.7x less wall.

Left for the 0.5x CPU and 0.25x memory target: the whole contig is fetched (64 MB here, 250 MB on
chr1) where GATK reads windows; every known site keeps its contig as a `String`; and dbSNP's INFO
is decoded on one thread before counting starts, about 4 s of the CPU.

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
