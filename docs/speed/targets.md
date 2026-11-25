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

**Result: refused.** Linux arm64 (Docker, `rust:1.97.1-bookworm`, glibc 2.36), the streaming base
`c/streaming-guard` (f0498420) since main holds whole inputs, the same build with and without the
line, three rounds alternated, GNU `time`. PrintReads, chr20:30-32 Mb of NA12878 30x:

| | user + sys | peak RSS |
|---|---:|---:|
| glibc | 13.0 to 13.5 s | 32 to 37 MB |
| mimalloc | 13.1 to 13.7 s | 81 to 102 MB |

The reads are the same; the headers differ only by the output path in `@PG`. CPU holds and peak RSS
rises 2.5x, so the rule above refuses it. HaplotypeCaller and Mutect2 on the same slice were killed
by the shared Docker VM before finishing on both allocators, so they were not measured there.

**The others, and the threaded callers.** The same PrintReads row with jemalloc (`tikv-jemallocator`
0.6) and with both allocators told to return memory at once (`MIMALLOC_PURGE_DELAY=0
MIMALLOC_ARENA_EAGER_COMMIT=0`, `_RJEM_MALLOC_CONF=dirty_decay_ms:0,muzzy_decay_ms:0`): mimalloc so
tuned holds glibc's 34 to 41 MB at the same CPU; jemalloc takes 74 to 110 MB either way, at no less
CPU. Then HaplotypeCaller and Mutect2 on #1474's branch (`c/hc-multicore`, e8ed8e2f, which streams
its output), 4 threads, its `#[global_allocator]` line against the same build without it:

| Row | allocator | user + sys | peak RSS | wall |
|---|---|---:|---:|---:|
| HaplotypeCaller, 1 Mb | glibc | 8.2 to 9.8 s | 99 to 103 MB | 6.3 to 7.4 s |
| | mimalloc | 7.8 to 8.6 s | 156 to 171 MB | 5.2 to 5.7 s |
| | mimalloc, tuned | 6.9 to 8.9 s | 116 to 124 MB | 4.7 to 5.8 s |
| HaplotypeCaller, 1 Mb, `-ERC GVCF` | glibc | 16.7 to 20.3 s | 103 MB | 11.7 to 12.1 s |
| | mimalloc | 13.2 to 15.2 s | 169 to 172 MB | 6.5 to 7.6 s |
| | mimalloc, tuned | 14.6 to 22.1 s | 121 to 135 MB | 6.9 to 10.7 s |
| Mutect2, 200 kb | glibc | 2.8 to 2.9 s | 154 to 157 MB | 2.1 to 2.2 s |
| | mimalloc | 2.3 to 2.4 s | 212 to 238 MB | 1.6 to 1.8 s |
| | mimalloc, tuned | 2.4 to 2.5 s | 183 to 188 MB | 1.7 to 1.8 s |

Every VCF, GVCF and stats file is the same across allocators (the `.idx` differs by its file
timestamp only). Under threads mimalloc buys 10 to 20% of CPU and a third of the wall clock for 45 to
65% more peak RSS, 20% when tuned: the trade the rule refuses. So glibc stays, and #1474's
`#[global_allocator]` line should not land. What would let it back is an allocator that keeps the
threaded gain without the resident growth, measured against this table.

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

**Result: promising, not landed.** Linux arm64 (Docker, `rust:1.97.1-bookworm`), main, trained on
the 139 fixture rows of `current.json`, three rounds alternated, user + sys from GNU `time`, on a
machine other jobs shared. Every output of the 139 rows (215 files) is the same as the non-PGO
binary's, apart from one `# Started on` timestamp that also differs between two non-PGO runs.

| Run | base | PGO |
|---|---:|---:|
| the 139 fixture rows, total | 1.20 s | 0.93 s |
| `ModelSegments`, 2000 samples | 10.6 to 14.6 s | 9.0 to 12.1 s |
| `VariantRecalibrator`, the fixture repeated 40x | 222 to 254 s | 232 to 271 s |

Peak RSS moved by a median of 3% down. The fixture rows are the training set, so their total is
optimistic; `ModelSegments` gains about 14% by the best of three. `VariantRecalibrator` on main is
the quadratic resource scan of entry 6, which no profile can shorten. The spread between rounds of
one binary is as large as the gain, so landing it waits on x86-64 numbers from a CI job. The release
build then needs a training step: three builds rather than one, and training inputs it can reach.

## 11. `VariantRecalibrator`'s resident set: each resource kept as what the lookup reads

**Ratio.** On the 40x input (`vrbig4`: 104,000 calls and 140,120 resource records, 25 MB of
text), with #1554's index, macOS: 1.82 s of CPU against GATK's 9.8-12.5 s (0.15-0.19), but
502 MB of peak RSS against a 1.48 GB cgroup peak for GATK (its page cache is the 25 MB of input):
0.34, above the programme's 0.25 bound for memory, the priority.

**Where the memory goes.** A heap snapshot at the plateau: 4.19 million live blocks averaging
73 bytes, about seventeen per record. Every resource is read whole and kept as
`VariantContext`s: contig, ID, alleles and their bases, filters, INFO keys and values, a genotype
with its sample name and allele copies. Clearing only the INFO and the ID, after each read,
brought the peak to 458 MB: the rest of the record is the bulk.

**The optimisation.** Each resource is reduced, right after it is read, to what the lookup reads
of a training record, none of which depends on the variant looking it up: its start and stop, its
reference bases and alternates' display strings (`isAlleleInList`), the class its type asks
(`checkVariationClass(evalVC, trainVC)`), whether it is unfiltered with more than one allele, and
whether it has no genotype or is polymorphic in its samples (`isValidVariant`). The records are
then dropped.

**Why it cannot move a byte.** Every one of those values is a pure function of the training record
alone, computed with the functions the lookup called, and the lookup combines them in the order
it did: `!filtered && alleles > 1 && sameClass && (trust || untyped || polymorphic)` with each
operand already known, and pure operands make short-circuiting invisible. The one message that
names a training record prints the variant's contig, which the index matched, and the record's
start, which is kept. The calls are untouched, and so are the recal file and the tranches made
from them.

**The suites that prove it.** `VariantRecalibrator`'s array; the index test, now comparing each
slim entry to the one built from the record the scan finds; locally, the 40x row's recal file and
tranches identical before and after.

**Result** (the 40x row, recal file and tranches identical, the fixture row identical too; CPU
user + sys, RAM peak RSS; the port on macOS, min and max of three runs; GATK 4.7.0.0 in its arm64
image, one run, its RAM the cgroup `memory.peak`):

| | CPU | RAM | CPU / GATK | RAM / GATK |
|---|---:|---:|---:|---:|
| #1554 without this entry | 1.82 s | 462-502 MiB | 0.19 | 0.40-0.43 |
| with it | 1.92-1.99 s | 420-435 MiB | 0.20 | 0.36-0.37 |
| GATK | 9.77 s | 1,163 MiB | | |

Still above the 0.25 bound for memory: what remains is the calls, which `multi_variant::open`
reads whole for every walker over several VCFs, about 3.6 million live blocks at the plateau.
Streaming them is the scale session's multi-variant work (`scale-fix/multi-variant-stream`).

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
