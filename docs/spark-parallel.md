# Spark, threads and partitions

What #72 and #956 can claim, and the measurements behind it. Everything about the reference was
measured in the x86-64 oracle image under emulation on a laptop (GATK 4.6.2.0, Spark 3.5.0). None of
it is a golden.

## The port's own threads

The only parallel path in the Spark tools is bwa-mem4's pool. It is reached by `BwaSpark`,
`BwaAndMarkDuplicatesPipelineSpark`, `ReadsPipelineSpark`, `PathSeqBwaSpark`, `PathSeqFilterSpark`
and `PathSeqPipelineSpark`, and `bwa_threads` ties its `-t` to `RAYON_NUM_THREADS`. No gatk-rs crate
uses rayon. `PileupSpark` reads the core count only for a task number in a stderr message.

`tools/coverage/thread_independence.py` runs every row of those six tools and four rows of each other
Spark tool at 1 and 4 threads, and compares exit, stdout and each file's digest. The `spark-threads`
CI job runs it in the container.

| Run on an arm64 host | Rows | Accepted | Different |
| --- | --- | --- | --- |
| All 39 Spark tools | 284 | 127 | 1, then 0 |

The one difference was in a metrics file of `PathSeqFilterSpark`. Taking out its `# Started on` clock
gave 20 of 20 rows equal; the file itself was not compared by eye. The refused rows are mostly ones
the arrays mean to refuse, and on this host also those needing ISA-L, so the CI job is the measurement
of record.

The arrays' BWA inputs hold few reads, so they barely fill a batch. The script adds one larger row,
`BwaSpark` on the 700 KB `svbp_reads.bam`, single-ended: accepted, and equal at 1 and 4 threads.

## The reference's master, split and reducers

`tools/coverage/spark_master_measure.py` ran 15 tools, each at `local[1]` and `local[4]`, a small
`--bam-partition-size` and `--num-reducers 3`, with repeats. Each comparison was with the same layout
at `local[1]`, because `local[4]` over one partition is not parallel. The split gave 11 partitions on
the 700 KB BAM, 2 for recalibration and 6 for `HaplotypeCallerSpark`.

The tools: `CountReadsSpark`, `CountBasesSpark`, `FlagStatSpark`, `PrintReadsSpark`, `SortSamSpark`,
`MarkDuplicatesSpark`, `MeanQualityByCycleSpark`, `QualityScoreDistributionSpark`,
`CollectInsertSizeMetricsSpark`, `CollectQualityYieldMetricsSpark`, `PileupSpark`,
`BaseRecalibratorSpark`, `BQSRPipelineSpark`, `HaplotypeCallerSpark`, `PrintVariantsSpark`.
`FindBadGenomicKmersSpark` was dropped: its six runs did not finish inside the 45 minutes of the
first pass.

| Compared | Result |
| --- | --- |
| `local[1]` again | 15 of 15 the same, apart from the `# Started on` clock |
| `local[4]` against `local[1]`, same layout | 15 of 15 the same, apart from that clock, at each layout |
| `local[4]` again | 15 of 15 the same, apart from that clock |
| A smaller split | `PrintReadsSpark` writes other BAM bytes, `.bai` and `.sbi`; the records are the same. `PileupSpark` writes more part files with the same text once joined. The other 13 are the same, apart from that clock and the command line |
| `--num-reducers 3` | `SortSamSpark`, `MarkDuplicatesSpark` and `BQSRPipelineSpark` write other BAM bytes and indexes; the records are the same. The other 12 are the same |

So for these 15 the master does not change the output on these small inputs. The split and the
reducers do change bytes, through the BGZF blocks of the merged parts, and they are arguments, so the
arrays carry them as rows. This is a sample, not a proof for the other 24 tools.

The exception is the panel.

## `CreateReadCountPanelOfNormals`

**Mode.** `RowMatrix.computeSVD(k, true, 1e-9)` chooses on `n` samples: LocalARPACK when `k < n / 3`,
LocalLAPACK otherwise, while `n < 100` or `k > n / 2` with `n <= 15000`; DistARPACK beyond that.
`solver_mode` has this, with the integer divisions.

**Libraries.** The image loads `Java11BLAS`, `F2jLAPACK` and `F2jARPACK` (netlib 3.0.3); no native
library loads. So the arithmetic is Java's and is reachable. `f2j_gesdd.rs` ports `dgesdd` up to 128
samples, with `dbdsdc`'s divide and conquer past 25 in `f2j_dbdsdc.rs`, checked bit for bit by the
`f2j-dgesdd` suite. Past 25 that path raises every scaled diagonal entry under the epsilon to it, so
a value is never under `rCond` and LocalLAPACK's count is `min(k, n)`, a flat sample included.
`f2j_arpack.rs` ports the ARPACK path, `symmetricEigs` over F2J's `dsaupd` and `dseupd`, with the
product Spark hands it: Breeze's `Java11BLAS.dgemvN` (four columns, `Math.fma`) for LocalARPACK, and
for DistARPACK `multiplyGramianMatrixBy`, whose `treeAggregate` has depth 2 and folds the hundred
slices into ten by index mod ten before the driver folds the ten. `dgetv0`'s seed is static in F2J
and starts at `1, 3, 5, 7` in a fresh JVM, which is the tool's one call. The
`create-read-count-panel-arpack` suite checks every eigenvalue's bits, the vectors' fingerprints and
the rank on eleven panels from 63 to 150 samples, two of them of rank under the request. Jacobi
remains only past what LocalLAPACK reaches.

**The rank.** It comes from values under `1e-9` of the first, which are rounding noise, so it depends
on how the Gramian's sum is associated. `treeAggregate` has depth 1 here, and the driver folds each
slice's result as it arrives.

| Fresh JVM runs, one panel | Distinct values | Eigensamples |
| --- | --- | --- |
| `local[4]`, nine samples | 11 of 12 runs | 7 (10 runs), 6 (2 runs) |
| `local[1]`, nine samples | 1 of 12 | 7 |
| `local[1]`, 12 samples of 3 profiles | 1 of 10 | 6 |
| `local[1]`, 25 of 4 | 1 of 6 | 8 |
| `local[1]`, 12 of 5 | 2 of 10 | 10 |
| `local[1]`, 8 of 3, 250 intervals | 2 of 6 | 6 |

So `local[1]` is the usual order and not a promise: 2 of 44 runs had other bits, and the count stayed
put in all 44. The cause is arrival order, not the JIT. Swapping one adjacent pair of slices in the
fold reproduces the rare 12-of-5 run's top and tenth values exactly, in a test of all 99 swaps, and
the in-order fold gives the usual ones. `Math.pow(x, 2)` equalled `x * x` on 200000 values in each JIT
tier, and `dgesdd` gave one answer across them.

The port folds in order, which is the usual run. The `nine` panel's last value is `1.137e-9` of the
first against a threshold of `1e-9`, so a rare order could move its count, though none has at
`local[1]`.

**The suite.** `create-read-count-panel-rank` is `golden-pending`, at `local[1]`, 19 rows. Seven have a
count lower than `min(requested, surviving)`. It holds the count and the values above `1e-6` of the
first to eleven digits, and not the bits, which a golden could not hold.

Against a local dump, not a golden, the port matches all 19 counts and, on the 15 rows where it takes
the reference's path, those values. Its full bits on those 15 equal the usual run wherever that was
measured; the one run that differed was a dump that caught the rare 12-of-5 order.

## A sort past the sample

`sortByKey` cuts a sort with `RangePartitioner`, and its bounds decide where each part, and so each
BGZF block of the merged BAM, ends. Up to `sampleSizePerPartition` keys per input partition
(`ceil(3 * 20 * reducers / input partitions)`) the sketch is every key; past it, each input partition
is a reservoir sample (`SamplingUtils.reservoirSampleAndCount`) drawn with an `XORShiftRandom` seeded
from `byteswap32(index ^ (shift << 16))`, Scala 2.13's MurmurHash3 under it. `shift` is the id of the
RDD `rdd.map(_._1)` that the partitioner builds, so it counts the RDDs the tool made before the sort.
The reference logs it, plus one, as the `MapPartitionsRDD` of the `sortByKey` result stage.

`spark_range_partitioner.rs` ports the sketch and `determineBounds`. The port knows `shift` where it
was measured in the image, and there it samples; elsewhere it refuses by name when the sample would be
drawn, as it did before for every tool.

| Tool | The sort | `shift` | Measured stable over |
| --- | --- | --- | --- |
| `SortSamSpark` | `writeReads` | 15 | both orders, `-L`, `--read-filter`, no default filters, a smaller split |
| `ApplyBQSRSpark` | `writeReads` | 16 | |
| `BQSRPipelineSpark` | `writeReads` | 20 | |
| `MarkDuplicatesSpark` | the queryname sort, then `writeReads` | 14, 37 | not wired: #1490 and #1469 hold the tool |

Against the arm64 image, `SortSamSpark` over 6652 reads at 3, 4 (by name) and 5 reducers (5 input
partitions), and both BQSR writers over 542 reads at 4 to 6 reducers, cut their parts at the same
records and write the same decompressed bytes. The compressed bytes there differ only by that image's
zlib. The `sort-sam-sampled*`, `apply-bqsr-sampled` and `bqsr-pipeline-sampled` runs of `spark-walkers`
hold this against the x86-64 oracle.

Still refused: an imbalanced input partition (`fraction * n > sampleSizePerPartition`, which Spark
samples again with a `BernoulliSampler`), and the sorts of `MarkDuplicatesSpark`,
`BwaAndMarkDuplicatesPipelineSpark`, `ReadsPipelineSpark`, `RevertSamSpark --sanitize`,
`CompareDuplicatesSpark`, the aligners and PathSeq, whose `shift` and input partitions are not
measured.

## The gaps left in `runners_spark.rs`

13 `PORT_LIMITATION` refusals: cluster masters, CRAM input, `PileupSpark --metadata`, an evidence
limit and an input larger than one partition in the SV tools, `ParallelCopyGCSDirectoryIntoHDFSSpark`,
a short `TableFeature` row, an invalid read attribute, an unknown read group in an output map, Disq's
CRAM writer, R for `--chart` and for `--produce-plot`, and a sampled sort where `shift` is not
measured (above).

3 `shortcut` comments: the Tribble index bin of 8000, the lost executor's message, and the whole-file
index where Disq merges parts.
