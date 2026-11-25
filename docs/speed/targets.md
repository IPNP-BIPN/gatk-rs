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

## Every tool at scale: the ranked table

Every tool of the inventory (`tools/coverage/arrays/*.t2.json`, the Picard arrays and
`MarkDuplicates`, 209 in all) is either a step of `tools/speed/scale.json` (164 steps, 159 tools)
or carries a written reason under its `not_measured` key (50 tools). `bench.py --check` fails on any
tool that is neither. Each step runs both sides in the pinned container under one memory limit
and one thread count, sides alternating (`timer.py`, wait4 rusage: CPU is user plus system, RSS
the peak). The bar is the milestone's: CPU at most 0.5x GATK's and peak RSS at most 0.25x.

**These rows are a local pass, not the measurement.** Apple M4 Max, linux/arm64 in
`gatk-rs-oracle-arm64:4.7.0.0`, 6 GiB, 4 threads, one run per side, GATK heap `-Xmx3584m`, the port
at `c/streaming-guard` (f0498420, the streaming reads stack that main does not have yet). The
`Speed at scale` workflow regenerates this table on x86-64 and arm64 runners (14 GiB, three runs per
side) with `bench.py --rank`.

Gaps in this pass: the local Docker daemon went down mid-run, which lost `MarkDuplicatesSpark` and
with it every step after it in the germline, joint and somatic chains ("no reference output");
`BwaMemIndexImageCreator` has no BWA library on arm64, so the BWA and SV chains cannot run here;
the CRAM steps ran before `--disable-sequence-dictionary-validation` was added (the CRAM carries
the full 3366-contig header, the reference only chr20); and five GATK command lines were wrong and
are fixed in `scale.json` (`CompareIntervalLists`, `LocalAssembler`, `CollectInsertSizeMetricsSpark`,
`CollectMultipleMetricsSpark`, `GetPileupSummaries`).

The fixes the table has led to so far, each with its own entry here, its profile and its proof of
unchanged bytes: #1553 (`CheckReferenceCompatibility`, a reads input read whole for its
dictionary), #1555 (`FlagStat`), #1556 (`UnmarkDuplicates`, with `sam_output::rewrite_to`), #1557
(`AddOriginalAlignmentTags`), #1558 (`VariantsToTable`), #1559 (`GatherVcfs`), #1560 (every
variant walker: a quadratic index scan), #1561 (`UpdateVCFSequenceDictionary`) and the
`RemoveNearbyIndels` change.

| step | host | GATK CPU | port CPU | CPU ratio | GATK peak RSS | port peak RSS | RSS ratio | reps | CPU ≤ 0.5x and RSS ≤ 0.25x |
|---|---|---:|---:|---:|---:|---:|---:|---:|:-:|
| CollectAllelicCountsSpark | aarch64 - | 249 s | 9 s | 0.04 | 4.04 GiB | 0.05 GiB | 0.01 | 1 | exit 2 |
| BQSRPipelineSpark | aarch64 - | 670 s | 12 s | 0.02 | 3.83 GiB | 0.08 GiB | 0.02 | 1 | exit 2 |
| SortSamSpark | aarch64 - | 146 s | 18 s | 0.12 | 3.75 GiB | 5.98 GiB | 1.59 | 1 | killed: out of memory |
| SortVcf | aarch64 - | 34 s | 4 s | 0.13 | 3.57 GiB | 5.98 GiB | 1.68 | 1 | killed: out of memory |
| BaseRecalibratorSpark | aarch64 - | 192 s | 11 s | 0.06 | 2.90 GiB | 0.08 GiB | 0.03 | 1 | exit 2 |
| PileupSpark | aarch64 - | 83 s | 26 s | 0.31 | 2.90 GiB | 5.98 GiB | 2.06 | 1 | killed: out of memory |
| PrintReadsSpark | aarch64 - | 71 s | 24 s | 0.33 | 2.22 GiB | 5.94 GiB | 2.68 | 1 | killed: out of memory |
| ExtractOriginalAlignmentRecordsByNameSpark | aarch64 - | 27 s | 21 s | 0.80 | 2.19 GiB | 4.66 GiB | 2.13 | 1 | killed: out of memory |
| CalcMetadataSpark | aarch64 - | 45 s | 16 s | 0.37 | 1.85 GiB | 0.05 GiB | 0.03 | 1 | exit 2 |
| CollectF1R2Counts | aarch64 - | 271 s | 18 s | 0.07 | 1.77 GiB | 4.76 GiB | 2.69 | 1 | killed: out of memory |
| FlagStatSpark | aarch64 - | 28 s | 19 s | 0.67 | 1.59 GiB | 5.98 GiB | 3.76 | 1 | killed: out of memory |
| QualityScoreDistributionSpark | aarch64 - | 34 s | 30 s | 0.90 | 1.58 GiB | 5.98 GiB | 3.79 | 1 | killed: out of memory |
| CollectQualityYieldMetricsSpark | aarch64 - | 28 s | 19 s | 0.68 | 1.49 GiB | 5.98 GiB | 4.00 | 1 | killed: out of memory |
| MeanQualityByCycleSpark | aarch64 - | 24 s | 19 s | 0.76 | 1.45 GiB | 5.98 GiB | 4.13 | 1 | killed: out of memory |
| CollectBaseDistributionByCycleSpark | aarch64 - | 31 s | 21 s | 0.68 | 1.29 GiB | 5.98 GiB | 4.63 | 1 | killed: out of memory |
| VcfFormatConverter | aarch64 - | 24 s | 5 s | 0.23 | 1.27 GiB | 3.54 GiB | 2.78 | 1 | killed: out of memory |
| RevertSamSpark | aarch64 - | 82 s | 18 s | 0.22 | 1.13 GiB | 5.98 GiB | 5.27 | 1 | killed: out of memory |
| CollectAllelicCounts | aarch64 - | 42 s | 9 s | 0.22 | 1.12 GiB | 0.05 GiB | 0.04 | 1 | exit 2 |
| CallableLoci | aarch64 - | 135 s | 15 s | 0.11 | 1.04 GiB | 4.01 GiB | 3.86 | 1 | killed: out of memory |
| FixVcfHeader | aarch64 - | 18 s | 8 s | 0.43 | 0.98 GiB | 3.82 GiB | 3.89 | 1 | killed: out of memory |
| LeftAlignAndTrimVariants | aarch64 - | 15 s | 7 s | 0.49 | 0.86 GiB | 4.51 GiB | 5.27 | 1 | killed: out of memory |
| DepthOfCoverage | aarch64 - | 688 s | 6 s | 0.01 | 0.77 GiB | 0.05 GiB | 0.07 | 1 | exit 3 |
| AddOriginalAlignmentTags | aarch64 - | 113 s | 17 s | 0.15 | 0.74 GiB | 5.98 GiB | 8.03 | 1 | killed: out of memory |
| FlagStat | aarch64 - | 20 s | 14 s | 0.67 | 0.64 GiB | 5.98 GiB | 9.34 | 1 | killed: out of memory |
| UpdateVcfSequenceDictionary | aarch64 - | 18 s | 4 s | 0.23 | 0.63 GiB | 4.08 GiB | 6.48 | 1 | killed: out of memory |
| UnmarkDuplicates | aarch64 - | 88 s | 17 s | 0.19 | 0.59 GiB | 5.98 GiB | 10.20 | 1 | killed: out of memory |
| ASEReadCounter | aarch64 - | 128 s | 17 s | 0.13 | 0.58 GiB | 3.62 GiB | 6.29 | 1 | killed: out of memory |
| SelectVariants.gnomad | aarch64 - | 13 s | 7 s | 0.58 | 0.56 GiB | 5.68 GiB | 10.18 | 1 | killed: out of memory |
| LearnReadOrientationModel | aarch64 - | 10 s | 0 s | 0.00 | 0.56 GiB | 0.04 GiB | 0.07 | 1 | exit 101 |
| SplitVcfs | aarch64 - | 13 s | 5 s | 0.40 | 0.55 GiB | 3.99 GiB | 7.23 | 1 | killed: out of memory |
| UpdateVCFSequenceDictionary | aarch64 - | 13 s | 4 s | 0.33 | 0.53 GiB | 4.13 GiB | 7.75 | 1 | killed: out of memory |
| RemoveNearbyIndels | aarch64 - | 15 s | 58 s | 3.77 | 0.53 GiB | 5.31 GiB | 10.04 | 1 | killed: out of memory |
| MergeVcfs | aarch64 - | 17 s | 7 s | 0.43 | 0.51 GiB | 4.51 GiB | 8.78 | 1 | killed: out of memory |
| VcfToIntervalList | aarch64 - | 7 s | 8 s | 1.06 | 0.51 GiB | 4.79 GiB | 9.40 | 1 | killed: out of memory |
| MakeSitesOnlyVcf | aarch64 - | 8 s | 5 s | 0.65 | 0.50 GiB | 3.80 GiB | 7.68 | 1 | killed: out of memory |
| Pileup | aarch64 - | 122 s | 17 s | 0.14 | 0.46 GiB | 5.15 GiB | 11.08 | 1 | killed: out of memory |
| SplitCRAM | aarch64 - | 4 s | 10 s | 2.79 | 0.45 GiB | 2.64 GiB | 5.81 | 1 | killed: out of memory |
| VariantsToTable | aarch64 - | 75 s | 32 s | 0.43 | 0.28 GiB | 5.12 GiB | 18.19 | 1 | killed: out of memory |
| PrintFileDiagnostics | aarch64 - | 8 s | 0 s | 0.02 | 0.28 GiB | 0.29 GiB | 1.04 | 1 | exit 3 |
| CheckReferenceCompatibility | aarch64 - | 2 s | 50 s | 20.82 | 0.28 GiB | 5.98 GiB | 21.43 | 1 | killed: out of memory |
| GatherVcfsCloud | aarch64 - | 3 s | 1 s | 0.33 | 0.27 GiB | 0.70 GiB | 2.60 | 1 | exit 3 |
| GatherVcfs | aarch64 - | 2 s | 5 s | 2.28 | 0.27 GiB | 4.26 GiB | 15.98 | 1 | killed: out of memory |
| ValidateVariants | aarch64 - | 9 s | 8 s | 0.92 | 1.19 GiB | 5.36 GiB | 4.49 | 1 | ❌ neither |
| IndexFeatureFile | aarch64 - | 6 s | 7 s | 1.13 | 0.67 GiB | 1.37 GiB | 2.05 | 1 | ❌ neither |
| CountVariants | aarch64 - | 7 s | 2 s | 0.29 | 0.61 GiB | 0.73 GiB | 1.19 | 1 | 🔴 CPU only |
| PrintDistantMates | aarch64 - | 36 s | 37 s | 1.05 | 1.81 GiB | 1.46 GiB | 0.81 | 1 | ❌ neither |
| CountVariantsSpark | aarch64 - | 10 s | 2 s | 0.17 | 1.28 GiB | 0.72 GiB | 0.57 | 1 | 🔴 CPU only |
| PreprocessIntervals | aarch64 - | 2 s | 153 s | 65.06 | 0.34 GiB | 0.14 GiB | 0.40 | 1 | ❌ neither |
| PrintBGZFBlockInformation | aarch64 - | 3 s | 0 s | 0.01 | 0.30 GiB | 0.07 GiB | 0.22 | 1 | ✅ met |
| PrintReadsHeader | aarch64 - | 2 s | 7 s | 3.10 | 0.29 GiB | 0.05 GiB | 0.18 | 1 | 🔴 RSS only |
| GetSampleName | aarch64 - | 2 s | 6 s | 2.87 | 0.30 GiB | 0.05 GiB | 0.16 | 1 | 🔴 RSS only |
| CountBases | aarch64 - | 24 s | 32 s | 1.32 | 0.53 GiB | 0.05 GiB | 0.09 | 1 | 🔴 RSS only |
| DumpTabixIndex | aarch64 - | 3 s | 0 s | 0.00 | 0.27 GiB | 0.02 GiB | 0.07 | 1 | ✅ met |
| PrintReads | aarch64 - | 78 s | 81 s | 1.03 | 0.73 GiB | 0.05 GiB | 0.07 | 1 | 🔴 RSS only |
| CountReads | aarch64 - | 24 s | 31 s | 1.26 | 0.81 GiB | 0.05 GiB | 0.06 | 1 | 🔴 RSS only |
| FilterIntervals | aarch64 - | 3 s | 4 s | 1.04 | 0.53 GiB | 0.03 GiB | 0.06 | 1 | 🔴 RSS only |
| CountReadsSpark | aarch64 - | 21 s | 33 s | 1.60 | 1.35 GiB | 0.05 GiB | 0.04 | 1 | 🔴 RSS only |
| AnnotateIntervals | aarch64 - | 4 s | 0 s | 0.11 | 0.47 GiB | 0.02 GiB | 0.04 | 1 | ✅ met |
| CountBasesSpark | aarch64 - | 30 s | 33 s | 1.10 | 1.61 GiB | 0.05 GiB | 0.03 | 1 | 🔴 RSS only |
| BaseRecalibrator.cram | aarch64 - | GATK did not finish (exit 2) |
| BwaMemIndexImageCreator | aarch64 - | GATK did not finish (exit 3) |
| CRAMIssue8768Detector | aarch64 - | GATK did not finish (exit 1) |
| CollectInsertSizeMetricsSpark | aarch64 - | GATK did not finish (exit 1) |
| CollectMultipleMetricsSpark | aarch64 - | GATK did not finish (exit 52) |
| CompareIntervalLists | aarch64 - | GATK did not finish (exit 3) |
| CountReads.cram | aarch64 - | GATK did not finish (exit 2) |
| FlagStat.cram | aarch64 - | GATK did not finish (exit 2) |
| GetPileupSummaries | aarch64 - | GATK did not finish (exit 137) |
| LocalAssembler | aarch64 - | GATK did not finish (exit 1) |
| PrintReads.cram | aarch64 - | GATK did not finish (exit 2) |
| AlleleFrequencyQC | ? | not run: no reference output from GenotypeGVCFs to run on |
| AnalyzeCovariates | ? | not run: no reference output from BaseRecalibrator to run on |
| AnnotateVcfWithBamDepth | ? | not run: no reference output from GenotypeGVCFs to run on |
| AnnotateVcfWithExpectedAlleleFraction | ? | not run: no reference output from GenotypeGVCFs, CalculateMixingFractions to run on |
| ApplyBQSR | ? | not run: no reference output from MarkDuplicatesSpark, BaseRecalibrator to run on |
| ApplyBQSRSpark | ? | not run: no reference output from BaseRecalibrator to run on |
| ApplyVQSR | ? | not run: no reference output from GenotypeGVCFs, VariantRecalibrator to run on |
| BaseRecalibrator | ? | not run: no reference output from MarkDuplicatesSpark to run on |
| BwaAndMarkDuplicatesPipelineSpark | aarch64 - | not run: no reference output from BwaMemIndexImageCreator to run on |
| BwaSpark | aarch64 - | not run: no reference output from BwaMemIndexImageCreator to run on |
| CalculateContamination | aarch64 - | not run: no reference output from GetPileupSummaries to run on |
| CalculateGenotypePosteriors | ? | not run: no reference output from GenotypeGVCFs to run on |
| CalculateMixingFractions | ? | not run: no reference output from GenotypeGVCFs to run on |
| CallCopyRatioSegments | aarch64 - | not run: no reference output from ModelSegments to run on |
| CollectReadCounts | aarch64 - | not run: no reference output from ApplyBQSR to run on |
| CollectSVEvidence | ? | not run: timer failed: failed to connect to the docker API at unix:///Users/benjamin/.docker/run/docker.sock; check if the path is correct and if the daemon is running: dial unix /Users/benjamin/.docker/run/docker.sock: connect: no such file or directory |
| CombineGVCFs | ? | not run: no reference output from HaplotypeCaller to run on |
| CombineSegmentBreakpoints | aarch64 - | not run: no reference output from ModelSegments, CallCopyRatioSegments to run on |
| CompareBaseQualities | ? | not run: no reference output from MarkDuplicatesSpark, ApplyBQSR to run on |
| CompareDuplicatesSpark | ? | not run: no reference output from MarkDuplicatesSpark to run on |
| Concordance | ? | not run: no reference output from GenotypeGVCFs, FilterMutectCalls to run on |
| CondenseDepthEvidence | ? | not run: no reference output from CollectSVEvidence to run on |
| ConvertCountsToDepthFile | aarch64 - | not run: no reference output from CollectReadCounts to run on |
| CountFalsePositives | ? | not run: no reference output from FilterMutectCalls to run on |
| DenoiseReadCounts | aarch64 - | not run: no reference output from CollectReadCounts to run on |
| DiscoverVariantsFromContigAlignmentsSAMSpark | ? | not run: no reference output from FindBreakpointEvidenceSpark to run on |
| EvaluateInfoFieldConcordance | ? | not run: no reference output from GenotypeGVCFs, ApplyVQSR to run on |
| ExampleMultiFeatureWalker | ? | not run: no reference output from CollectSVEvidence to run on |
| ExtractSVEvidenceSpark | ? | not run: timer failed: failed to connect to the docker API at unix:///Users/benjamin/.docker/run/docker.sock; check if the path is correct and if the daemon is running: dial unix /Users/benjamin/.docker/run/docker.sock: connect: no such file or directory |
| ExtractVariantAnnotations | ? | not run: no reference output from GenotypeGVCFs to run on |
| FastaAlternateReferenceMaker | ? | not run: no reference output from GenotypeGVCFs to run on |
| FilterAlignmentArtifacts | ? | not run: no reference output from ApplyBQSR, FilterMutectCalls, BwaMemIndexImageCreator to run on |
| FilterMutectCalls | ? | not run: no reference output from Mutect2 to run on |
| FilterVariantTranches | ? | not run: no reference output from ApplyVQSR to run on |
| FilterVcf | ? | not run: no reference output from GenotypeGVCFs to run on |
| FindBadGenomicKmersSpark | aarch64 - | not run: timer failed: time="2026-11-25T13:22:20+01:00" level=error msg="error waiting for container: unexpected EOF" |
| FindBreakpointEvidenceSpark | ? | not run: no reference output from BwaMemIndexImageCreator, FindBadGenomicKmersSpark to run on |
| GatherBQSRReports | ? | not run: no reference output from BaseRecalibrator to run on |
| GatherPileupSummaries | aarch64 - | not run: no reference output from GetPileupSummaries to run on |
| GatherTranches | ? | not run: no reference output from VariantRecalibrator to run on |
| GenomicsDBImport | ? | not run: no reference output from HaplotypeCaller to run on |
| GenotypeGVCFs | ? | not run: no reference output from HaplotypeCaller to run on |
| GnarlyGenotyper | ? | not run: no reference output from ReblockGVCF to run on |
| HaplotypeCaller | ? | not run: no reference output from ApplyBQSR to run on |
| HaplotypeCallerSpark | ? | not run: no reference output from ApplyBQSR to run on |
| MarkDuplicates | ? | not run: timer failed: failed to connect to the docker API at unix:///Users/benjamin/.docker/run/docker.sock; check if the path is correct and if the daemon is running: dial unix /Users/benjamin/.docker/run/docker.sock: connect: no such file or directory |
| MarkDuplicatesSpark | ? | not run: timer failed: failed to connect to the docker API at unix:///Users/benjamin/.docker/run/docker.sock; check if the path is correct and if the daemon is running: dial unix /Users/benjamin/.docker/run/docker.sock: connect: no such file or directory |
| MergeAnnotatedRegions | aarch64 - | not run: no reference output from ModelSegments to run on |
| MergeAnnotatedRegionsByAnnotation | aarch64 - | not run: no reference output from CallCopyRatioSegments to run on |
| MergeMutectStats | ? | not run: no reference output from Mutect2 to run on |
| ModelSegments | aarch64 - | not run: no reference output from DenoiseReadCounts to run on |
| Mutect2 | ? | not run: no reference output from ApplyBQSR to run on |
| PathSeqBuildKmers | ? | not run: timer failed: failed to connect to the docker API at unix:///Users/benjamin/.docker/run/docker.sock; check if the path is correct and if the daemon is running: dial unix /Users/benjamin/.docker/run/docker.sock: connect: no such file or directory |
| PathSeqFilterSpark | ? | not run: no reference output from PathSeqBuildKmers to run on |
| PrintReadCounts | aarch64 - | not run: no reference output from CollectReadCounts to run on |
| PrintSVEvidence | ? | not run: no reference output from CollectSVEvidence to run on |
| PrintVariantsSpark | aarch64 - | not run: timer failed:  |
| ReadsPipelineSpark | ? | not run: timer failed: time="2026-11-25T13:23:10+01:00" level=error msg="Error waiting for container: Canceled: grpc: the client connection is closing: context canceled" |
| ReblockGVCF | ? | not run: no reference output from HaplotypeCaller to run on |
| ReferenceBlockConcordance | ? | not run: no reference output from HaplotypeCaller, ReblockGVCF to run on |
| RenameSampleInVcf | ? | not run: no reference output from GenotypeGVCFs to run on |
| ScoreVariantAnnotations | ? | not run: no reference output from GenotypeGVCFs, TrainVariantAnnotationsModel to run on |
| SelectVariants | ? | not run: no reference output from ApplyVQSR to run on |
| SiteDepthtoBAF | ? | not run: no reference output from CollectSVEvidence to run on |
| StructuralVariantDiscoverer | ? | not run: no reference output from FindBreakpointEvidenceSpark to run on |
| StructuralVariationDiscoveryPipelineSpark | ? | not run: no reference output from BwaMemIndexImageCreator, FindBadGenomicKmersSpark to run on |
| SvDiscoverFromLocalAssemblyContigAlignmentsSpark | ? | not run: no reference output from FindBreakpointEvidenceSpark to run on |
| TagGermlineEvents | aarch64 - | not run: no reference output from CallCopyRatioSegments to run on |
| TrainVariantAnnotationsModel | ? | not run: no reference output from ExtractVariantAnnotations to run on |
| VCFComparator | ? | not run: no reference output from HaplotypeCaller, ReblockGVCF to run on |
| ValidateBasicSomaticShortMutations | ? | not run: no reference output from FilterMutectCalls to run on |
| VariantAnnotator | ? | not run: no reference output from GenotypeGVCFs to run on |
| VariantEval | ? | not run: no reference output from GenotypeGVCFs to run on |
| VariantFiltration | ? | not run: no reference output from GenotypeGVCFs to run on |
| VariantRecalibrator | ? | not run: no reference output from GenotypeGVCFs to run on |

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
