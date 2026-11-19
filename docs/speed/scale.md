# Speed and RAM at scale

Milestone S, S.7 (#1453). The baseline in [baseline.md](baseline.md) was measured on the
conformance fixtures: a few hundred reads, where the reference's number is mostly the JVM starting.
This one runs the best-practices chain on a real input and asks the question a user asks: on a 30x
genome, at the same number of threads, which side is faster, which uses less memory, and do they
still write the same bytes?

## What is measured

**Data** (`tools/speed/scale.json`, derived by `tools/speed/scale_data.py`). chr20 of NA12878 from
the 1000 Genomes high-coverage set: NYGC's 30x PCR-free 2x150 CRAM, cut to chr20 by range
requests and written as BAM. That gives 17.7 M reads (about 41x on chr20), 923 MiB. The rest is the
Broad's GRCh38 chr20, and chr20 of dbSNP 138, Mills, known indels, HapMap 3.3, Omni 2.5, 1000G
phase 1, af-only gnomAD and the 1000G panel of normals. Each file is the output of a stated
command in a stated image (samtools 1.24, bcftools 1.24, both by digest), and is pinned by its
sha256. CI derives the data again and refuses to run if any hash differs.

**Chain.** MarkDuplicatesSpark, BaseRecalibrator, ApplyBQSR, HaplotypeCaller `-ERC GVCF`,
GenotypeGVCFs, VariantRecalibrator (SNP), ApplyVQSR, SelectVariants, VariantFiltration, Mutect2
(tumour-only), FilterMutectCalls, PreprocessIntervals, CollectReadCounts, DenoiseReadCounts and
ModelSegments. Each step runs on the reference's outputs of the steps before it, so both sides of
every step get the same input, and a port failure upstream does not hide the steps below it.

**Rules.** Both sides run in the pinned oracle image (`gatk-rs-oracle:4.7.0.0`), one run each per
step, in one container with 4 CPUs (`--cpus`) and a 14 GiB memory limit with no swap. A 16 GiB
GitHub x86-64 runner keeps the rest. Threads are the same on both sides: 4 for the reference's
pair-HMM and Spark (`local[4]`), and `RAYON_NUM_THREADS=4` for the port. The reference's heap is
the JVM's default for that container, as `gatk` gives it. Wall clock, CPU (user + system) and peak
RSS come from `wait4` on the tool's own process. A side that hits the memory limit or the step's
timeout is killed and recorded that way, not retried. The two sides' files are compared one by one,
under the masks `run_array.py` applies. A divergence is a correctness bug, reported with its first
differing line.

## Predictions, written before the first measurement

Written on 2026-11-19, before any step had run at scale on either side. Today the port reads its
inputs into memory whole, and a streaming refactor (S.6) is running elsewhere. These predictions
describe the port as it is before that refactor.

1. **Memory is where the port loses.** On every step that reads the BAM (MarkDuplicatesSpark,
   BaseRecalibrator, ApplyBQSR, HaplotypeCaller, Mutect2, CollectReadCounts), the port's peak RSS
   is above the reference's, by roughly the size of the decoded reads (several GiB against the
   JVM's 2 to 4 GiB). At least one of MarkDuplicatesSpark and ApplyBQSR is killed at 14 GiB. The
   fixture ratio (port RSS 3 to 21% of the reference's) does not carry over to scale.
2. **The VCF steps win on both counts.** GenotypeGVCFs, ApplyVQSR, SelectVariants,
   VariantFiltration, FilterMutectCalls and the CNV steps after CollectReadCounts finish in 0.2 to
   0.6 of the reference's wall clock and use less memory: their inputs are megabytes, so there is
   little to hold.
3. **The pair-HMM steps are the open question.** The reference runs GKL's AVX pair-HMM on 4
   threads, and the port runs its own scalar, bit-exact one. HaplotypeCaller and Mutect2 come out
   between 1x and 3x the reference's wall clock. If either one comes out faster, the assembly and
   pileup work is the larger share, not the pair-HMM.
4. **BaseRecalibrator and ApplyBQSR**, which the reference runs on one thread, come out below 1x on
   wall clock if they finish at all.
5. **Bytes.** At least one step diverges at scale: 17.7 M reads reach paths a fixture of a few
   hundred reads does not. MarkDuplicatesSpark, whose covering array is at 0.963 and not 1.000, is
   the most likely one.

## The baseline

Pending: the first `Speed at scale` run on x86-64.
