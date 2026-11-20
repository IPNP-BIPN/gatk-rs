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

## The target

Per tool, at equal threads: **the port's CPU time (user + system) ≤ GATK's, and its peak RSS ≤
GATK's.** Both, with no trade of one for the other, and byte identity unchanged. The last column
of the table says where each tool stands:

- ✅ met: both are at or below GATK's.
- 🔴 one only: the port beats GATK on CPU or on RSS but not on both. **These rows are the perf
  backlog**, because one number has been bought with the other or never worked on.
- ❌ neither, or the port did not finish (killed, refused, or failed).

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

**Source.** `Speed at scale` run 44629287048, on commit `720773d1`, on GitHub's x86-64 runners.
The runner was an AMD EPYC 9V74 (Zen 4, which has AVX-512; that run did not record the flag, and
later runs do), 4 cores, 16 GiB. Only the prep group was measured. The other four groups failed
before any step ran, because the data cache was evicted in the hours between jobs (the workflow
now hands the data on as an artifact). The rest of the chain, and arm64, come with the next run.

| step | reference wall | port wall | port / ref | reference CPU | port CPU | reference peak RSS | port peak RSS | bytes equal | CPU and RSS ≤ GATK |
|---|---:|---:|---:|---:|---:|---:|---:|:-:|:-:|
| MarkDuplicatesSpark | 471 s | 14 s (killed: out of memory) | - | 1771 s | 22 s | 4.45 GiB | 13.98 GiB | - | ❌ port did not finish |
| MarkDuplicates | 88 s | 0 s (exit 2) | - | 135 s | 0 s | 4.20 GiB | 0.01 GiB | - | ❌ port did not finish |
| BaseRecalibrator | 208 s | 4 s (exit 2) | - | 219 s | 14 s | 0.52 GiB | 8.05 GiB | - | ❌ port did not finish |
| ApplyBQSR | 211 s | 15 s (killed: out of memory) | - | 325 s | 34 s | 0.46 GiB | 13.92 GiB | - | ❌ port did not finish |

Why each port side stopped:

- **MarkDuplicatesSpark** was killed at the 14 GiB limit after 14 s. GATK peaked at 4.45 GiB.
- **MarkDuplicates** (Picard's): the port does not carry it yet ("MarkDuplicates is a GATK tool
  that this port does not carry yet").
- **BaseRecalibrator** refused the bgzipped `--known-sites`, a bug, fixed by #1472. It had already
  reached 8.05 GiB, against GATK's 0.52 GiB.
- **ApplyBQSR** was killed at 14 GiB after 15 s. GATK peaked at 0.46 GiB.

Every step that reads the BAM loads it whole today. That is S.6's streaming refactor, not this
benchmark's to fix.

**The predictions, so far.**

1. Memory: confirmed, and worse than predicted. Both steps predicted to die at 14 GiB died there,
   and BaseRecalibrator was at 15 times GATK's peak before it stopped.
2. VCF steps, 3. pair-HMM steps: not yet measured at scale.
4. BaseRecalibrator and ApplyBQSR below 1x on wall clock: not testable; neither finished.
5. Bytes: no x86-64 comparison yet, because no port step finished. The local smoke run on a 2 Mb
   slice already found two divergences, VariantRecalibrator's `##contig` lines and
   CollectReadCounts' `@SQ` `M5`/`UR`.

**Every CPU.** The workflow runs the chain on x86-64 (`ubuntu-latest`) and on arm64
(`ubuntu-24.04-arm`), both sides native: on arm64 the reference is GATK without GKL's x86 native
code, as an ARM user runs it. Its last job compares every output's sha256 between the two CPUs,
the port's and GATK's. The port's bytes must be the same everywhere.
