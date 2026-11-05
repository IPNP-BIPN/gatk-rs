#!/bin/bash
# One pass of DetermineGermlineContigPloidy, GermlineCNVCaller and PostprocessGermlineCNVCalls on
# the fixture, as `gatk` runs them inside broadinstitute/gatk:4.6.2.0, with few iterations.
#
# Usage: measure_gcnv.sh <fixture directory> <output directory> [samples]
set -e
. /opt/miniconda/etc/profile.d/conda.sh
conda activate gatk
fixture=$1
out=$2
samples=${3:-6}
mkdir -p "$out"
cd "$out"
export PYTENSOR_FLAGS="base_compiledir=$out/.pytensor"
java_options="-Xmx2g -XX:ActiveProcessorCount=${JAVA_CORES:-2}"
inputs=""
for i in $(seq 0 $((samples - 1))); do inputs="$inputs -I $fixture/sample$i.counts.tsv"; done
echo "=== ploidy $(date +%s)"
gatk --java-options "$java_options" DetermineGermlineContigPloidy $inputs \
  -L "$fixture/intervals.interval_list" --interval-merging-rule OVERLAPPING_ONLY \
  --contig-ploidy-priors "$fixture/priors.tsv" --output ploidy --output-prefix p \
  --max-advi-iter-first-epoch 200 --max-advi-iter-subsequent-epochs 50 --num-thermal-advi-iters 100 \
  --max-training-epochs 3 --min-training-epochs 1 --max-calling-iters 3
echo "=== caller $(date +%s)"
gatk --java-options "$java_options" GermlineCNVCaller --run-mode COHORT $inputs \
  -L "$fixture/intervals.interval_list" --interval-merging-rule OVERLAPPING_ONLY \
  --contig-ploidy-calls ploidy/p-calls --output gcnv --output-prefix g \
  --max-advi-iter-first-epoch 200 --max-advi-iter-subsequent-epochs 50 --num-thermal-advi-iters 100 \
  --max-training-epochs 3 --min-training-epochs 1 --max-calling-iters 3
echo "=== postprocess $(date +%s)"
gatk --java-options "$java_options" PostprocessGermlineCNVCalls \
  --model-shard-path gcnv/g-model --calls-shard-path gcnv/g-calls \
  --contig-ploidy-calls ploidy/p-calls --sample-index 0 \
  --allosomal-contig chrX --autosomal-ref-copy-number 2 \
  --output-genotyped-intervals post.intervals.vcf.gz --output-genotyped-segments post.segments.vcf.gz \
  --output-denoised-copy-ratios post.denoised.tsv \
  --create-output-variant-index false --add-output-vcf-command-line false
echo "=== done $(date +%s)"
find . -path ./.pytensor -prune -o -type f -print | wc -l
