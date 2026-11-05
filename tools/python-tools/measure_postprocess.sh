#!/bin/bash
# PostprocessGermlineCNVCalls alone, on the shards an earlier measure_gcnv.sh left, to see whether
# what it does with fixed inputs is stable when what made the inputs is not.
#
# Usage: measure_postprocess.sh <directory measure_gcnv.sh wrote> <name for these outputs> [sample]
set -e
. /opt/miniconda/etc/profile.d/conda.sh
conda activate gatk
cd "$1"
name=$2
sample=${3:-0}
mkdir -p "$name"
export PYTENSOR_FLAGS="base_compiledir=/tmp/pytensor-$name${PYTENSOR_EXTRA:+,$PYTENSOR_EXTRA}"
gatk --java-options "-Xmx2g -XX:ActiveProcessorCount=${JAVA_CORES:-2}" PostprocessGermlineCNVCalls \
  --model-shard-path gcnv/g-model --calls-shard-path gcnv/g-calls \
  --contig-ploidy-calls ploidy/p-calls --sample-index "$sample" \
  --allosomal-contig chrX --autosomal-ref-copy-number 2 \
  --output-genotyped-intervals "$name/post.intervals.vcf" --output-genotyped-segments "$name/post.segments.vcf" \
  --output-denoised-copy-ratios "$name/post.denoised.tsv" \
  --create-output-variant-index false --add-output-vcf-command-line false
