#!/bin/bash
# DetermineGermlineContigPloidy on the gcnv-frozen fixture, traced for the port (ploidy_trace/).
# Inside broadinstitute/gatk:4.6.2.0 with tools/python-tools at /harness. Usage: ploidy_trace.sh <out>
set -e
. /harness/gcnv_frozen_env.sh
. /opt/miniconda/etc/profile.d/conda.sh
conda activate gatk
out=$1
mkdir -p "$out/fixture" "$out/trace"
python /harness/make_gcnv_fixture.py "$out/fixture" >&2
cd "$out"
export PYTENSOR_FLAGS="allow_gc=False,base_compiledir=$out/.pytensor${PYTENSOR_EXTRA:+,$PYTENSOR_EXTRA}"
export PYTHONPATH=/harness/ploidy_trace PLOIDY_TRACE="$out/trace"
inputs=""
for i in 0 1 2 3 4 5; do inputs="$inputs -I fixture/sample$i.counts.tsv"; done
gatk --java-options "-Xmx2g -XX:ActiveProcessorCount=2" DetermineGermlineContigPloidy $inputs \
  -L fixture/intervals.interval_list --interval-merging-rule OVERLAPPING_ONLY \
  --contig-ploidy-priors fixture/priors.tsv --output ploidy --output-prefix p \
  --max-advi-iter-first-epoch 200 --max-advi-iter-subsequent-epochs 50 --num-thermal-advi-iters 100 \
  --max-training-epochs 3 --min-training-epochs 1 --max-calling-iters 3
