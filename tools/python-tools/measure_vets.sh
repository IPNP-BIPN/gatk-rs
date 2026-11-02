#!/bin/bash
# One pass of ExtractVariantAnnotations, TrainVariantAnnotationsModel and ScoreVariantAnnotations
# on the fixture, as `gatk` runs them inside broadinstitute/gatk:4.6.2.0.
#
# Usage: measure_vets.sh <fixture directory> <output directory> [hyperparameters json]
set -e
. /opt/miniconda/etc/profile.d/conda.sh
conda activate gatk
fixture=$1
out=$2
hyper=${3:-}
mkdir -p "$out"
cd "$out"
java_options="-Xmx2g -XX:ActiveProcessorCount=${JAVA_CORES:-2}"
if [ ! -f "$fixture/calls.vcf.idx" ]; then
  gatk --java-options "$java_options" IndexFeatureFile -I "$fixture/calls.vcf"
  gatk --java-options "$java_options" IndexFeatureFile -I "$fixture/resource.vcf"
fi
annotations="-A QD -A MQ -A FS -A SOR"
gatk --java-options "$java_options" ExtractVariantAnnotations -V "$fixture/calls.vcf" $annotations \
  --resource:train,training=true,calibration=true "$fixture/resource.vcf" \
  --maximum-number-of-unlabeled-variants 400 -O extract \
  --add-output-vcf-command-line false --create-output-variant-index false --do-not-gzip-vcf-output
extra=""
[ -n "$hyper" ] && extra="--hyperparameters-json $hyper"
gatk --java-options "$java_options" TrainVariantAnnotationsModel \
  --annotations-hdf5 extract.annot.hdf5 \
  -O train $extra
gatk --java-options "$java_options" ScoreVariantAnnotations -V "$fixture/calls.vcf" $annotations \
  --resource:train,training=true,calibration=true "$fixture/resource.vcf" \
  --model-prefix train -O score \
  --add-output-vcf-command-line false --create-output-variant-index false --do-not-gzip-vcf-output
ls -la
python - <<'P'
import numpy as np, sys
c = np.core._multiarray_umath.__cpu_features__
print('cpu', ' '.join(k for k, v in c.items() if v))
P
