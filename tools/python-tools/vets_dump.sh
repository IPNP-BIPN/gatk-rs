#!/bin/bash
# The `vets-isolation-forest` suite's dump, inside broadinstitute/gatk:4.7.0.0: two passes over the
# fixture and one line per file they wrote. `dump_vets.py` says what the lines hold.
set -e
# numpy takes SVML loops for power and log on AVX-512, which the port does not follow; with them off
# the reference computes what it computes on an AVX2 CPU, whichever runner this lands on.
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR"
work=$(mktemp -d)
python3=/opt/miniconda/envs/gatk/bin/python
mkdir -p "$work/fixture"
$python3 /harness/make_vets_fixture.py "$work/fixture" >&2
printf '{"random_state": 1, "n_estimators": 7, "n_jobs": 2}' > "$work/fixture/seven.json"
SCORE_EXTRA="--snp-calibration-sensitivity-threshold 0.9 --indel-calibration-sensitivity-threshold 0.95" \
  bash /harness/measure_vets.sh "$work/fixture" "$work/default" >&2
bash /harness/measure_vets.sh "$work/fixture" "$work/seven" "$work/fixture/seven.json" >&2
$python3 /harness/dump_vets.py "$work/fixture" default "$work/default" seven "$work/seven"
