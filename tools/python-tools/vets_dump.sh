#!/bin/bash
# The `vets-isolation-forest` suite's dump, inside broadinstitute/gatk:4.6.2.0: two passes over the
# fixture and one line per file they wrote. `dump_vets.py` says what the lines hold.
set -e
work=$(mktemp -d)
python3=/opt/miniconda/envs/gatk/bin/python
$python3 /harness/make_vets_fixture.py "$work/fixture" >&2
printf '{"random_state": 1, "n_estimators": 7, "n_jobs": 2}' > "$work/fixture/seven.json"
SCORE_EXTRA="--snp-calibration-sensitivity-threshold 0.9 --indel-calibration-sensitivity-threshold 0.95" \
  bash /harness/measure_vets.sh "$work/fixture" "$work/default" >&2
bash /harness/measure_vets.sh "$work/fixture" "$work/seven" "$work/fixture/seven.json" >&2
$python3 /harness/dump_vets.py "$work/fixture" default "$work/default" seven "$work/seven"
