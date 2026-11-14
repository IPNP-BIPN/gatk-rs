#!/bin/bash
# The `gcnv-frozen` suite's dump, inside broadinstitute/gatk:4.7.0.0: one pass of the gCNV trio over
# the fixture under gcnv_frozen_env.sh, then the CPU it ran on, one line per file it wrote, and
# the shards' bytes so the port can replay PostprocessGermlineCNVCalls on them.
# What it is for: the same lines from runners of different CPU families, or the configuration fails.
set -e
. /harness/gcnv_frozen_env.sh
work=/tmp/gcnv-frozen
python3=/opt/miniconda/envs/gatk/bin/python
mkdir -p "$work/fixture"
$python3 /harness/make_gcnv_fixture.py "$work/fixture" >&2
bash /harness/measure_gcnv.sh "$work/fixture" "$work/run" >&2
# Every sample through PostprocessGermlineCNVCalls, plain VCFs, for the port's end-to-end test.
for i in 0 1 2 3 4 5; do bash /harness/measure_postprocess.sh "$work/run" "post-s$i" "$i" >&2; done
rm -rf /tmp/pytensor-post-s*
$python3 /harness/dump_gcnv.py "$work/run"
