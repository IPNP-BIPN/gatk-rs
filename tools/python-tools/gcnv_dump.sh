#!/bin/bash
# The `gcnv-frozen` suite's dump, inside broadinstitute/gatk:4.6.2.0: one pass of the gCNV trio over
# the fixture under gcnv_frozen_env.sh, then the CPU it ran on and one line per file it wrote.
# What it is for: the same lines from runners of different CPU families, or the configuration fails.
set -e
. /harness/gcnv_frozen_env.sh
work=/tmp/gcnv-frozen
python3=/opt/miniconda/envs/gatk/bin/python
mkdir -p "$work/fixture"
$python3 /harness/make_gcnv_fixture.py "$work/fixture" >&2
bash /harness/measure_gcnv.sh "$work/fixture" "$work/run" >&2
$python3 /harness/dump_gcnv.py "$work/run"
