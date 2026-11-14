#!/bin/bash
# probe_fp.py under the default environment, under the frozen one, and with each knob alone.
. /opt/miniconda/etc/profile.d/conda.sh
conda activate gatk
cd /tmp
# PYTENSOR_ADD is read by the child shell, so that `run name PYTENSOR_ADD=...` reaches the flags.
run() {
  echo "=== $1"; shift
  env "$@" bash -c 'PYTENSOR_FLAGS="base_compiledir=/tmp/pt-$$$PYTENSOR_ADD" python /harness/probe_fp.py pytensor' 2>&1 | grep -v WARN
}
run default
run glibc GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4
run mkl MKL_CBWR=COMPATIBLE MKL_NUM_THREADS=1
run numpy NPY_DISABLE_CPU_FEATURES="SSSE3 SSE41 POPCNT SSE42 AVX F16C FMA3 AVX2"
run cxx "PYTENSOR_ADD=,gcc__cxxflags=-march=x86-64 -mtune=generic -ffp-contract=off"
. /harness/gcnv_frozen_env.sh
run frozen "PYTENSOR_ADD=,$PYTENSOR_EXTRA"
