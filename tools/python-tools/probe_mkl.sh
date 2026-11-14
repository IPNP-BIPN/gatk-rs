#!/bin/bash
. /opt/miniconda/etc/profile.d/conda.sh
conda activate gatk
for cbwr in unset COMPATIBLE COMPATIBLE,STRICT AVX2 AVX2,STRICT; do
  for t in 1 4; do
    echo "=== MKL_CBWR=$cbwr MKL_NUM_THREADS=$t"
    if [ $cbwr = unset ]; then env MKL_NUM_THREADS=$t python /harness/probe_mkl.py
    else env MKL_CBWR=$cbwr MKL_NUM_THREADS=$t python /harness/probe_mkl.py; fi
    [ $cbwr = COMPATIBLE,STRICT ] && env MKL_CBWR=$cbwr MKL_NUM_THREADS=$t GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4 python /harness/probe_mkl.py | sed 's/^/glibc-masked /'
  done
done 2>&1
python -c "import mkl; print(mkl.get_version_string())"
