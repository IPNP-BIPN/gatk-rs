# Sourced: the configuration under which the gCNV reference is meant to compute the same bits on
# any x86-64 runner. docs/gcnv-reproducible-oracle.md says what each line removes and how it was checked.
# pytensor adds its -march=native equivalent (-mfma ... on this host) unless gcc__cxxflags names a
# -march (pytensor/link/c/cmodule.py, compile_args); baseline x86-64 has no FMA, and contraction is off.
export PYTENSOR_EXTRA="gcc__cxxflags=-march=x86-64 -mtune=generic -ffp-contract=off,openmp=false"
# The image's BLAS is MKL 2022.2.1 (libblas.so.3 -> libmkl_rt.so.2); MKL picks kernels per CPU.
# COMPATIBLE and one thread make it repeatable, but past about 40 terms its order of addition is
# MKL's own and unknown. The image also carries Ubuntu's reference BLAS (libblas3 3.10.0, with the
# cblas_ symbols), whose dgemm adds the products in order, one at a time: preloaded, it answers every
# BLAS call numpy and pytensor make. GCNV_BLAS=mkl keeps MKL.
export MKL_CBWR=COMPATIBLE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE OMP_NUM_THREADS=1
if [ "${GCNV_BLAS:-netlib}" = netlib ]; then
  export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
fi
# numpy 1.26.2: only its baseline (SSE SSE2 SSE3) loops.
export NPY_DISABLE_CPU_FEATURES="SSSE3 SSE41 POPCNT SSE42 AVX F16C FMA3 AVX2 AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR"
# glibc 2.35's libm takes its FMA/AVX2 variants of exp, log, pow, sin... through ifunc, and MKL calls
# that libm (libmkl_core imports exp, log, pow): masking AVX2 or FMA alone moves MKL's COMPATIBLE dgemm.
export GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4,-AVX512F
