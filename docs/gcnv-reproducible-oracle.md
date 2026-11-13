# A reproducible gCNV reference, and what porting against it would take

Part of #806, following [python-tools.md](python-tools.md). Reference: `broadinstitute/gatk:4.7.0.0`
(gcnvkernel 0.9, pymc 5.10.1, pytensor 2.18.3, numpy 1.26.2, MKL 2022.2.1, glibc 2.35, g++ 12.4.0).
The measurements were made under 4.6.2.0; 4.7.0.0 ships the same stack, its `gcnvkernel` sources
hash identically (sha256 of the sorted `*.py` sums) and the same reference BLAS sits at the same path.

**Everything here was measured under amd64 emulation (Rosetta) on an arm64 Mac, whose emulated CPU
reports AVX2 and FMA and no AVX-512. Whether the frozen configuration gives the same bits on another
CPU family is exactly what it cannot show; the `gcnv-frozen` suite is there for CI to show it.**

## Verdict

* **A frozen configuration exists, uses only the image's own files, and is repeatable here**
  (`tools/python-tools/gcnv_frozen_env.sh`, below). Under GATK's own flags the same emulated CPU
  gives different bytes in 81 of 132 files once those knobs move, so they are not cosmetic.
* **Porting against it is feasible in principle, and only with the reference BLAS.** With MKL kept,
  even in `MKL_CBWR=COMPATIBLE` on one thread, a dgemm whose inner dimension passes about 40 adds its
  products in an order no simple scheme reproduces (best fit 68 % of entries), and the caller's
  gradients have inner dimensions of intervals and samples. Preloading Ubuntu's reference BLAS, which
  is already in the image, makes every product a plain left-to-right sum.
* **Cost: on the order of 3 to 4 months of one person, roughly 15,000 lines of Rust**, for 46
  compiled functions (1,840 nodes) plus pymc's ADVI driver, three RNGs, glibc's libm and the I/O.
  Start with `PostprocessGermlineCNVCalls`' segments: 11 functions, 390 nodes, no RNG, one BLAS call.
* **No golden is committed.** The suite is `golden-pending`; a golden is worth something only once two
  runners of different CPU families have produced the same lines.

## The frozen configuration

`gcnv_frozen_env.sh`, sourced before `gatk`:

| Knob | Value | What it removes |
| --- | --- | --- |
| `PYTENSOR_FLAGS` (via `PYTENSOR_EXTRA`) | `gcc__cxxflags=-march=x86-64 -mtune=generic -ffp-contract=off,openmp=false` | The host's `-march=native` equivalent and `-mfma`, so no contraction in the 408 compiled modules |
| `LD_PRELOAD` | `/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0` (Ubuntu `libblas3` 3.10.0-2ubuntu1) | MKL, for every `cblas_*` numpy calls and every `dgemm_`/`dgemv_` pytensor's C calls |
| `MKL_CBWR`, `MKL_NUM_THREADS`, `MKL_DYNAMIC`, `OMP_NUM_THREADS` | `COMPATIBLE`, `1`, `FALSE`, `1` | MKL's per-CPU kernels and thread splits, for whatever still reaches MKL (and for `GCNV_BLAS=mkl`) |
| `NPY_DISABLE_CPU_FEATURES` | everything above SSE3 | numpy's dispatched loops (AVX-512 `exp`/`log`, SVML) |
| `GLIBC_TUNABLES` | `glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4,-AVX512F` | libm's FMA/AVX2 ifunc variants of `exp`, `log`, `pow`... |

GATK's scripts prepend their own defaults (`openmp=true,blas__ldflags=-lmkl_rt,optimizer=fast_run,...`)
and append the user's `PYTENSOR_FLAGS`, so ours win. `GCNV_BLAS=mkl` keeps MKL, for comparison.

### Each knob, checked to take effect

`probe_fp.sh` runs `probe_fp.py` (hashes of libm, numpy ufuncs, MKL products and a small pytensor
function) with each knob alone; `probe_mkl.sh` and `probe_mkl_order.py` look at MKL.

* **pytensor does drop native detection when a `-march` is given.** `compile_args` in
  `pytensor/link/c/cmodule.py` sets `detect_march = False` as soon as `gcc__cxxflags` holds a
  `-march=`. Default: `-march=westmere -mavx -mavx2 -mfma ... --param l2-cache-size=8192 -mtune=westmere`;
  frozen: `-march=x86-64 -mtune=generic -ffp-contract=off -DNPY_NO_DEPRECATED_API=... -m64 -fPIC`.
  The cache keys of the gCNV runs say the same (409 modules, 68 MB, default; 408, 66 MB, frozen).
  The flags alone change a three-op elementwise function (`u*v + exp(u/10)*v - log1p(v)`): the
  compiler was contracting.
* **The BLAS is MKL**: `libblas.so.3`, `liblapack.so.3`, `libcblas.so.3` all resolve to
  `libmkl_rt.so.2` (conda's `libblas 3.9.0-16_linux64_mkl`, `mkl 2022.2.1`). MKL's answer moves with
  the thread count (dgemm and ddot, even under `COMPATIBLE,STRICT`), with `MKL_CBWR`
  (unset/`AVX2`, `COMPATIBLE`, `SSE4_2`, `AVX2,STRICT` give four different dgemms), not with the
  alignment of operands or output (0 to 7 doubles off a 64-byte line: one hash each).
* **MKL depends on glibc's libm.** `libmkl_core.so.2` imports `exp`, `log`, `pow`; masking AVX2 or
  FMA alone in glibc moves `COMPATIBLE` dgemm, unrelated tunables do not. So the glibc mask is
  needed even with MKL pinned, and is one more reason to take MKL out.
* **The reference BLAS adds in order and has no FMA**: 0 `vfmadd`-family instructions in its
  disassembly, and with it preloaded `MKL_VERBOSE` prints nothing and every probed product
  (m,k,n from (6,3,3) to (6,180,5)) equals `sum over k from 0 up` on 100 % of entries.
* **glibc's tunable takes effect**: `ld.so --list-diagnostics` drops the `x86-64-v3` hwcaps subdir
  (`0x6` to `0x4`), and `exp`, `log`, `pow` from `libm.so.6` hash differently (`log1p`, `expm1` do not:
  they have no ifunc variant).
* **numpy's mask takes effect** (`__cpu_features__` keeps only MMX SSE SSE2 SSE3), and changes nothing
  on this CPU: without AVX-512 numpy's float64 `exp`/`log` already call libm.

### MKL's order of addition (why MKL has to go)

`probe_mkl_order.py`, under the frozen environment with MKL kept, share of entries a candidate order
reproduces:

| m, k, n | in order | FMA chain | best blocked scheme tried |
| --- | --- | --- | --- |
| 6, 3, 3 / 180, 5, 6 / 6, 5, 180 / 5, 6, 180 / 180, 6, 5 | **1.000** | 0.54 to 0.62 | |
| 6, 8 / 12 / 17 / 40, 5 | **1.000** | | |
| 5, 180, 6 | 0.200 | 0.200 | 0.30 (blocks of 4) |
| 6, 180, 5 | 0.600 | 0.100 | 0.53 |
| 180, 180, 6 | | | 0.22 (blocks of 16) |
| 6, 1000, 180 | | | 0.68 (blocks of 128) |

Blocks of 2 to 512 terms, each summed in order or over 2 or 4 interleaved accumulators, then added in
sequence or as a tree. Under GATK's flags (AVX2 kernels) the same products are an FMA chain on every
shape but the vector dot. Reproducing `libmkl_def`'s blocking would mean reading its disassembly, and
its blocking may follow libm results (above), so it is not pursued.

## Measurements of the trio

Fixture of `make_gcnv_fixture.py` (6 samples, 180 intervals, 3 contigs), `measure_gcnv.sh` (200/50 ADVI
iterations, 3 epochs), 132 output files, the pytensor cache left out. Each stage took 5 to 7 minutes
emulated on a shared machine.

| Run | Configuration | Docker CPUs |
| --- | --- | --- |
| `default1` | GATK's own | 4 |
| `mkl1` | frozen, MKL kept (before the preload was added) | 4 |
| `netlib1` | frozen | 4 |
| `netlib2` | frozen | 8, Java 4 |

| Pair | Same bytes | Differ |
| --- | --- | --- |
| `netlib1` / `netlib2` | **132** | **0** |
| `netlib1` / the suite's own run (`run_suite.py --suites gcnv-frozen`) | **132** | **0** |
| `default1` / `mkl1` | 52 | 80 |
| `default1` / `netlib1` | 51 | 81 |
| `mkl1` / `netlib1` | 71 | 61 |

* Between GATK's flags and the frozen ones, every posterior moves: the ploidy model and calls
  (`contig_ploidy.tsv`, `mu_/std_psi_*`, `mean_bias_j`), the whole caller model, every sample's
  `log_q_c_tc`, `log_c_emission_tc`, denoised copy ratios, the ELBO histories and `post.denoised.tsv`.
  The two Postprocess VCFs stay the same: their quantized calls and qualities absorb these last digits
  on this fixture, which is luck and not a property.
* MKL against the reference BLAS moves only the caller and what follows it; the ploidy model's dots
  are 3 by 3, where MKL adds in order too.
* `PostprocessGermlineCNVCalls` alone on `netlib1`'s shards, frozen, three times (4 CPUs, Java on 1,
  and with MKL kept): the three outputs have the same bytes every time, and the same content as the
  full run's (whose VCFs are gzipped). The segmenter's one `CGemv` is small enough that MKL adds in
  order there too.

Earlier runs under GATK's flags (python-tools.md) were stable across 1, 4 and 8 CPUs on this host.
The point of the frozen configuration is not stability on one CPU, it is that every choice that
followed the CPU is now a constant.

## What actually runs: inventory of the compiled graph

`graph_log/sitecustomize.py`, put on `PYTHONPATH` with `GCNV_GRAPH_LOG` set, wraps
`FunctionMaker.create` once GATK's script has set its flags, and logs every node of every compiled
function with how its thunk runs (`c` for a compiled module, `py` for `perform`). From `mkl1` (the
graph does not depend on the BLAS):

| Tool | Functions | Nodes | `perform` | Scan | Random draws | BLAS ops | Fused composites |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DetermineGermlineContigPloidy | 11 | 326 | 23 | 3 | 5 | 3 | 52 |
| GermlineCNVCaller | 24 | 1124 | 67 | 9 | 11 | 9 | 150 |
| PostprocessGermlineCNVCalls | 11 | 390 | 7 | 6 | 0 | 1 | 44 |
| All | 46 | 1840 | 97 | 18 | 17 | 13 | 246 |

By op, C: `DimShuffle` 318, `Elemwise{Composite}` 246, `ScalarFromTensor` 177, `Subtensor` 124,
`Sum` 121, `Add` 70, `Max` 62, `Elemwise{exp}` 41, `IncSubtensor` 40, `IsInf` 31, `Gemm` 7, `CGemv` 2,
`Dot22` 1, `Softmax` 1, `MaxAndArgmax` 1, and shape bookkeeping. Python `perform`: `Alloc` 22,
`Scan` 18, `log` 14, `DeepCopyOp` 12, `NormalRV` 9, `Cast` 8, `ExponentialRV` 4, `Dot` 3 (numpy, so
BLAS), `Eye` 2, `Nonzero` 2, `TruncatedNormalRV`, `UniformRV`, `LogNormalRV` 1 each.

What a port has to reproduce bit for bit, given the frozen configuration:

* **The rewritten graphs, not the model.** pytensor's `fast_run` fuses elementwise chains into
  `Composite`s whose C is the scalar expression in a fixed order; with `-ffp-contract=off` and no
  fast-math, g++ keeps that order, so a port follows the generated expression. Dumping the final
  graphs (and the C of each composite) is the first tool to write. Risk: rewrites that depend on
  static shapes would make the graph depend on the cohort.
* **Reductions**: `Sum`, `Max` are pytensor's C `CAReduce`, a sequential loop in its own nesting
  order; `logsumexp` is gcnvkernel's max-shifted one. numpy's pairwise sum appears only where Python
  code reduces (postprocess, `scipy.special.logsumexp`).
* **BLAS**: 13 nodes; with the reference BLAS each entry is a left-to-right sum (`ddot` unrolls by 5,
  still left to right). With MKL: not reproducible, see above.
* **libm**: glibc 2.35's SSE2 variants of `exp`, `log`, `log1p`, `pow` (softplus and the log densities
  go through them). The isolation-forest port reproduced the FMA variant of `pow`; the mask makes the
  generic C the target, which is easier.
* **Random numbers**: gcnvkernel's `approx.rng` (MT19937 `RandomState`, seed 1984, already ported for
  the isolation forest) draws `randint(2**30)` seeds; pymc's `reseed_rngs` turns each into
  `SeedSequence(seed).spawn(n)` and `Generator(PCG64)`; draws are numpy's ziggurat normal and
  exponential, uniform, lognormal, and one truncated normal through `scipy.stats.truncnorm`
  (scipy's `ndtri`/`log_ndtr`, cephes C).
* **Scan**: 18 loops (forward-backward, Viterbi, posterior sample means), each an inner compiled
  function; the order is the loop's, so no difficulty beyond the inner graphs.
* **The driver**: pymc's mean-field ADVI, gcnvkernel's Adamax with annealing, epochs, convergence
  tests on floats, posterior sampling, and the TSV writers' float formatting.

## Estimate

| Piece | Weeks |
| --- | --- |
| Graph dump (final fgraphs and composite C, per tool and model configuration) | 1 |
| Op semantics in Rust, about 45 op types, pytensor's loop orders | 3 to 4 |
| glibc 2.35 libm, SSE2 variants of `exp`, `log`, `log1p`, `pow` | 1 |
| PCG64, `SeedSequence.spawn`, numpy's ziggurat, scipy's truncated normal | 1 to 2 |
| pymc ADVI, gcnvkernel's Adamax, annealing, sampling, convergence | 3 to 4 |
| Postprocess: Viterbi, segment qualities, numpy pairwise sums | 1 to 2 |
| I/O (TSV formats, configs, shards) and the Java ends | 1 |
| Chasing first divergences through thousands of ADVI iterations | 2 to 4 |
| **Total** | **13 to 19, about 15,000 lines** |

One divergent last bit anywhere in an iteration propagates through every later one, so there is no
partial credit inside a tool; each tool passes or fails whole.

**Recommendation.** Do not port against GATK's own configuration: it has no CPU-independent target.
First let CI run `gcnv-frozen` on two CPU families; if the lines agree, freeze that golden and start
with the segments of `PostprocessGermlineCNVCalls` (11 functions, 390 nodes, 6 scans, one `CGemv`, no
random draws), which closes the refusal left by #549 and builds the graph-dump and op machinery the
two inference tools need. If the lines disagree, the diff of the dump says which file, and the probes
above say which knob to look at first.

**Done for the segments**: [gcnv-postprocess.md](gcnv-postprocess.md).

**A frozen reference is a different reference.** A port byte-identical to it is not byte-identical to
GATK as shipped on any given machine, since the shipped configuration differs between machines
already. Whether the port then accepts only the frozen configuration, and refuses otherwise, is a
decision for #806.

## Reproducing

All inside `broadinstitute/gatk:4.7.0.0`, with `tools/python-tools` at `/harness`:

```
bash /harness/probe_fp.sh                 # each knob alone, hashes
bash /harness/probe_mkl.sh                # MKL: threads, CBWR, alignment
. /harness/gcnv_frozen_env.sh; python /harness/probe_mkl_order.py
. /harness/gcnv_frozen_env.sh; bash /harness/measure_gcnv.sh <fixture> <out>
PYTHONPATH=/harness/graph_log GCNV_GRAPH_LOG=<file> bash /harness/measure_gcnv.sh <fixture> <out>
bash /harness/gcnv_dump.sh                # the gcnv-frozen suite's dump
```

and `python3 tools/conformance/run_suite.py --suites gcnv-frozen` from the repository, which writes the
candidate to `tools/conformance/pending/` and compares nothing.
