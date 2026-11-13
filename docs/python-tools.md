# The tools that run their arithmetic in Python

Issues #806 and #1394. Reference: `broadinstitute/gatk:4.7.0.0`, whose environment is 4.6.2.0's
package for package (Python 3.10.13, numpy 1.26.2, scikit-learn 1.3.2, h5py 3.10.0, dill 0.3.7,
pymc 5.10.1, pytensor 2.18.3). The pinned oracle image has no `python`, so every number here comes
from the publisher's image.

The measurements below were made under amd64 emulation on an arm64 Mac, whose emulated CPU reports
AVX2 and FMA and no AVX-512, against 4.6.2.0. CI has since re-derived them on real x86-64: suite
`vets-isolation-forest` is oracle-backed and held through the move to 4.7.0.0 unchanged. See
[What needs CI](#what-needs-ci).

## One line per tool

| Tool | Verdict | Measurement |
| --- | --- | --- |
| `TrainVariantAnnotationsModel` | **Byte-identical reachable, and landed**, for `PYTHON_IFOREST` with a seed | h5py score files identical; scorer identical but for its random script name |
| `ScoreVariantAnnotations` | **Byte-identical reachable, and landed**, for `PYTHON_IFOREST` | VCF identical; `.scores.hdf5` doubles identical (the reference's file holds the time) |
| `DetermineGermlineContigPloidy` | **Keeps the interpreter** | Stable, but the numbers come from C compiled per host CPU from a rewritten graph |
| `GermlineCNVCaller` | **Keeps the interpreter** | Same |
| `PostprocessGermlineCNVCalls` | **Keeps the interpreter for segments.** Intervals VCF and denoised table are Java and reachable, not landed | Stable on fixed shards; the segmenter is the same compiled C |

No tool is in the "only a stated numeric tolerance" class. A tolerance would need a reimplementation
of the inference to measure against, and none exists.

## Stability of the reference

### Variant annotations

2000 sites (740 SNP training, 234 indel training) from `make_vets_fixture.py`. Seven runs: two plain,
one CPU with MKL and OpenMP at 1, eight CPUs, another `PYTHONHASHSEED`, `n_jobs=4`, and `random_state=1`.

| File | Across the runs |
| --- | --- |
| `extract.vcf`, `score.vcf` | Same bytes, except that `random_state=1` moves `score.vcf`, as it should |
| `train.{snp,indel}.{training,calibration}Scores.hdf5` (h5py) | Same bytes, except under `random_state=1` |
| `*.annot.hdf5` | 12 to 24 bytes differ from offset 1988: the time Java's HDF5 library wrote |
| `score.scores.hdf5` | 1 to 2 bytes differ at offset 1972: the same |
| `train.*.scorer.pkl` | 17 bytes differ at offset 163: the script's name, `/tmp/isolation-forest.<random>.py`, in the lambda's code object |

A name one digit shorter moves the frames, so a byte diff can run to 505017 bytes while the opcodes
inside the frames are the same but for that name. `PYTHONHASHSEED` was not compared beyond this: the
port's scorer was compared against the plain, one-CPU, `n_jobs=4` and `random_state=1` runs.

### gCNV

Six samples, 180 intervals, three contigs, 200/50 ADVI iterations, three epochs, from
`make_gcnv_fixture.py`. Four runs: two on 4 CPUs, one on 1 CPU with MKL and OpenMP at 1, one on 8 CPUs
with both at 8. **All 132 output files have the same bytes in all four**, including the Postprocess
VCFs. Each stage took about 5 minutes emulated, Postprocess about 35 seconds.

`PostprocessGermlineCNVCalls` alone, four times on one run's shards at 4, 4, 1 and 8 CPUs: the three
outputs have the same bytes every time.

So stability is not what stops these tools. Those are few iterations on a small fixture, which is
evidence and not proof for a real cohort.

## What now matches, and how it was compared

`tools/python-tools/compare_port.py` against reference runs; the port ran with a stand-in `python`
that imports anything.

| Reference run | `score.vcf` | Four h5py score files | Two scorers |
| --- | --- | --- | --- |
| Default | identical, 238774 bytes | identical | identical but the name, 1299191 and 1235291 bytes |
| `random_state=1` | identical | identical | same, compared without frames (name of 19 digits) |
| `n_jobs=4` | identical | identical | same, compared without frames |
| One CPU | identical | identical | same, compared without frames |

The doubles of `score.scores.hdf5` were compared bit for bit for the default and `random_state=1` runs
(`crates/gatk-tools/tests/isolation_forest_measured.rs`). `2**y` agrees with numpy on all 2,000,000
points of `ys.bin`.

Not compared: the reference's `Score` loading a scorer the port wrote. The file differs from a
reference scorer only in digits of a file name, so this is expected to work and is untested. Hidden
by `%.4f`: a last-place difference in the Hipparchus distribution would not show in the VCF.

### What the port had to reproduce

* MT19937 and `sample_without_replacement` of numpy; the tree builder's `our_rand_r`, its splitter and
  float32 thresholds.
* **glibc's `pow`, FMA variant**, for `2**-depth/c`, taken from the disassembly; it differs from the
  SSE2 variant on 1458 of 2,000,000 points by one unit in the last place.
* **glibc's `log`** on the whole numbers 1 to 256, kept as a table. fdlibm differs at 3, 48, 74 and 196,
  which is what first failed the comparison.
* `pickle._Pickler` at protocol 4 with dill's reducers, including `_Framer`'s frame boundaries.
* Hipparchus' `EmpiricalDistribution` for the calibration sensitivity.

### What is still refused

Each ends with the port's own sentence.

* `PYTHON_SCRIPT`: a user's script is Python. There is nothing to port.
* `--unlabeled-annotations-hdf5` with `PYTHON_IFOREST`: the default script rejects the argument and the
  reference fails with exit status 2 and temporary names in its message.
* Hyperparameters whose effect was not measured: anything but `random_state`, `n_estimators`, `n_jobs`
  and the defaults of the others. Without a `random_state` a forest has no byte target.
* A scorer the default script at scikit-learn 1.3.2 did not write, or with `_max_samples` over 256.
* A `--double-format` other than `%.Nf`.
* `JAVA_BGMM` was already ported: the reference's own `NotImplementedException`, "BGMM module will be
  implemented in separate PR.". Nothing more is reachable.

`python` and the five packages must still be importable, because that refusal is the reference's
answer on a machine without them and the frozen `interpreter-tools` golden holds it. No Python runs
beyond those import checks.

## Why gCNV keeps the interpreter

**The reference's arithmetic is not in any file it ships.** At run time pytensor rewrites the model's
graph and compiles it: 409 C modules, 68 MB, in the measured run, with
`-march=westmere -mavx -mavx2 -mfma ... -mno-avx512f -fopenmp -fno-math-errno` read from the cache keys.
Those flags come from the host CPU. With `-mfma`, g++ 12.4.0 may contract products and sums, so the
rounding depends on the CPU and on the rewritten graph. A port would have to match the optimizer's
output and the compiler's contraction choices for every op, for each CPU family. The evidence is the
flags: a second CPU family was not available to measure the bytes moving.

### Survey (#1394)

What the model needs: `Normal`, `Lognormal`, `Exponential`, `Uniform`, `Deterministic`, `Potential`;
negative binomial and heavy-tail log densities; forward-backward and Viterbi through `scan` with
`logsumexp`; automatic differentiation of all that; mean-field ADVI with GATK's Adamax, annealing and
seeded sampling.

| Crate | What it offers | Against that |
| --- | --- | --- |
| `scirs2-stats` 0.4.0 | Mean-field ADVI with Adam over a log density you write | Closest match. Not Adamax, and no route to these bytes |
| `fugue-ppl` 0.1.0 | Monadic PPL listing VI among its methods | Early release; its own page warns about completeness |
| `modppl` 0.3.1 | Trace-based programmable inference | Research prototype |
| `opensrdk-probability` 0.8.9 | A `variational` module with a `stein` submodule | ADVI not confirmed |
| `linfa-clustering` 0.8.1 | A Gaussian mixture | Not variational; no PyMC equivalent |

This rests on two web searches and their snippets. No crate was built or run, and maintenance was not
checked. A library could only ever reach a tolerance, since none compiles the same graph.

**linfa against the reference's model was not run, because there is no reference GMM to compare with.**
`JAVA_BGMM` is a stub, `BayesianGaussianMixtureModeller` is unreachable from the command line, the
default model is an isolation forest, and gCNV has no mixture of that form. A ULP distance from
scikit-learn's mixture would describe a model GATK does not run.

### `PostprocessGermlineCNVCalls`

`GermlineCNVIntervalVariantComposer` writes the intervals VCF from the shards' posteriors in Java, and
`concatenateDenoisedCopyRatioFiles` writes the table. Neither calls Python, so both are reachable and
not landed. Between them `segment_gcnv_calls.py` runs the compiled forward-backward and Viterbi, and
`GermlineCNVSegmentVariantComposer` only formats what that returned, so the segments VCF keeps the
interpreter. The tool would still end in a refusal.

## What #806 can and cannot claim

Can claim, **once CI confirms**: `TrainVariantAnnotationsModel` and `ScoreVariantAnnotations` reproduce
the reference for the default backend with a seeded forest, the scorer modulo its random script name.

Cannot claim: anything for the gCNV trio beyond the checks before the import, and any byte target for
a scorer's name or for the time inside a Java-written HDF5 file.

## What needs CI

* ~~That any of this holds on real x86-64.~~ Done: `vets-isolation-forest` is oracle-backed, in
  `broadinstitute/gatk:4.7.0.0` pulled by the job that runs it.
* **The arguments past the interpreter.** `TrainVariantAnnotationsModel`'s covering array runs its
  reference in the Broad's image and the port behind a stand-in `python` (`$past_the_interpreter` in
  `tools/coverage/fixtures.json`), so its rows reach the fit instead of stopping at the import. That
  turned up one divergence, now fixed: the scorer names the script under `--tmp-dir`, which the port
  had written as `/tmp` whatever the argument said. `ScoreVariantAnnotations` does the same, over
  scorers the reference trains into the corpus before the rows (`run_array.train_scorers`), and
  turned up a second divergence, now fixed: the port wrote the genotypes under
  `--sites-only-vcf-output`. Its `PYTHON_SCRIPT` rows still differ, as the port's own refusal.
* **AVX-512.** On an AVX-512 CPU numpy takes other loops for `power` (SVML) and `log`. The port follows
  the path without AVX-512. A first golden made on an AVX-512 runner differed from the port on every
  row, so `vets_dump.sh` sets `NPY_DISABLE_CPU_FEATURES` to the AVX-512 family: the reference then
  computes what it computes on an AVX2 CPU on any runner, and the dump's comment line shows it.
* **Stability on real hardware**, for both tool families.
* **That gCNV's bytes move across CPU families**, which is inferred from compile flags and not measured.

## Reproducing

`tools/python-tools/`: `make_vets_fixture.py` and `make_gcnv_fixture.py` write the fixtures,
`measure_vets.sh`, `measure_gcnv.sh` and `measure_postprocess.sh` run inside the image,
`compare_runs.py` compares runs, and `compare_port.py` compares the port with a run.
