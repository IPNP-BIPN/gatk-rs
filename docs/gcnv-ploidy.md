# DetermineGermlineContigPloidy against the frozen reference

Part of #549, outside milestone X. The target is GATK 4.7.0.0 in `broadinstitute/gatk:4.7.0.0`
**under `tools/python-tools/gcnv_frozen_env.sh`** ([gcnv-reproducible-oracle.md](gcnv-reproducible-oracle.md)),
as for [PostprocessGermlineCNVCalls](gcnv-postprocess.md). Cohort mode only.

## What runs where

| Step | Reference | CLI | Library |
| --- | --- | --- | --- |
| Import check, then `doWork`: validation, intervals, counts summed per contig, the script's arguments | Java | `runners_gcnv.rs` | |
| `cohort_determine_ploidy_and_depth.py` / `case_determine_ploidy_and_depth.py`: the model, ADVI, sampler, caller, writers | gcnvkernel 0.9, pymc 5.10.1, pytensor 2.18 | the same scripts, run by `python` | `gatk_tools::gcnv_ploidy`, cohort mode |

The CLI runs GATK's own scripts as GATK does, so every mode and argument reaches gcnvkernel; its
covering array is measured past the interpreter under the frozen environment. The native port is a
library, kept and checked against the frozen golden; the CLI does not call it.

## What the numbers are

Six compiled pytensor functions: two initial-point functions, the ADVI step, the log-emission
sampler (a `Scan` over the draws), the caller, and `softplus(rho)` for the writers. The port is the
C pytensor generated for each, op for op, read from the graph dump of
`tools/python-tools/ploidy_trace.sh` (which also records every call's inputs and outputs):

* constants as pytensor folded them: Exponential scales kept in float32 (`1/f32(1e-4)` is
  `10000.000252621256`), `log(0.25) - log(0.75)` for the caller's admixing, `1/f32(√2)` in the
  Hellinger distance, `1 - 0.9 = 0.09999999999999998` in Adamax;
* `CAReduce` sums in index order from zero; the contig exclusion `dot` with its inner index
  ascending (with three contigs every order gives the same bits);
* glibc 2.35's SSE2 `exp`, `log`, `log1p`, `expm1` (`glibc_libm`) and its fdlibm `lgamma`
  (positive arguments only, checked bit for bit against the image), pytensor's own `_psi`;
* the RNGs: `RandomState(1984)` hands out three seeds; the step draws its noise from
  `SeedSequence(seed).spawn(1)[0]`, the sampler from `spawn(4)[0]`, both PCG64 with numpy 1.26's
  ziggurat normal (`numpy_random`, `numpy_ziggurat`);
* the numpy driver around them: the running mean of the log emission, `np.median` of the relative
  updates, `np.mean` (pairwise) of the caller's Hellinger distances;
* Python's `repr` of floats, `write_ndarray_to_tsv`, `json.dump(indent=1)`.

## Out of the library's reach

* `--model` (case mode, `case_determine_ploidy_and_depth.py`);
* `--disable-annealing`, `--disable-caller`, `--disable-sampler`;
* `--mean-bias-standard-deviation` other than 1 (the truncated normal's constant is folded for 1);
* any run whose total ADVI iterations reach `--convergence-snr-averaging-window`: the ELBO
  convergence tracker (a recursive least squares through MKL's matrix inverse) is not ported, and
  below that window it can never stop ADVI. GATK's defaults (1000 + 99 × 1000 iterations, window
  5000) are refused for this reason.

## What is proven, and where

* `crates/gatk-cli/tests/gcnv_ploidy_end_to_end.rs` makes the fixture (`make_gcnv_fixture.py`),
  runs `measure_gcnv.sh`'s command line through the CLI with a stand-in `python` that keeps the
  files and arguments the CLI hands the script, feeds those to `gatk_tools::gcnv_ploidy`, and
  compares all 38 `ploidy/` files of the `gcnv-frozen` golden byte for byte.
* `crates/gatk-tools/tests/gcnv_ploidy_trace.rs` replays a trace of the reference
  (`GCNV_PLOIDY_TRACE=<ploidy_trace.sh output>/trace`) call by call: every step's loss, `mu` and
  `rho`, every draw, every caller update. It is a debugging aid and skips without the variable.

Not shown: more than three contigs (where the `dot`'s order starts to matter), other fixtures, or
GATK's default iteration counts.
