# PostprocessGermlineCNVCalls against the frozen reference

Part of #549. The target is GATK 4.7.0.0 in `broadinstitute/gatk:4.7.0.0` **under
`tools/python-tools/gcnv_frozen_env.sh`** ([gcnv-reproducible-oracle.md](gcnv-reproducible-oracle.md)):
pytensor's C at baseline x86-64 with `-ffp-contract=off` and no OpenMP, Ubuntu's reference BLAS over
MKL, numpy's baseline loops, glibc's libm without FMA or AVX2. Under GATK's own configuration the
bytes depend on the CPU, so there is nothing CPU-independent to be identical to.

## The three outputs

| Output | Reference | Port |
| --- | --- | --- |
| Intervals VCF | `GermlineCNVIntervalVariantComposer`, Java | `runners_postprocess_gcnv.rs` |
| Denoised copy ratios | `concatenateDenoisedCopyRatioFiles`, Java | `runners_postprocess_gcnv.rs` |
| Segments VCF | `segment_gcnv_calls.py` (gcnvkernel's Viterbi segmentation), then `GermlineCNVSegmentVariantComposer` | `gatk_tools::gcnv_segments`, then `runners_postprocess_gcnv.rs` |

`python` must be on the path and `gcnvkernel` must import, as for the reference; without them the port
gives the reference's refusal, which the `interpreter-tools` golden holds. No Python runs beyond that
check.

## What the segmenter's numbers are

Four compiled pytensor functions (HMM specs, forward-backward, Viterbi, constrained path log
probability) and the qualities numpy and scipy compute around them. What is reproduced is the graph
`fast_run` leaves and the C pytensor generates for it, read from the module cache:

* every `logsumexp` is rewritten twice over, so it is `m2 + log(sum(e)) + m1` with `m1 = max(x)`,
  `m2 = max(x - m1)`, `e = isinf(m2) ? exp(m2) : exp(x - m1 - m2)`, added left to right;
* `CAReduce` sums run along the axis in order; the one `CGemv` has an inner dimension of two;
* `exp`, `log`, `log1p`, `expm1` are glibc 2.35's SSE2 builds (`gatk_tools::glibc_libm`, with
  `__log_data` read out of the image's `libm.so.6`); none is correctly rounded;
* the shards are parsed by pandas 2.1.3's `precise_xstrtod`, which is not correctly rounded either
  (it differs from the nearest double on 41 % of `log_q_tau_tk`);
* numpy 1.26.2's pairwise sum, `scipy.special.logsumexp`, `np.logaddexp`, Python's `'{:.2f}'`, then
  Java's `FastMath.round` and htsjdk's QUAL formatting;
* gcnvkernel's own off-by-one in the end quality (`beta[end + 1]`) is kept.

## What is proven, and where

Locally, against the reference's postprocessing under the frozen configuration, on three sets of
shards for the six samples of the fixture `make_gcnv_fixture.py` writes:

| Shards | Intermediates (`gcnv_segments_trace.rs`) | Final files (`gcnv_postprocess_end_to_end.rs`) |
| --- | --- | --- |
| `netlib1`: the trio run under `gcnv_frozen_env.sh` | 18 HMM specs, 18 forward-backward, 18 Viterbi paths, 202 constrained log probabilities, 412 qualities: every bit | 18 of 18 files byte-identical |
| `default1`: the trio run under GATK's own flags, postprocessed frozen | the same counts, every bit | 18 of 18 |
| `frozen1`: an earlier frozen run whose shards differ from `netlib1`'s | the same counts, every bit | 18 of 18 |

`glibc_libm`'s four functions match numpy in the image on 1.2 million points each.

**On real x86-64**: the `gcnv-frozen` golden is CI run 42048939031's (Xeon 6973P-C). The trio's 132
files have the same bytes there as in run 42040149376 (Xeon Platinum 8573C) and in the emulated
AVX2-only run of [gcnv-reproducible-oracle.md](gcnv-reproducible-oracle.md), and `gcnv_postprocess_end_to_end.rs` replays its shards and gets the 18
Postprocess files byte for byte. All three CPUs are x86-64 with or without AVX-512; an AMD runner has
not been seen yet.

Not ported: `--input-intervals-vcf` and `--clustered-breakpoints` (the joint-calling path) and
`--reference`, which the port refuses in its own words.

## Reproducing

Inside `broadinstitute/gatk:4.7.0.0` with `tools/python-tools` at `/harness`, after
`. /harness/gcnv_frozen_env.sh`:

```
python /harness/pp_trace.py <segment_gcnv_calls.py> <trace.tsv> <graph.txt> \
  --ploidy_calls_path <run>/ploidy/p-calls --model_shards <run>/gcnv/g-model \
  --calls_shards <run>/gcnv/g-calls --output_path <dir> --sample_index <i>
python /harness/probe_gcnv_libm.py <dir>
```

`segment_gcnv_calls.py` is the one GATK extracts from its jar. Then, from the repository:

```
GCNV_PP_TRACE=<dir holding s<i>/trace.tsv> GCNV_PP_RUN=<run> \
  cargo test --release -p gatk-tools --test gcnv_segments_trace
GCNV_LIBM_VECTORS=<dir> cargo test --release -p gatk-tools --lib glibc_libm
GATK_RS_GCNV_DUMP=<gcnv_dump.sh output> \
  cargo test --release -p gatk-cli --test gcnv_postprocess_end_to_end
```
