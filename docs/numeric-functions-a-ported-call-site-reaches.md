# Which numeric functions a ported call site reaches, and what each one is worth

htsjdk-rs decision 0023 replaced "the jmath corpus reaches 100%" — unreachable by construction,
because its columns are `java.lang.Math` and the remaining divergent functions can only be made
exact by transcribing GPL2 source — with a rule:

> every function a ported call site reaches is exact, and every one that cannot be is named at the
> call site.

The rule was written down and the list was not. This is the list, produced by walking the call
sites rather than the library: what a ported path actually reaches, and what is known about each.

## Reached through `jmath`

Counts are occurrences in `crates/*/src` on 2026-11-25, comments included.

| function | occurrences | status |
|---|---|---|
| `math::log` | 52 | **exact**. Correctly rounded in the reference, so rounding the true result suffices (decision 0006) |
| `math::log10` | 96 | **exact**, for the same reason |
| `math::log1p`, `math::expm1` | 6, 3 | **exact on every host**. `Math.log1p` and `Math.expm1` delegate to `StrictMath` and have no intrinsic, so FDLIBM is the reference everywhere (htsjdk-rs decision 0051, IPNP-BIPN/htsjdk-rs#301) |
| `math::sqrt` | 2 | **exact**. IEEE 754 mandates the rounding, so every implementation already agrees |
| `math::round`, `fast_math::round` | 4 | **exact**, oracle-backed |
| `strict_math::exp` | 1, in `NaturalLogUtils` | **exact** against `StrictMath.exp`; the call site's own reference is `Math.exp`, which is **1 ulp** away (decision 0025) |
| `gamma::*` | 1 | **exact**, oracle-backed (`gamma-erf-normal`) |
| `normal::*` | 1 | **exact**, oracle-backed |
| `percentile::*` | 1 | **exact**, oracle-backed |
| `saddle_point::hypergeometric_log_probability` | 1 | **exact**, oracle-backed |

## Reached through the host libm

These are the ones the rule exists for. Since 2026-11-25 there are exactly two ways a ported call
site reaches the host libm, and both are named functions, so `git grep java_math_` lists every
such site but one, named below the table:

| named function | the Java it stands for | sites | status |
|---|---|---|---|
| `gatk_engine::natural_log_utils::java_math_exp` | `Math.exp` | 11 | **exact on aarch64**, where it is `strict_math::exp` (0 of 2,000,000 differ, below). **On x86-64 it is the host `exp`**, ~1 ulp from the HotSpot `dexp` stub on 0.0289% of the corpus; the stub is GPL2 and cannot be transcribed (decision 0014), and fdlibm was measured to be the worse stand-in on three suites (below) |
| `gatk_engine::math_utils::java_math_pow` | `Math.pow` | 28, plus `pow10`'s 34 | **the host `pow` on every target**. On x86-64 `Math.pow` is the `dpow` stub (99.9378% agreement on the corpus); fdlibm is 1 ulp from it at worst (decision 0027) and was **measured to be worse here**: switching `pow10` broke a passing byte-identity claim. On aarch64 `Math.pow` *is* fdlibm (0 of 4,000,000 pairs differ, decision 0051), but the goldens are x86-64 Java's and the host is the closer stand-in for them there too, so the helper does not switch |

What goes through them, by Java method:

- `java_math_exp`: `MathUtils.normalDistribution` (activity profile), `NaturalLogUtils.normalizeFromLogToLinearSpace` for LROM's E-step, FilterMutectCalls' posteriors and the clustering EM (measured, below), `NaturalLogUtils.logSumExp` in ModelSegments, `MathUtils.log1mexp` in DRAGstr, `BetaBinomialDistribution.probability`, `SomaticClusteringModel`'s printed weights, `GCBiasCorrector`'s decay weights, `PairHMMLikelihoodCalculationEngine`'s PCR indel model.
- `java_math_pow`: `MathUtils.pow10` and every `Math.pow(10.0, x)` it stands for (QualityUtils, BAQ, the genotype likelihood tables, the GQ of `GenotypeLikelihoods.getGQLog10FromLikelihoods` in subsetting, posteriors and GVCF merging, `log10sumLog10`, `normalizeFromLog10`), `AllelePseudoDepth.calculateWeights` (twice, by `weightDecay`), the PairHMM's `log10OneMinusPow10` and Jacobian table, `FisherExactTest`, `AddFlowSNVQuality`, `FlowBasedRead`, `CallCopyRatioSegments`' `Math.pow(2, x)`, DRAGstr's `gpMin` and `MathUtils.doubles`, and VQSR's model.

Why the helpers keep the host `exp` and `pow` on x86-64 bit for bit: `java_math_pow` is
`#[inline(always)]`, so each call site is the same IR as the `powf` it replaced, and LLVM's
library-call rewrites (`pow(2, x)` to `exp2`, `pow(10, x)` to `exp10`) apply exactly as when the
goldens were measured; the sites that used `black_box` still do.

One `Math.exp` site keeps the host `exp` on every target, by name:
`allele_frequency_calculator.rs`, `log1mexp_with_host_exp`, for `MathUtils.log10OneMinusPow10`.
A no-variant posterior of about 1e-19 passes `exp(a)` through `log1p` unchanged, so the `exp`
reaches the output as it is: fdlibm's answer is two ulp from the x86-64 golden's and the host's is
identical, on aarch64 as well, so `java_math_exp`'s aarch64 fdlibm turns this suite red locally
(measured 2026-11-25).

### Host calls that remain and are exact

`arithmetic.json` still counts these, because the audit counts every host transcendental. Each is
exact against what its Java computes, measured on `gatk-rs-oracle:4.7.0.0` (`-Xbatch` and `-Xint`,
2026-11-25):

| call site | Java | why it is exact |
|---|---|---|
| `allele_frequency_qc.rs:33,39`, `call_copy_ratio_segments.rs:153`, `series_stats.rs:182` | `Math.pow(x, 2.)` or `x * x` | `powi(2)` is `x * x`; `Math.pow(x, 2.0) == x * x` on 2,000,000 of 2,000,000 inputs |
| `family_priors.rs:318` | `Math.pow(deNovoPrior, mvCount)`, `mvCount` 1 or 2 | `Math.pow(x, 1.0) == x` and `Math.pow(x, 2.0) == x * x` on 2,000,000 inputs each |
| `mann_whitney.rs:167` | `Math.pow(count, 3)` | integer counts, exactly representable: equal on every count 0 to 200,000 |
| `variant_eval_engine.rs:664,814` | `Math.pow(10.0, -6)`, `Math.pow(10, n)` | equal to `powi`'s repeated multiplication (and `1 / 10^n` for negative `n`) for every `n` in -22..=22 |
| `pair_hmm.rs:19`, `sv_contig_alignments.rs:799` | `Math.pow(2, 1020)`, `Math.pow(2, 971)` | a power of two is exact: `Math.pow(2, k)` equals `2^k` for every `k` in -1022..=1023 |
| `gatk-barclay/src/lib.rs:576` | `Float.parseFloat` of a hex literal, not a `Math` call | scaling by `2^exponent` is exact while the exponent is in the normal range; outside it the parse is a different question from this list's |
| `allele_pseudo_depth.rs:470` | none | a unit test, not a call site |

Nothing else in `crates/*/src` calls a host transcendental outside a comment.

### `java_math_exp`: exact on aarch64, host libm on x86-64 (2026-11-23)

Java's `Math.exp` is not one function. On x86-64 HotSpot replaces it with its `dexp` stub; on
aarch64 there is no such stub, so it is fdlibm, bit for bit `StrictMath.exp`. The rows above go
through `natural_log_utils::java_math_exp`, which picks at compile time (`target_arch`, never a
runtime CPU feature):

| target | `java_math_exp` is | inputs differing from that target's Java 17 `Math.exp`, of 2,000,000 x <= 0 |
|---|---|---:|
| aarch64 | `jmath::strict_math::exp` | **0** (glibc's `exp` there: 192,457) |
| x86-64 | the host libm's `exp`, unchanged | 4,860 (glibc, FMA variant; 4,494 without FMA; 2026-11-20 sweep, other inputs) |

Exact on x86-64 would need a transcription of the HotSpot stub, which decision 0014 withdrew for
licence reasons. The measurements below are all x86-64 ones and still hold there. The x86-64
goldens are x86-64 Java's output, so on aarch64 two suites differ from them exactly where x86-64
and aarch64 Java do, and their tests name those rows rather than skip: `activity_profile.rs` (the
kernels at sigma 17 for sizes 25 and 50, and the band-pass profile built on one) and
`mutect_engine_arithmetic.rs` (the posterior at `100.0,-1.0`). Re-running both dumps on arm64 Java
17 and pointing the tests at that output passes every row.

### The activity profile's `exp`, measured (2026-08-19)

The row said the choice was worth revisiting once the call site had a byte claim. It has one:
`MathUtils.normalDistribution` is the only `exp` a band-pass profile reaches, and the
`activityprofile` suite pins **twenty kernels, 266 values, as raw bits** against the oracle.

So the question is answerable rather than arguable, and the answer is that the current choice is
correct:

| what the kernel is built with | values matching the oracle bit-for-bit |
|---|---|
| the host `exp` (what the port calls) | **266 of 266** |
| `jmath::strict_math::exp` | 256 of 266 — ten wrong, by one ulp, in two of the twenty kernels |

The Java is `Math.exp`, not `StrictMath.exp` (`MathUtils.normalDistribution`, line 949), which is
why the exact-against-`StrictMath` function is the *worse* choice here: it is faithful to a
different reference. `Math.exp` is the HotSpot intrinsic, decision 0014 withdrew its transcription,
and on this call site's inputs the host libm happens to agree with it everywhere.

This is the same shape as the three `powf` rows, and it makes the shape a rule rather than a
coincidence: **where the reference is `Math.*`, a bounded `StrictMath.*` port is not an improvement
unless it is measured to be one.** All three measurements so far went the other way.

The claim is not that the host `exp` equals `Math.exp` everywhere -- it does not, 99.9711% is the
rate. It is that on the points this call site reaches it does, and the suite is what holds that: a
future swap moves kernel bits and turns the suite red, which is the guard this row wanted.

### The read orientation model's `exp`, measured (2026-10-01)

`LearnReadOrientationModelEngine.computeResponsibilities` normalises twelve log posteriors with
`NaturalLogUtils.normalizeFromLogToLinearSpace`, so every E-step of every context goes through
`Math.exp` once per non-impossible state and once more per term of the log-sum. Its covering array
pins whole prior tables as text.

| what the E-step is built with | covering-array rows matching the oracle |
|---|---|
| `jmath::strict_math::exp` (what `NaturalLogUtils` calls) | 6 of 8 (two rows, one prior each, one ulp out) |
| the host `exp` | **8 of 8**, and 12 of 12 on the widened array |

The reference was then run on the two failing rows with `-XX:+UnlockDiagnosticVMOptions
-XX:-UseLibmIntrinsic`, which sends `Math.exp` to the shared runtime rather than the intrinsic: its
priors equal the fdlibm build's byte for byte. So nothing else in the EM differs, and the choice
between the two stand-ins is the whole question. As for the activity profile, the host `exp` is
the closer one on the points this call site reaches. Every other caller of
`normalize_from_log_to_linear_space` keeps fdlibm, because none of them has been measured.

### FilterMutectCalls' `exp`, measured (2026-11-19)

FilterMutectCalls learns two passes of state (the strand artifact EM, `SomaticClusteringModel`) and
prints the learned priors in its filtering stats to the last digit, so an ulp in pass 0 is still
visible at the end. A harness in the filters' package (reflection on `Mutect2FilteringEngine`, every
filter's answer and every learned value as raw bits, per pass) was compared against a traced port
on `mutect.vcf`, at contamination 0.1 with the orientation priors, and at 0.0 without them:

| what the port is built with | trace lines matching the oracle bit-for-bit (0.1 + priors / 0.0) |
|---|---|
| fdlibm in `posterior_probability_of_error` and in the clustering EM (and the strand fix below not yet in) | not counted, from pass 0 on / 10 of 314 differ |
| host `exp` in `posterior_probability_of_error`, fdlibm in the EM, strand fix in | 66 of 350 differ / 0 of 314: pass 0's answers all match, the EM's clusters and weights still differ from pass 0's learning on |
| host `exp` in both | **350 of 350 / 314 of 314** |

(The 0.0 column also needed a separate fix that is not a stand-in: the port computed the strand
filter's accumulated E-steps with the initial prior and shape instead of the learned ones.) With
both, the FilterMutectCalls t=2 covering array goes from 24/25 to 25/25: row 16's `Ln prior of SNV`
was the one row that printed the drift. As before, this is a claim about the points these call
sites reach, not that the host `exp` equals `Math.exp`.

## What this changes

**The gap is two named functions, not a list of call sites.** Every ported call site that reaches
a function no permissive source can make exact goes through `java_math_exp` or `java_math_pow`, and
both name the measured exception in their documentation. Everything else is exact: `log`, `log10`,
`log1p`, `expm1` and `sqrt` on every host, the `powi` sites by measurement, the commons-math and
saddle-point ports by their oracle suites.

**`BinomialDistribution` and SVD, the other two names H.5 used to carry, are ports rather than
libm calls.** SVD is commons-math3 3.5's `SingularValueDecomposition`
(`gatk-engine/src/singular_value_decomposition.rs`, reached through `KernelSegmenter`) and the
binomial is `jmath::binomial`; whatever they reach of `Math` is in the tables above.
