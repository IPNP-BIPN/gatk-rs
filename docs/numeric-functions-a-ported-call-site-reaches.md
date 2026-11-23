# Which numeric functions a ported call site reaches, and what each one is worth

htsjdk-rs decision 0023 replaced "the jmath corpus reaches 100%" — unreachable by construction,
because its columns are `java.lang.Math` and the remaining divergent functions can only be made
exact by transcribing GPL2 source — with a rule:

> every function a ported call site reaches is exact, and every one that cannot be is named at the
> call site.

The rule was written down and the list was not. This is the list, produced by walking the call
sites rather than the library: what a ported path actually reaches, and what is known about each.

## Reached through `jmath`

| function | call sites | status |
|---|---|---|
| `math::log` | 14 | **exact**. Correctly rounded in the reference, so rounding the true result suffices (decision 0006) |
| `math::sqrt` | 2 | **exact**. IEEE 754 mandates the rounding, so every implementation already agrees |
| `math::round`, `fast_math::round` | 4 | **exact**, oracle-backed |
| `strict_math::exp` | 1, in `NaturalLogUtils` | **exact** against `StrictMath.exp`; the call site's own reference is `Math.exp`, which is **1 ulp** away (decision 0025) |
| `gamma::*` | 1 | **exact**, oracle-backed (`gamma-erf-normal`) |
| `normal::*` | 1 | **exact**, oracle-backed |
| `percentile::*` | 1 | **exact**, oracle-backed |
| `saddle_point::hypergeometric_log_probability` | 1 | **exact**, oracle-backed |

## Reached through the host libm, bypassing `jmath`

These are the ones the rule exists for. Each is a place where a ported output depends on whatever
libm the machine ships.

| call site | function | status |
|---|---|---|
| `gatk-engine/src/math_utils.rs:79` — `pow10` | `powf` | **99.9378% against `Math.pow`**, unbounded per call. A bounded alternative now exists (`strict_math::pow`, 1 ulp, decision 0027) and was **measured to be worse here**: switching broke a passing byte-identity claim, because on the points these suites reach the libm is simply closer |
| `gatk-engine/src/fisher_exact.rs:109` — `pow10` | `powf` | as above |
| `gatk-annotation/src/allele_pseudo_depth.rs:339,341` — `calculateWeights` | `powf` ×2 | as above. Which of the two runs depends on `weightDecay`, so the exposure varies with an argument |
| `gatk-engine/src/natural_log_utils.rs:172,174` — `log1mexp` | `ln_1p`, `exp_m1` | **not exact**, and the module says so at the call site. Neither has a ported equivalent, and nothing in G1 reaches this branch |
| `gatk-engine/src/activity_profile.rs` — `normal_distribution` | `exp` | `java_math_exp`: **exact on aarch64**, where Java's `Math.exp` is fdlibm and the helper is `strict_math::exp` (0 of 2,000,000 inputs differ, 2026-11-23). **x86-64 still uses the host (glibc) `exp`**, ~0.24% of inputs 1 ulp from the HotSpot intrinsic, which was **measured to be the closer stand-in there**: the `activityprofile` suite compares 266 kernel values as raw bits and the host `exp` matches every one, where `strict_math::exp` moves 10 of them by an ulp. See below |
| `gatk-engine/src/natural_log_utils.rs` — `normalize_from_log_to_linear_space_host_exp` | `exp` | `java_math_exp`: **exact on aarch64** (fdlibm, as above); **x86-64 still uses the host `exp`**, for `LearnReadOrientationModel`'s E-step only. Built on `strict_math::exp`, two of the eight covering-array rows printed one prior an ulp from the oracle's; built on the host `exp`, every row of the array matches. The reference run with `-XX:-UseLibmIntrinsic` matches the fdlibm build exactly, which pins the difference on the `Math.exp` intrinsic. See below |
| `gatk-engine/src/mutect_engine.rs` — `posterior_probability_of_error` | `exp` | `java_math_exp` (exact on aarch64, host `exp` on x86-64), through `normalize_from_log_to_linear_space_host_exp`. Its callers are FilterMutectCalls' filters only (contamination, normal artifact, slippage, the clustering model's sequencing-error probability). Measured with the filters' answers as raw bits: see below |
| `gatk-engine/src/somatic_clustering_model.rs` — `performEMIteration`, `backgroundProbGivenSomatic` | `exp` | `java_math_exp` (exact on aarch64, host `exp` on x86-64), as above. The model's `log_sum_exp` calls keep fdlibm: no measured row depended on them |
| `gatk-engine/src/mann_whitney.rs:167` | `powi(3)` | **exact**: an integer power of a `f64` is repeated multiplication, and the reference computes it the same way |

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

Nothing in the code, and two things in what is claimed.

**The gap is nine call sites, not "the corpus".** Every one is named above with the number attached,
so a future byte claim over any of them knows what it is standing on.

**`BinomialDistribution` and SVD, the other two names in H.5, are not reached by anything ported.**
Both are Mutect2-family — read orientation, panel of normals, the somatic filters — which is
Milestone G3. Porting them now would be work ahead of its consumer, and the honest entry for them
is "waits for G3" rather than "remaining".
