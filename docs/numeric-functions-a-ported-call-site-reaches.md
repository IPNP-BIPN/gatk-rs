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
| `gatk-engine/src/activity_profile.rs:98` | `exp` | the **host** `exp`, not `strict_math::exp`, and now **measured to be the right one**: the `activityprofile` suite compares 266 kernel values as raw bits and the host `exp` matches every one, where `strict_math::exp` moves 10 of them by an ulp. See below |
| `gatk-engine/src/natural_log_utils.rs` — `normalize_from_log_to_linear_space_host_exp` | `exp` | the **host** `exp`, for `LearnReadOrientationModel`'s E-step only. Built on `strict_math::exp`, two of the eight covering-array rows printed one prior an ulp from the oracle's; built on the host `exp`, every row of the array matches. The reference run with `-XX:-UseLibmIntrinsic` matches the fdlibm build exactly, which pins the difference on the `Math.exp` intrinsic. See below |
| `gatk-engine/src/mann_whitney.rs:167` | `powi(3)` | **exact**: an integer power of a `f64` is repeated multiplication, and the reference computes it the same way |

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

## What this changes

Nothing in the code, and two things in what is claimed.

**The gap is seven call sites, not "the corpus".** Every one is named above with the number attached,
so a future byte claim over any of them knows what it is standing on.

**`BinomialDistribution` and SVD, the other two names in H.5, are not reached by anything ported.**
Both are Mutect2-family — read orientation, panel of normals, the somatic filters — which is
Milestone G3. Porting them now would be work ahead of its consumer, and the honest entry for them
is "waits for G3" rather than "remaining".
