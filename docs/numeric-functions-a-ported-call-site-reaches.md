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
| `gatk-engine/src/activity_profile.rs:98` | `exp` | the **host** `exp`, not `strict_math::exp`. 99.9711% against `Math.exp`. Worth revisiting when the activity profile gets its byte claim |
| `gatk-engine/src/mann_whitney.rs:167` | `powi(3)` | **exact**: an integer power of a `f64` is repeated multiplication, and the reference computes it the same way |

## What this changes

Nothing in the code, and two things in what is claimed.

**The gap is six call sites, not "the corpus".** Every one is named above with the number attached,
so a future byte claim over any of them knows what it is standing on.

**`BinomialDistribution` and SVD, the other two names in H.5, are not reached by anything ported.**
Both are Mutect2-family — read orientation, panel of normals, the somatic filters — which is
Milestone G3. Porting them now would be work ahead of its consumer, and the honest entry for them
is "waits for G3" rather than "remaining".
