# Costs already on the record

Milestone S, S.4 (#111). Places where the port pays for something, collected before the baseline
so the baseline had something to check itself against, and then checked. None of them is a bug.
Each entry says what the cost is, why it was accepted, whether it can be removed, and what the
baseline (S.3, #110; [baseline.md](baseline.md)) says it is worth.

The baseline's verdict, in one line: one of these costs is most of the port's time where the port
is slow, and none of the others shows at all on the rows the claim covers.

| Cost | Removable? | Measured |
|---|---|---|
| `jmath` transcendentals in software | the **result** no, the **algorithm** yes | **confirmed**, the dominant cost |
| fixed reduction order | no | demoted: absent from every profile |
| the exact decimal expansion | no | demoted: no tool near the slow end |
| cloning where the reference aliases (`AllelePseudoDepth`, `normalize_log`) | yes | demoted: no measurable cost |
| gzip-decoding goldens per test | yes, test-only | demoted: 11 ms in an 87 s run |

## `jmath` transcendentals in software

**The cost.** The reference's `Math.log` is a hardware-assisted intrinsic; `jmath::math::log`
returns the same correctly rounded value by summing a 22-term atanh series in double-double, 22
double-double divisions per call (`crates/jmath/src/log.rs` in htsjdk-rs). `strict_math::exp` and
`pow` are FDLIBM in software where HotSpot has intrinsics.

**Why it was accepted.** It is the claim: decision 0006 of htsjdk-rs measured `Math.log` to be
correctly rounded, so the port must produce the correctly rounded value, and the host libm does
not promise it (S.2's ratchet refuses new host transcendentals for that reason).

**Removable?** The result is not negotiable; the way it is computed is. A correctly rounded
function has exactly one answer, so any algorithm that provably returns it cannot move a byte.

**Measured: confirmed.** A native profile of the baseline row of `ModelSegments`, the slowest tool
(18.4 times a warm JVM), puts 94% of samples in that series (`jmath::log::ln_dd`); the second,
`VariantRecalibrator` (4.6), is dominated by the same function. It is S.5's first target (#112).

## Fixed reduction order

**The cost.** `MathUtils.sum`, `sumArrayFunction` and their ports accumulate in index order, one
term at a time, so no reduction is vectorised or reassociated.

**Removable?** No. It is the claim, and S.2 found a corpus that sees it: reversing the reduction
over reads moved `somatic-likelihoods` by two ulp. It is on the slow-on-purpose list.

**Measured: demoted.** No profile taken for the baseline shows a reduction loop; the two outliers
are spent in the logarithm. On the rows the claim covers, a few hundred terms summed in order
costs nothing anyone could measure.

## The exact decimal expansion

**The cost.** `gatk-annotation/src/decimal_format.rs` expands a double to its exact decimal value,
up to 1,074 digits, so `DecimalFormat`'s rounding is decided on the value and not on a shortest
representation. The naive version took a conformance run from 0.12 s to 6 s; it now runs only
where `may_be_a_short_decimal` says the value could be a short decimal that ties, a pessimistic
bound on how often, not on how much.

**Removable?** No: it is decision 0017 of htsjdk-rs, and on the slow-on-purpose list.

**Measured: demoted.** The tools that write formatted annotations sit at the table's median
(`VariantAnnotator` 0.058 steady, `VariantFiltration` 0.061). The gate bounds how often the
expansion runs, and the baseline says that is often enough to be invisible.

## Cloning where the reference aliases

**The cost.** `AllelePseudoDepth::compose_input_likelihood_matrix` clones each emitted allele's row
per call where the reference wraps the matrix in a `SubsettedLikelihoodMatrix`;
`natural_log_utils::normalize_log` returns a new vector where the reference can normalise in
place.

**Why it was accepted.** Rust will not hand out the aliasing the reference relies on, and the
reference's aliasing is itself the source of the memo bug the `allele-pseudo-depth` suite
reproduces.

**Removable?** Yes, in principle: a borrowed view would do for the natural-log branch, and an
in-place normalisation for a caller that owns its array.

**Measured: demoted.** Peak RSS is 3 to 21% of the reference's on every tool, the largest being
`BaseRecalibrator` (58 MiB against 272 MiB), so the copies do not show in memory. In time, the
tools on these paths are below one: `LearnReadOrientationModel`, the slowest of them at 0.64
steady, runs in 38 ms natively, too short for a profile to attribute anything to an allocation.
`AllelePseudoDepth` is not reachable from any measured tool's row, so it has no number of its own;
its suite measures what it produces, and nothing measures what it costs. Not worth removing
before a larger input says otherwise.

## Gzip-decoding goldens per test

**The cost.** Every golden-backed test decompresses its golden.

**Removable?** Yes, and test-only: nothing a user runs pays it.

**Measured: demoted.** The workspace holds 306 goldens, 1.0 MB compressed and 9.2 MB decoded;
decoding all of them takes 11 ms, against 87 s for `cargo test --release` on the same machine.
