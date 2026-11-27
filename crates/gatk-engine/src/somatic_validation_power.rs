//! `PowerCalculationUtils`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup` (GATK 4.7.0.0;
//! `PowerCalculationUtils.java` is unchanged since 4.6.2.0).
//!
//! How likely a validation pileup of a given depth is to show a variant the discovery pileup
//! called, and how many reads count as showing it.
//!
//! # The power is one minus a beta-binomial tail
//!
//! The distribution's shapes are the discovery counts PLUS ONE, `alt + 1` and
//! `total - alt + 1`, so for non-negative counts they are at least one and the constructor's
//! `alpha` and `beta` refusals cannot fire. The tail is
//! `1 - cumulativeProbability(minCount - 1)`, and that cumulative probability is a plain
//! uncompensated loop over `probability(i)`, each of which takes the exponential of
//! `logProbability(i)`:
//!
//! ```java
//! double result = 0.0;
//! for (int i = 0; i <= k; i ++) {
//!     result += probability(i);
//! }
//! return result;
//! ```
//!
//! # The minimum count is a quantile floored at two
//!
//! `max(inverseCumulativeProbability(0.99), 2)` over a binomial of the validation depth at the
//! noise ratio. The floor is what makes a pileup of zero reads and a pileup of 317 clean ones both
//! answer two.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, PowerError>` is "a `T` or one of the refusals below": `Ok(..)` is the return and
//!   `Err(..)` the throw. `.map_err(PowerError::Distribution)` wraps a lower layer's error in this
//!   file's type, and the `?` after it propagates that error to the caller at once, as an
//!   uncaught exception would.
//! - `&ReadPileup<'a>` is a borrowed, read-only pileup. The `'a` is a lifetime: it records that
//!   the pileup's elements borrow reads owned elsewhere, and lets the elements copied out of it
//!   keep borrowing them. It has no Java counterpart and does not change what is computed.
//! - `.iter().filter(..).count()` is Java's `stream().filter(..).count()`.
//! - `i32` is Java's `int`, `i64` its `long`, `f64` its `double`; `x as f64` widens exactly for
//!   the counts here.

use crate::beta_binomial::{BetaBinomialDistribution, BetaBinomialError};
use crate::pileup::PileupElement;
use crate::read_pileup::ReadPileup;
use crate::variant_context_utils::{
    choose_allele_for_read, does_read_contain_allele, Allele, PileupAlleleError, Trilean,
};

/// `PowerCalculationUtils.P_VALUE_FOR_NOISE`.
pub const P_VALUE_FOR_NOISE: f64 = 0.99;

/// `PowerCalculationUtils.MINIMUM_NUM_READS_FOR_SIGNAL_COUNT`.
pub const MINIMUM_NUM_READS_FOR_SIGNAL_COUNT: i32 = 2;

/// What the power calculation refuses. Each variant is one throw site in the reference, kept as
/// a value so the caller can report the class and message the reference would have printed.
#[derive(Debug, Clone, PartialEq)]
pub enum PowerError {
    /// `ParamUtils.isPositiveOrZero` on the validation depth.
    NegativeTotalCount,
    /// `ParamUtils.inRange` on the noise ratio.
    RatioOutOfRange,
    /// The distribution refused.
    Distribution(BetaBinomialError),
    /// The quantile refused.
    Quantile(jmath::binomial::BinomialError),
    /// `ParamUtils.isPositiveOrZero` on the base-quality cutoff, or the allele choice underneath.
    Allele(PileupAlleleError),
}

impl PowerError {
    /// The exception class the reference throws. Every `ParamUtils` check throws
    /// `IllegalArgumentException`; only the allele layer has others of its own.
    pub fn java_class(&self) -> &'static str {
        match self {
            PowerError::Allele(error) => error.java_class(),
            _ => "java.lang.IllegalArgumentException",
        }
    }

    /// The message, including the reference's own doubled word in the second:
    /// `"Cannot have have a ratio that is outside of 0.0 - 1.0."` is copied from
    /// `PowerCalculationUtils.java` line 50 as it stands.
    ///
    /// The two wrapped layers print their `Debug` form, which is this port's own text and not a
    /// Java message.
    pub fn message(&self) -> String {
        match self {
            PowerError::NegativeTotalCount => "Cannot have a negative total count.".to_string(),
            PowerError::RatioOutOfRange => {
                "Cannot have have a ratio that is outside of 0.0 - 1.0.".to_string()
            }
            PowerError::Distribution(error) => format!("{error:?}"),
            PowerError::Quantile(error) => format!("{error:?}"),
            PowerError::Allele(error) => error.message(),
        }
    }
}

/// `calculatePower(validationTumorTotalCount, discoveryTumorAltCount, discoveryTumorTotalCount,
/// minCountForSignal)`: one minus the beta-binomial probability of fewer than
/// `minimum_count_for_signal` alternate reads.
///
/// ```java
/// final BetaBinomialDistribution betaBinomialDistribution =
///         new BetaBinomialDistribution(null, discoveryTumorAltCount + 1,
///             discoveryTumorTotalCount - discoveryTumorAltCount + 1, validationTumorTotalCount);
/// return (1 - betaBinomialDistribution.cumulativeProbability(minCountForSignal - 1));
/// ```
///
/// (`PowerCalculationUtils.java` lines 36 to 38, the constructor call wrapped at a space.)
///
/// The reference adds one in `int` and widens the sum to `double` at the call; this widens first
/// and adds `1.0`. For any count a pileup can hold the two give the same `double`, because both
/// the count and the count plus one are exact in a `double`. The `null` generator has no
/// counterpart: nothing here samples.
pub fn calculate_power(
    validation_total_count: i32,
    discovery_alt_count: i32,
    discovery_total_count: i32,
    minimum_count_for_signal: i32,
) -> Result<f64, PowerError> {
    let distribution = BetaBinomialDistribution::new(
        f64::from(discovery_alt_count) + 1.0,
        f64::from(discovery_total_count - discovery_alt_count) + 1.0,
        validation_total_count,
    )
    .map_err(PowerError::Distribution)?;
    let cumulative = distribution
        .cumulative_probability(minimum_count_for_signal - 1)
        .map_err(PowerError::Distribution)?;
    Ok(1.0 - cumulative)
}

/// `calculateMinCountForSignal`: the binomial quantile at 0.99, floored at two.
///
/// ```java
/// final BinomialDistribution binomialDistribution =
///     new BinomialDistribution(validationTumorTotalCount, maxSignalRatioInNormal);
/// return Math.max(binomialDistribution.inverseCumulativeProbability(P_VALUE_FOR_NOISE),
///     MINIMUM_NUM_READS_FOR_SIGNAL_COUNT);
/// ```
///
/// (`PowerCalculationUtils.java` lines 52 and 53, each wrapped at a space.)
///
/// The distribution is commons-math3's, reached through `jmath::binomial`. The two `ParamUtils`
/// checks before it run in the reference's order, so a negative depth is reported even when the
/// ratio is also out of range.
pub fn calculate_min_count_for_signal(
    validation_total_count: i32,
    max_signal_ratio_in_normal: f64,
) -> Result<i32, PowerError> {
    // `ParamUtils.isPositiveOrZero(validationTumorTotalCount, ...)`, the first of the two checks.
    if validation_total_count < 0 {
        return Err(PowerError::NegativeTotalCount);
    }
    // `ParamUtils.inRange`, inclusive at both ends. `(0.0..=1.0)` is the closed interval and
    // `.contains` is false for NaN, so a NaN ratio is refused here too.
    if !(0.0..=1.0).contains(&max_signal_ratio_in_normal) {
        return Err(PowerError::RatioOutOfRange);
    }
    let quantile = jmath::binomial::inverse_cumulative_probability(
        validation_total_count,
        max_signal_ratio_in_normal,
        P_VALUE_FOR_NOISE,
    )
    .map_err(PowerError::Quantile)?;
    // `Math.max` on two `int`s.
    Ok(quantile.max(MINIMUM_NUM_READS_FOR_SIGNAL_COUNT))
}

/// `retrievePileupElements`: not a deletion, and at or above the cutoff.
///
/// ```java
/// return Utils.stream(readPileup.iterator())
///             .filter(pe -> !pe.isDeletion())
///             .filter(pe -> pe.getQual() >= minBaseQualityCutoff)
///             .collect(Collectors.toList());
/// ```
///
/// The quality here is the pileup element's own base quality, which for a deletion would be the
/// constant 16 (`return isDeletion() ? DELETION_QUAL : read.getBaseQuality(offset);`, with
/// `DELETION_QUAL = (byte) 16` in `PileupElement.java`), so the deletion filter runs first and
/// the constant never decides anything.
///
/// Java compares a signed `byte` widened to `int`; `qual()` is an unsigned byte here, and the
/// two readings differ only for a stored quality above 127.
///
/// `.cloned()` copies each element out of the pileup, the way the reference collects references
/// to the same elements into a new list. The copies still borrow the same reads.
fn elements_passing_quality<'a>(
    pileup: &ReadPileup<'a>,
    minimum_base_quality: i32,
) -> Vec<PileupElement<'a>> {
    pileup
        .elements
        .iter()
        .filter(|element| !element.is_deletion())
        .filter(|element| i32::from(element.qual()) >= minimum_base_quality)
        .cloned()
        .collect()
}

/// `calculateMaxAltRatio`: the fraction of the pileup that is not the reference allele.
///
/// ```java
/// return numReference + numAlternate == 0 ? 0.0 :
///     (double) numAlternate / ((double) numReference + (double) numAlternate);
/// ```
///
/// (`PowerCalculationUtils.java` line 84, wrapped at a space.)
///
/// THE TWO FILTERS ARE NOT COMPLEMENTS. An element is alternate when it does not contain the
/// reference allele OR it precedes an indel; it is reference when it does contain the reference
/// allele AND precedes neither. An element that answers UNKNOWN, which is a read ending inside the
/// allele, is in neither count, so a pileup of nothing but such reads has a denominator of zero and
/// the ratio is the literal 0.0 rather than a NaN. The caller depends on that: a NaN ratio is
/// what makes `calculateBasicValidationResult` return `null` (no row), and the zero guard is why
/// that branch is never taken from here.
pub fn calculate_max_alt_ratio(
    pileup: &ReadPileup<'_>,
    reference: &Allele,
    minimum_base_quality: i32,
) -> Result<f64, PowerError> {
    // `ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Cannot have a negative
    // minBaseQualityCutoff.")`, which comes first, before the pileup is looked at.
    if minimum_base_quality < 0 {
        return Err(PowerError::Allele(
            PileupAlleleError::NegativeMinimumBaseQualityRatio,
        ));
    }
    let passing = elements_passing_quality(pileup, minimum_base_quality);
    // The closures `|element| ...` are Java's lambdas `pe -> ...`.
    let alternate = passing
        .iter()
        .filter(|element| {
            does_read_contain_allele(element, reference) == Trilean::False
                || element.is_before_deletion_start()
                || element.is_before_insertion()
        })
        .count();
    let reference_count = passing
        .iter()
        .filter(|element| {
            does_read_contain_allele(element, reference) == Trilean::True
                && !element.is_before_deletion_start()
                && !element.is_before_insertion()
        })
        .count();
    if reference_count + alternate == 0 {
        return Ok(0.0);
    }
    // Each count widened to a `double` separately and then added, as the reference's three
    // `(double)` casts do. `as f64` from a count is exact at any pileup size.
    Ok(alternate as f64 / (reference_count as f64 + alternate as f64))
}

/// `calculateNumReadsSupportingAllele`: the elements whose chosen allele is this alternate.
///
/// The choice is made against a list of exactly one alternate, so an element supporting a different
/// alternate is not counted here even when the pileup was built for several:
///
/// ```java
/// .filter(pe -> altAllele.equals(GATKVariantContextUtils.chooseAlleleForRead(pe,
///     referenceAllele, Collections.singletonList(altAllele), minBaseQualityCutoff)))
/// ```
///
/// (`PowerCalculationUtils.java` line 114, wrapped at a space.)
///
/// `altAllele.equals(null)` is false, which is the `None` branch here. `Allele.equals` compares
/// the reference flag as well as the bases, and so does the `==` on [`Allele`].
pub fn calculate_num_reads_supporting_allele(
    pileup: &ReadPileup<'_>,
    reference: &Allele,
    alternate: &Allele,
    minimum_base_quality: i32,
) -> Result<i64, PowerError> {
    if minimum_base_quality < 0 {
        return Err(PowerError::Allele(
            PileupAlleleError::NegativeMinimumBaseQualityRatio,
        ));
    }
    // `Collections.singletonList(altAllele)`: a list of exactly one alternate.
    let alternates = [alternate.clone()];
    let mut count = 0;
    for element in elements_passing_quality(pileup, minimum_base_quality) {
        let chosen = choose_allele_for_read(&element, reference, &alternates, minimum_base_quality)
            .map_err(PowerError::Allele)?;
        // `chosen.as_ref()` views the `Option<Allele>` as an `Option<&Allele>`, so it can be
        // compared with `Some(alternate)` without copying the allele.
        if chosen.as_ref() == Some(alternate) {
            count += 1;
        }
    }
    Ok(count)
}
