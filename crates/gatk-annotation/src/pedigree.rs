//! `PossibleDeNovo` and `TransmittedSingleton`, and the `MendelianViolation` under them, ported
//! from GATK 4.7.0.0 (`PossibleDeNovo.java`, `TransmittedSingleton.java`,
//! `PedigreeAnnotation.java` and `MendelianViolation.java` are unchanged since 4.6.2.0).
//!
//! The two annotations that need a pedigree: which children carry a variant neither parent has, and
//! which parents passed a singleton on.
//!
//! # `TransmittedSingleton` tests the child's depth three times and calls it three things
//!
//! ```java
//! final boolean childIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
//! final boolean momIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
//! final boolean dadIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
//! ```
//!
//! All three read `getChildID()`. The parents' depths are never looked at, so a trio whose child is
//! deep enough passes the depth test however shallow the parents are, and the annotation's own
//! documented caveat ("high depth (>20)" for all three samples) is not what the code does. The
//! golden's `singleton shallow-parents` row (`pedigree.txt.gz`) separates the two readings: both
//! parents have a `DP` of 1, the child 30, and the row is `transmittedSingleton=[mom]`.
//!
//! # A missing `DP` is `-1`, not zero, and the comparison is against a threshold of zero
//!
//! `Genotype.getDP()` returns `-1` when the field is absent. `PossibleDeNovo`'s depth threshold
//! defaults to **zero**, so `-1 >= 0` is false and a trio with no `DP` at all fails the
//! high-confidence branch and falls through to the low-confidence one, which does not test depth.
//! The golden's `denovo no-depth` row is `loConfDeNovo=[kid]`.
//!
//! # `isViolation` compares **allele objects**, and a half-called parent takes the full test
//!
//! ```java
//! final Allele childRef = gChild.getAlleles().get(0);
//! return !(gMom.getAlleles().contains(childRef)
//!         && gDad.getAlleles().contains(gChild.getAlleles().get(1)) ||
//!     gMom.getAlleles().contains(gChild.getAlleles().get(1))
//!         && gDad.getAlleles().contains(childRef));
//! ```
//!
//! (Two lines in `MendelianViolation.java`, wrapped here at the `&&`s.) The final test indexes the
//! child's alleles at 0 and 1, so it is diploid-only, and it asks each parent to *contain* one of
//! them rather than to be able to transmit it. A parent that is `0/.` is `MIXED`, which
//! `isCalled()` counts as called, so the uncalled-parent branches before this test fire only for a
//! parent that is **entirely** uncalled, and a `0/.` parent reaches the containment test, where it
//! contains the reference. The golden's `violation mother-half-called` row (a half-called mother,
//! a hom-ref father and a het child) is `true`.
//!
//! # The de novo allele-frequency cutoff is a `Math.max` of a count and a fraction
//!
//! ```java
//! final double percentNumberOfSamplesCutoff = vc.getNSamples()*percentOfSamplesCutoff;
//! final double AFcutoff = Math.max(flatNumberOfSamplesCutoff, percentNumberOfSamplesCutoff);
//! ```
//!
//! Four samples, or a thousandth of the cohort, whichever is larger. So the flat four wins for
//! every cohort below four thousand samples, and the fraction only starts to matter above that.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<&Genotype>` is "a borrowed genotype, or nothing": `vc.genotype(id)` answers `None`
//!   where Java's `vc.getGenotype(id)` answers `null`. `.unwrap_or(-1)` supplies the value for the
//!   nothing case, which is how `getDP()` and `getGQ()` turn an absent field into `-1`.
//! - `.expect("...")` unwraps an `Option` and aborts the program with that message if it is empty.
//!   It is used only where an earlier check (`contextHasTrioGQs`) has already proven presence.
//! - `matches!(value, A | B)` is true when `value` is one of the listed variants.
//! - `let (Some(a), Some(b), Some(c)) = (x, y, z) else { ... };` unwraps three options at once or
//!   takes the early exit.
//! - `x as f64` converts an integer to a `double`; Java does the same implicitly when an `int` is
//!   compared with a `double`.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

/// `GATKVCFConstants.HI_CONF_DENOVO_KEY`.
pub const HI_CONF_DENOVO_KEY: &str = "hiConfDeNovo";
/// `GATKVCFConstants.LO_CONF_DENOVO_KEY`.
pub const LO_CONF_DENOVO_KEY: &str = "loConfDeNovo";
/// `GATKVCFConstants.TRANSMITTED_SINGLETON`.
pub const TRANSMITTED_SINGLETON_KEY: &str = "transmittedSingleton";
/// `GATKVCFConstants.NON_TRANSMITTED_SINGLETON`.
pub const NON_TRANSMITTED_SINGLETON_KEY: &str = "nonTransmittedSingleton";
/// `VCFConstants.ALLELE_COUNT_KEY`.
pub const ALLELE_COUNT_KEY: &str = "AC";

/// `PossibleDeNovo.hi_GQ_threshold`.
const HI_GQ_THRESHOLD: i32 = 20;
/// `PossibleDeNovo.lo_GQ_threshold`.
const LO_GQ_THRESHOLD: i32 = 10;
/// `PossibleDeNovo.percentOfSamplesCutoff`.
const PERCENT_OF_SAMPLES_CUTOFF: f64 = 0.001;
/// `PossibleDeNovo.flatNumberOfSamplesCutoff`, an `int` 4 in the reference. Its only use is
/// `Math.max(int, double)`, which widens it to `double`, so it is stored as one.
const FLAT_NUMBER_OF_SAMPLES_CUTOFF: f64 = 4.0;
/// `TransmittedSingleton`'s thresholds: `HI_GQ_THRESHOLD`, `HI_DP_THRESHOLD` and
/// `CALL_RATE_THRESHOLD`, prefixed here because `PossibleDeNovo` has a GQ threshold of its own.
const TS_HI_GQ_THRESHOLD: i32 = 20;
const TS_HI_DP_THRESHOLD: i32 = 20;
const CALL_RATE_THRESHOLD: f64 = 0.90;

/// `htsjdk.variant.variantcontext.GenotypeType`, which `MendelianViolation` uses as the three
/// axes of its inheritance table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeType {
    /// No alleles at all.
    Unavailable,
    /// Every allele is a no-call.
    NoCall,
    HomRef,
    Het,
    HomVar,
    /// Some alleles called and some not.
    Mixed,
}

/// `Genotype.determineType()` (htsjdk 5.0.0).
///
/// **How**: one pass that notes whether any allele is a no-call, remembers the first called
/// allele, and notes whether a later called allele differs from it; then the verdict in the
/// reference's order.
///
/// Note the ordering: a genotype with **any** no-call and at least one called allele is `MIXED`,
/// which `isCalled()` treats as called and `isHet()` does not.
pub fn genotype_type(genotype: &Genotype) -> GenotypeType {
    if genotype.alleles.is_empty() {
        return GenotypeType::Unavailable;
    }
    let mut saw_no_call = false;
    let mut saw_multiple = false;
    let mut first: Option<&Allele> = None;
    for allele in &genotype.alleles {
        if allele.is_no_call() {
            saw_no_call = true;
        } else if first.is_none() {
            first = Some(allele);
        // `!allele.equals(firstCallAllele)`: `==` and `!=` on alleles are value equality, Java's
        // `equals`, not reference identity.
        } else if Some(allele) != first {
            saw_multiple = true;
        }
    }
    if saw_no_call {
        return match first {
            None => GenotypeType::NoCall,
            Some(_) => GenotypeType::Mixed,
        };
    }
    let Some(first) = first else {
        // The reference throws `IllegalStateException` here ("BUG: there are no alleles present
        // in this genotype but the alleles list is not null"). Unreachable: a non-empty list with
        // no no-call has a called allele.
        return GenotypeType::Unavailable;
    };
    if saw_multiple {
        GenotypeType::Het
    } else if first.is_reference() {
        GenotypeType::HomRef
    } else {
        GenotypeType::HomVar
    }
}

/// `Genotype.isCalled()`: everything but `NO_CALL` and `UNAVAILABLE`, so `MIXED` counts as called.
pub fn is_called(genotype: &Genotype) -> bool {
    !matches!(
        genotype_type(genotype),
        GenotypeType::NoCall | GenotypeType::Unavailable
    )
}

/// `Genotype.isHomRef()`.
pub fn is_hom_ref(genotype: &Genotype) -> bool {
    genotype_type(genotype) == GenotypeType::HomRef
}

/// `Genotype.isHet()`. A `MIXED` genotype such as `0/.` is not het.
pub fn is_het(genotype: &Genotype) -> bool {
    genotype_type(genotype) == GenotypeType::Het
}

/// `Genotype.isHomVar()`.
pub fn is_hom_var(genotype: &Genotype) -> bool {
    genotype_type(genotype) == GenotypeType::HomVar
}

/// `Genotype.isNoCall()`: every allele a no-call. An empty genotype is `UNAVAILABLE`, not this.
pub fn is_no_call(genotype: &Genotype) -> bool {
    genotype_type(genotype) == GenotypeType::NoCall
}

/// `Genotype.getDP()`, which is **-1** when the field is absent rather than zero.
pub fn depth(genotype: &Genotype) -> i32 {
    genotype.dp.unwrap_or(-1)
}

/// `Genotype.getGQ()`, likewise -1 when absent.
pub fn genotype_quality(genotype: &Genotype) -> i32 {
    genotype.gq.unwrap_or(-1)
}

/// `Trio`, reduced to the three sample identifiers the annotations use.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trio {
    pub family_id: String,
    pub maternal_id: String,
    pub paternal_id: String,
    pub child_id: String,
}

/// `PedigreeAnnotation.contextHasTrioGQs`: all three present **and** carrying a `GQ`.
///
/// ```java
/// return   (!mom.isEmpty() && vc.hasGenotype(mom) && vc.getGenotype(mom).hasGQ())
/// ```
///
/// and the same for `dad` and `kid`. `.all(..)` stops at the first member that fails, as `&&`
/// does; `.map(..).unwrap_or(false)` turns a missing genotype into a failure.
pub fn context_has_trio_gqs(vc: &VariantContext, trio: &Trio) -> bool {
    [&trio.maternal_id, &trio.paternal_id, &trio.child_id]
        .iter()
        .all(|id| {
            !id.is_empty()
                && vc
                    .genotype(id)
                    .map(|genotype| genotype.gq.is_some())
                    .unwrap_or(false)
        })
}

/// `MendelianViolation.isViolation(gMom, gDad, gChild)`, the static three-genotype form.
///
/// **How**: the reference's five steps in order, because the order decides which fires first: a
/// no-call child is never a violation; an all-hom-ref trio is not; a wholly uncalled mother, then
/// father, is judged from the other parent's homozygosity alone; otherwise the containment test
/// quoted in the module note.
///
/// **What is not reproduced**: a child with fewer than two alleles. The reference indexes
/// `getAlleles().get(1)` and throws `IndexOutOfBoundsException`; here it is "not a violation".
pub fn is_violation(mother: &Genotype, father: &Genotype, child: &Genotype) -> bool {
    if is_no_call(child) {
        return false;
    }
    if is_hom_ref(mother) && is_hom_ref(father) && is_hom_ref(child) {
        return false;
    }
    if !is_called(mother) {
        return (is_hom_ref(father) && is_hom_var(child))
            || (is_hom_var(father) && is_hom_ref(child));
    }
    if !is_called(father) {
        return (is_hom_ref(mother) && is_hom_var(child))
            || (is_hom_var(mother) && is_hom_ref(child));
    }
    // Diploid-only: indexes the child's alleles at 0 and 1. `.first()` and `.get(1)` answer `None`
    // past the end where Java's `get` throws.
    let Some(child_first) = child.alleles.first() else {
        return false;
    };
    let Some(child_second) = child.alleles.get(1) else {
        return false;
    };
    !(mother.alleles.contains(child_first) && father.alleles.contains(child_second)
        || mother.alleles.contains(child_second) && father.alleles.contains(child_first))
}

/// `MendelianViolation.isViolation(mother, father, child, vc)` and `getParentsRefRefChildHet`
/// together, which is how `PossibleDeNovo` uses them:
///
/// ```java
/// mendelianViolation.isViolation(trio.getMother(), trio.getFather(), trio.getChild(), vc) &&
/// mendelianViolation.getParentsRefRefChildHet() > 0) {
/// ```
///
/// The instance form recomputes both from scratch on every call, and `PossibleDeNovo` requires
/// **both** that there was a violation and that the ref/ref/het cell of the inheritance table is
/// non-zero. Since the table is reset per call and holds one family, the second condition says the
/// violation was specifically a ref/ref parent pair with a het child, and it is the stricter of the
/// two.
///
/// **The `complete_trios_only` argument**: every `PossibleDeNovo` constructor builds
/// `new MendelianViolation(minGenotypeQualityP)`, which is `this(minGenotypeQualityP, true,
/// false)`, so the reference runs with `completeTriosOnly` **false**, while [`possible_de_novo`]
/// passes `true`. The answer is the same either way: the only extra trios `true` rejects are
/// those with one parent uncalled, and such a trio can never satisfy the final ref/ref/het test,
/// which needs both parents `HOM_REF`.
pub fn is_de_novo_violation(
    vc: &VariantContext,
    trio: &Trio,
    min_genotype_quality: f64,
    complete_trios_only: bool,
) -> bool {
    let (Some(mother), Some(father), Some(child)) = (
        vc.genotype(&trio.maternal_id),
        vc.genotype(&trio.paternal_id),
        vc.genotype(&trio.child_id),
    ) else {
        // `abortOnSampleNotFound` is true for the constructor `PossibleDeNovo` uses, so this is an
        // `IllegalArgumentException`. It is unreachable here because `contextHasTrioGQs` ran first.
        return false;
    };
    // `updateViolations`' `noCall++` branches: either one ends the call with no violation.
    if complete_trios_only && (!is_called(mother) || !is_called(father) || !is_called(child)) {
        return false;
    }
    if (!is_called(mother) && !is_called(father)) || !is_called(child) {
        return false;
    }
    // The `lowQual++` branch. With a minimum of zero, which is the only value `PossibleDeNovo`
    // uses, it never fires, so an absent `GQ` (read as -1) is not rejected here.
    if min_genotype_quality > 0.0
        && ((genotype_quality(mother) as f64) < min_genotype_quality
            || (genotype_quality(father) as f64) < min_genotype_quality
            || (genotype_quality(child) as f64) < min_genotype_quality)
    {
        return false;
    }
    if !is_violation(mother, father, child) {
        return false;
    }
    // `getParentsRefRefChildHet() > 0`: the one cell of the inheritance table this call filled.
    is_hom_ref(mother) && is_hom_ref(father) && is_het(child)
}

/// `VariantContext.getCalledChrCount(allele)`: how many copies of one allele the genotypes carry.
///
/// The reference counts only unfiltered genotypes:
///
/// ```java
/// if (!g.isFiltered()) {
///     n += g.countAllele(a);
/// }
/// ```
///
/// The port counts every genotype, so a genotype whose `isFiltered()` is true and that carries
/// the allele is counted here and not by the reference.
pub fn called_chromosome_count_for(vc: &VariantContext, allele: &Allele) -> usize {
    vc.genotypes
        .iter()
        .map(|genotype| genotype.alleles.iter().filter(|a| *a == allele).count())
        .sum()
}

/// `PossibleDeNovo.annotate`: the high- and low-confidence child lists.
///
/// **What**: for each trio at a biallelic site whose three members have a `GQ` and whose child is
/// a het with two hom-ref parents, the child goes in the high-confidence list if every threshold
/// below holds, else in the low-confidence list if the weaker ones do. Both lists are then
/// dropped when the alternate allele is common (the cutoff in the module note).
///
/// `parent_gq_threshold` and `depth_threshold` are the `--denovo-parent-gq-threshold` and
/// `--denovo-depth-threshold` arguments, defaulting to 20 and 0 in the reference. The returned
/// pairs are in the order of the reference's `LinkedHashMap` puts.
pub fn possible_de_novo(
    vc: &VariantContext,
    trios: &[Trio],
    parent_gq_threshold: i32,
    depth_threshold: i32,
) -> Vec<(String, Vec<String>)> {
    if trios.is_empty() {
        return Vec::new();
    }
    let mut high = Vec::new();
    let mut low = Vec::new();
    for trio in trios {
        // `vc.isBiallelic()` is `getNAlleles() == 2`. The reference's `&&` chain is split into
        // two `continue`s; the order of the tests is kept.
        if !(vc.alleles.len() == 2 && context_has_trio_gqs(vc, trio)) {
            continue;
        }
        if !is_de_novo_violation(vc, trio, 0.0, true) {
            continue;
        }
        // `contextHasTrioGQs` just proved all three present, so these cannot fail.
        let child = vc.genotype(&trio.child_id).expect("a child");
        let mother = vc.genotype(&trio.maternal_id).expect("a mother");
        let father = vc.genotype(&trio.paternal_id).expect("a father");
        let (child_gq, mom_gq, dad_gq) = (
            genotype_quality(child),
            genotype_quality(mother),
            genotype_quality(father),
        );
        // A tuple on each side assigns three variables at once.
        let (child_dp, mom_dp, dad_dp) = (depth(child), depth(mother), depth(father));

        if child_gq >= HI_GQ_THRESHOLD
            && mom_gq >= parent_gq_threshold
            && dad_gq >= parent_gq_threshold
            && child_dp >= depth_threshold
            && mom_dp >= depth_threshold
            && dad_dp >= depth_threshold
        {
            high.push(trio.child_id.clone());
        // Strictly positive parent GQs: a GQ of zero, or an absent one read as -1, fails.
        } else if child_gq >= LO_GQ_THRESHOLD && mom_gq > 0 && dad_gq > 0 {
            low.push(trio.child_id.clone());
        }
    }

    // `vc.getNSamples()` is the number of genotypes. `java_max` is `Math.max(double, double)`.
    let percent_cutoff = vc.genotypes.len() as f64 * PERCENT_OF_SAMPLES_CUTOFF;
    let af_cutoff =
        gatk_engine::math_utils::java_max(FLAT_NUMBER_OF_SAMPLES_CUTOFF, percent_cutoff);
    // "we assume we're biallelic above so use the first alt": read outside the loop, so a site
    // with no alternate at all is an IndexOutOfBounds here even though the loop skipped it. The
    // port answers an empty result instead.
    let Some(first_alternate) = vc.alternate_alleles().first() else {
        return Vec::new();
    };
    let de_novo_allele_count = called_chromosome_count_for(vc, first_alternate) as f64;

    let mut out = Vec::new();
    if !high.is_empty() && de_novo_allele_count < af_cutoff {
        out.push((HI_CONF_DENOVO_KEY.to_string(), high));
    }
    if !low.is_empty() && de_novo_allele_count < af_cutoff {
        out.push((LO_CONF_DENOVO_KEY.to_string(), low));
    }
    out
}

/// `TransmittedSingleton.annotate`: the parents that transmitted, and those that did not.
///
/// `allele_count` is the site's `AC` attribute, defaulting to zero, which is what gates the two
/// branches: exactly two copies for a transmitted singleton and exactly one for a non-transmitted
/// one. It is the caller's reading of `vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0)`.
///
/// **The call-rate gate** counts genotypes whose `GQ` is above 20, not called genotypes, though
/// the reference's documentation says "Call Rate". A genotype with no `GQ` reads as -1 and so
/// counts as a failed call:
///
/// ```java
/// long highQualCalls = vc.getGenotypes().stream()
///         .filter(gt -> gt.getGQ() > HI_GQ_THRESHOLD).count();
/// ```
///
/// (One line in `TransmittedSingleton.java`, wrapped here.)
pub fn transmitted_singleton(
    vc: &VariantContext,
    trios: &[Trio],
    allele_count: i32,
) -> Vec<(String, Vec<String>)> {
    if vc.alleles.len() != 2 || trios.is_empty() {
        return Vec::new();
    }
    // Strictly greater than the threshold here, where the per-trio GQ tests below are
    // greater-or-equal.
    let high_quality_calls = vc
        .genotypes
        .iter()
        .filter(|genotype| genotype_quality(genotype) > TS_HI_GQ_THRESHOLD)
        .count();
    // `(double) highQualCalls / vc.getNSamples()`: a `double` division.
    let call_rate = high_quality_calls as f64 / vc.genotypes.len() as f64;
    if call_rate < CALL_RATE_THRESHOLD {
        return Vec::new();
    }
    let mut transmitted = Vec::new();
    let mut non_transmitted = Vec::new();
    for trio in trios {
        if !context_has_trio_gqs(vc, trio) {
            continue;
        }
        let mother = vc.genotype(&trio.maternal_id).expect("a mother");
        let father = vc.genotype(&trio.paternal_id).expect("a father");
        let child = vc.genotype(&trio.child_id).expect("a child");

        let one_parent_has_allele =
            (is_het(mother) && is_hom_ref(father)) || (is_hom_ref(mother) && is_het(father));
        // The het parent's identifier, or the father's whenever the mother is not het. The
        // non-transmitted branch does not test `one_parent_has_allele`, so there a trio whose
        // mother is not het reports the father whatever his genotype.
        let matching_parent = if is_het(mother) {
            &trio.maternal_id
        } else {
            &trio.paternal_id
        };

        let mom_high_gq = genotype_quality(mother) >= TS_HI_GQ_THRESHOLD;
        let dad_high_gq = genotype_quality(father) >= TS_HI_GQ_THRESHOLD;
        let child_high_gq_het = is_het(child) && genotype_quality(child) >= TS_HI_GQ_THRESHOLD;
        let child_high_gq_hom_ref =
            is_hom_ref(child) && genotype_quality(child) >= TS_HI_GQ_THRESHOLD;

        // All three read the **child's** depth in the reference. See the module note.
        let child_deep = depth(child) >= TS_HI_DP_THRESHOLD;
        let mom_deep = depth(child) >= TS_HI_DP_THRESHOLD;
        let dad_deep = depth(child) >= TS_HI_DP_THRESHOLD;

        // Two independent `if`s, as in the reference, not an `if`/`else`: a trio is tested for
        // both lists.
        if child_deep
            && mom_deep
            && dad_deep
            && allele_count == 2
            && child_high_gq_het
            && one_parent_has_allele
            && mom_high_gq
            && dad_high_gq
        {
            transmitted.push(matching_parent.clone());
        }
        if child_deep
            && mom_deep
            && dad_deep
            && allele_count == 1
            && child_high_gq_hom_ref
            && mom_high_gq
            && dad_high_gq
        {
            non_transmitted.push(matching_parent.clone());
        }
    }
    let mut out = Vec::new();
    if !transmitted.is_empty() {
        out.push((TRANSMITTED_SINGLETON_KEY.to_string(), transmitted));
    }
    if !non_transmitted.is_empty() {
        out.push((NON_TRANSMITTED_SINGLETON_KEY.to_string(), non_transmitted));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genotype(name: &str, alleles: Vec<Allele>, gq: Option<i32>, dp: Option<i32>) -> Genotype {
        let mut g = Genotype::new(name, alleles);
        g.gq = gq;
        g.dp = dp;
        g
    }

    fn reference() -> Allele {
        Allele::from_str("A", true).expect("an allele")
    }

    fn alternate() -> Allele {
        Allele::from_str("C", false).expect("an allele")
    }

    #[test]
    fn a_half_called_genotype_is_mixed_and_counts_as_called() {
        let g = genotype(
            "s",
            vec![reference(), Allele::no_call()],
            Some(30),
            Some(20),
        );
        assert_eq!(genotype_type(&g), GenotypeType::Mixed);
        assert!(is_called(&g));
        assert!(!is_het(&g));
    }

    #[test]
    fn two_hom_ref_parents_with_a_het_child_are_a_violation() {
        let mother = genotype("m", vec![reference(), reference()], Some(50), Some(30));
        let father = genotype("f", vec![reference(), reference()], Some(50), Some(30));
        let child = genotype("c", vec![reference(), alternate()], Some(50), Some(30));
        assert!(is_violation(&mother, &father, &child));
    }

    #[test]
    fn a_missing_depth_is_minus_one() {
        let g = genotype("s", vec![reference(), reference()], Some(50), None);
        assert_eq!(depth(&g), -1);
    }
}
