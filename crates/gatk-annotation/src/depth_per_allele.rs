//! `DepthPerAlleleBySample`, `AlleleFraction` and `DepthPerSampleHC`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.7.0.0; all three files are
//! unchanged since 4.6.2.0).
//!
//! `AD`, `AF` and the HaplotypeCaller's `DP`: three genotype fields, all counted off a
//! **marginalised** likelihood matrix.
//!
//! # The allele order the marginalisation uses is a `HashMap`'s
//!
//! ```java
//! final Map<Allele, List<Allele>> alleleSubset =
//!         alleles.stream().collect(Collectors.toMap(a -> a, Arrays::asList));
//! final AlleleLikelihoods<EVIDENCE, Allele> subsettedLikelihoods =
//!         likelihoods.marginalize(alleleSubset);
//! ```
//!
//! `Collectors.toMap` builds a `HashMap`, and `marginalize` takes its **key set** as the new
//! allele array. So the new matrix's allele order is `HashMap` iteration order over
//! `Allele.hashCode`, which is `Arrays.hashCode(bases) * 31 + Boolean.hashCode(isRef)`.
//!
//! `searchBestAllele` breaks a tie by keeping the first index, so a read that supports two
//! alleles equally gets whichever the hash put first as its best allele. That does **not** reach
//! `AD`: an exact tie has a confidence of 0, `isInformative` is false, and the count skips the read
//! rather than crediting the first allele. In `tests/data/depth_per_allele.txt.gz` the
//! `triallelic` case has four reads, the fourth tied between `C` and `G` at -1 (its `marginal`
//! rows, in hash order `C,A*,G`), and `ad triallelic` is `1,1,1` with `dp` 3, not `1,2,1`.
//!
//! `AD` reads the counts back **by allele**, so the order can only reach it through which reads are
//! informative. With the default priorities that is order-independent for finite likelihoods; the
//! exception is a NaN, which is informative or not depending on whether it sits at index 0 (see
//! the module note in `gatk_engine::allele_likelihoods`, row `best with-nan`). The order is still
//! reproduced, by [`gatk_engine::java_hash::hash_map_order`], so the marginalised matrix itself
//! matches the golden's `marginal` rows allele for allele.
//!
//! # `AD` is indexed by the variant, and after marginalising the matrix holds nothing else
//!
//! ```java
//! counts[0] = alleleCounts.get(vc.getReference()); //first one in AD is always ref
//! ```
//!
//! The counts map is keyed on `vc.getAlleles()`. `containsAll` has already refused a matrix that
//! lacks one of them, and the marginalisation keeps exactly those alleles, so every informative
//! read lands in a bucket that is read. A matrix allele the variant lacks is not counted anywhere:
//! the marginalisation drops it and a read that supported it is judged on what is left. In
//! `depth_per_allele.txt.gz`, `matrix-has-extra-allele` gives its third read -1 on `T` and -10 on
//! `A` and `C`; the `marginal` rows show that read at -10 on both survivors, a tie, and `ad` is
//! `1,1`.
//!
//! # `AF` prefers an existing `AD` and recomputes only when there is none
//!
//! ```java
//! if (g.hasAD()) { ... } else if (likelihoods != null) { ... }
//! ```
//!
//! And it drops the first entry, so `AF` has one fewer number than `AD` and the two cannot be
//! compared position by position.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a `T` or nothing". `None` stands for the reference's early `return` that
//!   leaves the genotype without the key. `x?` on an `Option` returns `None` from the whole
//!   function when `x` is `None`, and otherwise unwraps it; `.ok()?` first turns a `Result` (a
//!   value or an error) into an `Option`, so an error also becomes "no key";
//! - `let Some(x) = y else { continue };` unwraps `y` or skips to the next loop iteration;
//! - `&[Allele]` is a read-only borrowed view of a list the caller keeps;
//! - `|a| ...` is a closure, Java's `a -> ...`; `.iter().map(..).collect()` is Java's stream, map,
//!   collect, and `.iter_mut().find(..)` finds the first element and lends it out writable;
//! - `wrapping_mul` and `wrapping_add` are Java's `int` arithmetic, which wraps on overflow where
//!   plain `*` and `+` would panic in a debug build;
//! - a body's last expression, without a semicolon, is its value.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::java_hash::{byte_array_hash_code, hash_map_order};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

/// `GATKVCFConstants.ALLELE_FRACTION_KEY`.
pub const ALLELE_FRACTION_KEY: &str = "AF";

/// `Allele.hashCode`: `Arrays.hashCode(bases) * 31 + Boolean.hashCode(isRef)`.
///
/// `Boolean.hashCode` is 1231 for true and 1237 for false, which are the constants the JDK
/// documents and not an implementation detail.
pub fn allele_hash_code(allele: &Allele) -> i32 {
    let bases = byte_array_hash_code(allele.display_string().as_bytes());
    bases
        .wrapping_mul(31)
        .wrapping_add(if allele.is_reference() { 1231 } else { 1237 })
}

/// The order `Collectors.toMap(a -> a, ...)` then `keySet().toArray()` produces.
///
/// The insertion order is the stream's, which is the `LinkedHashSet` of `vc.getAlleles()`, so the
/// variant's order goes in and the hash order comes out.
pub fn marginalisation_order(alleles: &[Allele]) -> Vec<Allele> {
    let entries: Vec<(Allele, i32)> = alleles
        .iter()
        .map(|allele| (allele.clone(), allele_hash_code(allele)))
        .collect();
    // A bucket past the treeify threshold is refused rather than guessed; no variant has that
    // many alleles, so the fallback is the insertion order and is unreachable in practice.
    hash_map_order(&entries).unwrap_or_else(|_| alleles.to_vec())
}

/// `DepthPerAlleleBySample.annotateWithLikelihoods`: the `AD` array.
///
/// `None` is the reference's early return: a null genotype, an uncalled one, or no likelihoods.
pub fn allele_depths(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    sample: &str,
    called: bool,
) -> Option<Vec<i32>> {
    let likelihoods = likelihoods?;
    if !called {
        return None;
    }
    // `Utils.validateArg(likelihoods.alleles().containsAll(alleles), ...)`: the matrix must hold
    // every allele the variant declares, or the annotation throws rather than counting what it can.
    // Here the throw is a `None`, so the genotype gets no `AD` and the run goes on.
    if !vc
        .alleles
        .iter()
        .all(|allele| likelihoods.index_of_allele(allele).is_some())
    {
        return None;
    }

    // `likelihoods.marginalize(alleleSubset)`, with each new allele standing for itself alone and
    // the new allele array in `HashMap` key order; see the module note.
    let order = marginalisation_order(&vc.alleles);
    let new_to_old: Vec<(Allele, Vec<Allele>)> = order
        .iter()
        .map(|allele| (allele.clone(), vec![allele.clone()]))
        .collect();
    let marginalised = likelihoods.marginalize(&new_to_old).ok()?;

    let sample_index = marginalised.index_of_sample(sample)?;
    let mut counts: Vec<(Allele, i32)> = vc
        .alleles
        .iter()
        .map(|allele| (allele.clone(), 0))
        .collect();
    // `bestAllelesBreakingTies(sample).stream().filter(ba -> ba.isInformative())`, then one count
    // per read for its best allele. `None` asks for the default priorities (reference first).
    for best in marginalised.best_alleles_breaking_ties_for_sample(sample_index, None) {
        if !best.is_informative() {
            continue;
        }
        // Neither skip below can be taken: an informative read always has a best allele, and the
        // marginalised matrix holds only the variant's alleles, which are all keys here.
        let Some(allele) = best.allele else { continue };
        if let Some(slot) = counts.iter_mut().find(|(a, _)| *a == allele) {
            slot.1 += 1;
        }
    }

    // "first one in AD is always ref", then the alternates in the variant's order.
    let mut out = Vec::with_capacity(counts.len());
    let reference = vc.alleles.iter().find(|a| a.is_reference())?;
    out.push(
        counts
            .iter()
            .find(|(a, _)| a == reference)
            .map(|(_, c)| *c)
            .unwrap_or(0),
    );
    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        out.push(
            counts
                .iter()
                .find(|(a, _)| a == allele)
                .map(|(_, c)| *c)
                .unwrap_or(0),
        );
    }
    Some(out)
}

/// `MathUtils.normalizeSumToOne`: each value divided by the sum.
///
/// The reference's one guard, `Utils.validateArg(sum >= 0.0, ...)`, refuses a negative sum. This
/// function has no such refusal, which matters only for a genotype `AD` whose entries sum below
/// zero; one counted by [`allele_depths`] never does.
///
/// There is no guard for a zero sum: an all-zero `AD` divides zero by zero, so a site with no
/// informative reads gets an `AF` of `NaN` rather than of zero. The golden's `empty` and
/// `uninformative` rows are both that, and a consumer reading `AF` as a number has to handle it.
fn normalize_sum_to_one(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let sum: f64 = values.iter().sum();
    values.iter().map(|value| value / sum).collect()
}

/// `AlleleFraction.annotate`: the normalised depths with the **reference dropped**.
///
/// `genotype_ad` is the genotype's own `AD` when it has one, which wins over the likelihoods.
///
/// The reference's first test, `g.hasExtendedAttribute(getKeyNames().get(0))` (an `AF` already on
/// the genotype is never overwritten), is not here: the caller in `gatk-tools`'
/// `variant_annotator_engine.rs` makes it before calling.
pub fn allele_fractions(
    vc: &VariantContext,
    genotype_ad: Option<&[i32]>,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    sample: &str,
    called: bool,
) -> Option<Vec<f64>> {
    if !called {
        return None;
    }
    let depths: Vec<i32> = match genotype_ad {
        Some(ad) => ad.to_vec(),
        None => allele_depths(vc, likelihoods, sample, called)?,
    };
    // `Arrays.stream(AD).mapToDouble(x -> x)`: `as f64` widens each `int` exactly.
    let doubles: Vec<f64> = depths.iter().map(|d| *d as f64).collect();
    let fractions = normalize_sum_to_one(&doubles);
    // `Arrays.copyOfRange(all, 1, all.length)`: the reference's entry is dropped, so AF is one
    // shorter than AD and the two do not line up by index.
    Some(fractions[1..].to_vec())
}

/// `DepthPerSampleHC`: the `DP` the HaplotypeCaller writes, which is **informative reads only**.
///
/// Not the pileup depth: a read whose likelihoods are within the informative threshold of each
/// other is not counted, so `DP` here can be lower than the `DP` an INFO-level `Coverage` reports
/// over the same site.
pub fn informative_depth(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    sample: &str,
    called: bool,
) -> Option<i32> {
    let likelihoods = likelihoods?;
    if !called {
        return None;
    }
    let sample_index = likelihoods.index_of_sample(sample)?;
    // `gb.DP(0)`: a sample with no reads gets a written `DP=0`, not an absent key.
    if likelihoods.sample_evidence_count(sample_index) == 0 {
        return Some(0);
    }
    // The variant's alleles must all be in the matrix, which is then marginalized to them: a read
    // tied between an allele the site dropped and one it kept is informative once the dropped one
    // is gone.
    if !vc
        .alleles
        .iter()
        .all(|allele| likelihoods.index_of_allele(allele).is_some())
    {
        return None;
    }
    let new_to_old: Vec<(Allele, Vec<Allele>)> = marginalisation_order(&vc.alleles)
        .iter()
        .map(|allele| (allele.clone(), vec![allele.clone()]))
        .collect();
    let marginalised = likelihoods.marginalize(&new_to_old).ok()?;
    let sample_index = marginalised.index_of_sample(sample)?;
    Some(
        marginalised
            .best_alleles_breaking_ties_for_sample(sample_index, None)
            .iter()
            .filter(|best| best.is_informative())
            .count() as i32,
    )
}

/// Unused parameter kept so the signature matches the reference's.
pub fn reference_context_is_ignored(_reference: Option<&ReferenceContext>) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_allele_hash_is_the_bases_and_the_reference_flag() {
        let a_ref = Allele::from_str("A", true).unwrap();
        let a_alt = Allele::from_str("A", false).unwrap();
        // Same bases, different flag, different hash: the flag is worth 1231 against 1237.
        assert_eq!(allele_hash_code(&a_ref) - allele_hash_code(&a_alt), -6);
    }

    #[test]
    fn the_marginalisation_order_is_not_the_variants_order() {
        let alleles = vec![
            Allele::from_str("A", true).unwrap(),
            Allele::from_str("C", false).unwrap(),
            Allele::from_str("G", false).unwrap(),
            Allele::from_str("T", false).unwrap(),
        ];
        let order = marginalisation_order(&alleles);
        assert_eq!(order.len(), 4);
        // Every allele survives; whether the order moved is what the golden measures.
        for allele in &alleles {
            assert!(order.contains(allele));
        }
    }
}
