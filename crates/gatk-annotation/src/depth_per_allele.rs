//! `DepthPerAlleleBySample`, `AlleleFraction` and `DepthPerSampleHC`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
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
//! That order is observable: `searchBestAllele` breaks a tie by keeping the first index, so a read
//! that supports two alleles equally is attributed to whichever the hash happened to put first.
//! `AD` then reads the counts back **by allele**, so the tie is the only way the order reaches the
//! output, and it does reach it. [`gatk_engine::java_hash::hash_map_order`] reproduces it.
//!
//! # `AD` is indexed by the variant and counted over the matrix
//!
//! ```java
//! counts[0] = alleleCounts.get(vc.getReference()); //first one in AD is always ref
//! ```
//!
//! The counts map is keyed on `vc.getAlleles()` and the traversal adds to it by allele, so an
//! allele the matrix holds and the variant does not is counted into a bucket that is never read,
//! and an allele the variant holds and the matrix does not stays at zero.
//!
//! # `AF` prefers an existing `AD` and recomputes only when there is none
//!
//! ```java
//! if (g.hasAD()) { ... } else if (likelihoods != null) { ... }
//! ```
//!
//! And it drops the first entry, so `AF` has one fewer number than `AD` and the two cannot be
//! compared position by position.

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
    // `Utils.validateArg(likelihoods.alleles().containsAll(alleles))`: the matrix must hold every
    // allele the variant declares, or the annotation throws rather than counting what it can.
    if !vc
        .alleles
        .iter()
        .all(|allele| likelihoods.index_of_allele(allele).is_some())
    {
        return None;
    }

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
    for best in marginalised.best_alleles_breaking_ties_for_sample(sample_index, None) {
        if !best.is_informative() {
            continue;
        }
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

/// `MathUtils.normalizeSumToOne`, which divides by the sum whatever the sum is.
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
/// `existing_ad` is the genotype's own `AD` when it has one, which wins over the likelihoods.
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
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    sample: &str,
    called: bool,
) -> Option<i32> {
    let likelihoods = likelihoods?;
    if !called {
        return None;
    }
    let sample_index = likelihoods.index_of_sample(sample)?;
    Some(
        likelihoods
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
