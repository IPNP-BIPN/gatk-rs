//! `OrientationBiasReadCounts` and `FragmentDepthPerAlleleBySample`, ported from GATK 4.6.2.0.
//!
//! `F1R2`, `F2R1` and `FAD`: the two annotations that count **fragments** rather than reads. Both
//! are `JumboGenotypeAnnotation`s, which is the interface that receives a fragment-indexed
//! likelihood matrix alongside the read-indexed one.
//!
//! # A read pair votes once, and the pair's orientation is read off its **first** read
//!
//! ```java
//! // read orientation is a property of fragments; that is, if one read is F1 the other read is R2,
//! // hence both reads of the fragment are F1R2.  Thus we can arbitrarily take the first read of each fragment
//! ```
//!
//! "Arbitrarily" is doing real work: the first read is the first in the sample's own order, so which
//! of a pair is consulted depends on the order the reads were added to the matrix, and every test on
//! that fragment (usability, base quality, orientation) is applied to that read alone. A pair whose
//! second read has a low base quality at the site is counted; the same pair added the other way
//! round is not.
//!
//! # `F1R2` and `F2R1` are **not** complementary
//!
//! `isF2R1` is `read.isReverseStrand() == read.isFirstOfPair()`, and an unpaired read has
//! `isFirstOfPair()` false. So a forward unpaired read lands in `F2R1` and a reverse one in `F1R2`,
//! which is the opposite of what the names suggest. Both keys are always written, one entry per
//! allele of the **variant**, so the two arrays have the same length and their sum is not the depth:
//! a fragment that fails any filter is in neither.
//!
//! # The counts are keyed on the **matrix's** alleles and read back by the variant's
//!
//! ```java
//! final int[] f1r2 = vc.getAlleles().stream().mapToInt(a -> f1r2Counts.get(a).intValue()).toArray();
//! ```
//!
//! The maps are built from `fragmentLikelihoods.alleles()`, so an allele the variant declares and
//! the matrix does not is a `NullPointerException` here rather than a zero. The reverse, an allele
//! the matrix holds and the variant does not, is counted into a bucket nothing reads.
//!
//! # `FAD` is `AD` over fragments, and reuses the same marginalisation
//!
//! `FragmentDepthPerAlleleBySample` calls `DepthPerAlleleBySample.annotateWithLikelihoods` with the
//! fragment matrix, so it inherits the `HashMap` allele order that [`crate::depth_per_allele`]
//! records, and differs from `AD` only in what a "one" counts.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::fragment::Fragment;
use gatk_engine::read::mapping_quality;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

use crate::per_allele::BaseQuality;

/// `GATKVCFConstants.F1R2_KEY`.
pub const F1R2_KEY: &str = "F1R2";
/// `GATKVCFConstants.F2R1_KEY`.
pub const F2R1_KEY: &str = "F2R1";
/// `GATKVCFConstants.FRAGMENT_ALLELE_DEPTHS`.
pub const FRAGMENT_ALLELE_DEPTHS_KEY: &str = "FAD";

/// `OrientationBiasReadCounts.MINIMUM_BASE_QUALITY`.
const MINIMUM_BASE_QUALITY: i32 = 20;
/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// What these annotations refuse.
#[derive(Debug, Clone, PartialEq)]
pub enum FragmentCountError {
    /// The `NullPointerException` from reading a count for an allele the matrix does not hold.
    AlleleNotInMatrix { allele: String },
    /// `Utils.validateArg(fragmentLikelihoods.alleles().containsAll(alleles))`.
    AllelesNotASubset,
}

/// `OrientationBiasReadCounts.isUsableRead`.
fn is_usable_read(read: &BamRecord) -> bool {
    let quality = mapping_quality(read);
    quality != 0 && quality != MAPPING_QUALITY_UNAVAILABLE
}

/// `OrientationBiasReadCounts.annotate`: the `(F1R2, F2R1)` pair of per-allele counts.
///
/// `None` is the reference's early return: a null genotype or a null fragment matrix, which writes
/// neither key.
pub fn orientation_bias_counts(
    vc: &VariantContext,
    sample: &str,
    fragment_likelihoods: Option<&AlleleLikelihoods<Fragment>>,
) -> Result<Option<(Vec<i32>, Vec<i32>)>, FragmentCountError> {
    let Some(likelihoods) = fragment_likelihoods else {
        return Ok(None);
    };
    let Some(sample_index) = likelihoods.index_of_sample(sample) else {
        // `bestAllelesBreakingTies(sampleName)` on an unknown sample is an
        // `IllegalArgumentException` upstream; no caller reaches it.
        return Ok(None);
    };

    let matrix_alleles: Vec<Allele> = (0..likelihoods.number_of_alleles())
        .filter_map(|index| likelihoods.get_allele(index).cloned())
        .collect();
    let mut f1r2: Vec<(Allele, i32)> = matrix_alleles
        .iter()
        .map(|allele| (allele.clone(), 0))
        .collect();
    let mut f2r1 = f1r2.clone();

    for best in likelihoods.best_alleles_breaking_ties_for_sample(sample_index, None) {
        if !best.is_informative() {
            continue;
        }
        let Some(fragment) = likelihoods
            .sample_evidence(sample_index)
            .and_then(|fragments| fragments.get(best.evidence_index))
        else {
            continue;
        };
        // "arbitrarily take the first read of each fragment": every test below is on this one read.
        let Some(read) = fragment.reads.first() else {
            continue;
        };
        if !is_usable_read(read) {
            continue;
        }
        // `.orElse(0)`, so a read with no base at the site scores zero and fails the threshold
        // rather than being skipped for a different reason.
        let quality = BaseQuality::base_quality(read, vc).unwrap_or(0);
        if quality < MINIMUM_BASE_QUALITY {
            continue;
        }
        let Some(allele) = best.allele.as_ref() else {
            continue;
        };
        let table = if Fragment::is_f2r1(read) {
            &mut f2r1
        } else {
            &mut f1r2
        };
        if let Some(slot) = table.iter_mut().find(|(a, _)| a == allele) {
            slot.1 += 1;
        }
    }

    // Read back by the **variant's** alleles: a missing one is a NullPointerException.
    let mut f1r2_out = Vec::with_capacity(vc.alleles.len());
    let mut f2r1_out = Vec::with_capacity(vc.alleles.len());
    for allele in &vc.alleles {
        let Some((_, first)) = f1r2.iter().find(|(a, _)| a == allele) else {
            return Err(FragmentCountError::AlleleNotInMatrix {
                allele: allele.display_string(),
            });
        };
        let Some((_, second)) = f2r1.iter().find(|(a, _)| a == allele) else {
            return Err(FragmentCountError::AlleleNotInMatrix {
                allele: allele.display_string(),
            });
        };
        f1r2_out.push(*first);
        f2r1_out.push(*second);
    }
    Ok(Some((f1r2_out, f2r1_out)))
}

/// `FragmentDepthPerAlleleBySample.annotate`: `FAD`.
///
/// `None` is the reference's early return: no genotype, an uncalled one, or no fragment matrix.
/// The subset check is an `IllegalArgumentException` rather than a quiet skip, modelled as an error.
pub fn fragment_allele_depths(
    vc: &VariantContext,
    sample: &str,
    called: bool,
    fragment_likelihoods: Option<&AlleleLikelihoods<Fragment>>,
) -> Result<Option<Vec<i32>>, FragmentCountError> {
    let Some(likelihoods) = fragment_likelihoods else {
        return Ok(None);
    };
    if !called {
        return Ok(None);
    }
    if !vc
        .alleles
        .iter()
        .all(|allele| likelihoods.index_of_allele(allele).is_some())
    {
        return Err(FragmentCountError::AllelesNotASubset);
    }

    // The same marginalisation `AD` uses, over the same `HashMap` allele order.
    let order = crate::depth_per_allele::marginalisation_order(&vc.alleles);
    let new_to_old: Vec<(Allele, Vec<Allele>)> = order
        .iter()
        .map(|allele| (allele.clone(), vec![allele.clone()]))
        .collect();
    let Ok(marginalised) = likelihoods.marginalize(&new_to_old) else {
        return Ok(None);
    };
    let Some(sample_index) = marginalised.index_of_sample(sample) else {
        return Ok(None);
    };

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

    let reference = vc.alleles.iter().find(|a| a.is_reference());
    let mut out = Vec::with_capacity(counts.len());
    if let Some(reference) = reference {
        out.push(
            counts
                .iter()
                .find(|(a, _)| a == reference)
                .map(|(_, count)| *count)
                .unwrap_or(0),
        );
    }
    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        out.push(
            counts
                .iter()
                .find(|(a, _)| a == allele)
                .map(|(_, count)| *count)
                .unwrap_or(0),
        );
    }
    Ok(Some(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_two_orientation_keys_are_not_complementary_for_an_unpaired_read() {
        let mut forward = BamRecord {
            read_name: "r".to_string(),
            flags: 0,
            ..Default::default()
        };
        assert!(Fragment::is_f2r1(&forward));
        forward.flags = 0x10;
        assert!(!Fragment::is_f2r1(&forward));
    }
}
