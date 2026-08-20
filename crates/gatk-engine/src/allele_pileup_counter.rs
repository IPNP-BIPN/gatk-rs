//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.AllelePileupCounter`
//! (GATK 4.6.2.0).
//!
//! How many reads of a pileup carry each allele, when the alleles are known before the counting.
//!
//! # What it drops, and what it does not
//!
//! One filter on the read: a mapping quality of zero or the unavailable 255. Nothing else is
//! filtered out. A base under the quality cutoff is not dropped from the pileup; it is chosen as no
//! allele by [`crate::variant_context_utils::choose_allele_for_read`] and therefore increments
//! nothing, which looks the same in the counts and is not the same thing at all: the cutoff belongs
//! to the choice, and the mapping quality belongs here.
//!
//! # The map is not grown by counting
//!
//! Only the alleles the counter was built with have entries, and the increment checks for the key
//! rather than inserting. An alternate the pileup carries but the caller did not ask about is
//! counted nowhere.

use crate::read_pileup::ReadPileup;
use crate::variant_context_utils::{choose_allele_for_read, Allele, PileupAlleleError};

/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
pub const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// What building a counter refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum CounterError {
    /// `UserException.BadInput`. UNREACHABLE through htsjdk, which refuses to tag a symbolic
    /// allele as the reference before the counter is ever built, and the golden carries htsjdk's
    /// refusal instead. Ported because the branch is there.
    SymbolicReference,
    /// `Utils.validateArg` on a reference allele that is not marked as the reference.
    NonReferenceReference(Allele),
    /// `Utils.validateArg` on an alternate that is marked as the reference.
    ReferenceAlternate(Vec<Allele>),
    /// `ParamUtils.isPositiveOrZero` on the cutoff.
    NegativeMinimumBaseQuality,
}

impl CounterError {
    pub fn java_class(&self) -> &'static str {
        match self {
            CounterError::SymbolicReference => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            _ => "java.lang.IllegalArgumentException",
        }
    }

    /// The message, with alleles written the way `Allele.toString` writes them: the bases, and a
    /// trailing `*` when the allele is the reference.
    pub fn message(&self) -> String {
        match self {
            CounterError::SymbolicReference => {
                "A symbolic reference allele was specified.".to_string()
            }
            CounterError::NonReferenceReference(allele) => {
                format!(
                    "Reference allele was non-reference: {}",
                    java_string(allele)
                )
            }
            CounterError::ReferenceAlternate(alleles) => format!(
                "One or more alternate alleles were reference: {}",
                alleles
                    .iter()
                    .map(java_string)
                    .collect::<Vec<String>>()
                    .join(", ")
            ),
            CounterError::NegativeMinimumBaseQuality => {
                "Minimum base quality must be positive or zero.".to_string()
            }
        }
    }
}

/// `Allele.toString()`: the bases, then a `*` for the reference.
fn java_string(allele: &Allele) -> String {
    let mut text = String::from_utf8_lossy(&allele.bases).into_owned();
    if allele.is_reference {
        text.push('*');
    }
    text
}

/// `AllelePileupCounter`.
#[derive(Debug, Clone)]
pub struct AllelePileupCounter {
    reference: Allele,
    alternates: Vec<Allele>,
    minimum_base_quality: i32,
    /// The alternates in the order they were given, then the reference, which is the order the
    /// reference's own `forEach` plus `put` fills a map in. A `HashMap` has no order of its own, so
    /// a caller that needs one sorts; this keeps insertion order so that nothing here is arbitrary.
    counts: Vec<(Allele, i32)>,
}

impl AllelePileupCounter {
    /// The constructor, with its four checks in the reference's order.
    pub fn new(
        reference: &Allele,
        alternates: &[Allele],
        minimum_base_quality: i32,
    ) -> Result<AllelePileupCounter, CounterError> {
        if reference.is_symbolic() {
            return Err(CounterError::SymbolicReference);
        }
        if !reference.is_reference {
            return Err(CounterError::NonReferenceReference(reference.clone()));
        }
        if alternates.iter().any(|allele| allele.is_reference) {
            return Err(CounterError::ReferenceAlternate(alternates.to_vec()));
        }
        if minimum_base_quality < 0 {
            return Err(CounterError::NegativeMinimumBaseQuality);
        }
        let mut counts: Vec<(Allele, i32)> = alternates
            .iter()
            .map(|allele| (allele.clone(), 0))
            .collect();
        counts.push((reference.clone(), 0));
        Ok(AllelePileupCounter {
            reference: reference.clone(),
            alternates: alternates.to_vec(),
            minimum_base_quality,
            counts,
        })
    }

    /// The constructor that counts a pileup straight away.
    pub fn with_pileup(
        reference: &Allele,
        alternates: &[Allele],
        minimum_base_quality: i32,
        pileup: &ReadPileup<'_>,
    ) -> Result<AllelePileupCounter, CounterError> {
        let mut counter = AllelePileupCounter::new(reference, alternates, minimum_base_quality)?;
        counter.add_pileup(pileup).map_err(|_| {
            // The only refusal underneath is the negative cutoff, which the constructor has
            // already rejected, so this cannot be reached.
            CounterError::NegativeMinimumBaseQuality
        })?;
        Ok(counter)
    }

    /// `addPileup`. A pileup the reference is never given is a null there and nothing here.
    pub fn add_pileup(&mut self, pileup: &ReadPileup<'_>) -> Result<(), PileupAlleleError> {
        if self.reference.is_symbolic() {
            return Ok(());
        }
        for element in &pileup.elements {
            if !is_usable_read(element.mapping_qual()) {
                continue;
            }
            let chosen = choose_allele_for_read(
                element,
                &self.reference,
                &self.alternates,
                self.minimum_base_quality,
            )?;
            let Some(allele) = chosen else {
                continue;
            };
            // `containsKey` and then `increment`: an allele the map does not hold is not added.
            if let Some(entry) = self.counts.iter_mut().find(|(key, _)| *key == allele) {
                entry.1 += 1;
            }
        }
        Ok(())
    }

    /// `getCountMap`, as pairs in insertion order.
    pub fn count_map(&self) -> &[(Allele, i32)] {
        &self.counts
    }

    /// The count of one allele, or `None` when the map has no such key.
    pub fn count(&self, allele: &Allele) -> Option<i32> {
        self.counts
            .iter()
            .find(|(key, _)| key == allele)
            .map(|(_, count)| *count)
    }
}

/// `isUsableRead`: a mapping quality of zero or the unavailable 255 is not usable.
pub fn is_usable_read(mapping_quality: u8) -> bool {
    mapping_quality != 0 && mapping_quality != MAPPING_QUALITY_UNAVAILABLE
}
