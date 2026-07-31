//! `UniqueAltReadCount`, `BaseQualityHistogram` and `ReferenceBases`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! Three annotations that read the likelihood matrix or the reference window and nothing else, and
//! whose output shape is more interesting than their arithmetic.
//!
//! # `AS_UNIQ_ALT_READ_COUNT` counts **distinct fragments**, not reads
//!
//! ```java
//! .map(ba -> new ImmutablePair<>(ba.evidence.getStart(), ba.evidence.getFragmentLength()))
//! .collect(Collectors.groupingBy(x -> x, Collectors.counting()));
//! return duplicateReadMap.size();
//! ```
//!
//! Two reads with the same start and the same template length are one entry, so a hundred PCR
//! duplicates of one fragment count once. The pair is `(start, fragmentLength)` and not the read
//! name, so two genuinely distinct fragments that happen to share both are also one.
//!
//! And the value is a **string** joined with `|`, because it is an allele-specific raw annotation.
//! The reference's own comment on the bracket-stripping regex is *"Who actually wants brackets at
//! the ends of their string? Who???"*.
//!
//! # `BQHIST` is a flat list, allele-major inside quality-major
//!
//! ```text
//! [q1, count(ref, q1), count(alt, q1), q2, count(ref, q2), count(alt, q2), ...]
//! ```
//!
//! One entry per distinct quality seen anywhere, then one count per allele of the **matrix**, in
//! matrix order. So the list's length depends on how many alleles the matrix holds rather than on
//! how many the variant declares, and reading it requires knowing that number.
//!
//! # `REF_BASES` pads with `N` and can be off-centre
//!
//! ```java
//! final int basesToDiscardInFront = max(vc.getStart() - ref.getWindow().getStart() - 10, 0);
//! ```
//!
//! Twenty-one bases centred on the variant, taken out of whatever window the reference context
//! happens to carry. If the window starts too late the discard clamps to zero and the string is
//! *not* centred on the variant; if it ends too early the string is padded on the right with `N`
//! to twenty-one characters. Either way the length is fixed and the centring is not.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::read::mapping_quality;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use crate::per_allele::BaseQuality;

/// `GATKVCFConstants.AS_UNIQUE_ALT_READ_SET_COUNT_KEY`.
pub const AS_UNIQUE_ALT_READ_SET_COUNT_KEY: &str = "AS_UNIQ_ALT_READ_COUNT";
/// `GATKVCFConstants.BASE_QUAL_HISTOGRAM_KEY`.
pub const BASE_QUAL_HISTOGRAM_KEY: &str = "BQHIST";
/// `GATKVCFConstants.REFERENCE_BASES_KEY`.
pub const REFERENCE_BASES_KEY: &str = "REF_BASES";

/// `AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM`.
const ALLELE_SPECIFIC_RAW_DELIM: &str = "|";

/// `ReferenceBases.NUM_BASES_ON_EITHER_SIDE`.
const NUM_BASES_ON_EITHER_SIDE: i64 = 10;
const REFERENCE_CONTEXT_LENGTH: usize = (2 * NUM_BASES_ON_EITHER_SIDE + 1) as usize;

/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// `UniqueAltReadCount`: one count per **alternate** allele, joined with `|`.
pub struct UniqueAltReadCount;

impl InfoFieldAnnotation for UniqueAltReadCount {
    fn key_names(&self) -> Vec<&'static str> {
        vec![AS_UNIQUE_ALT_READ_SET_COUNT_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // No null guard: the reference dereferences the likelihoods straight away, so a null
        // matrix is a NullPointerException there. Here it is an empty answer, and the suite has
        // no row for it because the dump cannot produce one without crashing.
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        let counts: Vec<String> = vc
            .alleles
            .iter()
            .filter(|allele| !allele.is_reference())
            .map(|alt| {
                let mut fragments: Vec<(i32, i32)> = Vec::new();
                for best in likelihoods.best_alleles_breaking_ties(None) {
                    if best.allele.as_ref() != Some(alt) || !best.is_informative() {
                        continue;
                    }
                    let Some(read) = likelihoods
                        .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
                        .and_then(|reads| reads.get(best.evidence_index))
                    else {
                        continue;
                    };
                    // `(getStart(), getFragmentLength())`, which is the pair the grouping keys on.
                    let key = (read.alignment_start, read.inferred_insert_size);
                    if !fragments.contains(&key) {
                        fragments.push(key);
                    }
                }
                fragments.len().to_string()
            })
            .collect();
        vec![(
            AS_UNIQUE_ALT_READ_SET_COUNT_KEY.to_string(),
            // `StringUtils.join(list, "|")` after the brackets are stripped: a String, and it is
            // one string even for a single alternate allele.
            AnnotationValue::Str(counts.join(ALLELE_SPECIFIC_RAW_DELIM)),
        )]
    }
}

/// `BaseQualityHistogram`: `BQHIST`, a flat list of quality then one count per matrix allele.
pub struct BaseQualityHistogram;

impl InfoFieldAnnotation for BaseQualityHistogram {
    fn key_names(&self) -> Vec<&'static str> {
        vec![BASE_QUAL_HISTOGRAM_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        // One multiset per allele of the **matrix**, so the row length follows the matrix.
        let allele_count = likelihoods.number_of_alleles();
        let mut per_allele: Vec<Vec<i32>> = vec![Vec::new(); allele_count];

        for best in likelihoods.best_alleles_breaking_ties(None) {
            if !best.is_informative() {
                continue;
            }
            let Some(read) = likelihoods
                .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            if !is_usable_read(read) {
                continue;
            }
            let Some(allele) = best.allele.as_ref() else {
                continue;
            };
            let Some(index) = likelihoods.index_of_allele(allele) else {
                continue;
            };
            if let Some(quality) = BaseQuality::base_quality(read, vc) {
                per_allele[index].push(quality);
            }
        }

        // `distinct().sorted()` over every allele's values, in matrix order.
        let mut qualities: Vec<i32> = per_allele.iter().flatten().copied().collect();
        qualities.sort_unstable();
        qualities.dedup();

        let mut output: Vec<AnnotationValue> = Vec::new();
        for quality in qualities {
            output.push(AnnotationValue::Int(quality));
            for values in &per_allele {
                let count = values.iter().filter(|value| **value == quality).count();
                output.push(AnnotationValue::Int(count as i32));
            }
        }
        vec![(
            BASE_QUAL_HISTOGRAM_KEY.to_string(),
            // `ImmutableMap.of(KEY, List<Integer>)`: a list, not a joined string.
            AnnotationValue::List(output),
        )]
    }
}

/// `OrientationBiasReadCounts.isUsableRead`, which `BaseQualityHistogram` reuses.
fn is_usable_read(read: &BamRecord) -> bool {
    let quality = mapping_quality(read);
    quality != 0 && quality != MAPPING_QUALITY_UNAVAILABLE
}

/// `ReferenceBases`: `REF_BASES`, twenty-one bases around the variant, padded with `N`.
pub struct ReferenceBases;

impl ReferenceBases {
    /// `ReferenceBases.annotate(ref, vc)`, over the window's bases as the caller already has them.
    ///
    /// The bases are passed in rather than fetched, because the reference context's own fetch
    /// needs a file source and this computation does not care where they came from.
    pub fn local_bases(window_start: i64, window_bases: &[u8], vc: &VariantContext) -> String {
        let discard = (vc.start - window_start - NUM_BASES_ON_EITHER_SIDE).max(0) as usize;
        let all: String = String::from_utf8_lossy(window_bases).into_owned();
        // `substring` on a start past the end throws; the reference never guards it, and neither
        // does this, because the window always covers the variant when a walker built it.
        let end = (discard + REFERENCE_CONTEXT_LENGTH).min(all.len());
        let mut local = all[discard.min(all.len())..end].to_string();
        if local.len() < REFERENCE_CONTEXT_LENGTH {
            // Padded on the **right** only, so a variant near the start of a contig is annotated
            // with a string that is not centred on it and says nothing about the fact.
            local.push_str(&"N".repeat(REFERENCE_CONTEXT_LENGTH - local.len()));
        }
        local
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn site(start: i64) -> VariantContext {
        let mut vc = VariantContext::new(
            "chr1",
            start,
            vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("C", false).unwrap(),
            ],
        );
        vc.stop = start;
        vc
    }

    #[test]
    fn the_window_is_centred_when_it_can_be_and_padded_when_it_cannot() {
        let bases: Vec<u8> = (0..30).map(|i| b"ACGT"[i % 4]).collect();
        // A variant ten bases into the window: centred, twenty-one bases.
        let centred = ReferenceBases::local_bases(100, &bases, &site(110));
        assert_eq!(centred.len(), 21);
        assert!(!centred.contains('N'));
        // A variant at the window's start: the discard clamps to zero, so it is not centred.
        let off_centre = ReferenceBases::local_bases(100, &bases, &site(100));
        assert_eq!(off_centre.len(), 21);
        assert_eq!(&off_centre[..4], "ACGT");
        // A window too short on the right: padded with N.
        let padded = ReferenceBases::local_bases(100, &bases[..12], &site(105));
        assert_eq!(padded.len(), 21);
        assert!(padded.ends_with('N'));
    }
}
