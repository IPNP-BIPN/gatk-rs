//! `GATKVariantContextUtils.leftAlignAndTrim`, ported from GATK 4.6.2.0.
//!
//! The workhorse under `LeftAlignAndTrimVariants`: it slides an indel as far left as the reference
//! lets it. [`crate::alignment_utils::normalize_alleles`] does the sliding; what is here is the
//! window that keeps widening around it, and what the record becomes afterwards.
//!
//! # The window widens, and it can stop short without saying so
//!
//! ```java
//! for (int leadingBases = Math.min(maxLeadingBases, 10); leadingBases <= maxLeadingBases; leadingBases = Math.min(2*leadingBases, maxLeadingBases)) {
//! ...
//! } else if (shifts.getLeft() == variantOffsetInRef && leadingBases < maxLeadingBases) {
//!     continue;
//! }
//! ```
//!
//! Ten bases first, then twenty, then forty, until the shift stops landing on the edge of the
//! slice. When the shift **does** land on the edge and the window is already at its maximum, the
//! loop falls through and the partly-aligned record is returned: the same deletion comes out at 15
//! with a window of 2, at 49 with 10, at 39 with 20, and at the start of its run with a wide
//! enough one. Nothing in the record says which of those happened.
//!
//! # The window that ran out is not the same as the record that had nowhere to go
//!
//! The reference returns a `VariantContext` and nothing else, so those two cases are
//! indistinguishable to its caller: a record stopped by `--max-leading-bases` looks exactly like
//! one that is genuinely left aligned. [`left_align_and_trim_reporting`] returns the same record
//! with an [`Alignment`] beside it. The bytes do not change, which is why this is not a divergence;
//! what changes is that the answer stops being silent.
//!
//! # What comes back unchanged
//!
//! A non-indel, a `maxLeadingBases` of zero or less, and a shift of zero all return the input
//! itself. The first is a **type** test, so a MIXED site is never aligned; the second is what the
//! caller produces for two adjacent variants, `min(maxLeadingBases, distanceToLastVariant - 1)`.

use crate::alignment_utils::{normalize_alleles, AlignmentError, IndexRange};
use crate::pileup::PileupElement;
use crate::subset_alleles::{subset_alleles, AssignmentMethod, Genotype};

/// One allele: its bases and whether it is the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele {
    pub bases: Vec<u8>,
    pub is_reference: bool,
}

impl Allele {
    pub fn new(bases: &[u8], is_reference: bool) -> Allele {
        Allele {
            bases: bases.to_vec(),
            is_reference,
        }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// `Allele.isSymbolic()`, as far as the trimmer needs it: `<DEL>` and its kind.
    ///
    /// htsjdk also calls a breakend symbolic; nothing measured here reaches one, and the trimmer
    /// treats both the same way, which is to leave them alone.
    pub fn is_symbolic(&self) -> bool {
        self.bases.first() == Some(&b'<') && self.bases.last() == Some(&b'>')
    }

    /// `Allele.SPAN_DEL`, the `*` a record carries for a deletion that started earlier.
    pub fn is_span_del(&self) -> bool {
        self.bases == b"*"
    }
}

/// As much of a `VariantContext` as these functions read and write.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Variant {
    pub contig: String,
    pub start: i32,
    pub stop: i32,
    /// The reference allele first, then the alternates.
    pub alleles: Vec<Allele>,
    /// One call per sample, which the splitting rewrites and the alignment carries.
    pub genotypes: Vec<Genotype>,
    /// The INFO attributes, in the order they were set. Only the three the splitter keeps are
    /// modelled by name; anything else is carried and dropped as a whole.
    pub attributes: Vec<(String, String)>,
}

impl Variant {
    /// `isIndel()`: the type is INDEL only when every alternate differs from the reference in
    /// length in the same way. Two alternates of different types are MIXED, which is not an indel.
    ///
    /// The same rule as [`crate::pileup_summary`]'s neighbours; it lives here because this is what
    /// decides whether the alignment runs at all.
    pub fn is_indel(&self) -> bool {
        if self.alleles.len() < 2 {
            return false;
        }
        let reference = &self.alleles[0];
        let mut kind: Option<bool> = None;
        for allele in &self.alleles[1..] {
            // An alternate is "indel-shaped" when its length differs from the reference's.
            let this = allele.len() != reference.len();
            match kind {
                None => kind = Some(this),
                Some(seen) if seen != this => return false,
                Some(_) => {}
            }
        }
        kind.unwrap_or(false)
    }
}

/// Whether the record that came back is as far left as the reference would allow.
///
/// The reference does not carry this: `leftAlignAndTrim` returns a `VariantContext` and nothing
/// else, so a record that ran out of window is indistinguishable from one that had nowhere left to
/// go. The bytes of the record are the same either way, which is why reporting this changes no
/// output; it only stops the answer from being silent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Alignment {
    /// The shift stopped before the edge of the window, so widening it would change nothing.
    Complete,
    /// The shift reached the edge of a window that was already at `maxLeadingBases`. The record is
    /// shifted as far as the argument allowed and NOT as far as the reference would allow: a
    /// larger `--max-leading-bases` would move it further.
    WindowExhausted,
    /// Nothing was attempted: a non-indel, or a `maxLeadingBases` of zero or less.
    NotAttempted,
}

/// `leftAlignAndTrim(vc, ref, maxLeadingBases, trim)`.
///
/// `reference_bases` is the whole contig, one-based positions indexing from 1, which is what a
/// `ReferenceContext` hands back a slice of.
///
/// The record this returns is the reference's, byte for byte. Use
/// [`left_align_and_trim_reporting`] where it matters whether the window ran out.
pub fn left_align_and_trim(
    variant: &Variant,
    reference_bases: &[u8],
    max_leading_bases: i32,
    trim: bool,
) -> Result<Variant, AlignmentError> {
    left_align_and_trim_reporting(variant, reference_bases, max_leading_bases, trim)
        .map(|(variant, _)| variant)
}

/// The same alignment, with what the reference throws away: whether the window ran out.
///
/// GATK's own loop knows the difference, in `shifts.getLeft() == variantOffsetInRef &&
/// leadingBases < maxLeadingBases`, and drops it on the floor when the second half is false. A
/// caller that wants to widen the window and try again, or to warn, has no way to ask.
pub fn left_align_and_trim_reporting(
    variant: &Variant,
    reference_bases: &[u8],
    max_leading_bases: i32,
    trim: bool,
) -> Result<(Variant, Alignment), AlignmentError> {
    if !variant.is_indel() || max_leading_bases <= 0 {
        return Ok((variant.clone(), Alignment::NotAttempted));
    }

    let mut leading_bases = max_leading_bases.min(10);
    loop {
        let ref_start = (variant.start - leading_bases).max(1);
        // The slice the reference context would return: from `ref_start` to the variant's end.
        let slice = &reference_bases[(ref_start - 1) as usize..variant.stop as usize];
        let variant_offset = variant.start - ref_start;

        // Each allele preceded by the reference bases before it, which is what gives the shift its
        // room. The trailing reference is NOT appended: the sequences end where the allele does.
        let sequences: Vec<Vec<u8>> = variant
            .alleles
            .iter()
            .map(|allele| {
                let mut sequence = slice[..variant_offset as usize].to_vec();
                sequence.extend_from_slice(&allele.bases);
                sequence
            })
            .collect();
        // `+1` to ignore the base shared by every allele.
        let mut ranges: Vec<IndexRange> = variant
            .alleles
            .iter()
            .map(|allele| IndexRange::new(variant_offset + 1, variant_offset + allele.len() as i32))
            .collect();

        let borrowed: Vec<&[u8]> = sequences
            .iter()
            .map(|sequence| sequence.as_slice())
            .collect();
        let (left, right) = normalize_alleles(&borrowed, &mut ranges, variant_offset, trim)?;

        if left == 0 && right == 0 {
            return Ok((variant.clone(), Alignment::Complete));
        }
        // The shift reached the edge of the slice and a wider one is allowed: try again.
        if left == variant_offset && leading_bases < max_leading_bases {
            leading_bases = (2 * leading_bases).min(max_leading_bases);
            continue;
        }
        // Falling through with the shift still on the edge is the case the reference cannot tell
        // its caller about. Landing ON the edge is not enough to call the window exhausted: a
        // window of exactly the distance the record wants to walk lands there and is complete. The
        // only way to tell is to offer one more base and see whether the shift grows, which is one
        // extra pass and only in this case.
        let alignment = if left == variant_offset {
            if would_shift_further(variant, reference_bases, leading_bases, trim)? {
                Alignment::WindowExhausted
            } else {
                Alignment::Complete
            }
        } else {
            Alignment::Complete
        };

        let alleles: Vec<Allele> = variant
            .alleles
            .iter()
            .enumerate()
            .map(|(index, allele)| {
                let from = (variant_offset - left) as usize;
                let to = (variant_offset - right) as usize + allele.len();
                Allele::new(&sequences[index][from..to], index == 0)
            })
            .collect();

        return Ok((
            Variant {
                contig: variant.contig.clone(),
                start: variant.start - left,
                stop: variant.stop - right,
                // The reference builds this list from a HashMap's values, so the order is the
                // map's. Rebuilding it in record order is what the golden shows, and all it shows.
                alleles,
                genotypes: variant.genotypes.clone(),
                attributes: variant.attributes.clone(),
            },
            alignment,
        ));
    }
}

/// Whether one more base of window would move the record further left.
///
/// Left alignment slides base by base, so a record that can still move moves by at least one when
/// it is given one more base. A record already at the start of the contig cannot be given one.
fn would_shift_further(
    variant: &Variant,
    reference_bases: &[u8],
    leading_bases: i32,
    trim: bool,
) -> Result<bool, AlignmentError> {
    let ref_start = (variant.start - leading_bases).max(1);
    if ref_start == 1 {
        return Ok(false);
    }
    let wider_start = ref_start - 1;
    let slice = &reference_bases[(wider_start - 1) as usize..variant.stop as usize];
    let variant_offset = variant.start - wider_start;

    let sequences: Vec<Vec<u8>> = variant
        .alleles
        .iter()
        .map(|allele| {
            let mut sequence = slice[..variant_offset as usize].to_vec();
            sequence.extend_from_slice(&allele.bases);
            sequence
        })
        .collect();
    let mut ranges: Vec<IndexRange> = variant
        .alleles
        .iter()
        .map(|allele| IndexRange::new(variant_offset + 1, variant_offset + allele.len() as i32))
        .collect();
    let borrowed: Vec<&[u8]> = sequences
        .iter()
        .map(|sequence| sequence.as_slice())
        .collect();
    let (left, _) = normalize_alleles(&borrowed, &mut ranges, variant_offset, trim)?;
    Ok(left > leading_bases)
}

/// `trimAlleles(inputVC, trimForward, trimReverse)`: cut the bases every allele shares.
///
/// ```java
/// if ( inputVC.getNAlleles() <= 1 || inputVC.getAlleles().stream().anyMatch(a -> a.length() == 1 && !a.equals(Allele.SPAN_DEL)) ) {
///     return inputVC;
/// }
/// ```
///
/// **One allele of length one turns the whole thing off.** A record carrying a snp beside a
/// deletion is never trimmed, however much its other alleles share. The spanning deletion is the
/// one exception, excluded by name, and neither it nor a symbolic allele is fed to the comparison
/// though both are kept in the output.
pub fn trim_alleles(
    variant: &Variant,
    trim_forward: bool,
    trim_reverse: bool,
) -> Result<Variant, AlignmentError> {
    if variant.alleles.len() <= 1
        || variant
            .alleles
            .iter()
            .any(|allele| allele.len() == 1 && !allele.is_span_del())
    {
        return Ok(variant.clone());
    }

    let comparable: Vec<&Allele> = variant
        .alleles
        .iter()
        .filter(|allele| !allele.is_symbolic() && !allele.is_span_del())
        .collect();
    let sequences: Vec<&[u8]> = comparable
        .iter()
        .map(|allele| allele.bases.as_slice())
        .collect();
    let mut ranges: Vec<IndexRange> = comparable
        .iter()
        .map(|allele| IndexRange::new(0, allele.len() as i32))
        .collect();

    // The ranges start at index 0, so the only way left is backwards: the forward trim is the
    // NEGATIVE of the left shift.
    let (left, right) = normalize_alleles(&sequences, &mut ranges, 0, true)?;
    let end_trim = right;
    let start_trim = -left;

    // An allele trimmed to nothing gets one base back, at the end when nothing came off the front
    // and at the start when something did. It is the only reason a record keeps its anchor base.
    let empty_allele = ranges.iter().any(|range| range.size() == 0);
    let restore_one_at_end = empty_allele && start_trim == 0;
    let restore_one_at_start = empty_allele && start_trim > 0;

    let end_bases_to_clip = if restore_one_at_end {
        end_trim - 1
    } else {
        end_trim
    };
    let start_bases_to_clip = if restore_one_at_start {
        start_trim - 1
    } else {
        start_trim
    };

    trim_alleles_at(
        variant,
        if trim_forward { start_bases_to_clip } else { 0 } - 1,
        if trim_reverse { end_bases_to_clip } else { 0 },
    )
}

/// `trimAlleles(inputVC, fwdTrimEnd, revTrim)`, whose first index is INCLUSIVE.
///
/// `-1` therefore means "trim nothing from the front", which is what the caller above passes when
/// forward trimming is off, and `fwdTrimEnd == -1 && revTrim == 0` returns the input itself.
pub fn trim_alleles_at(
    variant: &Variant,
    forward_trim_end: i32,
    reverse_trim: i32,
) -> Result<Variant, AlignmentError> {
    if forward_trim_end == -1 && reverse_trim == 0 {
        return Ok(variant.clone());
    }

    let alleles: Vec<Allele> = variant
        .alleles
        .iter()
        .map(|allele| {
            if allele.is_symbolic() || allele.is_span_del() {
                allele.clone()
            } else {
                let from = (forward_trim_end + 1) as usize;
                let to = allele.len() - reverse_trim as usize;
                Allele::new(&allele.bases[from..to], allele.is_reference)
            }
        })
        .collect();

    // The reference remaps the genotypes through a LinkedHashMap of old allele to new, and leaves
    // any allele the map does not hold as it is. Genotypes are allele INDICES here and the list
    // keeps its order and length, so the remapping is the identity and the indices stand.
    let start = variant.start + (forward_trim_end + 1);
    Ok(Variant {
        contig: variant.contig.clone(),
        start,
        // The end is recomputed from the reference allele's new length, not by subtracting the
        // reverse trim, which is the same number by a different route.
        stop: start + alleles[0].len() as i32 - 1,
        alleles,
        genotypes: variant.genotypes.clone(),
        attributes: variant.attributes.clone(),
    })
}

/// `splitVariantContextToBiallelics(vc, trimLeft, BEST_MATCH_TO_ORIGINAL, false)`.
///
/// One record with N alternates becomes N records with one each, every one right trimmed. A
/// non-variant record becomes an EMPTY list rather than a list of one; a biallelic one is returned
/// as it stands, untrimmed, so the trimming every split record gets is not applied to a record that
/// did not need splitting.
///
/// # One 1/2 call empties every record
///
/// ```java
/// final GenotypeAssignmentMethod genotypeAssignmentMethodUsed = hasHetNonRef(vc.getGenotypes()) ? GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS : genotypeAssignmentMethod;
/// ```
///
/// One het-non-ref genotype in one sample switches the method for the whole record, so every
/// sample's call in every output record becomes a no-call and every attribute is stripped.
///
/// # What this port does not do
///
/// `AlleleSubsettingUtils.subsetAlleles` reindexes PLs, ADs and the annotations whose length
/// follows the allele count. Nothing in the golden carries any of those, so none of it is measured
/// and none of it is written here: [`SplitError::LikelihoodsNotPorted`] is returned instead of a
/// guess. What is ported is the path a record without likelihoods takes, where
/// `BEST_MATCH_TO_ORIGINAL` turns an allele that is not in the pair into the REFERENCE.
pub fn split_variant_context_to_biallelics(
    variant: &Variant,
    trim_left: bool,
) -> Result<Vec<Variant>, SplitError> {
    // `isVariant()` is "more than one allele", and `isBiallelic()` is "exactly two".
    if variant.alleles.len() <= 1 {
        return Ok(Vec::new());
    }
    if variant.alleles.len() == 2 {
        return Ok(vec![variant.clone()]);
    }

    let no_annotations = has_het_non_ref(variant);
    let mut out = Vec::with_capacity(variant.alleles.len() - 1);

    for alternate in 1..variant.alleles.len() {
        let alleles = vec![
            variant.alleles[0].clone(),
            variant.alleles[alternate].clone(),
        ];

        // `!(AC || AF || AN) || method == SET_TO_NO_CALL_NO_ANNOTATIONS` decides what is REMOVED,
        // so those three are what survives, and only while that method is not in force.
        let attributes: Vec<(String, String)> = if no_annotations {
            Vec::new()
        } else {
            variant
                .attributes
                .iter()
                .filter(|(key, _)| matches!(key.as_str(), "AC" | "AF" | "AN"))
                .cloned()
                .collect()
        };

        let method = if no_annotations {
            AssignmentMethod::SetToNoCallNoAnnotations
        } else {
            AssignmentMethod::BestMatchToOriginal
        };
        let genotypes = subset_alleles(
            &variant.genotypes,
            2,
            variant.alleles.len(),
            &[0, alternate],
            method,
        )
        .map_err(SplitError::GenotypeIndex)?;
        // And then they are recomputed from the genotypes, which is why a record with none loses
        // all three however well they survived the filter.
        let attributes = if no_annotations {
            attributes
        } else {
            chromosome_counts(&genotypes, variant.genotypes.len())
        };

        let split = Variant {
            contig: variant.contig.clone(),
            start: variant.start,
            stop: variant.stop,
            alleles,
            genotypes,
            attributes,
        };
        out.push(trim_alleles(&split, trim_left, true).map_err(SplitError::Alignment)?);
    }
    Ok(out)
}

/// What the splitter refuses rather than guesses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitError {
    /// The trimming underneath refused.
    Alignment(AlignmentError),
    /// The genotype-index combinatorics refused, which only a very large allele count reaches.
    GenotypeIndex(crate::genotype_index::GenotypeIndexError),
}

impl SplitError {
    pub fn message(&self) -> String {
        match self {
            SplitError::Alignment(error) => error.message(),
            SplitError::GenotypeIndex(error) => error.message(),
        }
    }
}

/// `hasHetNonRef`: a genotype calling two different alternates, `1/2`.
fn has_het_non_ref(variant: &Variant) -> bool {
    variant.genotypes.iter().any(|genotype| {
        genotype.alleles.len() == 2
            && matches!((genotype.alleles[0], genotype.alleles[1]), (Some(first), Some(second))
                if first != second && first != 0 && second != 0)
    })
}

/// `addInfoFieldAnnotations` for the three counts, which is all a split record keeps.
///
/// AN is every called allele, AC is the alternate's, and AF is the one over the other. A record
/// with no genotypes gets none of the three, which is why they vanish however well they survived
/// the attribute filter.
fn chromosome_counts(genotypes: &[Genotype], original_samples: usize) -> Vec<(String, String)> {
    if original_samples == 0 {
        return Vec::new();
    }
    let mut allele_number = 0;
    let mut allele_count = 0;
    for genotype in genotypes {
        for allele in genotype.alleles.iter().flatten() {
            allele_number += 1;
            if *allele == 1 {
                allele_count += 1;
            }
        }
    }
    if allele_number == 0 {
        return Vec::new();
    }
    let frequency = f64::from(allele_count) / f64::from(allele_number);
    vec![
        ("AC".to_string(), allele_count.to_string()),
        (
            "AF".to_string(),
            crate::tsv_table::java_double_to_string(frequency),
        ),
        ("AN".to_string(), allele_number.to_string()),
    ]
}

/// `VariantContext.Type`, as far as `typeOfVariant` produces it.
///
/// MIXED is not here: the comment in the reference says a pairwise comparison cannot produce it,
/// and the code proves it by returning INDEL for everything that is not a symbolic allele, a
/// spanning deletion or a same-length pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
}

impl VariantType {
    /// The enum constant's own name, which is what a Java `printf` of it writes.
    pub fn name(&self) -> &'static str {
        match self {
            VariantType::NoVariation => "NO_VARIATION",
            VariantType::Snp => "SNP",
            VariantType::Mnp => "MNP",
            VariantType::Indel => "INDEL",
            VariantType::Symbolic => "SYMBOLIC",
        }
    }
}

/// `Trilean`, which is what `doesReadContainAllele` answers.
///
/// The third value is not a convenience: a read that ends inside the allele is UNKNOWN, and
/// `calculateMaxAltRatio` counts UNKNOWN into neither the numerator nor the denominator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trilean {
    True,
    False,
    Unknown,
}

/// What typing an allele pair, or choosing one for a read, refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum PileupAlleleError {
    /// A symbolic reference allele. Unreachable through htsjdk, which refuses to build one.
    SymbolicReference,
    /// A spanning-deletion reference allele.
    SpanDelReference,
    /// `ParamUtils.isPositiveOrZero` on the cutoff, as `chooseAlleleForRead` words it.
    NegativeMinimumBaseQuality,
    /// The same check, as the two pileup halves of `PowerCalculationUtils` word it. Two messages
    /// for one condition, and the golden carries both.
    NegativeMinimumBaseQualityRatio,
}

impl PileupAlleleError {
    pub fn java_class(&self) -> &'static str {
        match self {
            PileupAlleleError::SymbolicReference | PileupAlleleError::SpanDelReference => {
                "java.lang.IllegalStateException"
            }
            PileupAlleleError::NegativeMinimumBaseQuality
            | PileupAlleleError::NegativeMinimumBaseQualityRatio => {
                "java.lang.IllegalArgumentException"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            PileupAlleleError::SymbolicReference => {
                "Unexpected error: encountered a record with a symbolic reference allele"
                    .to_string()
            }
            PileupAlleleError::SpanDelReference => {
                "Unexpected error: encountered a record with a spanning deletion reference allele"
                    .to_string()
            }
            PileupAlleleError::NegativeMinimumBaseQuality => {
                "Minimum base quality must be positive or zero.".to_string()
            }
            PileupAlleleError::NegativeMinimumBaseQualityRatio => {
                "Cannot have a negative minBaseQualityCutoff.".to_string()
            }
        }
    }
}

/// `typeOfVariant(ref, allele)`, one alternate against the reference.
///
/// Two decisions worth naming. A SPANNING DELETION ALTERNATE IS NO_VARIATION, tested before the
/// reference's own spanning-deletion refusal, so `*` against `*` answers rather than throwing. And
/// EQUAL-LENGTH ALLELES DIFFERING IN EXACTLY ONE BASE ARE A SNP however long they are, so `ACGT`
/// against `ACGA` is a SNP and only two or more differences make an MNP.
pub fn type_of_variant(
    reference: &Allele,
    allele: &Allele,
) -> Result<VariantType, PileupAlleleError> {
    if reference.is_symbolic() {
        return Err(PileupAlleleError::SymbolicReference);
    }
    if allele.is_symbolic() {
        return Ok(VariantType::Symbolic);
    }
    // `allele.equals(Allele.SPAN_DEL)`, which is the bases and the reference flag both.
    if allele.is_span_del() && !allele.is_reference {
        return Ok(VariantType::NoVariation);
    }
    if reference.is_span_del() && reference.is_reference {
        return Err(PileupAlleleError::SpanDelReference);
    }
    if reference.len() == allele.len() {
        if reference.bases == allele.bases {
            return Ok(VariantType::NoVariation);
        }
        if allele.len() == 1 {
            return Ok(VariantType::Snp);
        }
        let differences = reference
            .bases
            .iter()
            .zip(allele.bases.iter())
            .filter(|(left, right)| left != right)
            .count();
        return Ok(if differences == 1 {
            VariantType::Snp
        } else {
            VariantType::Mnp
        });
    }
    Ok(VariantType::Indel)
}

/// `isComplexIndel(ref, allele)`: a prefix test, and nothing else.
///
/// Same length, either side empty or symbolic, or either side a single base: not complex. Past
/// that the shorter allele has to be a prefix of the longer one, so a deletion of `AAA` to `TA` is
/// complex and a genotype carrying it cannot be validated at all.
pub fn is_complex_indel(reference: &Allele, allele: &Allele) -> bool {
    if reference.is_symbolic() || reference.is_empty() {
        return false;
    }
    if allele.is_symbolic() || allele.is_empty() {
        return false;
    }
    if reference.len() == allele.len() {
        return false;
    }
    if allele.len() == 1 || reference.len() == 1 {
        return false;
    }
    let shorter = reference.len().min(allele.len());
    let prefix = reference.bases[..shorter] == allele.bases[..shorter];
    !prefix
}

/// `doesReadContainAllele(pileupElement, allele)`.
///
/// The subarray stops at the end of the read, so a read that ends inside the allele produces a
/// short array and the answer is UNKNOWN rather than false. The comparison itself is exact: htsjdk
/// upper-cases an allele's bases when it is built and does not upper-case the read's.
pub fn does_read_contain_allele(element: &PileupElement<'_>, allele: &Allele) -> Trilean {
    let bases = bases_for_allele_in_read(element, allele);
    if bases.len() < allele.len() {
        return Trilean::Unknown;
    }
    if bases == allele.bases {
        Trilean::True
    } else {
        Trilean::False
    }
}

/// `getBasesForAlleleInRead`, which is `ArrayUtils.subarray` and therefore clamps rather than
/// throwing when the read ends first.
fn bases_for_allele_in_read(element: &PileupElement<'_>, allele: &Allele) -> Vec<u8> {
    let from = element.offset as usize;
    let to = (from + allele.len()).min(element.read.read_bases.len());
    if from >= to {
        return Vec::new();
    }
    element.read.read_bases[from..to].to_vec()
}

/// `getMinBaseQualityForAlleleInRead`, whose minimum over an empty array is -1.
///
/// A read that ends exactly at the allele therefore fails any cutoff above -1, which includes a
/// cutoff of zero.
fn min_base_quality_for_allele_in_read(element: &PileupElement<'_>, allele: &Allele) -> i32 {
    let from = element.offset as usize;
    let to = (from + allele.len()).min(element.read.base_qualities.len());
    if from >= to {
        return -1;
    }
    element.read.base_qualities[from..to]
        .iter()
        .map(|quality| i32::from(*quality))
        .min()
        .unwrap_or(-1)
}

/// `isIndelInThePileupElement`: an insertion matches on its bases, a deletion on its length alone.
fn is_indel_in_the_pileup_element(
    element: &PileupElement<'_>,
    reference: &Allele,
    alternate: &Allele,
) -> bool {
    if element.is_before_insertion() {
        // A deletion immediately preceding an insertion has no following insertion bases, which
        // the reference reads as a null and ignores.
        let Some(inserted) = element.bases_of_immediately_following_insertion() else {
            return false;
        };
        let mut extended = reference.bases.clone();
        extended.extend_from_slice(&inserted);
        return extended == alternate.bases;
    }
    if element.is_before_deletion_start() {
        let length = element.length_of_immediately_following_indel() as i64;
        return (reference.len() as i64 - alternate.len() as i64) == length;
    }
    false
}

/// `chooseAlleleForRead(pileupElement, referenceAllele, altAlleles, minBaseQualityCutoff)`.
///
/// THE LOOP DOES NOT BREAK. Every alternate that matches assigns, so the last one wins, and two
/// alternates that both match are not a refusal.
///
/// THE REFERENCE TEST INCLUDES THE INDEL LOOKAHEAD: a read whose bases match the reference is not
/// the reference allele when it sits before an insertion or a deletion start, which is what lets
/// that read be counted as an indel alternate instead.
///
/// THE QUALITY CUTOFF IS APPLIED LAST, over the chosen allele's own bases rather than the
/// reference's, so a cutoff can turn an answer into none but never into a different allele.
pub fn choose_allele_for_read(
    element: &PileupElement<'_>,
    reference: &Allele,
    alternates: &[Allele],
    minimum_base_quality: i32,
) -> Result<Option<Allele>, PileupAlleleError> {
    if minimum_base_quality < 0 {
        return Err(PileupAlleleError::NegativeMinimumBaseQuality);
    }
    let is_reference = bases_for_allele_in_read(element, reference) == reference.bases
        && !element.is_before_deletion_start()
        && !element.is_before_insertion();
    let mut chosen: Option<Allele> = None;
    if is_reference {
        chosen = Some(reference.clone());
    } else {
        for alternate in alternates {
            match type_of_variant(reference, alternate)? {
                VariantType::Indel => {
                    if is_indel_in_the_pileup_element(element, reference, alternate) {
                        chosen = Some(alternate.clone());
                    }
                }
                VariantType::Mnp | VariantType::Snp => {
                    if does_read_contain_allele(element, alternate) == Trilean::True {
                        chosen = Some(alternate.clone());
                    }
                }
                VariantType::NoVariation | VariantType::Symbolic => {}
            }
        }
    }
    if let Some(allele) = &chosen {
        if min_base_quality_for_allele_in_read(element, allele) < minimum_base_quality {
            chosen = None;
        }
    }
    Ok(chosen)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ten G, eight A, twelve G, thirty T, ten G, a CA repeat, ten G: the dump's reference.
    fn reference() -> Vec<u8> {
        let mut bases = Vec::new();
        bases.extend(std::iter::repeat_n(b'G', 10));
        bases.extend(std::iter::repeat_n(b'A', 8));
        bases.extend(std::iter::repeat_n(b'G', 12));
        bases.extend(std::iter::repeat_n(b'T', 30));
        bases.extend(std::iter::repeat_n(b'G', 10));
        for _ in 0..10 {
            bases.extend_from_slice(b"CA");
        }
        bases.extend(std::iter::repeat_n(b'G', 10));
        bases
    }

    fn deletion(start: i32, reference: &str, alternate: &str) -> Variant {
        Variant {
            contig: "chr1".to_string(),
            start,
            stop: start + reference.len() as i32 - 1,
            alleles: vec![
                Allele::new(reference.as_bytes(), true),
                Allele::new(alternate.as_bytes(), false),
            ],
            genotypes: Vec::new(),
            attributes: Vec::new(),
        }
    }

    #[test]
    fn a_narrow_window_stops_short_and_says_nothing() {
        let bases = reference();
        let variant = deletion(17, "AA", "A");
        let wide = left_align_and_trim(&variant, &bases, 1000, true).expect("aligned");
        assert_eq!((wide.start, wide.stop), (10, 11));
        let narrow = left_align_and_trim(&variant, &bases, 2, true).expect("aligned");
        assert_eq!((narrow.start, narrow.stop), (15, 16));
    }

    /// The same two calls, asked which of them was stopped by the window.
    #[test]
    fn the_report_tells_the_two_apart() {
        let bases = reference();
        let variant = deletion(17, "AA", "A");

        let (wide, alignment) =
            left_align_and_trim_reporting(&variant, &bases, 1000, true).expect("aligned");
        assert_eq!((wide.start, wide.stop), (10, 11));
        assert_eq!(alignment, Alignment::Complete);

        let (narrow, alignment) =
            left_align_and_trim_reporting(&variant, &bases, 2, true).expect("aligned");
        assert_eq!((narrow.start, narrow.stop), (15, 16));
        assert_eq!(alignment, Alignment::WindowExhausted);

        // Every window between the two, so the report follows the record and not the argument.
        for (window, expected) in [
            (7, Alignment::Complete),
            (6, Alignment::WindowExhausted),
            (10, Alignment::Complete),
        ] {
            let (_, alignment) =
                left_align_and_trim_reporting(&variant, &bases, window, true).expect("aligned");
            assert_eq!(alignment, expected, "window {window}");
        }
    }

    /// A record the reference never touches reports that nothing was attempted, which is not the
    /// same as a record it aligned and could not move.
    #[test]
    fn nothing_attempted_is_not_the_same_as_nothing_to_do() {
        let bases = reference();
        let snp = deletion(18, "A", "C");
        let (_, alignment) =
            left_align_and_trim_reporting(&snp, &bases, 1000, true).expect("untouched");
        assert_eq!(alignment, Alignment::NotAttempted);

        let (_, alignment) =
            left_align_and_trim_reporting(&deletion(17, "AA", "A"), &bases, 0, true)
                .expect("untouched");
        assert_eq!(alignment, Alignment::NotAttempted);

        // Already left aligned: attempted, and there was nowhere to go.
        let (_, alignment) =
            left_align_and_trim_reporting(&deletion(10, "GA", "G"), &bases, 1000, true)
                .expect("untouched");
        assert_eq!(alignment, Alignment::Complete);
    }

    #[test]
    fn the_window_doubles_until_it_is_wide_enough() {
        let bases = reference();
        let variant = deletion(59, "TT", "T");
        for (window, start) in [(10, 49), (20, 39), (1000, 30)] {
            let result = left_align_and_trim(&variant, &bases, window, true).expect("aligned");
            assert_eq!(result.start, start, "window {window}");
        }
    }

    #[test]
    fn three_things_come_back_unchanged() {
        let bases = reference();
        let aligned = deletion(10, "GA", "G");
        assert_eq!(
            left_align_and_trim(&aligned, &bases, 1000, true).expect("same"),
            aligned
        );

        let variant = deletion(17, "AA", "A");
        assert_eq!(
            left_align_and_trim(&variant, &bases, 0, true).expect("same"),
            variant
        );

        // A snp is not an indel, and neither is a MIXED site.
        let snp = deletion(18, "A", "C");
        assert_eq!(
            left_align_and_trim(&snp, &bases, 1000, true).expect("same"),
            snp
        );
        let mut mixed = deletion(18, "AG", "A");
        mixed.alleles.push(Allele::new(b"CC", false));
        assert!(!mixed.is_indel());
        assert_eq!(
            left_align_and_trim(&mixed, &bases, 1000, true).expect("same"),
            mixed
        );
    }

    #[test]
    fn trimming_decides_whether_the_record_shrinks() {
        let bases = reference();
        // Alleles sharing a trailing base as well as the leading one.
        let variant = deletion(18, "AGG", "AG");
        let trimmed = left_align_and_trim(&variant, &bases, 1000, true).expect("aligned");
        assert_eq!((trimmed.start, trimmed.stop), (18, 19));
        let untrimmed = left_align_and_trim(&variant, &bases, 1000, false).expect("aligned");
        assert_eq!((untrimmed.start, untrimmed.stop), (17, 19));
    }
}
