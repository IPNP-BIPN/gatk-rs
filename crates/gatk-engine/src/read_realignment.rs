//! Ported from `org.broadinstitute.hellbender.utils.read.AlignmentUtils` (`createReadAlignedToRef`,
//! `applyCigarToCigar`, `readStartOnReferenceHaplotype`, `appendClippedElementsFromCigarToCigar`)
//! and `AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype` (GATK 4.7.0.0;
//! `AlignmentUtils.java` and `AssemblyBasedCallerUtils.java` are unchanged since 4.6.2.0): every
//! read realigned to the reference through the haplotype that explains it best.
//!
//! A read is aligned to its haplotype by Smith-Waterman (soft-clipping overhangs, its own soft
//! clips hard-clipped first), the haplotype's CIGAR against the reference (padded by 1000 matches
//! on the right) is cut from where the read starts, the read-to-haplotype CIGAR is projected
//! through it, the indels are left-aligned against the reference haplotype, and the read's
//! original clips are put back at both ends. A read the aligner cannot place keeps its alignment.
//!
//! The best haplotype breaks likelihood ties by `HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`: the
//! reference first, then the haplotype with the fewest CIGAR elements.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, RealignmentError>` is "a value or the exception"; `?` after a call hands the error
//!   to the caller, Java's `throws` at one call site. The `impl From<..>` blocks below are what
//!   let `?` turn a lower layer's error into this one without a wrapper at every call.
//! - `.ok_or_else(|| ..)` turns an `Option` (a value or nothing, Java's nullable) into a
//!   `Result`, so a missing value becomes the error the closure builds.
//! - `&dyn SmithWatermanAligner` is any aligner behind an interface, as Java's
//!   `SmithWatermanAligner` parameter is.
//! - `impl IntoIterator<Item = CigarElement>` accepts anything that yields cigar elements, Java's
//!   `Iterable<CigarElement>`.
//! - `x as i32` converts between Rust's unsigned `usize`/`u32` lengths and Java's `int`.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::alignment_utils::{left_align_indels, trim_cigar_by_bases, AlignmentError};
use crate::allele_likelihoods::AlleleLikelihoods;
use crate::cigar_builder::{CigarBuilder, CigarError};
use crate::clipping::{hard_clip_soft_clipped_bases, ClipError};
use crate::haplotype::Haplotype;
use crate::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};

/// `AlignmentUtils.HAPLOTYPE_TAG`.
pub const HAPLOTYPE_TAG: &[u8; 2] = b"HC";

/// What the realignment refuses: the layer it came from, or the reference's own refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RealignmentError {
    Alignment(AlignmentError),
    Cigar(CigarError),
    Clip(ClipError),
    /// `IllegalArgumentException` / `IllegalStateException` / `GATKException`.
    Illegal(String),
}

impl From<AlignmentError> for RealignmentError {
    fn from(e: AlignmentError) -> Self {
        RealignmentError::Alignment(e)
    }
}

impl From<CigarError> for RealignmentError {
    fn from(e: CigarError) -> Self {
        RealignmentError::Cigar(e)
    }
}

/// `new CigarBuilder().addAll(elements).make()`. The no-argument constructor is
/// `this(true)` (remove deletions at the ends) and `make()` is `make(false)` (refuse an empty
/// cigar), which are the two literals here.
fn build(elements: impl IntoIterator<Item = CigarElement>) -> Result<Cigar, CigarError> {
    let mut builder = CigarBuilder::new(true);
    for element in elements {
        builder.add(element)?;
    }
    builder.make(false)
}

/// `HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`: one for the reference, less one for every CIGAR
/// element past the first.
///
/// ```java
/// final int referenceTerm = (h.isReference() ? 1 : 0);
/// final int cigarTerm = cigar == null ? 0 : (1 - cigar.numCigarElements());
/// return (double) referenceTerm + cigarTerm;
/// ```
///
/// Summed as integers and then widened here, summed in `double` upstream; for integers this small
/// both are exact, so the two agree.
pub fn haplotype_alignment_tiebreaking_priority(haplotype: &Haplotype) -> f64 {
    let reference_term = if haplotype.is_reference() { 1 } else { 0 };
    let cigar_term = haplotype.cigar().map_or(0, |c| 1 - c.elements.len() as i64);
    (reference_term + cigar_term) as f64
}

/// `realignReadsToTheirBestHaplotype` followed by `changeEvidence`: each sample's reads replaced
/// by their realignment, in place.
///
/// The reference returns a `Map<GATKRead, GATKRead>` from original to realigned, and the caller
/// hands it to `changeEvidence`. Here the pairs are collected first and applied after the loop,
/// as that map is: `original` is read out of `likelihoods` while the loop runs, and Rust refuses
/// to change `likelihoods` while that borrow is live. The likelihoods themselves are untouched;
/// only the evidence is replaced.
pub fn realign_reads_to_their_best_haplotype(
    likelihoods: &mut AlleleLikelihoods<BamRecord, Haplotype>,
    ref_haplotype: &Haplotype,
    padded_reference_start: i32,
    aligner: &dyn SmithWatermanAligner,
    parameters: &SwParameters,
) -> Result<(), RealignmentError> {
    let priorities: Vec<f64> = (0..likelihoods.number_of_alleles())
        .map(|a| {
            haplotype_alignment_tiebreaking_priority(likelihoods.get_allele(a).expect("an allele"))
        })
        .collect();
    // The reference passes the priority function; here it is evaluated once per allele, in allele
    // order, and the values handed over (see `best_alleles_breaking_ties`).
    let best = likelihoods.best_alleles_breaking_ties(Some(&priorities));
    let mut replacements: Vec<(usize, usize, BamRecord)> = Vec::with_capacity(best.len());
    for allele in best {
        let sample = likelihoods
            .index_of_sample(&allele.sample)
            .expect("a sample of the matrix");
        let original =
            &likelihoods.sample_evidence(sample).expect("a sample")[allele.evidence_index];
        let haplotype = allele
            .allele
            .as_ref()
            .ok_or_else(|| RealignmentError::Illegal("no best haplotype for a read".to_string()))?;
        let realigned = create_read_aligned_to_ref(
            original,
            haplotype,
            ref_haplotype,
            padded_reference_start,
            allele.is_informative(),
            aligner,
            parameters,
        )?;
        replacements.push((sample, allele.evidence_index, realigned));
    }
    for (sample, index, read) in replacements {
        likelihoods.replace_evidence(sample, index, read);
    }
    Ok(())
}

/// `createReadAlignedToRef(originalRead, haplotype, refHaplotype, referenceStart, isInformative,
/// aligner, readToHaplotypeSWParameters)`.
///
/// The steps, each the reference's in its order:
///
/// 1. hard-clip the read's soft clips and align what is left to the haplotype, soft-clipping
///    overhangs; an offset of `-1` means the aligner failed, and the read is returned unchanged
///    ("sw can fail (reasons not clear) so if it happens just don't realign the read");
/// 2. tag an informative read with `HC`, the haplotype's `hashCode()`;
/// 3. find where the read starts on the reference, through the haplotype's cigar padded with
///    `1000M` on the right;
/// 4. cut that padded cigar at the read's start, project the read-to-haplotype cigar through it,
///    and left-align the indels against the reference haplotype;
/// 5. move the start right by any leading deletion the left-alignment removed, and put the
///    original clips back at both ends.
///
/// The final length check is the reference's `GATKException`; its message here stops after the
/// first clause, where the reference goes on to print the read and the intermediate cigars.
pub fn create_read_aligned_to_ref(
    original: &BamRecord,
    haplotype: &Haplotype,
    ref_haplotype: &Haplotype,
    reference_start: i32,
    is_informative: bool,
    aligner: &dyn SmithWatermanAligner,
    parameters: &SwParameters,
) -> Result<BamRecord, RealignmentError> {
    if reference_start < 1 {
        return Err(RealignmentError::Illegal(format!(
            "reference start much be >= 1 but got {reference_start}"
        )));
    }
    let haplotype_cigar = haplotype
        .cigar()
        .ok_or_else(|| RealignmentError::Illegal("haplotype without a cigar".to_string()))?;
    // `ReadClipper.hardClipSoftClippedBases(originalRead)`: no header, no extra bases.
    let minus_soft_clips =
        hard_clip_soft_clipped_bases(original, None, 0).map_err(RealignmentError::Clip)?;
    let soft_clipped_bases = original.read_bases.len() - minus_soft_clips.read_bases.len();
    let alignment = aligner.align(
        &haplotype.bases(),
        &minus_soft_clips.read_bases,
        parameters,
        SwOverhangStrategy::SoftClip,
    )?;
    if alignment.alignment_offset == -1 {
        // A read that cannot be aligned keeps its original alignment.
        return Ok(original.clone());
    }
    let sw_cigar = build(alignment.cigar.elements.iter().copied())?;
    // `originalRead.copy()`: the tag and the new position go on a copy, never on the input.
    let mut copied = original.clone();
    // "only informative reads are given the haplotype tag to enhance visualization". The value is
    // `Haplotype.hashCode()`, which is `Arrays.hashCode(getBases())`.
    if is_informative {
        copied.tags.insert(
            Tag::new(HAPLOTYPE_TAG),
            TagValue::Int(i64::from(haplotype.java_hash_code())),
        );
    }
    // `getConsolidatedPaddedCigar(1000)`, which is `new CigarBuilder().addAll(getCigar())` then
    // `.add(new CigarElement(padSize, CigarOperator.M)).make()`.
    // `std::iter::once` is a one-element iterator, chained after the haplotype's own elements.
    let padded_haplotype_cigar = build(haplotype_cigar.elements.iter().copied().chain(
        std::iter::once(CigarElement {
            length: 1000,
            op: Op::M,
        }),
    ))?;
    let read_start_on_ref_haplotype =
        read_start_on_reference_haplotype(&padded_haplotype_cigar, alignment.alignment_offset)?;
    let read_start_on_reference =
        reference_start + haplotype.alignment_start_hap_wrt_ref() + read_start_on_ref_haplotype;
    // The haplotype-to-reference cigar with everything before the read's start removed, and
    // everything after the read's end kept.
    let haplotype_to_ref = trim_cigar_by_bases(
        &padded_haplotype_cigar,
        alignment.alignment_offset,
        padded_haplotype_cigar.read_length() as i32 - 1,
    )?
    .cigar;
    let read_to_ref = apply_cigar_to_cigar(&sw_cigar, &haplotype_to_ref)?;
    let left_aligned = left_align_indels(
        &read_to_ref,
        &ref_haplotype.bases(),
        &minus_soft_clips.read_bases,
        read_start_on_ref_haplotype,
    )?;
    // "it's possible that left-alignment shifted a deletion to the beginning of a read and removed
    // it, shifting the first aligned base to the right".
    copied.alignment_start =
        read_start_on_reference + left_aligned.leading_deletion_bases_removed as i32;
    let new_cigar =
        append_clipped_elements_from_cigar_to_cigar(&left_aligned.cigar, &original.cigar);
    copied.cigar = new_cigar;
    if left_aligned.cigar.read_length() as usize + soft_clipped_bases != copied.read_bases.len() {
        return Err(RealignmentError::Illegal(format!(
            "Cigar {} with read length {} != read length {}",
            left_aligned.cigar.to_text(),
            left_aligned.cigar.read_length(),
            copied.read_bases.len()
        )));
    }
    Ok(copied)
}

/// `readStartOnReferenceHaplotype(haplotypeVsRefCigar, readStartOnHaplotype)`: how many reference
/// bases precede the read's first base, given its offset on the haplotype.
///
/// Walk the haplotype-versus-reference cigar until enough haplotype bases are consumed to reach
/// the read. If the element that got there also consumes reference, the overshoot is taken back;
/// if it does not (an insertion), the read starts where the insertion is anchored:
///
/// ```java
/// final int excess = element.getOperator().consumesReferenceBases()
///         ? haplotypeBasesConsumed - readStartOnHaplotype : 0;
/// return refBasesConsumedBeforeReadStart - excess;
/// ```
///
/// (The first statement is one line in the reference.)
fn read_start_on_reference_haplotype(
    haplotype_vs_ref: &Cigar,
    read_start_on_haplotype: i32,
) -> Result<i32, RealignmentError> {
    if read_start_on_haplotype == 0 {
        return Ok(0);
    }
    let mut ref_consumed = 0i32;
    let mut haplotype_consumed = 0i32;
    for element in &haplotype_vs_ref.elements {
        if element.op.consumes_reference_bases() {
            ref_consumed += element.length as i32;
        }
        if element.op.consumes_read_bases() {
            haplotype_consumed += element.length as i32;
        }
        if haplotype_consumed >= read_start_on_haplotype {
            let excess = if element.op.consumes_reference_bases() {
                haplotype_consumed - read_start_on_haplotype
            } else {
                0
            };
            return Ok(ref_consumed - excess);
        }
    }
    Err(RealignmentError::Illegal(
        "Cigar doesn't reach the read start".to_string(),
    ))
}

/// The `CigarPairTransform` classes: `M` covers `M`, `=` and `X`; `I` covers `I` and `S`.
///
/// This is `getCigarSet`. Anything else (`N`, `H`, `P`) has no class, which makes
/// [`apply_cigar_to_cigar`] refuse, as the reference's `getTransformer` throws
/// `IllegalStateException` when no transform matches.
fn class(op: Op) -> Option<char> {
    match op {
        Op::M | Op::Eq | Op::X => Some('M'),
        Op::I | Op::S => Some('I'),
        Op::D => Some('D'),
        _ => None,
    }
}

/// `applyCigarToCigar(firstToSecond, secondToThird)`: the projection, one element at a time.
///
/// Given an alignment of 1 to 2 and of 2 to 3, the alignment of 1 to 3. Each step reads the
/// current operator on both sides, emits one base of the combined operator (or nothing), and
/// advances either side by one base of its element. The table is the reference's
/// `cigarPairTransformers`, in its order; the `match` below lists the same nine rows:
///
/// ```text
/// 1->2  2->3   1->3  advance12 advance23
///  M     M      M        1        1
///  M     I      I        1        1
///  M     D      D        0        1
///  D     M      D        1        1
///  D     D      D        0        1
///  D     I      -        1        1
///  I     M      I        1        0
///  I     D      I        1        0
///  I     I      I        1        0
/// ```
///
/// The loop stops as soon as either cigar runs out, as the reference's `while` does, so a longer
/// second cigar (the padded haplotype here) is simply not read past the read's end. The builder
/// merges the one-base elements back into runs.
pub fn apply_cigar_to_cigar(
    first_to_second: &Cigar,
    second_to_third: &Cigar,
) -> Result<Cigar, RealignmentError> {
    let mut builder = CigarBuilder::new(true);
    let (n12, n23) = (
        first_to_second.elements.len(),
        second_to_third.elements.len(),
    );
    // `cigar12I, cigar23I` (which element) and `elt12I, elt23I` (how far into it).
    let (mut c12, mut c23, mut e12, mut e23) = (0usize, 0usize, 0u32, 0u32);
    while c12 < n12 && c23 < n23 {
        let elt12 = first_to_second.elements[c12];
        let elt23 = second_to_third.elements[c23];
        // Matching on the pair is the reference's first-match search through the list; the rows
        // do not overlap, so order cannot change which one applies. `None` for `op13` is the
        // `D`/`I` row's "skip no ops".
        let (op13, advance12, advance23) = match (class(elt12.op), class(elt23.op)) {
            (Some('M'), Some('M')) => (Some(Op::M), 1, 1),
            (Some('M'), Some('I')) => (Some(Op::I), 1, 1),
            (Some('M'), Some('D')) => (Some(Op::D), 0, 1),
            (Some('D'), Some('M')) => (Some(Op::D), 1, 1),
            (Some('D'), Some('D')) => (Some(Op::D), 0, 1),
            (Some('D'), Some('I')) => (None, 1, 1),
            (Some('I'), Some('M')) => (Some(Op::I), 1, 0),
            (Some('I'), Some('D')) => (Some(Op::I), 1, 0),
            (Some('I'), Some('I')) => (Some(Op::I), 1, 0),
            _ => {
                return Err(RealignmentError::Illegal(format!(
                    "No transformer for operators {:?} and {:?}",
                    elt12.op, elt23.op
                )))
            }
        };
        if let Some(op) = op13 {
            builder.add(CigarElement { length: 1, op })?;
        }
        e12 += advance12;
        e23 += advance23;
        if e12 == elt12.length {
            c12 += 1;
            e12 = 0;
        }
        if e23 == elt23.length {
            c23 += 1;
            e23 = 0;
        }
    }
    Ok(builder.make(false)?)
}

/// `appendClippedElementsFromCigarToCigar`: the original's leading and trailing clips around the
/// new CIGAR.
///
/// Both walks stop when `first` meets `last`, so in a cigar of clips alone the last element is
/// never copied. The trailing clips are collected from the end inwards and then
/// reversed, which the reference does "to preserve the original soft/hardclip ordering in mixed
/// clipping cases where softclips precede hardclips". An empty original cigar fails in both:
/// upstream `getFirstCigarElement()` answers `null` and `.getOperator()` on it throws, and here
/// `elements.len() - 1` underflows and the indexing panics.
fn append_clipped_elements_from_cigar_to_cigar(new: &Cigar, original: &Cigar) -> Cigar {
    let elements = &original.elements;
    let mut first = 0usize;
    let mut last = elements.len() - 1;
    let is_clip = |op: Op| matches!(op, Op::S | Op::H);
    let mut out = Vec::new();
    while is_clip(elements[first].op) && first != last {
        out.push(elements[first]);
        first += 1;
    }
    out.extend(new.elements.iter().copied());
    let mut tail = Vec::new();
    while is_clip(elements[last].op) && first != last {
        tail.push(elements[last]);
        last -= 1;
    }
    tail.reverse();
    out.extend(tail);
    Cigar::new(out)
}
