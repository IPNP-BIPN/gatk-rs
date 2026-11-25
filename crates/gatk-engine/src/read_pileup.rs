//! Ported from `org.broadinstitute.hellbender.utils.pileup.ReadPileup` and the part of
//! `org.broadinstitute.hellbender.utils.BaseUtils` it reaches (GATK 4.7.0.0; `ReadPileup.java`
//! and `BaseUtils.java` are unchanged since 4.6.2.0).
//!
//! Everything one locus looks like to a tool. [`crate::alignment_state::AlignmentStateMachine`]
//! decides where each read stops and [`crate::pileup::PileupElement`] decides what each stop looks
//! like; this collects them and answers the questions a caller asks of a position.
//!
//! Four of those answers are decisions rather than aggregation:
//!
//!  * **`getBaseCounts` counts `*` as an `A`.** `BaseUtils.baseIndexMap` maps the wildcard
//!    character to `Base.A.ordinal()`, with a comment saying so, so a read carrying `*` in its
//!    sequence inflates the A count rather than being skipped like an `N`;
//!  * **it skips deletions but not every non-base.** A deletion is excluded by an explicit test,
//!    and anything else that maps to `-1`, `N` included, is excluded by the index check. The two
//!    exclusions have different reasons and a port that merged them would still agree;
//!  * **`splitBySample` throws on a read with no sample** when `unknownSampleName` is null, and
//!    the exception names the first such read rather than counting them;
//!  * **the overlap fix truncates twice.** Agreeing bases sum and cap at 93; disagreeing ones
//!    multiply the winner by 0.8 in *double* arithmetic and then cast to `byte`, which truncates
//!    toward zero, so Q30 becomes Q24 and Q31 becomes Q24 as well.
//!
//! # Reading this file without reading Rust
//!
//! - `keep: impl Fn(&PileupElement<'a>) -> bool` is a parameter that takes a closure, Java's
//!   `Predicate<PileupElement>`. `|e| keep(e)` is a lambda.
//! - `.iter().filter(..).cloned().collect()` is Java's `stream().filter(..).collect(toList())`;
//!   `.cloned()` copies each element out, and a `PileupElement` copy shares the read rather than
//!   duplicating it, because the element only holds a borrowed `&BamRecord`.
//! - `Option<String>` is a nullable `String`. Where the reference returns a `Set<String>` that may
//!   contain `null`, this returns a `Vec<Option<String>>` with `None` in that place.
//! - `Result<T, String>` is "a value or an error message": `Err(..)` is where the reference
//!   throws.
//! - `x as i8` keeps the low eight bits, so it wraps exactly like Java's `(byte)` cast of an
//!   `int`. `f as u8` from a double truncates toward zero, as Java's `(byte)` does for a
//!   quality up to 127 (Java's `byte` is signed, so the two read a larger quality differently).

use crate::pileup::PileupElement;
use crate::read_group;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `QualityUtils.MAX_SAM_QUAL_SCORE`, which is `SAMUtils.MAX_PHRED_SCORE`.
pub const MAX_SAM_QUAL_SCORE: u8 = 93;

/// `ReadPileup.SAMTOOLS_OVERLAP_LOW_CONFIDENCE`.
pub const SAMTOOLS_OVERLAP_LOW_CONFIDENCE: f64 = 0.8;

/// `BaseUtils.simpleBaseToBaseIndex`, which now lives in [`crate::base_utils`] because the
/// covariates need it too. Kept here as a re-export so the pileup code reads as it did.
pub use crate::base_utils::simple_base_to_base_index;

/// `ReadPileup`: a locus and the elements at it, in the order they were assembled.
///
/// The reference holds a `Locatable loc`; only its contig and start are ever read, so those two
/// are stored directly.
pub struct ReadPileup<'a> {
    pub contig: String,
    /// `loc.getStart()`, 1-based.
    pub start: i32,
    pub elements: Vec<PileupElement<'a>>,
}

impl<'a> ReadPileup<'a> {
    /// `ReadPileup(Locatable, List<PileupElement>)`.
    pub fn new(contig: &str, start: i32, elements: Vec<PileupElement<'a>>) -> Self {
        ReadPileup {
            contig: contig.to_string(),
            start,
            elements,
        }
    }

    /// `size`.
    pub fn size(&self) -> usize {
        self.elements.len()
    }

    /// `isEmpty`.
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// `makeFilteredPileup`, which shares the underlying reads rather than copying them.
    pub fn filtered(&self, keep: impl Fn(&PileupElement<'a>) -> bool) -> ReadPileup<'a> {
        ReadPileup {
            contig: self.contig.clone(),
            start: self.start,
            elements: self.elements.iter().filter(|e| keep(e)).cloned().collect(),
        }
    }

    /// `getBaseCounts`: A, C, G, T, deletions skipped and unmapped bases dropped.
    ///
    /// "Unmapped" here means a byte `simpleBaseToBaseIndex` answers `-1` for: see the module note
    /// on `*`, which is counted as an `A`. `[i32; 4]` is Java's `new int[4]`, a fixed-size array.
    pub fn base_counts(&self) -> [i32; 4] {
        let mut counts = [0i32; 4];
        for element in &self.elements {
            if element.is_deletion() {
                continue;
            }
            let index = simple_base_to_base_index(element.base());
            if index != -1 {
                counts[index as usize] += 1;
            }
        }
        counts
    }

    /// `getBases`: one byte per element, `D` for a deletion.
    pub fn bases(&self) -> Vec<u8> {
        self.elements.iter().map(|e| e.base()).collect()
    }

    /// `getBaseQuals`: one byte per element, 16 for a deletion.
    pub fn quals(&self) -> Vec<u8> {
        self.elements.iter().map(|e| e.qual()).collect()
    }

    /// `getOffsets`.
    pub fn offsets(&self) -> Vec<i32> {
        self.elements.iter().map(|e| e.offset).collect()
    }

    /// `getNumberOfElements`: how many elements pass `keep`.
    pub fn number_of_elements(&self, keep: impl Fn(&PileupElement<'a>) -> bool) -> usize {
        self.elements.iter().filter(|e| keep(e)).count()
    }

    /// `sortedIterator`: by the **read's** start, and by nothing else.
    ///
    /// Java's `Stream.sorted` is stable, so elements from reads starting at the same position keep
    /// the order they were assembled in. Sorting by anything more, the offset for instance, would
    /// reorder exactly those ties. Rust's `sort_by_key` is stable too.
    ///
    /// ```java
    /// .sorted((l, r) -> Integer.compare(l.getRead().getStart(), r.getRead().getStart()))
    /// ```
    pub fn sorted(&self) -> Vec<PileupElement<'a>> {
        let mut sorted = self.elements.clone();
        sorted.sort_by_key(|e| e.read.alignment_start);
        sorted
    }

    /// `getReadGroupIDs`, as a set. The reference returns a `Set`, so a read with no read group
    /// contributes a null entry rather than being dropped.
    ///
    /// The `Vec` is in first-seen order, where the reference's `Collectors.toSet()` has no defined
    /// order; `contains` before `push` is what makes it a set.
    pub fn read_group_ids(&self) -> Vec<Option<String>> {
        let mut ids: Vec<Option<String>> = Vec::new();
        for element in &self.elements {
            let id = read_group_of(element.read);
            if !ids.contains(&id) {
                ids.push(id);
            }
        }
        ids
    }

    /// `getSamples`: `ReadUtils.getSampleName`, which is null when the read has no read group *or*
    /// when the read group declares no sample.
    pub fn samples(&self, header: &SamHeader) -> Vec<Option<String>> {
        let mut samples: Vec<Option<String>> = Vec::new();
        for element in &self.elements {
            let sample = sample_name(element.read, header);
            if !samples.contains(&sample) {
                samples.push(sample);
            }
        }
        samples
    }

    /// `getPileupForSample`: the elements whose read's sample equals `sample`, with `None`
    /// matching the reads that have no sample (`Objects.equals` of two nulls is true).
    pub fn pileup_for_sample(&self, sample: Option<&str>, header: &SamHeader) -> ReadPileup<'a> {
        let wanted = sample.map(|s| s.to_string());
        self.filtered(|e| sample_name(e.read, header) == wanted)
    }

    /// `splitBySample`.
    ///
    /// Returns `Err` where the reference throws `UserException.ReadMissingReadGroup`: a read with
    /// no sample, and no `unknownSampleName` to file it under. The error carries that read's name,
    /// as the reference's message does.
    ///
    /// The reference returns a `HashMap`; this returns the pairs in the first-seen order of
    /// [`Self::samples`], so a caller that needs the reference's iteration order has to impose it.
    pub fn split_by_sample(
        &self,
        header: &SamHeader,
        unknown_sample_name: Option<&str>,
    ) -> Result<Vec<(String, ReadPileup<'a>)>, String> {
        let mut split = Vec::new();
        for sample in self.samples(header) {
            // `sample.as_deref()` borrows the `Option<String>` as an `Option<&str>`.
            let pileup = self.pileup_for_sample(sample.as_deref(), header);
            // A `match` on the pair (sample, unknown name): the compiler checks that the arms
            // cover every combination, which is the reference's `if (sample != null) ... else`.
            match (&sample, unknown_sample_name) {
                (Some(name), _) => split.push((name.clone(), pileup)),
                (None, Some(unknown)) => split.push((unknown.to_string(), pileup)),
                (None, None) => {
                    // `pileupBySample.iterator().next().getRead()`: the first read of that
                    // sample's pileup, which cannot be empty because the sample came from it.
                    let name = pileup
                        .elements
                        .first()
                        .map(|e| e.read.read_name.clone())
                        .unwrap_or_default();
                    return Err(format!("Read {name} is missing a read group"));
                }
            }
        }
        Ok(split)
    }

    /// `getPileupString`: the samtools-like line a tool prints.
    ///
    /// ```java
    /// return String.format("%s %s %c %s %s",
    /// ```
    pub fn pileup_string(&self, reference_base: char) -> String {
        format!(
            "{} {} {} {} {}",
            self.contig,
            self.start,
            reference_base,
            String::from_utf8_lossy(&self.bases()),
            self.quals_string(),
        )
    }

    /// `getQualsString`: the qualities as printable SAM characters, each one `33 + q`.
    fn quals_string(&self) -> String {
        self.quals().iter().map(|q| (q + 33) as char).collect()
    }
}

/// `GATKRead.getReadGroup`: the `RG` tag's value, or `None` when it is absent or not a string.
fn read_group_of(read: &BamRecord) -> Option<String> {
    match read.tags.get(htsjdk_bam::tag::Tag::new(b"RG")) {
        Some(htsjdk_bam::tag::TagValue::Str(id)) => Some(id.clone()),
        _ => None,
    }
}

/// `ReadUtils.getSampleName`: the read group's `SM`, or null if either the group or the field is
/// absent. The two nulls are indistinguishable downstream, which is why `splitBySample` cannot say
/// which of the two happened.
pub fn sample_name(read: &BamRecord, header: &SamHeader) -> Option<String> {
    read_group::attribute(read, header, "SM").map(|s| s.to_string())
}

/// `ReadPileup.fixPairOverlappingQualities`, as a pure function over the two qualities.
///
/// Returns the new `(first, second)` qualities. The reference writes them back into the reads'
/// quality arrays, which this leaves to the caller so that the arithmetic can be compared on its
/// own.
///
/// Both branches truncate, and differently. Agreeing bases *sum*, and the sum is capped at 93 when
/// it exceeds it **or when it goes negative**, because the reference stores it in a `byte` first
/// and then tests the stored value: two Q80 bases sum to 160, which is -96 as a signed byte, and
/// the negative test is what catches it. Disagreeing bases multiply the winner by 0.8 in double
/// arithmetic and cast back, truncating toward zero.
///
/// The reference also skips the whole fix when either element is a deletion
/// (`if (!secondElement.isDeletion() && !firstElement.isDeletion())`). That test is not here,
/// because this function sees only bases and qualities; the caller has to make it.
pub fn fix_pair_overlapping_qualities(
    first_base: u8,
    first_qual: u8,
    second_base: u8,
    second_qual: u8,
) -> (u8, u8) {
    if first_base == second_base {
        let sum = first_qual as i32 + second_qual as i32;
        // The reference's `(byte)` store happens before the test, so this reproduces the store.
        let stored = sum as i8;
        let capped = if stored < 0 || stored > MAX_SAM_QUAL_SCORE as i8 {
            MAX_SAM_QUAL_SCORE
        } else {
            stored as u8
        };
        (capped, 0)
    } else if first_qual >= second_qual {
        (
            (SAMTOOLS_OVERLAP_LOW_CONFIDENCE * first_qual as f64) as u8,
            0,
        )
    } else {
        (
            0,
            (SAMTOOLS_OVERLAP_LOW_CONFIDENCE * second_qual as f64) as u8,
        )
    }
}

/// `ReadPileup(Locatable, Iterable<GATKRead>)`: the constructor a `VariantWalker` uses.
///
/// Three things happen before an element exists, and all three are the reference's:
///
///  * the reads are filtered by `PASSES_VENDOR_QUALITY_CHECK` **and** `NOT_DUPLICATE`, and by
///    nothing else. A failing-vendor-quality read never reaches a pileup even when the tool
///    disabled its own filters;
///  * each surviving read is stepped forward until its genome position reaches the locus, and
///    contributes only if it lands *exactly* on it. A read whose deletion spans the locus does
///    land on it, because the state machine stops once per deleted reference base;
///  * a read that steps off its right edge first contributes nothing. That is silent in the
///    reference only while the off-edge position is not the locus: stepping off advances the genome
///    offset one past the read's last aligned base, so for a read ending exactly one base before
///    the locus the reference's test passes and `makePileupElement` throws
///    `IllegalStateException` ("Cannot make a pileup element from an edge alignment state"). The
///    edge checks below drop that read instead.
///
/// ```java
/// while ( asm.stepForwardOnGenome() != null && asm.getGenomePosition() < loc.getStart()) { }
/// return asm.getGenomePosition() == loc.getStart() ? asm.makePileupElement() : null;
/// ```
pub fn pileup_from_reads<'a>(
    contig: &str,
    locus: i32,
    reads: &'a [BamRecord],
    passes_vendor_quality_check: impl Fn(&BamRecord) -> bool,
    not_duplicate: impl Fn(&BamRecord) -> bool,
) -> ReadPileup<'a> {
    let mut elements = Vec::new();
    for read in reads {
        if !passes_vendor_quality_check(read) || !not_duplicate(read) {
            continue;
        }
        let mut machine = crate::alignment_state::AlignmentStateMachine::new(read);
        // `while (asm.stepForwardOnGenome() != null && asm.getGenomePosition() < loc.getStart())`.
        // The loop leaves the machine wherever it stopped, and the test after it is what decides.
        // `while let Ok(Some(_)) = ..` keeps looping only while the step returns a position.
        // Off the right edge, or a malformed cigar, both end the loop. A malformed cigar is
        // where the reference throws `UserException.MalformedRead`; here the read is dropped.
        while let Ok(Some(_)) = machine.step_forward_on_genome() {
            if machine.genome_position() >= locus {
                break;
            }
        }
        if machine.genome_position() == locus && !machine.is_left_edge() && !machine.is_right_edge()
        {
            if let Some(element) = PileupElement::from_state(read, &machine) {
                elements.push(element);
            }
        }
    }
    ReadPileup::new(contig, locus, elements)
}
