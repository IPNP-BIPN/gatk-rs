//! `CheckPileup`, ported from `org.broadinstitute.hellbender.tools.walkers.qc.CheckPileup`
//! (GATK 4.6.2.0).
//!
//! A `LocusWalker` that compares GATK's own pileup against a samtools mpileup file, locus by locus,
//! and prints every disagreement. The truth file is read by [`gatk_engine::sam_pileup`].
//!
//! # The pileup it reports is not the qualities the reads carry
//!
//! `fixOverlaps` runs by default and does what samtools does to an overlapping pair: where two mates
//! cover one locus with the same base, one quality becomes their sum and the other becomes **zero**.
//! Two mates carrying `I` and `J` come out `r!`. `--ignore-overlaps` turns it off, and the same run
//! then agrees with a truth file that recorded the reads' own qualities.
//!
//! # The comparison stops at the first failure
//!
//! Size, then location, then bases, then qualities, and the message names the one that failed. The
//! bases are compared case-insensitively and the qualities are not.
//!
//! The message prints the qualities as **raw Phred bytes cast to characters**, where the report line
//! above it prints them as SAM qualities: the same two qualities appear as `IJ` in one and `()` in
//! the other.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::locus_iterator::AlignmentContext;
use gatk_engine::read_pileup::{fix_pair_overlapping_qualities, ReadPileup};
use gatk_engine::sam_pileup::SamPileupFeature;

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CheckPileupError {
    /// A locus the truth file does not cover.
    NoTruthData { location: String, bases: String },
    /// A locus the two disagree on, with the comparison that failed.
    Mismatch(String),
}

impl CheckPileupError {
    /// The message `UserException.BadInput` carries, without its prefix.
    pub fn message(&self) -> String {
        match self {
            CheckPileupError::NoTruthData { location, bases } => format!(
                "No pileup data available at {location} given GATK's output of {bases} -- this walker requires samtools mpileup data over all bases"
            ),
            CheckPileupError::Mismatch(difference) => format!(
                "The input pileup doesn't match the GATK's internal pileup: {difference}"
            ),
        }
    }
}

/// The tool's own arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CheckPileupArguments {
    /// `--ignore-overlaps`.
    pub ignore_overlaps: bool,
    /// `--continue-after-error`.
    pub continue_after_error: bool,
}

/// What one run produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckPileupResult {
    /// The file the tool wrote, which a failing run leaves behind too.
    pub report: String,
    /// `onTraversalSuccess`, which is returned rather than written to the file.
    pub summary: String,
}

/// `pileupDiff(a, b)`: the first comparison that fails, or nothing.
///
/// The order is the reference's, and so is the asymmetry: the bases are upper-cased on both sides
/// and the qualities are not touched. The qualities are a parameter rather than the pileup's own,
/// because the overlap fixing happens before this is called and changes them.
pub fn pileup_diff(
    contig: &str,
    position: i32,
    bases: &[u8],
    quals: &[u8],
    truth: &SamPileupFeature,
) -> Option<String> {
    if bases.len() != truth.size() {
        return Some(format!(
            "Sizes not equal: {} vs. {}",
            bases.len(),
            truth.size()
        ));
    }
    if contig != truth.contig || position != truth.position {
        return Some(format!(
            "Locations not equal: {contig}:{position}-{position} vs. {}",
            truth.contig
        ));
    }
    let ours = String::from_utf8_lossy(bases).to_uppercase();
    let theirs = truth.bases_string().to_uppercase();
    if ours != theirs {
        return Some(format!(
            "Bases not equal: {} vs. {}",
            String::from_utf8_lossy(bases),
            truth.bases_string()
        ));
    }
    // The qualities are compared as the characters the raw Phred bytes cast to, which is what
    // `new String(byte[])` does and is not the SAM encoding the report line uses.
    let our_quals = raw_chars(quals);
    let their_quals = raw_chars(&truth.base_quals());
    if our_quals != their_quals {
        return Some(format!("Quals not equal: {our_quals} vs. {their_quals}"));
    }
    None
}

/// `new String(byte[])`: each byte as the character of that code point.
fn raw_chars(quals: &[u8]) -> String {
    quals.iter().map(|&q| q as char).collect()
}

/// `ReadPileup.fixOverlaps`, as the pair of replacement qualities per element.
///
/// Returns the qualities the pileup should report, in element order. It is a function of the
/// pileup rather than a mutation of the reads, because the reads themselves are not changed by this
/// tool and a caller that wrote the new qualities back would change what the next locus sees.
pub fn fixed_overlap_quals(pileup: &ReadPileup<'_>) -> Vec<u8> {
    let elements = pileup.sorted();
    let mut quals: Vec<u8> = elements.iter().map(|element| element.qual()).collect();

    // `FragmentCollection.create`: the first element of a name is held until its mate arrives, and
    // only a read whose mate could overlap is held at all.
    let mut waiting: Vec<(String, usize)> = Vec::new();
    for (index, element) in elements.iter().enumerate() {
        let read = element.read;
        let could_overlap = gatk_engine::read::is_paired(read)
            && !gatk_engine::read::mate_is_unmapped(read)
            && read.mate_alignment_start != 0
            && read.mate_alignment_start <= gatk_engine::read_utils::end(read);
        if !could_overlap {
            continue;
        }
        match waiting.iter().position(|(name, _)| name == &read.read_name) {
            Some(position) => {
                let (_, first) = waiting.remove(position);
                let (left, right) = fix_pair_overlapping_qualities(
                    elements[first].base(),
                    quals[first],
                    element.base(),
                    quals[index],
                );
                quals[first] = left;
                quals[index] = right;
            }
            None => waiting.push((read.read_name.clone(), index)),
        }
    }
    quals
}

/// `apply`: one locus compared, as the line it prints and the refusal it raises.
///
/// Returns the line to write, if any, and the error the run would throw. Both can be present: the
/// reference writes the line and then throws, so a failing run leaves the line that explains it.
pub fn apply(
    context: &AlignmentContext<'_>,
    reference_base: u8,
    truth: Option<&SamPileupFeature>,
    arguments: &CheckPileupArguments,
) -> (Option<String>, Option<CheckPileupError>) {
    let quals = if arguments.ignore_overlaps {
        context.pileup.quals()
    } else {
        fixed_overlap_quals(&context.pileup)
    };
    let bases = context.pileup.bases();
    let printed = pileup_string(
        &context.contig,
        context.position,
        reference_base,
        &bases,
        &quals,
    );

    let Some(truth) = truth else {
        let line = format!("No truth pileup data available at {printed}\n");
        let location = SimpleInterval {
            contig: context.contig.clone(),
            start: context.position,
            end: context.position,
        };
        return (
            Some(line),
            Some(CheckPileupError::NoTruthData {
                location: format!("{}:{}-{}", location.contig, location.start, location.end),
                bases: String::from_utf8_lossy(&bases).into_owned(),
            }),
        );
    };

    let difference = pileup_diff(&context.contig, context.position, &bases, &quals, truth);
    match difference {
        None => (None, None),
        Some(difference) => (
            Some(format!("{printed} vs. {}\n", truth_string(truth))),
            Some(CheckPileupError::Mismatch(difference)),
        ),
    }
}

/// `ReadPileup.getPileupString(ref)`, with the qualities the caller decided on.
fn pileup_string(
    contig: &str,
    position: i32,
    reference_base: u8,
    bases: &[u8],
    quals: &[u8],
) -> String {
    format!(
        "{} {} {} {} {}",
        contig,
        position,
        reference_base as char,
        String::from_utf8_lossy(bases),
        quals.iter().map(|&q| (q + 33) as char).collect::<String>(),
    )
}

/// `SAMPileupFeature.getPileupString()`.
fn truth_string(truth: &SamPileupFeature) -> String {
    format!(
        "{} {} {} {} {}",
        truth.contig,
        truth.position,
        truth.reference_base as char,
        truth.bases_string(),
        truth
            .base_quals()
            .iter()
            .map(|&q| (q + 33) as char)
            .collect::<String>(),
    )
}

/// `onTraversalSuccess`, which is returned rather than written to the file.
pub fn summary(loci: i64, bases: i64) -> String {
    format!("Validated {loci} sites covered by {bases} bases\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use gatk_engine::sam_pileup::{SamPileupElement, SamPileupFeature};

    fn truth(bases: &str, quals: &[u8]) -> SamPileupFeature {
        SamPileupFeature {
            contig: "chr1".to_string(),
            position: 21,
            reference_base: b'C',
            elements: bases
                .bytes()
                .zip(quals.iter())
                .map(|(base, &qual)| SamPileupElement { base, qual })
                .collect(),
        }
    }

    #[test]
    fn the_summary_is_the_references_sentence() {
        assert_eq!(summary(4, 8), "Validated 4 sites covered by 8 bases\n");
    }

    #[test]
    fn the_qualities_are_compared_as_raw_bytes() {
        // 40 and 41 print as `IJ` in a pileup string and as `()` in the message.
        assert_eq!(raw_chars(&[40, 41]), "()");
        let feature = truth("CC", &[40, 40]);
        assert_eq!(raw_chars(&feature.base_quals()), "((");
    }
}
