//! `SplitNCigarReads`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.rnaseq.SplitNCigarReads` (GATK 4.6.2.0), with the
//! two transformers it installs: `NDNCigarReadTransformer` and `MappingQualityReadTransformer`.
//!
//! The tool [`gatk_engine::overhang_fixing_manager`] was built for. Every `N` in a read's cigar is a
//! splice, and a read spanning k of them becomes k+1 reads: each keeps **all** the bases and soft
//! clips everything outside its own section, so `3M2N5M` comes out as `3M5S` at the read's start and
//! `3S5M` five bases later, both carrying the same bases and qualities.
//!
//! # Two passes over the same reads
//!
//! It is a `MultiplePassReadWalker`. The first pass splits and clips but writes nothing: it records
//! which reads the overhang clipper moved. The second pass splits again, this time with the manager
//! writing, so a read whose **mate** was moved has its mate position and `MC` tag repaired from what
//! the first pass learned.
//!
//! # Three edges of the split that are not obvious
//!
//! ```java
//! while (read.getCigarElement(cigarFirstIndex).getOperator().equals(CigarOperator.D)) cigarFirstIndex++;
//! while (read.getCigarElement(cigarSecondIndex-1).getOperator().equals(CigarOperator.D)) cigarSecondIndex--;
//! if (cigarFirstIndex > cigarSecondIndex) throw new IllegalArgumentException(...);
//! ```
//!
//!  * a section that would begin or end on a **deletion** is trimmed back first, so `3M1D2N5M` and
//!    `3M2N1D5M` both come out the same shape as `3M2N5M`;
//!  * a cigar that **ends** in `N` emits one piece and loses the `N`: `8M2N` becomes `8M`;
//!  * a cigar that **begins** with `N` is passed through untouched, `2N8M` and all, because the
//!    leading `N` produces no section and the tool returns the read it was handed.
//!
//! `N-D-N` is where the throw lands: the middle section is empty, and the message names the cigar.
//! `--refactor-cigar-string` merges the motif into one `N` before the read reaches the filters, and
//! it merges **one motif per two elements skipped**, so `3M2N1D2N1D2N2M` comes out `3M5N1D2N2M`
//! rather than fully merged.
//!
//! # The mapping quality transform is on by default
//!
//! 255 becomes 60, and only 255: it exists so that STAR's "uniquely mapped" reads are not dropped by
//! `HaplotypeCaller`'s mapping quality filter. `--skip-mapping-quality-transform` leaves it alone.
//! The transform runs **after** the read filters, and this tool's only filter is `ALLOW_ALL_READS`,
//! so a read whose cigar disagrees with its bases is split and written like any other.

use gatk_engine::cigar_utils;
use gatk_engine::clipping;
use gatk_engine::overhang_fixing_manager::{
    OverhangArguments, OverhangError, OverhangFixingManager, ReferenceQuery,
};
use gatk_engine::read;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK SplitNCigarReads";

/// `MATE_CIGAR_TAG`.
pub const MATE_CIGAR_TAG: &[u8; 2] = b"MC";

/// The mapping quality the transform rewrites, and what it rewrites it to.
pub const FROM_QUALITY: u8 = 255;
pub const TO_QUALITY: u8 = 60;

/// This tool's own arguments, whose defaults are the reference's.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SplitArguments {
    /// `--refactor-cigar-string`.
    pub refactor_ndn_cigar_reads: bool,
    /// `--skip-mapping-quality-transform`.
    pub skip_mq_transform: bool,
    /// `--process-secondary-alignments`, which the manager also reads.
    pub process_secondary_alignments: bool,
    /// The manager's own arguments: the queue size, the overhang tolerances and whether to fix.
    pub overhang: OverhangArguments,
}

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitError {
    /// `splitReadBasedOnCigar`: a section with nothing in it, which `N-D-N` produces.
    CannotSplit(String),
    /// Anything the manager or the clipper refused.
    Overhang(OverhangError),
    /// A clip the clipper refused.
    Clip(clipping::ClipError),
}

impl SplitError {
    /// The message the reference carries.
    pub fn message(&self) -> String {
        match self {
            SplitError::CannotSplit(cigar) => format!(
                "Cannot split this read (might be an empty section between Ns, for example 1N1D1N): {cigar}"
            ),
            SplitError::Overhang(error) => error.message(),
            SplitError::Clip(error) => format!("{error:?}"),
        }
    }
}

/// `NDNCigarReadTransformer.refactorNDNtoN`.
///
/// One `N-D-N` becomes one `N` of the three lengths added. The loop then skips **two** elements, so
/// a run of motifs sharing an element is only partly merged: `3M2N1D2N1D2N2M` comes out
/// `3M5N1D2N2M`, where a second pass would have merged the rest.
pub fn refactor_ndn_to_n(cigar: &Cigar) -> Cigar {
    let elements = &cigar.elements;
    let mut refactored: Vec<CigarElement> = Vec::with_capacity(elements.len());
    let mut i = 0;
    while i < elements.len() {
        let element = elements[i];
        // `thereAreAtLeast2MoreElements`: `index < cigarLength - 2`, so a trailing N is left alone.
        if element.op == Op::N && i + 2 < elements.len() {
            let next = elements[i + 1];
            let next_next = elements[i + 2];
            if next.op == Op::D && next_next.op == Op::N {
                refactored.push(CigarElement {
                    length: element.length + next.length + next_next.length,
                    op: Op::N,
                });
                i += 3;
                continue;
            }
        }
        refactored.push(element);
        i += 1;
    }
    Cigar::new(refactored)
}

/// `MappingQualityReadTransformer(255, 60)`.
pub fn transform_mapping_quality(record: &mut BamRecord, from: u8, to: u8) {
    if record.mapping_quality == from {
        record.mapping_quality = to;
    }
}

/// The pieces one read is split into, and the splices it declares along the way.
pub type SplitFamily = (Vec<BamRecord>, Vec<(i32, i32)>);

/// `splitNCigarRead(read, manager, emitReads, header, secondaryAlignments)`.
///
/// With `emit` set the family is handed to the manager and the read itself is returned; without it
/// the manager is never touched and the **first** piece is returned, which is how the tool predicts
/// what a mate's cigar will become.
pub fn split_n_cigar_read(
    record: &BamRecord,
    header: &SamHeader,
    secondary_alignments: bool,
) -> Result<SplitFamily, SplitError> {
    let elements = record.cigar.elements.clone();
    let mut split_reads: Vec<BamRecord> = Vec::with_capacity(2);

    // A secondary alignment is passed to the manager whole, and its mate information is still
    // repaired: the split is what it is spared, not the traversal.
    if !secondary_alignments && read::is_secondary_alignment(record) {
        return Ok((vec![record.clone()], Vec::new()));
    }

    let mut section_has_match = false;
    let mut first_cigar_index = 0usize;
    let mut splices: Vec<(i32, i32)> = Vec::new();

    for (i, element) in elements.iter().enumerate() {
        let op = element.op;
        if matches!(op, Op::M | Op::Eq | Op::X | Op::I | Op::D) {
            section_has_match = true;
        }
        if op == Op::N {
            if section_has_match {
                let (piece, splice) =
                    split_read_based_on_cigar(record, first_cigar_index, i, header, true)?;
                if let Some(splice) = splice {
                    splices.push(splice);
                }
                split_reads.push(piece);
            }
            first_cigar_index = i + 1;
            section_has_match = false;
        }
    }

    // No N at all, or nothing before the first one: the read is handed over as it arrived.
    if split_reads.is_empty() {
        return Ok((vec![record.clone()], splices));
    }
    // The last section, which a cigar ending in N does not have.
    if first_cigar_index < elements.len() && section_has_match {
        let (piece, _) =
            split_read_based_on_cigar(record, first_cigar_index, elements.len(), header, false)?;
        split_reads.push(piece);
    }

    Ok((split_reads, splices))
}

/// `splitReadBasedOnCigar`: one section of the read, with everything else soft clipped.
///
/// Returns the piece and, when asked for, the splice it sits beside: the `N` element at
/// `cigar_end_index`, whose **untrimmed** end is used on purpose, so deletions at the end of the
/// section do not move the splice.
pub fn split_read_based_on_cigar(
    record: &BamRecord,
    cigar_start_index: usize,
    cigar_end_index: usize,
    header: &SamHeader,
    want_splice: bool,
) -> Result<(BamRecord, Option<(i32, i32)>), SplitError> {
    let elements = &record.cigar.elements;
    let mut first = cigar_start_index;
    let mut second = cigar_end_index;

    // A section that begins or ends on a deletion is trimmed back before anything is measured.
    while first < elements.len() && elements[first].op == Op::D {
        first += 1;
    }
    while second > 0 && elements[second - 1].op == Op::D {
        second -= 1;
    }
    if first > second {
        return Err(SplitError::CannotSplit(record.cigar.to_text()));
    }

    let start_ref_index = gatk_engine::read_utils::unclipped_start(record)
        + cigar_utils::count_ref_bases_and_clips(&elements[0..first]);
    let stop_ref_index =
        start_ref_index + cigar_utils::count_ref_bases_and_clips(&elements[first..second]) - 1;

    let splice = if want_splice && cigar_end_index < elements.len() {
        let splice_start = start_ref_index
            + cigar_utils::count_ref_bases_and_clips(&elements[first..cigar_end_index]);
        let splice_end = splice_start + elements[cigar_end_index].length as i32 - 1;
        Some((splice_start, splice_end))
    } else {
        None
    };

    let piece = clipping::soft_clip_to_region_including_clipped_bases(
        record,
        Some(header),
        start_ref_index,
        stop_ref_index,
    )
    .map_err(SplitError::Clip)?;
    Ok((piece, splice))
}

/// `SplitNCigarReads`: the whole tool, both passes.
///
/// `reference` answers the manager's one query per splice. The output is what the manager handed to
/// the writer, in that order, which is **not** traversal order: the writer is created with
/// `presorted = false` because the manager makes no such promise.
pub fn split_n_cigar_reads(
    source: &ReadsDataSource,
    arguments: &SplitArguments,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    reference: ReferenceQuery<'_>,
) -> Result<(Vec<u8>, Option<Vec<u8>>), SplitToolError> {
    let header = source.header().clone();
    let raw = crate::read_walker::traverse(source, &options.intervals, &|_| true)
        .map_err(SplitToolError::Reads)?;

    // `makePreReadFilterTransformer` runs before the filters, `makePostReadFilterTransformer`
    // after, and this tool's only filter is ALLOW_ALL_READS unless the caller replaced it.
    let mut records = Vec::with_capacity(raw.len());
    for mut record in raw {
        if arguments.refactor_ndn_cigar_reads {
            record.cigar = refactor_ndn_to_n(&record.cigar);
        }
        if !filter(&record) {
            continue;
        }
        if !arguments.skip_mq_transform {
            transform_mapping_quality(&mut record, FROM_QUALITY, TO_QUALITY);
        }
        records.push(record);
    }

    let overhang = OverhangArguments {
        process_secondary_reads: arguments.process_secondary_alignments,
        ..arguments.overhang.clone()
    };
    let mut manager = OverhangFixingManager::new(&header, overhang);

    // The first pass records which reads the clipper moved; the second writes them out.
    for pass in 0..2 {
        if pass == 1 {
            manager
                .activate_writing()
                .map_err(|error| SplitToolError::Split(SplitError::Overhang(error)))?;
        }
        for record in &records {
            let mut record = record.clone();
            // The MC tag is rewritten to what the mate's cigar will become, worked out by running
            // the split over an artificial read carrying that cigar and nothing else.
            if let Some(TagValue::Str(mate_cigar)) =
                record.tags.get(Tag::new(MATE_CIGAR_TAG)).cloned()
            {
                let predicted = predicted_mate_cigar(&mate_cigar, &header, arguments)?;
                record
                    .tags
                    .insert(Tag::new(MATE_CIGAR_TAG), TagValue::Str(predicted));
            }
            manager.set_predicted_mate_information(&mut record);

            let (family, splices) =
                split_n_cigar_read(&record, &header, arguments.process_secondary_alignments)
                    .map_err(SplitToolError::Split)?;
            // The reference records each splice as it builds the piece beside it, before the family
            // is queued, so every splice this read declares is already held when its own family is
            // checked against them.
            for (start, end) in splices {
                let contig = contig_of(&record, &header).unwrap_or_default().to_string();
                manager
                    .add_splice_position(&contig, start, end, reference)
                    .map_err(|error| SplitToolError::Split(SplitError::Overhang(error)))?;
            }
            manager
                .add_read_group(&family)
                .map_err(|error| SplitToolError::Split(SplitError::Overhang(error)))?;
        }
        manager
            .flush()
            .map_err(|error| SplitToolError::Split(SplitError::Overhang(error)))?;
    }

    let mut written = std::mem::take(&mut manager.written);
    // `createSAMWriter(OUTPUT, false)`: the writer is told the reads are NOT sorted, so htsjdk's
    // sorting collection puts them in coordinate order before they reach the file. The manager
    // hands over whole families at once, so without this the second piece of a family comes out
    // beside the first instead of where its own coordinate puts it.
    written.sort_by(htsjdk_bam::coordinate::compare);
    // TWO @PG records, not one. `onTraversalStart` keeps `getHeaderForSAMWriter()` in a field and
    // then calls `createSAMWriter`, which asks for the header again; each call appends a program
    // record and makes its ID unique, so the output header carries `GATK SplitNCigarReads` and
    // `GATK SplitNCigarReads.1` with the same command line. No other tool ported so far does this,
    // and it is worth two lines here rather than a note somewhere else: it is 1266 bytes of the
    // output file.
    let once = header_for_sam_writer(source.header(), TOOL_NAME, options);
    let out_header = header_for_sam_writer(&once, TOOL_NAME, options);
    write_records(&out_header, &written, options.create_output_bam_index)
        .map_err(SplitToolError::Reads)
}

/// What the tool as a whole can fail with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitToolError {
    Reads(ReadsError),
    Split(SplitError),
}

/// The cigar a mate's `MC` tag will hold after the mate is itself split.
///
/// The reference builds an artificial read carrying that cigar, splits it with the manager left out
/// of it, and takes the **first** piece's cigar.
fn predicted_mate_cigar(
    mate_cigar: &str,
    header: &SamHeader,
    arguments: &SplitArguments,
) -> Result<String, SplitToolError> {
    let mut artificial = artificial_read(mate_cigar);
    if arguments.refactor_ndn_cigar_reads {
        artificial.cigar = refactor_ndn_to_n(&artificial.cigar);
    }
    let (pieces, _) =
        split_n_cigar_read(&artificial, header, arguments.process_secondary_alignments)
            .map_err(SplitToolError::Split)?;
    Ok(pieces[0].cigar.to_text())
}

/// `ArtificialReadUtils.createArtificialRead(header, cigar)`: a read at the first contig's first
/// base whose bases and qualities are as long as the cigar says.
fn artificial_read(cigar_text: &str) -> BamRecord {
    let cigar = parse_cigar(cigar_text);
    let length = cigar.read_length() as usize;
    BamRecord {
        read_name: "read".to_string(),
        reference_index: 0,
        alignment_start: 1,
        mapping_quality: 0,
        cigar,
        read_bases: vec![b'A'; length],
        base_qualities: vec![30; length],
        ..Default::default()
    }
}

/// `TextCigarCodec.decode`, for the one place a cigar arrives as text: an `MC` tag.
fn parse_cigar(text: &str) -> Cigar {
    if text == "*" {
        return Cigar::new(Vec::new());
    }
    let mut elements = Vec::new();
    let mut length = 0u32;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
            continue;
        }
        let op = match byte {
            b'M' => Op::M,
            b'I' => Op::I,
            b'D' => Op::D,
            b'N' => Op::N,
            b'S' => Op::S,
            b'H' => Op::H,
            b'P' => Op::P,
            b'=' => Op::Eq,
            b'X' => Op::X,
            _ => continue,
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    Cigar::new(elements)
}

fn contig_of<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a str> {
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.as_str())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(text: &str) -> Cigar {
        parse_cigar(text)
    }

    #[test]
    fn one_motif_is_merged_per_two_elements_skipped() {
        assert_eq!(refactor_ndn_to_n(&cigar("3M2N1D2N4M")).to_text(), "3M5N4M");
        // The second motif shares an element with the first, and is left where it is.
        assert_eq!(
            refactor_ndn_to_n(&cigar("3M2N1D2N1D2N2M")).to_text(),
            "3M5N1D2N2M"
        );
    }

    #[test]
    fn a_trailing_n_has_too_few_elements_after_it_to_merge() {
        assert_eq!(refactor_ndn_to_n(&cigar("3M2N")).to_text(), "3M2N");
        assert_eq!(refactor_ndn_to_n(&cigar("2N3M")).to_text(), "2N3M");
        // N-I-N is not the motif.
        assert_eq!(
            refactor_ndn_to_n(&cigar("3M2N1I2N4M")).to_text(),
            "3M2N1I2N4M"
        );
    }

    #[test]
    fn only_the_one_mapping_quality_is_rewritten() {
        let mut record = BamRecord {
            mapping_quality: 255,
            ..Default::default()
        };
        transform_mapping_quality(&mut record, FROM_QUALITY, TO_QUALITY);
        assert_eq!(record.mapping_quality, 60);

        let mut untouched = BamRecord {
            mapping_quality: 254,
            ..Default::default()
        };
        transform_mapping_quality(&mut untouched, FROM_QUALITY, TO_QUALITY);
        assert_eq!(untouched.mapping_quality, 254);
    }
}
