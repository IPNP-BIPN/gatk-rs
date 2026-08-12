//! Ported from `org.broadinstitute.hellbender.tools.walkers.ReadAnonymizer` (GATK 4.6.2.0).
//!
//! The ninth whole tool of the record-transform archetype, and the first whose transform makes the
//! read a **different length**. Every base that disagrees with the reference is replaced by the
//! reference base, and the cigar is rewritten to say so.
//!
//! # A deletion lengthens the read and an insertion shortens it
//!
//! ```java
//! case X:
//! case D:
//!     for (int i = 0; i < cigarElement.getLength(); ++i) {
//!         newReadBases.add(refBases[refIndex + i]);
//!         newBaseQualities.add((byte)refQual);
//!     }
//!     iterCigarOp = useSimpleCigar ? CigarOperator.M : CigarOperator.EQ;
//! ...
//! case I:
//!     iterCigarOp = currentNewCigarOp;
//!     iterCigarOpCount = 0;
//! ```
//!
//! A deletion puts the reference bases the read did not have **into** it, so `4M2D6M` over ten bases
//! comes out twelve. An insertion drops its bases and contributes a cigar element of length **zero**
//! under whatever operator was already accumulating, so it merges into the previous element rather
//! than ending it.
//!
//! # And every M, X and D becomes one operator
//!
//! `=` by default, `M` under `--use-simple-cigar`. Consecutive elements of different kinds therefore
//! collapse into one, which is why `4M2D6M` comes out as a single `12=`.
//!
//! The last element is added **unconditionally** after the loop, so a read whose cigar ends in an
//! insertion emits whatever the accumulator held at that point.
//!
//! # Its filters are a sixth pattern, and the first without Wellformed
//!
//! Seven filters, and `WellformedReadFilter` is not among them: valid alignment start, valid
//! alignment end, read length equals cigar length, sequence is stored, matching bases and quals,
//! mapped, and alignment agrees with the header. Six of the seven are what `WellformedReadFilter`
//! itself is made of, listed separately.

use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK ReadAnonymizer";

/// `--ref-base-quality`, whose default is sixty and whose declared range is `[0, 60]`.
pub const DEFAULT_REF_BASE_QUALITY: u8 = 60;

/// This tool's own arguments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnonymizerArguments {
    /// `--ref-base-quality`. The parser refuses anything above sixty before the tool runs, so the
    /// port takes a `u8` and the range check belongs to the caller.
    pub ref_base_quality: u8,
    /// `--use-simple-cigar`: `M` instead of `=`.
    pub use_simple_cigar: bool,
}

impl Default for AnonymizerArguments {
    fn default() -> Self {
        AnonymizerArguments {
            ref_base_quality: DEFAULT_REF_BASE_QUALITY,
            use_simple_cigar: false,
        }
    }
}

/// `anonymizeRead(read, referenceContext)`.
///
/// `reference_bases` is the read's own span, which is what `referenceContext.getBases(readInterval)`
/// hands over.
pub fn anonymize_read(
    record: &BamRecord,
    reference_bases: &[u8],
    arguments: &AnonymizerArguments,
) -> BamRecord {
    let replacement = if arguments.use_simple_cigar {
        Op::M
    } else {
        Op::Eq
    };

    let mut new_bases: Vec<u8> = Vec::new();
    let mut new_quals: Vec<u8> = Vec::new();
    let mut new_cigar: Vec<CigarElement> = Vec::new();

    let mut read_index = 0usize;
    let mut ref_index = 0usize;
    // `null` until the first element sets it, which is what makes a leading insertion contribute
    // nothing at all.
    let mut current_op: Option<Op> = None;
    let mut current_count: u32 = 0;

    for element in &record.cigar.elements {
        let length = element.length as usize;
        let (iter_op, iter_count): (Option<Op>, u32) = match element.op {
            // Nothing to copy and nothing to replace.
            Op::H | Op::N | Op::P => (Some(element.op), element.length),
            // Kept as they are, bases and qualities both.
            Op::S | Op::Eq => {
                for i in 0..length {
                    new_bases.push(record.read_bases[read_index + i]);
                    new_quals.push(record.base_qualities[read_index + i]);
                }
                (Some(element.op), element.length)
            }
            Op::M => {
                for i in 0..length {
                    new_bases.push(reference_bases[ref_index + i]);
                    // A base that already agreed keeps its own quality; one that was replaced takes
                    // the reference quality.
                    if record.read_bases[read_index + i] == reference_bases[ref_index + i] {
                        new_quals.push(record.base_qualities[read_index + i]);
                    } else {
                        new_quals.push(arguments.ref_base_quality);
                    }
                }
                (Some(replacement), element.length)
            }
            // A deletion's reference bases are put INTO the read, so it grows.
            Op::X | Op::D => {
                for i in 0..length {
                    new_bases.push(reference_bases[ref_index + i]);
                    new_quals.push(arguments.ref_base_quality);
                }
                (Some(replacement), element.length)
            }
            // Dropped, and contributing a zero-length element under the operator already in hand.
            Op::I => (current_op, 0),
        };

        if iter_op == current_op {
            current_count += iter_count;
        } else {
            if let Some(op) = current_op {
                new_cigar.push(CigarElement {
                    length: current_count,
                    op,
                });
            }
            current_op = iter_op;
            current_count = iter_count;
        }

        if element.op.consumes_reference_bases() {
            ref_index += length;
        }
        if element.op.consumes_read_bases() {
            read_index += length;
        }
    }

    // Added unconditionally, which is where a trailing insertion's accumulator lands.
    if let Some(op) = current_op {
        new_cigar.push(CigarElement {
            length: current_count,
            op,
        });
    }

    let mut out = record.clone();
    out.cigar = Cigar::new(new_cigar);
    out.read_bases = new_bases;
    out.base_qualities = new_quals;

    // `clearAttributes()` then `setReadGroup(readGroup)`: everything else goes.
    let read_group = out.tags.get(htsjdk_bam::tag::Tag::new(b"RG")).cloned();
    out.tags = htsjdk_bam::tag::Tags::new();
    if let Some(value) = read_group {
        out.tags.insert(htsjdk_bam::tag::Tag::new(b"RG"), value);
    }
    out
}

/// `ReadAnonymizer`: every read that survives the filters, anonymised and written back out.
///
/// `contig_bases` is the whole reference contig; each read's own span is taken from it, as
/// `ReferenceContext.getBases(readInterval)` does.
pub fn read_anonymizer(
    source: &ReadsDataSource,
    contig_bases: &[u8],
    arguments: &AnonymizerArguments,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    let records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    let anonymised: Vec<BamRecord> = records
        .iter()
        .map(|record| {
            let start = gatk_engine::read_utils::start(record).max(1) as usize - 1;
            let end = (gatk_engine::read_utils::end(record) as usize).min(contig_bases.len());
            anonymize_read(record, &contig_bases[start..end], arguments)
        })
        .collect();
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    write_records(&header, &anonymised, options.create_output_bam_index)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(cigar: &[(u32, Op)], bases: &[u8]) -> BamRecord {
        BamRecord {
            read_name: "r".to_string(),
            reference_index: 0,
            alignment_start: 1,
            cigar: Cigar::new(
                cigar
                    .iter()
                    .map(|(length, op)| CigarElement {
                        length: *length,
                        op: *op,
                    })
                    .collect(),
            ),
            read_bases: bases.to_vec(),
            base_qualities: (0..bases.len()).map(|i| (10 + i * 3) as u8).collect(),
            mapping_quality: 60,
            ..BamRecord::default()
        }
    }

    #[test]
    fn a_deletion_lengthens_the_read() {
        let record = read(&[(4, Op::M), (2, Op::D), (6, Op::M)], b"ACGTACGTAC");
        let out = anonymize_read(&record, b"ACGTACGTACGT", &AnonymizerArguments::default());
        assert_eq!(out.read_bases.len(), 12);
        assert_eq!(out.base_qualities.len(), 12);
        // And the three elements collapse into one, because M and D both become `=`.
        assert_eq!(out.cigar.elements.len(), 1);
        assert_eq!(out.cigar.elements[0].length, 12);
        assert_eq!(out.cigar.elements[0].op, Op::Eq);
    }

    #[test]
    fn an_insertion_shortens_it_and_merges_into_the_element_before() {
        let record = read(&[(4, Op::M), (2, Op::I), (4, Op::M)], b"ACGTACGTAC");
        let out = anonymize_read(&record, b"ACGTACGT", &AnonymizerArguments::default());
        assert_eq!(out.read_bases.len(), 8);
        assert_eq!(out.cigar.elements.len(), 1);
        assert_eq!(out.cigar.elements[0].length, 8);
    }

    #[test]
    fn a_matching_base_keeps_its_quality_and_a_replaced_one_does_not() {
        let record = read(&[(4, Op::M)], b"ACGT");
        // Two agree, two do not.
        let out = anonymize_read(&record, b"ACTT", &AnonymizerArguments::default());
        assert_eq!(out.read_bases, b"ACTT".to_vec());
        assert_eq!(out.base_qualities, vec![10, 13, 60, 19]);
    }

    #[test]
    fn every_attribute_but_the_read_group_is_cleared() {
        let mut record = read(&[(4, Op::M)], b"ACGT");
        record.tags.insert(
            htsjdk_bam::tag::Tag::new(b"RG"),
            htsjdk_bam::tag::TagValue::Str("rg1".to_string()),
        );
        record.tags.insert(
            htsjdk_bam::tag::Tag::new(b"NM"),
            htsjdk_bam::tag::TagValue::Int(3),
        );
        let out = anonymize_read(&record, b"ACGT", &AnonymizerArguments::default());
        assert_eq!(out.tags.len(), 1);
        assert!(out.tags.get(htsjdk_bam::tag::Tag::new(b"RG")).is_some());
    }
}
