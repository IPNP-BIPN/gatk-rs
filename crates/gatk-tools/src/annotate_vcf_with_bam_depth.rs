//! `AnnotateVcfWithBamDepth`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithBamDepth`
//! (GATK 4.6.2.0).
//!
//! Each record is written back out with one INFO field added, counting the reads of a separate bam
//! that cover it. The count is not a pileup: it is five conditions on the read's own coordinates
//! and flags.
//!
//! ```java
//! if (!read.failsVendorQualityCheck() && !read.isDuplicate() && !read.isUnmapped()
//!         && read.getEnd() > read.getStart() && new SimpleInterval(read).contains(vc)) {
//!     depth.increment();
//! }
//! ```
//!
//! # A read one base long is never counted
//!
//! `read.getEnd() > read.getStart()` is a strict inequality, so a `1M` read sitting exactly on the
//! variant contributes nothing. Nothing in the tool's documentation says a one-base read is not
//! coverage; the golden's `BAM_DEPTH=0` at a site such a read covers is what says it.
//!
//! # The read must contain the variant's whole span
//!
//! `contains`, not `overlaps`. A four-base deletion is counted only by reads covering all four
//! bases, and a record carrying `END` is asked about that whole block: the golden's `<DEL>` block
//! from 80 to 110 is `BAM_DEPTH=0` even though reads cover its start.
//!
//! # The annotation is written even when it is zero, and it overwrites
//!
//! Every record gets `BAM_DEPTH`, and `VariantContextBuilder.attribute` replaces whatever the input
//! carried: the golden's record with `BAM_DEPTH=99` comes out as `BAM_DEPTH=0`.
//!
//! # The header ends up with two BAM_DEPTH lines
//!
//! ```java
//! final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
//! headerLines.add(new VCFInfoHeaderLine(POOLED_BAM_DEPTH_ANNOTATION_NAME, 1, Integer, "pooled bam depth"));
//! ```
//!
//! The set is a set of **lines**, not of IDs, so an input that already declares `BAM_DEPTH` with a
//! different description keeps its line and the tool's is added beside it. The output VCF therefore
//! declares the same INFO ID twice, sorted by description, which is what the golden shows and what
//! a port deduplicating by ID would silently repair.

use htsjdk_vcf::variant::VariantContext;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK AnnotateVcfWithBamDepth";

/// `POOLED_BAM_DEPTH_ANNOTATION_NAME`.
pub const BAM_DEPTH: &str = "BAM_DEPTH";

/// The INFO line the tool adds, whatever the input already declared.
pub const BAM_DEPTH_HEADER_LINE: &str =
    "##INFO=<ID=BAM_DEPTH,Number=1,Type=Integer,Description=\"pooled bam depth\">";

/// As much of a read as the five conditions look at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Read<'a> {
    pub contig: &'a str,
    /// 1-based, inclusive.
    pub start: i32,
    /// 1-based, inclusive, from the cigar's reference length.
    pub end: i32,
    pub flags: u16,
}

const DUPLICATE: u16 = 0x400;
const VENDOR_QUALITY_CHECK_FAILED: u16 = 0x200;
const UNMAPPED: u16 = 0x4;

impl Read<'_> {
    fn is_duplicate(&self) -> bool {
        self.flags & DUPLICATE != 0
    }

    fn fails_vendor_quality_check(&self) -> bool {
        self.flags & VENDOR_QUALITY_CHECK_FAILED != 0
    }

    fn is_unmapped(&self) -> bool {
        self.flags & UNMAPPED != 0
    }
}

/// `apply`'s condition, read by read.
///
/// The order is the reference's, which matters only for a read that would fail more than one of
/// them: the answer is the same either way, and keeping it makes the two readable side by side.
pub fn counts_towards_depth(read: &Read, variant: &VariantContext) -> bool {
    !read.fails_vendor_quality_check()
        && !read.is_duplicate()
        && !read.is_unmapped()
        // A one-base read is excluded by this strict inequality.
        && read.end > read.start
        // Containment of the variant's whole span, not overlap.
        && read.contig == variant.contig
        && read.start <= variant.start as i32
        && variant.stop as i32 <= read.end
}

/// The `BAM_DEPTH` of one record.
pub fn bam_depth(reads: &[Read], variant: &VariantContext) -> i32 {
    reads
        .iter()
        .filter(|read| counts_towards_depth(read, variant))
        .count() as i32
}

/// `new VariantContextBuilder(vc).attribute(BAM_DEPTH, depth)`: the annotation replaces whatever
/// the record carried, and is written even when it is zero.
pub fn annotate(variant: &VariantContext, depth: i32) -> VariantContext {
    let mut annotated = variant.clone();
    annotated.attributes.retain(|(key, _)| key != BAM_DEPTH);
    annotated.attributes.push((
        BAM_DEPTH.to_string(),
        htsjdk_vcf::variant::Value::Int(depth as i64),
    ));
    annotated
}

/// The metadata lines of the output header, as the `HashSet` and the writer leave them.
///
/// The input's lines are kept as they are and the tool's `BAM_DEPTH` line is added beside them, so
/// an input already declaring `BAM_DEPTH` produces **two** lines with that ID. Only an exactly
/// identical line collapses, the set being a set of lines.
pub fn header_lines(input_lines: &[String]) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    for line in input_lines {
        if !lines.contains(line) {
            lines.push(line.clone());
        }
    }
    let added = BAM_DEPTH_HEADER_LINE.to_string();
    if !lines.contains(&added) {
        lines.push(added);
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn variant(start: i64, reference: &str, alternate: &str, end: Option<i64>) -> VariantContext {
        let alleles = vec![
            Allele::create(reference.as_bytes(), true).expect("a reference"),
            Allele::create(alternate.as_bytes(), false).expect("an alternate"),
        ];
        let mut context = VariantContext::new("chr1", start, alleles);
        context.stop = end.unwrap_or(start + reference.len() as i64 - 1);
        context
    }

    fn read(start: i32, end: i32, flags: u16) -> Read<'static> {
        Read {
            contig: "chr1",
            start,
            end,
            flags,
        }
    }

    #[test]
    fn a_read_one_base_long_is_never_counted() {
        let site = variant(40, "A", "C", None);
        assert!(!counts_towards_depth(&read(40, 40, 0), &site));
        assert!(counts_towards_depth(&read(40, 41, 0), &site));
    }

    #[test]
    fn the_read_must_contain_the_whole_span() {
        // The golden's four-base deletion at 60, and the two reads over it.
        let deletion = variant(60, "ACCC", "A", None);
        assert!(!counts_towards_depth(&read(58, 62, 0), &deletion));
        assert!(counts_towards_depth(&read(58, 77, 0), &deletion));

        // The block carrying END, whose span is what containment is asked about.
        let block = variant(80, "A", "<DEL>", Some(110));
        assert!(!counts_towards_depth(&read(78, 97, 0), &block));
    }

    #[test]
    fn duplicates_vendor_failed_and_unmapped_reads_are_excluded() {
        let site = variant(20, "A", "C", None);
        assert!(counts_towards_depth(&read(15, 34, 0), &site));
        assert!(!counts_towards_depth(&read(15, 34, DUPLICATE), &site));
        assert!(!counts_towards_depth(
            &read(15, 34, VENDOR_QUALITY_CHECK_FAILED),
            &site
        ));
        assert!(!counts_towards_depth(&read(15, 34, UNMAPPED), &site));
    }

    #[test]
    fn the_depth_of_the_goldens_first_site_is_two_of_four_reads() {
        let site = variant(20, "A", "C", None);
        let reads = [
            read(15, 34, 0),
            read(15, 34, 0),
            read(15, 34, DUPLICATE),
            read(15, 34, VENDOR_QUALITY_CHECK_FAILED),
        ];
        assert_eq!(bam_depth(&reads, &site), 2);
    }

    #[test]
    fn the_annotation_overwrites_and_is_written_at_zero() {
        let mut carried = variant(160, "A", "C", None);
        carried
            .attributes
            .push((BAM_DEPTH.to_string(), htsjdk_vcf::variant::Value::Int(99)));
        let annotated = annotate(&carried, 0);
        assert_eq!(
            annotated.attributes,
            vec![(BAM_DEPTH.to_string(), htsjdk_vcf::variant::Value::Int(0))]
        );
    }

    #[test]
    fn an_input_declaring_bam_depth_leaves_two_lines_in_the_header() {
        let input = vec![
            "##INFO=<ID=BAM_DEPTH,Number=1,Type=Integer,Description=\"was already here\">"
                .to_string(),
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the block\">".to_string(),
        ];
        let lines = header_lines(&input);
        assert_eq!(
            lines
                .iter()
                .filter(|line| line.contains("ID=BAM_DEPTH"))
                .count(),
            2,
            "the set is a set of lines, not of IDs"
        );
    }
}
