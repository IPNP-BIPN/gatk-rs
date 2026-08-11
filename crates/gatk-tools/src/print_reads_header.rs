//! Ported from `org.broadinstitute.hellbender.tools.PrintReadsHeader` (GATK 4.6.2.0).
//!
//! The eleventh whole tool of the record-transform archetype, the second that is not a walker, and
//! the smallest: forty-six lines of which four do anything.
//!
//! ```java
//! public void traverse() {
//!     final SAMFileHeader bamHeader = getHeaderForReads();
//!     try ( final OutputStreamWriter outputWriter = new OutputStreamWriter(outputFile.getOutputStream()) ) {
//!         final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
//!         codec.encode(outputWriter, bamHeader);
//!     } catch (IOException e ) { ... }
//! }
//! ```
//!
//! # The header it prints is not the header it was given
//!
//! `encode(writer, header)` is the two-argument overload, which passes
//! `keepExistingVersionNumber = false`. `writeHDLine(false)` then builds a **fresh**
//! `SAMFileHeader`, copies every attribute except `VN` into it, and lets the constructor's own `VN`
//! stand. Two consequences, and this port reproduces both:
//!
//!  * the version becomes [`htsjdk_bam::header::CURRENT_VERSION`], whatever the file said;
//!  * `VN` leads the `@HD` line, whatever position it held in the file, because the fresh header
//!    set it before any attribute was copied in.
//!
//! Neither is observable on a file htsjdk itself wrote, and the reason is the finding this tool's
//! measurement turned up: `SAMFileWriterImpl.writeHeader(SAMFileHeader)` goes through the *same*
//! two-argument overload, so **no BAM htsjdk produced can carry a non-current `VN` in the first
//! place**. The three-argument call that would keep it lives in
//! `BAMFileWriter.writeHeader(BinaryCodec, SAMFileHeader)` and is reachable only from the standalone
//! block-copy reheader.
//!
//! That makes it htsjdk-rs's problem rather than this tool's: `BamWriter` calls
//! [`htsjdk_bam::header::SamHeader::encode`], whose own comment claims the writer passes `true`, so
//! the port keeps a version htsjdk replaces on any input whose `VN` is not current. Filed as
//! htsjdk-rs#164. This module does the right thing locally and does not paper over the other one.
//!
//! # Nothing is appended
//!
//! `getHeaderForReads()`, not `getHeaderForSAMWriter()`, so the `@PG` chain is the file's own and no
//! record is added to it. Every other tool of this archetype does the opposite, which is the only
//! reason this one needs a module rather than a line.
//!
//! # The output is text, in the platform's charset
//!
//! `new OutputStreamWriter(stream)` with no charset argument. The pinned container reports UTF-8,
//! which is what the golden's bytes show and what this port writes. A container with a different
//! default would produce different bytes from the same header, which is why the dump records the
//! charset's name beside the bytes rather than leaving it to be inferred.

use htsjdk_bam::header::{Attributes, SamHeader, CURRENT_VERSION};

use gatk_engine::reads::ReadsDataSource;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK PrintReadsHeader";

/// `SAMTextHeaderCodec.encode(writer, header)`: the two-argument overload, which does **not** keep
/// the header's version.
///
/// `writeHDLine(false)` builds a fresh `SAMFileHeader` and copies every attribute except `VN` into
/// it, so `VN` is the current version and leads the line. The rest of the header is untouched:
/// [`SamHeader::encode`] is the `true` behaviour and is correct for everything below `@HD`.
pub fn encode_without_existing_version(header: &SamHeader) -> String {
    let mut rebuilt = header.clone();
    let mut attributes = Attributes::new();
    // The constructor's own VN, set before anything is copied in, which is what puts it first.
    attributes.set("VN", CURRENT_VERSION);
    for (key, value) in header.attributes.iter() {
        if key != "VN" {
            attributes.set(key, value);
        }
    }
    rebuilt.attributes = attributes;
    rebuilt.encode()
}

/// `PrintReadsHeader`: the reads header, as text.
///
/// The bytes rather than the string, because the reference writes through an `OutputStreamWriter`
/// with no charset and the golden compares what landed on disk. UTF-8 is the pinned container's
/// default and Rust's only encoding, so the two agree by construction here and would not on a
/// container that defaulted to anything else.
pub fn print_reads_header(source: &ReadsDataSource) -> Vec<u8> {
    encode_without_existing_version(source.header()).into_bytes()
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::{ProgramRecord, SequenceRecord};

    fn header() -> SamHeader {
        let mut header = SamHeader::default();
        header.set_sort_order("coordinate");
        header.sequences.push(SequenceRecord::new("chr1", 100));
        header.programs.push(ProgramRecord::new("upstream"));
        header.add_comment("a comment");
        header
    }

    #[test]
    fn the_version_is_replaced_rather_than_kept() {
        let mut old = header();
        old.attributes.set("VN", "1.5");
        assert!(
            old.encode().starts_with("@HD\tVN:1.5"),
            "the `true` encoding keeps it"
        );
        assert!(
            encode_without_existing_version(&old).starts_with("@HD\tVN:1.6"),
            "and this one does not"
        );
    }

    #[test]
    fn and_the_version_leads_the_line_whatever_position_it_held() {
        // SO first and VN last, which is an order the `true` encoding keeps and the reference's
        // fresh header does not. Built rather than mutated from the default, whose constructor
        // would have put VN in front already.
        let mut attributes = Attributes::new();
        attributes.set("SO", "coordinate");
        attributes.set("VN", "1.5");
        let moved = SamHeader {
            attributes,
            ..SamHeader::default()
        };
        assert_eq!(
            moved.encode().lines().next(),
            Some("@HD\tSO:coordinate\tVN:1.5")
        );
        assert_eq!(
            encode_without_existing_version(&moved).lines().next(),
            Some("@HD\tVN:1.6\tSO:coordinate")
        );
    }

    #[test]
    fn everything_below_the_hd_line_is_untouched() {
        let header = header();
        let ours = encode_without_existing_version(&header);
        let theirs = header.encode();
        assert_eq!(
            ours.lines().skip(1).collect::<Vec<_>>(),
            theirs.lines().skip(1).collect::<Vec<_>>()
        );
    }

    #[test]
    fn no_program_record_is_appended() {
        let header = header();
        let text = encode_without_existing_version(&header);
        assert_eq!(
            text.lines().filter(|line| line.starts_with("@PG")).count(),
            1,
            "the file's own, and nothing else"
        );
        assert!(!text.contains(TOOL_NAME));
    }
}
