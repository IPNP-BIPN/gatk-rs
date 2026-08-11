//! Ported from `org.broadinstitute.hellbender.tools.ConvertHeaderlessHadoopBamShardToBam` and the
//! two `SparkUtils` methods it calls (GATK 4.6.2.0).
//!
//! The thirteenth whole tool of the record-transform archetype, and the first that is not a
//! `GATKTool` at all: it extends `CommandLineProgram` and implements `doWork()`. No traversal, no
//! read filter, no `@PG` record, no engine.
//!
//! ```java
//! try ( FileOutputStream outStream = new FileOutputStream(destination) ) {
//!     writeBAMHeaderToStream(header, outStream);
//!     FileUtils.copyFile(bamShard, outStream);
//!     outStream.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
//! }
//! ```
//!
//! # The shard is copied, not read
//!
//! `FileUtils.copyFile` between a header block and a terminator, so the output's data blocks are the
//! shard's **own bytes**: not decompressed, not re-deflated, not re-blocked. A port that read the
//! records and wrote them back would produce a valid BAM with different bytes, and every read-level
//! assertion would still pass. That is why the golden measures a byte search and a three-part
//! layout rather than a round trip.
//!
//! # This is the one path where the version is kept
//!
//! `writeBAMHeaderToStream` calls `new SAMTextHeaderCodec().encode(stringWriter, header, true)`,
//! which is [`htsjdk_bam::header::SamHeader::encode`] exactly. [`crate::print_reads_header`] takes
//! the other branch, and htsjdk-rs#164 is about the ordinary BAM writer taking the other branch too.
//! Here `true` is correct, so this module calls [`htsjdk_bam::writer::write_bam_header_block`],
//! which is the port of this very method.
//!
//! # The header block carries no terminator
//!
//! `blockCompressedOutputStream.flush()` rather than `close()`, so the empty gzip block appears
//! exactly once, at the very end, after the copied shard. Getting that wrong produces a file that
//! most readers stop at halfway through.
//!
//! # A donor that is not a BAM is not refused
//!
//! `SamReaderFactory.makeDefault().validationStringency(SILENT)` on a headerless shard yields an
//! **empty** `SAMFileHeader`, and the tool writes a BAM whose whole header is `@HD VN:1.6` with the
//! shard appended. Garbage in, valid BAM out. That is the reference's behaviour and this port does
//! not improve on it: the caller decides what header to hand over.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::writer::write_bam_header_block;

use gatk_engine::reads::ReadsError;

/// `GATKTool.getToolName()` for this tool, which never reaches an output because there is no `@PG`.
pub const TOOL_NAME: &str = "GATK ConvertHeaderlessHadoopBamShardToBam";

/// `BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK`: the 28 bytes that end a BAM.
pub const EMPTY_GZIP_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// `SparkUtils.convertHeaderlessHadoopBamShardToBam`: header block, shard verbatim, terminator.
///
/// `shard` is the raw file, already BGZF-compressed and carrying no header block and no terminator
/// of its own. It is appended without being looked at, which is the whole point.
pub fn convert_headerless_shard(shard: &[u8], header: &SamHeader) -> Result<Vec<u8>, ReadsError> {
    let mut out = write_bam_header_block(header).map_err(|e| ReadsError::Io(e.to_string()))?;
    out.extend_from_slice(shard);
    out.extend_from_slice(&EMPTY_GZIP_BLOCK);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::SequenceRecord;

    fn header() -> SamHeader {
        let mut header = SamHeader::default();
        header.set_sort_order("coordinate");
        header.sequences.push(SequenceRecord::new("chr1", 1000));
        header
    }

    #[test]
    fn the_shard_appears_verbatim_between_the_header_and_the_terminator() {
        let shard: Vec<u8> = (0u8..91).collect();
        let out = convert_headerless_shard(&shard, &header()).unwrap();
        let header_block = write_bam_header_block(&header()).unwrap();

        assert_eq!(&out[..header_block.len()], &header_block[..]);
        assert_eq!(
            &out[header_block.len()..header_block.len() + shard.len()],
            &shard[..],
            "the shard is copied, not re-encoded"
        );
        assert_eq!(&out[out.len() - 28..], &EMPTY_GZIP_BLOCK);
        assert_eq!(out.len(), header_block.len() + shard.len() + 28);
    }

    #[test]
    fn an_empty_shard_leaves_a_header_and_a_terminator() {
        let out = convert_headerless_shard(&[], &header()).unwrap();
        let header_block = write_bam_header_block(&header()).unwrap();
        assert_eq!(out.len(), header_block.len() + 28);
    }

    #[test]
    fn the_terminator_appears_exactly_once() {
        let shard: Vec<u8> = (0u8..91).collect();
        let out = convert_headerless_shard(&shard, &header()).unwrap();
        let occurrences = out
            .windows(EMPTY_GZIP_BLOCK.len())
            .filter(|w| *w == EMPTY_GZIP_BLOCK)
            .count();
        assert_eq!(
            occurrences, 1,
            "the header block is flushed rather than closed"
        );
    }

    #[test]
    fn an_empty_header_is_written_rather_than_refused() {
        // What the reference produces when the donor is not a BAM: `@HD VN:1.6` and nothing else.
        let empty = SamHeader::default();
        assert_eq!(empty.encode(), "@HD\tVN:1.6\n");
        let out = convert_headerless_shard(&[1, 2, 3], &empty).unwrap();
        assert!(out.starts_with(b"\x1f\x8b"), "still a BGZF file");
    }
}
