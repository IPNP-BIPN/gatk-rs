//! `AddCommentsToBam.doWork`, over [`gatk_tools::add_comments_to_bam`].
//!
//! `BamFileIoUtils.reheaderBamFile` asserts the INPUT writable before it copies a block, because
//! the copy is written for the in-place case too: a corpus mounted read-only refuses every BAM.
//!
//! ```java
//! IOUtil.assertFileIsReadable(inputFile);
//! IOUtil.assertFileIsWritable(inputFile);
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `?` returns the error to the caller at once, Java's `throws`; `.map_err(f)` rewrites the
//!   error first, and `|error| ...` is a closure, Java's lambda;
//! - `&[u8]` is a borrowed, read-only view of a byte array; `bam[a..b]` is the half-open slice
//!   from `a` up to but not including `b`, and `bam[..end]` starts at zero;
//! - `(u32, htsjdk_bgzf::Deflater)` is a pair, the compression level and the deflater; a
//!   parameter written `(level, deflater): (...)` unpacks it on the way in;
//! - `while let Some(block) = ...` loops until the reader has no block left.

use gatk_barclay::Parser;
use gatk_tools::add_comments_to_bam::{add_comments, is_refused_by_name, CommentError};
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};

use super::{absolute, assert_readable, assert_writable, thrown};
use crate::runners::{argument, arguments, flag, is_coordinate_sorted, output_compression};

/// `AddCommentsToBam.doWork` and the `BamFileIoUtils.reheaderBamFile` it ends with.
///
/// The refusals come in the reference's order: INPUT readable, OUTPUT writable, a `.sam` name,
/// a comment with a newline, then `reheaderBamFile`'s own assertions on the INPUT. Picard tests
/// the name with `INPUT.getAbsolutePath().endsWith(".sam")`, which is why the absolute path is
/// what `is_refused_by_name` is handed.
pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    let comments = arguments(parser, "COMMENT");
    assert_readable(&input)?;
    assert_writable(&output)?;
    // A closure, so every place that raises a `CommentError` name its class the same way.
    let refuse = |error: CommentError| thrown(error.java_class(), error.message());
    if is_refused_by_name(&absolute(&input)) {
        return Err(refuse(CommentError::SamNotSupported));
    }
    // `if (comment.contains("\n"))`, checked over every comment before any header is written.
    if comments.iter().any(|comment| comment.contains('\n')) {
        return Err(refuse(CommentError::ContainsNewline));
    }
    // `reheaderBamFile`'s own two assertions, the second of them on the INPUT.
    assert_readable(&input)?;
    assert_writable(&input)?;
    let bam = std::fs::read(&input)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{input}: {error}")))?;
    // The library answers every refusal and the block copy; its two recompressed stretches are
    // htsjdk's defaults, and `reheaderBamFile` writes them with the run's level and deflater.
    add_comments(&bam, &absolute(&input), &comments).map_err(refuse)?;
    let written = reheader(&bam, &comments, output_compression(parser))
        .map_err(|error| Thrown::non_user(PORT_FAILURE, error))?;
    std::fs::write(&output, &written)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    // shortcut: CREATE_MD5_FILE and CREATE_INDEX are only reachable on a writable input, which the
    // read-only corpus never is; the index is the output's whole-file one, named `<output>.bai`.
    if flag(parser, "CREATE_MD5_FILE") {
        std::fs::write(
            format!("{output}.md5"),
            gatk_tools::gather_bam_files::md5_file(&written),
        )
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}.md5: {error}")))?;
    }
    if flag(parser, "CREATE_INDEX") && is_coordinate_sorted(&written) {
        if let Ok(bai) = htsjdk_bam::build_index::build_bam_index(&written) {
            std::fs::write(format!("{output}.bai"), bai).map_err(|error| {
                Thrown::non_user(PORT_FAILURE, format!("{output}.bai: {error}"))
            })?;
        }
    }
    Ok(0)
}

/// `BamFileIoUtils.reheaderBamFile`: the new header compressed as the run compresses, then the
/// input's records block-copied.
fn reheader(
    bam: &[u8],
    comments: &[String],
    compression: (u32, htsjdk_bgzf::Deflater),
) -> Result<Vec<u8>, String> {
    // `getFileHeader(INPUT)` then `samFileHeader.addComment(comment)` for each, in order.
    let plain = htsjdk_bgzf::read::decompress_all(bam).map_err(|e| format!("{e:?}"))?;
    let reader = htsjdk_bam::reader::BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
    let mut header = reader.header.text.clone();
    for comment in comments {
        header.add_comment(comment);
    }
    // `BAMFileWriter.writeHeader(outputStream, samFileHeader)`: the header's bytes are built
    // once, then decompressed and recompressed so they carry the run's level and deflater.
    let block =
        htsjdk_bam::writer::write_bam_header_block(&header).map_err(|e| format!("{e:?}"))?;
    let header_plain = htsjdk_bgzf::read::decompress_all(&block).map_err(|e| format!("{e:?}"))?;
    let mut out = compress(&header_plain, compression)?;
    // `blockCopyBamFile(inputFile, outputStream, true, false)`: skip the header, keep the
    // terminator.
    out.extend_from_slice(&block_copy(bam, true, false, compression)?);
    Ok(out)
}

/// One flushed run of BGZF blocks, compressed with the run's level and deflater.
///
/// What `new BlockCompressedOutputStream(outputStream, (Path) null)` then `blockOut.flush()`
/// writes: the data's blocks and no terminator, because the underlying stream stays open.
pub(super) fn compress(
    data: &[u8],
    (level, deflater): (u32, htsjdk_bgzf::Deflater),
) -> Result<Vec<u8>, String> {
    let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
    std::io::Write::write_all(&mut writer, data).map_err(|e| format!("{e:?}"))?;
    writer
        .into_inner_without_terminator()
        .map_err(|e| format!("{e:?}"))
}

/// `BamFileIoUtils.blockCopyBamFile` past the termination check: without its header, the rest of
/// the block the first record starts in is recompressed with the run's setting (htsjdk-rs's own
/// `block_copy` uses htsjdk's defaults), and every later block is copied untouched.
///
/// ```java
/// final long skipLast =
/// ((term == BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK) && skipTerminator)
///         ? BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length
///         : 0;
/// ```
pub(super) fn block_copy(
    bam: &[u8],
    skip_header: bool,
    skip_terminator: bool,
    compression: (u32, htsjdk_bgzf::Deflater),
) -> Result<Vec<u8>, String> {
    let terminated =
        htsjdk_bgzf::check_termination(bam) == htsjdk_bgzf::FileTermination::HasTerminatorBlock;
    let end = if terminated && skip_terminator {
        bam.len() - htsjdk_bgzf::EMPTY_GZIP_BLOCK.len()
    } else {
        bam.len()
    };
    if !skip_header {
        return Ok(bam[..end].to_vec());
    }
    let plain = htsjdk_bgzf::read::decompress_all(bam).map_err(|e| format!("{e:?}"))?;
    let first_record = htsjdk_bam::reader::BamReader::new(&plain)
        .map_err(|e| format!("{e:?}"))?
        .records_start();
    let mut out = Vec::new();
    // The block the first record starts in; a record at a block's very end starts the next one.
    // `first_record` is an offset into the decompressed stream, `start` the decompressed offset
    // where the current block begins.
    let mut blocks = htsjdk_bgzf::BgzfReader::new(bam);
    let mut start = 0;
    while let Some(block) = blocks.next_block().map_err(|e| format!("{e:?}"))? {
        let block_end = start + block.data.len();
        if first_record < block_end {
            out.extend_from_slice(&compress(&block.data[first_record - start..], compression)?);
            // `in.seek(pos)` at the next block's file offset, then the raw bytes to `end`.
            // `as usize` widens the `u64` file offset to an index.
            let copy_from = block.block_address as usize + block.block_compressed_size;
            if copy_from < end {
                out.extend_from_slice(&bam[copy_from..end]);
            }
            return Ok(out);
        }
        start = block_end;
    }
    // No record at all: what follows the header is only the terminator, if it is kept.
    if !skip_terminator && terminated {
        out.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    }
    Ok(out)
}
