//! `GatherVcfs.doWork`, over [`gatk_tools::picard_gather_vcfs`].
//!
//! Everything after the dictionary check runs inside a `try` that catches any `RuntimeException`:
//! the tool logs it, deletes the output and returns 1, which GATK turns into status 4. What it
//! logs is a Picard `Log` line, and its timestamp is the run's own.
//!
//! When every input and the output carry a block-compressed extension the records are not read
//! at all: the first file's blocks are copied, each later file's header is cut out of the block it
//! ends in, the rest of that block is recompressed as one new block, and the remaining blocks are
//! copied as they are. No index is written on that path, whatever `CREATE_INDEX` says.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};
use gatk_tools::picard_gather_vcfs::{
    assert_same_dictionary, checked_order, dictionary, gather, Arguments, Failure, Input,
};

use super::{
    absolute, assert_readable, assert_writable, thrown, unroll_paths, vcf_prefix, vcf_text,
    write_vcf,
};
use crate::runners::{argument, arguments, flag, output_compression};

pub(super) fn run(parser: &Parser, stderr: &mut String) -> Result<i32, Thrown> {
    stderr.push_str(&log("INFO", "Checking inputs."));
    let paths = unroll_paths(&arguments(parser, "INPUT"))?;
    for path in &paths {
        assert_readable(path)?;
    }
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_writable(&output)?;

    let block_copy = paths.iter().chain(std::iter::once(&output)).all(|path| {
        let lower = path.to_lowercase();
        !lower.ends_with(".bcf")
            && [".gz", ".gzip", ".bgz", ".bgzf"]
                .iter()
                .any(|extension| lower.ends_with(extension))
    });
    // Block copying reads each file's first record, which `VCFFileReader` decodes with the one
    // after it; the records past those are copied as bytes and never decoded.
    let texts: Vec<String> = paths
        .iter()
        .map(|path| {
            if block_copy {
                vcf_prefix(path, 2)
            } else {
                vcf_text(path)
            }
        })
        .collect::<Result<_, _>>()?;
    let inputs: Vec<Input> = paths
        .iter()
        .zip(&texts)
        .map(|(path, text)| Input { path, text })
        .collect();
    let args = Arguments {
        comments: arguments(parser, "COMMENT"),
        create_index: flag(parser, "CREATE_INDEX"),
        reorder_input_by_first_variant: flag(parser, "REORDER_INPUT_BY_FIRST_VARIANT"),
    };

    let gathered = if block_copy {
        checked_order(&inputs, &args).and_then(|order| {
            block_copied(parser, &order, &output).map_err(|error| Failure::Exit {
                class: error.exception.to_string(),
                message: error.message.unwrap_or_default(),
            })
        })
    } else {
        gather(&inputs, &args).and_then(|text| {
            write_vcf(parser, &output, &text, flag(parser, "CREATE_INDEX"))
                .map(|()| None)
                .map_err(|error| Failure::Exit {
                    class: error.exception.to_string(),
                    message: error.message.unwrap_or_default(),
                })
        })
    };
    match gathered {
        Ok(Some(bytes)) => {
            std::fs::write(&output, bytes)
                .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
            Ok(0)
        }
        Ok(None) => Ok(0),
        // The dictionary check before the `try`, and an `AssertionError`, which is no
        // `RuntimeException` and escapes it.
        Err(Failure::Thrown { class, message }) => {
            if class == "java.lang.AssertionError" {
                // `log.error("File #1: ...")` and `"File #2: ..."` before the rethrow.
                // shortcut: the culprit is looked for in the order given, which is the gathered
                // order unless REORDER_INPUT_BY_FIRST_VARIANT moved the files.
                let headers: Vec<_> = texts
                    .iter()
                    .map(|text| {
                        htsjdk_vcf::reader::read_vcf(text)
                            .ok()
                            .map(|file| file.header)
                    })
                    .collect();
                let first = headers[0].as_ref().and_then(dictionary);
                let culprit = headers.iter().position(|header| {
                    match (&first, header.as_ref().and_then(dictionary)) {
                        (Some(first), Some(found)) => {
                            assert_same_dictionary(first, &found).is_err()
                        }
                        _ => false,
                    }
                });
                stderr.push_str(&log("ERROR", &format!("File #1: {}", paths[0])));
                if let Some(culprit) = culprit {
                    stderr.push_str(&log("ERROR", &format!("File #2: {}", paths[culprit])));
                }
            }
            Err(thrown(&class, message))
        }
        Err(Failure::Exit { class, message }) => {
            stderr.push_str(&log(
                "ERROR",
                &format!("There was a problem with gathering the INPUT.{class}: {message}"),
            ));
            let _ = std::fs::remove_file(&output);
            Ok(1)
        }
    }
}

fn log(level: &str, message: &str) -> String {
    super::log_line(level, "GatherVcfs", message)
}

/// `gatherWithBlockCopying`. The output is written to `<output>.partial` as it is gathered and
/// renamed onto `output` at the end, so a refusal leaves none: `Ok(None)`. When that file cannot
/// be created the bytes are gathered in memory and returned for the caller to write, as they
/// always were, so a refusal from the inputs still comes before the output's own failure.
fn block_copied(
    parser: &Parser,
    paths: &[String],
    output: &str,
) -> Result<Option<Vec<u8>>, Thrown> {
    let partial = PartialFile(format!("{output}.partial"));
    let Ok(file) = std::fs::File::create(&partial.0) else {
        let mut out = Vec::new();
        gather_blocks(parser, paths, &mut out)?;
        return Ok(Some(out));
    };
    let failed =
        |error: std::io::Error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}"));
    let mut out = std::io::BufWriter::new(file);
    gather_blocks(parser, paths, &mut out)?;
    out.into_inner()
        .map_err(|error| failed(error.into_error()))?;
    std::fs::rename(&partial.0, output).map_err(failed)?;
    Ok(None)
}

/// A file removed when it is dropped, which after the rename onto the output is a no-op.
struct PartialFile(String);

impl Drop for PartialFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// The gathered bytes into `out`, each input read a block at a time.
fn gather_blocks(
    parser: &Parser,
    paths: &[String],
    out: &mut dyn std::io::Write,
) -> Result<(), Thrown> {
    use gatk_tools::check_terminator_block::check;
    use htsjdk_bgzf::termination::FileTermination;
    use std::io::{Read, Seek};

    let unreadable = |path: &str, error: std::io::Error| {
        Thrown::non_user(PORT_FAILURE, format!("{path}: {error}"))
    };
    let unwritable = |error: std::io::Error| Thrown::non_user(PORT_FAILURE, format!("{error}"));
    for (index, path) in paths.iter().enumerate() {
        let mut file = std::fs::File::open(path).map_err(|error| unreadable(path, error))?;
        let length = file
            .metadata()
            .map_err(|error| unreadable(path, error))?
            .len() as usize;
        // `check_termination` looks at no more than the last block's worth of the file.
        let tail = length.min(htsjdk_bgzf::MAX_COMPRESSED_BLOCK_SIZE);
        let mut end_bytes = vec![0; tail];
        file.seek(std::io::SeekFrom::Start((length - tail) as u64))
            .and_then(|_| file.read_exact(&mut end_bytes))
            .and_then(|()| file.seek(std::io::SeekFrom::Start(0)).map(|_| ()))
            .map_err(|error| unreadable(path, error))?;
        let (termination, _) = check(&end_bytes);
        if termination == FileTermination::Defective {
            return Err(thrown(
                "picard.PicardException",
                format!(
                    "{} does not have a valid GZIP block at the end of the file.",
                    absolute(path)
                ),
            ));
        }
        let mut file = std::io::BufReader::new(file);
        let mut position = 0;
        if index > 0 {
            let mut last_byte_newline = true;
            let mut block = Vec::new();
            while let Some(size) = next_block(&mut file, position, length, &mut block)
                .map_err(|error| unreadable(path, error))?
            {
                let block = htsjdk_bgzf::read::decompress_all(&block).map_err(|error| {
                    Thrown::non_user(PORT_FAILURE, format!("{path}: {error:?}"))
                })?;
                position += size;
                if block.is_empty() {
                    break;
                }
                let mut first_record = None;
                for (at, &byte) in block.iter().enumerate() {
                    let newline = byte == b'\n' || byte == b'\r';
                    if last_byte_newline && !newline && byte != b'#' {
                        first_record = Some(at);
                        break;
                    }
                    last_byte_newline = newline;
                }
                if let Some(at) = first_record {
                    // `new BlockCompressedOutputStream(out, null)`, written and flushed, never
                    // closed: the default level and deflater, and no terminator.
                    let (level, deflater) = output_compression(parser);
                    let mut writer =
                        htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
                    std::io::Write::write_all(&mut writer, &block[at..]).map_err(unwritable)?;
                    out.write_all(&writer.into_inner_without_terminator().map_err(unwritable)?)
                        .map_err(unwritable)?;
                    break;
                }
            }
        }
        let skip = if termination == FileTermination::HasTerminatorBlock {
            htsjdk_bgzf::EMPTY_GZIP_BLOCK.len()
        } else {
            0
        };
        let end = length.saturating_sub(skip);
        if position < end {
            let copied = std::io::copy(&mut (&mut file).take((end - position) as u64), out)
                .map_err(unwritable)?;
            if copied != (end - position) as u64 {
                return Err(unreadable(
                    path,
                    std::io::Error::from(std::io::ErrorKind::UnexpectedEof),
                ));
            }
        }
    }
    out.write_all(&htsjdk_bgzf::EMPTY_GZIP_BLOCK)
        .map_err(unwritable)
}

/// The BGZF block at `at` of a file `length` bytes long into `block`, and its length, by its
/// `BSIZE`: `None` when the file has no whole block there, which leaves the reader at `at`.
fn next_block(
    file: &mut std::io::BufReader<std::fs::File>,
    at: usize,
    length: usize,
    block: &mut Vec<u8>,
) -> std::io::Result<Option<usize>> {
    if at + 18 > length {
        return Ok(None);
    }
    block.resize(18, 0);
    std::io::Read::read_exact(file, block)?;
    let size = usize::from(u16::from_le_bytes([block[16], block[17]])) + 1;
    if at + size > length {
        file.seek_relative(-18)?;
        return Ok(None);
    }
    block.resize(size, 0);
    std::io::Read::read_exact(file, &mut block[18..])?;
    Ok(Some(size))
}
