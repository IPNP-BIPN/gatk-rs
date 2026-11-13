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
    absolute, assert_readable, assert_writable, thrown, unroll_paths, vcf_text, write_vcf,
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

    let texts: Vec<String> = paths
        .iter()
        .map(|path| vcf_text(path))
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

    let block_copy = paths.iter().chain(std::iter::once(&output)).all(|path| {
        let lower = path.to_lowercase();
        !lower.ends_with(".bcf")
            && [".gz", ".gzip", ".bgz", ".bgzf"]
                .iter()
                .any(|extension| lower.ends_with(extension))
    });
    let gathered = if block_copy {
        checked_order(&inputs, &args).and_then(|order| {
            block_copied(parser, &order).map_err(|error| Failure::Exit {
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

/// `gatherWithBlockCopying`: the output's bytes, which the caller writes. `Some` always; the
/// `Option` is the conventional path's, which writes its own file.
fn block_copied(parser: &Parser, paths: &[String]) -> Result<Option<Vec<u8>>, Thrown> {
    use gatk_tools::check_terminator_block::check;
    use htsjdk_bgzf::termination::FileTermination;

    let mut out = Vec::new();
    for (index, path) in paths.iter().enumerate() {
        let bytes = std::fs::read(path)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
        let (termination, _) = check(&bytes);
        if termination == FileTermination::Defective {
            return Err(thrown(
                "picard.PicardException",
                format!(
                    "{} does not have a valid GZIP block at the end of the file.",
                    absolute(path)
                ),
            ));
        }
        let mut position = 0;
        if index > 0 {
            let mut last_byte_newline = true;
            while let Some(size) = block_size(&bytes, position) {
                let block = htsjdk_bgzf::read::decompress_all(&bytes[position..position + size])
                    .map_err(|error| {
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
                    std::io::Write::write_all(&mut writer, &block[at..])
                        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
                    out.extend(
                        writer
                            .into_inner_without_terminator()
                            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?,
                    );
                    break;
                }
            }
        }
        let skip = if termination == FileTermination::HasTerminatorBlock {
            htsjdk_bgzf::EMPTY_GZIP_BLOCK.len()
        } else {
            0
        };
        let end = bytes.len().saturating_sub(skip);
        if position < end {
            out.extend_from_slice(&bytes[position..end]);
        }
    }
    out.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    Ok(Some(out))
}

/// The length of the BGZF block at `at`, from its `BSIZE`.
fn block_size(bytes: &[u8], at: usize) -> Option<usize> {
    let header = bytes.get(at..at + 18)?;
    let size = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
    (at + size <= bytes.len()).then_some(size)
}
