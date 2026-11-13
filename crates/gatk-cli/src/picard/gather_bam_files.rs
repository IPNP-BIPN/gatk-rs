//! `GatherBamFiles.doWork`.
//!
//! Two paths, chosen by the inputs' NAMES (`BamFileIoUtils.isBamFile` reads the extension): every
//! input a `.bam` is `gatherWithBlockCopying`, the first file's header and every file's blocks;
//! anything else is `gatherNormally`, every record written again under the first file's header,
//! in the format the OUTPUT's name asks for.
//!
//! The SAM/BAM reading and writing here is shared with `SetNmMdAndUqTags`.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_FAILURE, PORT_LIMITATION};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use super::add_comments_to_bam::block_copy;
use super::{absolute, assert_readable, assert_writable, thrown};
use crate::runners::{argument, arguments, flag, is_coordinate_sorted, output_compression};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let inputs = unroll(&arguments(parser, "INPUT"))?;
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    for input in &inputs {
        assert_readable(input)?;
    }
    assert_writable(&output)?;
    let create_index = flag(parser, "CREATE_INDEX");
    let create_md5 = flag(parser, "CREATE_MD5_FILE");

    if !inputs.iter().all(|input| input.ends_with(".bam")) {
        let (header, _) = read_reads(&inputs[0])?;
        let mut records = Vec::new();
        for input in &inputs {
            records.extend(read_reads(input)?.1);
        }
        return write_reads(parser, &output, &header, &records, create_index, create_md5)
            .map(|()| 0);
    }

    let compression = output_compression(parser);
    let mut out = Vec::new();
    for (index, input) in inputs.iter().enumerate() {
        let bam = std::fs::read(input)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{input}: {error}")))?;
        if htsjdk_bgzf::check_termination(&bam) == htsjdk_bgzf::FileTermination::Defective {
            return Err(thrown(
                "htsjdk.samtools.SAMException",
                format!(
                    "file://{} does not have a valid GZIP block at the end of the file.",
                    absolute(input)
                ),
            ));
        }
        out.extend_from_slice(
            &block_copy(&bam, index > 0, true, compression)
                .map_err(|error| Thrown::non_user(PORT_FAILURE, error))?,
        );
    }
    out.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    write(&output, &out)?;
    if create_md5 {
        write(
            &format!("{output}.md5"),
            gatk_tools::gather_bam_files::md5_file(&out).as_bytes(),
        )?;
    }
    // `StreamInflatingIndexingOutputStream`, which indexes what went through it, and only under a
    // coordinate header: its `BAMIndexer` refuses anything else in a thread nobody joins.
    // shortcut: inputs whose concatenation goes back to an earlier reference kill that thread
    // mid-index, leaving only the references it finished; this writes the whole index. Truncate it
    // the same way if an array ever gathers out-of-order shards with CREATE_INDEX.
    if create_index && is_coordinate_sorted(&out) {
        if let Ok(bai) = htsjdk_bam::build_index::build_bam_index(&out) {
            write(&index_path(&output), &bai)?;
        }
    }
    Ok(0)
}

/// `IOUtil.unrollFiles(inputs, ".bam", ".sam", ".cram")`: an input of any other name is a list
/// of paths, one per line.
fn unroll(inputs: &[String]) -> Result<Vec<String>, Thrown> {
    let mut out = Vec::new();
    for input in inputs {
        if [".bam", ".sam", ".cram"]
            .iter()
            .any(|ext| input.ends_with(ext))
        {
            out.push(input.clone());
            continue;
        }
        let text = std::fs::read_to_string(input).map_err(|_| {
            thrown(
                "java.lang.IllegalArgumentException",
                format!("had trouble reading from file://{}", absolute(input)),
            )
        })?;
        let listed: Vec<String> = text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_string)
            .collect();
        out.extend(unroll(&listed)?);
    }
    Ok(out)
}

/// `SamReaderFactory.open` over a BAM or a SAM text file: the header and every record.
pub(super) fn read_reads(path: &str) -> Result<(SamHeader, Vec<BamRecord>), Thrown> {
    let malformed = |error: String| thrown("htsjdk.samtools.SAMFormatException", error);
    let raw = std::fs::read(path)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
    if path.ends_with(".cram") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "A CRAM input to a Picard tool is not read by this port yet. This message is the port's own.",
        ));
    }
    if !raw.starts_with(&[0x1f, 0x8b]) {
        return htsjdk_bam::sam_file::read_sam(&String::from_utf8_lossy(&raw))
            .map_err(|error| malformed(format!("{error:?}")));
    }
    let plain = htsjdk_bgzf::read::decompress_all(&raw).map_err(|e| malformed(format!("{e:?}")))?;
    let reader =
        htsjdk_bam::reader::BamReader::new(&plain).map_err(|e| malformed(format!("{e:?}")))?;
    let header = reader.header.text.clone();
    let records = reader
        .map(|record| record.map_err(|e| malformed(format!("{e:?}"))))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((header, records))
}

/// `SAMFileWriterFactory.makeWriter(header, true, output, reference)`: SAM text for a `.sam`
/// name, a BAM for anything but a `.cram`, its `.bai` only when the header says `coordinate`.
pub(super) fn write_reads(
    parser: &Parser,
    output: &str,
    header: &SamHeader,
    records: &[BamRecord],
    create_index: bool,
    create_md5: bool,
) -> Result<(), Thrown> {
    // `SAMFileWriterImpl.assertPresorted`, which refuses at the first record out of the header's
    // order. The SAM text writer has no file name to print, the BAM one prints the output's URI.
    let name = if output.ends_with(".sam") {
        "null".to_string()
    } else {
        format!("file://{}", absolute(output))
    };
    assert_presorted(header, records, &name)?;
    let bytes = if output.ends_with(".sam") {
        htsjdk_bam::sam_file::write_sam(header, records)
            .ok_or_else(|| Thrown::non_user(PORT_FAILURE, "records failed to encode as SAM"))?
            .into_bytes()
    } else if output.ends_with(".cram") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "A CRAM output from a Picard tool is not written by this port yet. This message is the port's own.",
        ));
    } else {
        let (level, deflater) = output_compression(parser);
        let index = create_index && header.attributes.get("SO") == Some("coordinate");
        let (bam, bai) =
            gatk_tools::sam_output::write_records_with(header, records, index, level, deflater)
                .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
        if let Some(bai) = bai {
            write(&index_path(output), &bai)?;
        }
        bam
    };
    write(output, &bytes)?;
    if create_md5 {
        write(
            &format!("{output}.md5"),
            gatk_tools::gather_bam_files::md5_file(&bytes).as_bytes(),
        )?;
    }
    Ok(())
}

/// `SAMSortOrderChecker` over the records, for the two orders it has a comparator for.
fn assert_presorted(header: &SamHeader, records: &[BamRecord], name: &str) -> Result<(), Thrown> {
    use std::cmp::Ordering::Greater;
    let order = header.attributes.get("SO").unwrap_or("unsorted");
    if order != "coordinate" && order != "queryname" {
        return Ok(());
    }
    let sort_key = |record: &BamRecord| {
        if order == "queryname" {
            return record.read_name.clone();
        }
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
            .map_or("*", |sequence| sequence.name.as_str());
        format!("{contig}:{}", record.alignment_start)
    };
    for pair in records.windows(2) {
        let out_of_order = if order == "coordinate" {
            htsjdk_bam::coordinate::file_order_compare(&pair[0], &pair[1]) == Greater
        } else {
            htsjdk_bam::query_name::compare_read_names(&pair[0].read_name, &pair[1].read_name)
                == Greater
        };
        if out_of_order {
            return Err(thrown(
                "java.lang.IllegalArgumentException",
                format!(
                    "Alignments added out of order in SAMFileWriterImpl.addAlignment for {name}. \
                     Sort order is {order}. Offending records are at [{}] and [{}]",
                    // `isSorted` has already made the offending record the previous one, so the
                    // reference prints that record's key twice.
                    sort_key(&pair[1]),
                    sort_key(&pair[1])
                ),
            ));
        }
    }
    Ok(())
}

/// `IOUtil.basename(output) + ".bai"` beside the output: its last extension replaced.
fn index_path(output: &str) -> String {
    std::path::Path::new(output)
        .with_extension("bai")
        .display()
        .to_string()
}

fn write(path: &str, bytes: &[u8]) -> Result<(), Thrown> {
    std::fs::write(path, bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))
}
