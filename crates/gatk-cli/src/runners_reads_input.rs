//! `SamReaderFactory.open` for the tools that open their reads themselves rather than through
//! [`super::read_walker_startup`]: a BAM is opened as it is, and SAM text and CRAM are decoded a
//! record at a time into a scratch BAM, which the source then reads like any other. Nothing holds
//! the file: SAM text is read a line at a time and CRAM a container at a time.

use super::*;
use std::io::{BufRead, Read};
use std::path::Path;

/// The reads at `path` as a [`ReadsDataSource`], `index` being what the caller found beside them.
///
/// For a BAM, `index` is its `.bai`, opened as [`ReadsDataSource::open`] does; for a CRAM, a
/// `.crai` beside it, and the source is given the `.bai` built while transcoding instead, so an
/// interval query answers the same. SAM text has no index. Text that does not parse as SAM, and
/// anything else that is not CRAM, goes to the BAM open, so its refusal is the one it was.
pub(crate) fn open_reads(
    parser: &Parser,
    path: &Path,
    index: Option<&Path>,
) -> Result<ReadsDataSource, Thrown> {
    let (source, pending) = open_reads_deferred(parser, path, index)?;
    pending.map_or(Ok(source), Err)
}

/// [`open_reads`] with a CRAM's decode refusal held back rather than raised: htsjdk decodes a
/// container when the traversal reaches it, so a tool's checks before its traversal come first.
/// The source holds the records before the refusal; the caller raises it where it traverses.
pub(crate) fn open_reads_deferred(
    parser: &Parser,
    path: &Path,
    index: Option<&Path>,
) -> Result<(ReadsDataSource, Option<Thrown>), Thrown> {
    open(path, index, || {
        argument(parser, "reference").map(Some).ok_or_else(|| {
            Thrown::user(
                "A reference file is required when using CRAM files.  A reference is specified \
                 with the -R command line argument.",
            )
        })
    })
}

/// [`open_reads`] for a tool that opens its reads with a bare `new ReadsPathDataSource(path)`,
/// which is given no reference: a CRAM decodes until a slice needs one, and is refused there.
pub(crate) fn open_reads_without_reference(path: &Path) -> Result<ReadsDataSource, Thrown> {
    let (source, pending) = open(path, None, || Ok(None))?;
    pending.map_or(Ok(source), Err)
}

fn open(
    path: &Path,
    index: Option<&Path>,
    cram_reference: impl FnOnce() -> Result<Option<String>, Thrown>,
) -> Result<(ReadsDataSource, Option<Thrown>), Thrown> {
    let as_bam = || {
        match index {
            Some(index) => ReadsDataSource::open(path, index),
            None => ReadsDataSource::open_unindexed(path),
        }
        .map(|source| (source, None))
        .map_err(|error| Thrown::user(format!("{error:?}")))
    };
    let mut magic = [0u8; 4];
    let read = std::fs::File::open(path)
        .and_then(|mut file| file.read(&mut magic))
        .unwrap_or(0);
    let magic = &magic[..read];
    if magic.is_empty()
        || magic == gatk_tools::read_walker_refusal::BAM_MAGIC
        || gatk_tools::read_walker_refusal::is_block_compressed(magic)
    {
        return as_bam();
    }
    let transcoded = if magic == b"CRAM" {
        transcode_cram(path, cram_reference()?.as_deref(), index.is_some())?
    } else {
        match transcode_sam(path)? {
            Some(transcoded) => transcoded,
            None => return as_bam(),
        }
    };
    let (scratch, header, bai, pending) = transcoded;
    let opened = ReadsDataSource::open_unindexed(&scratch.0);
    drop(scratch);
    let mut source =
        opened.map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
    if let Some(bai) = bai {
        source
            .load_index_bytes(&bai)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
    }
    // The file's own header, not the one the BAM writer stamps a sort order onto.
    Ok((source.with_header(header), pending))
}

/// Whether `path` is SAM text: not a BAM, a CRAM or block compressed. htsjdk cannot query one,
/// and it has no index whatever sits beside it.
pub(crate) fn is_sam_text(path: &Path) -> bool {
    let mut magic = [0u8; 4];
    let read = std::fs::File::open(path)
        .and_then(|mut file| file.read(&mut magic))
        .unwrap_or(0);
    let magic = &magic[..read];
    !magic.is_empty()
        && magic != gatk_tools::read_walker_refusal::BAM_MAGIC
        && magic != b"CRAM"
        && !gatk_tools::read_walker_refusal::is_block_compressed(magic)
}

/// `ReadsPathDataSource.isQueryableByInterval` for one file: an index found beside it, which SAM
/// text never has.
pub(crate) fn is_indexed(path: &Path) -> bool {
    !is_sam_text(path) && htsjdk_bam::sam_files::find_index(path).is_some()
}

/// `GATKTool.validateSequenceDictionaries` for reads, a reference and a master dictionary, which
/// the argument turns off wholesale. A CRAM input asks more of both: the master and the reference
/// must hold every contig the reads do.
pub(crate) fn validate_reads_dictionaries(
    parser: &Parser,
    master: Option<&SamHeader>,
    reference: Option<&SamHeader>,
    reads: &SamHeader,
    cram: bool,
) -> Result<(), Thrown> {
    if flag(parser, "disable-sequence-dictionary-validation") {
        return Ok(());
    }
    let check = |name1: &str, dict1: &SamHeader, name2: &str, dict2: &SamHeader, superset: bool| {
        gatk_tools::sequence_dictionary::validate(
            name1,
            &dict1.sequences,
            name2,
            &dict2.sequences,
            superset,
            false,
        )
        .map_err(|refusal| Thrown {
            failure: Failure::User,
            exception: refusal.java_class(),
            message: Some(refusal.message()),
        })
    };
    if let Some(master) = master {
        check("master sequence dictionary", master, "reads", reads, cram)?;
        if let Some(reference) = reference {
            check(
                "master sequence dictionary",
                master,
                "reference",
                reference,
                cram,
            )?;
        }
    }
    if let Some(reference) = reference {
        check("reference", reference, "reads", reads, cram)?;
    }
    Ok(())
}

/// A scratch BAM in the temporary directory, removed when dropped. The source has opened it by
/// then, and holds what it needs of it.
struct Scratch(std::path::PathBuf);

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// The scratch BAM, the file's header, the `.bai` asked for, and a refusal held back.
type Transcoded = (Scratch, SamHeader, Option<Vec<u8>>, Option<Thrown>);

/// The writer of the scratch BAM. Level one: it is read once, right away, so time beats size.
fn scratch_writer(
    header: &SamHeader,
    indexed: bool,
) -> Result<
    (
        Scratch,
        htsjdk_bam::writer::BamWriter<std::io::BufWriter<std::fs::File>>,
    ),
    Thrown,
> {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let failure = |error: std::io::Error| Thrown::non_user(PORT_FAILURE, error.to_string());
    let scratch = Scratch(std::env::temp_dir().join(format!(
        ".gatk-rs-reads.{}.{}.bam",
        std::process::id(),
        NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    )));
    let file = std::io::BufWriter::with_capacity(
        1 << 20,
        std::fs::File::create(&scratch.0).map_err(failure)?,
    );
    let writer = htsjdk_bam::writer::BamWriter::with_compression(
        file,
        header,
        1,
        htsjdk_bgzf::Deflater::Gkl,
    )
    .map_err(failure)?;
    Ok((scratch, if indexed { writer.with_index() } else { writer }))
}

/// Close the scratch BAM, with its `.bai` when one was asked for.
fn finish(
    writer: htsjdk_bam::writer::BamWriter<std::io::BufWriter<std::fs::File>>,
    indexed: bool,
) -> Result<Option<Vec<u8>>, Thrown> {
    use std::io::Write;
    let failure = |error: std::io::Error| Thrown::non_user(PORT_FAILURE, error.to_string());
    let (mut file, bai) = if indexed {
        let (file, bai) = writer.finish_with_index().map_err(failure)?;
        (file, Some(bai))
    } else {
        (writer.finish().map_err(failure)?, None)
    };
    file.flush().map_err(failure)?;
    Ok(bai)
}

fn write_record(
    writer: &mut htsjdk_bam::writer::BamWriter<std::io::BufWriter<std::fs::File>>,
    record: &htsjdk_bam::record::BamRecord,
) -> Result<(), Thrown> {
    writer
        .write(record)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))
}

/// A CRAM decoded against `reference` one container at a time.
fn transcode_cram(
    path: &Path,
    reference: Option<&str>,
    indexed: bool,
) -> Result<Transcoded, Thrown> {
    let refusal = |error: htsjdk_cram::reader::CramError| {
        Thrown::non_user(error.java_exception(), error.message())
    };
    let file = std::fs::File::open(path).map_err(|error| Thrown::user(format!("{error}")))?;
    let mut reader = htsjdk_cram::reader::CramReader::from_read(std::io::BufReader::with_capacity(
        1 << 20,
        file,
    ))
    .map_err(refusal)?;
    let mut fasta: Box<dyn htsjdk_cram::reader::ReferenceSource> = match reference {
        Some(reference) => Box::new(htsjdk_cram::reader::FastaReference::open(Path::new(
            reference,
        ))),
        // `CRAMLazyReferenceSource`, which throws for the first sequence a slice needs.
        None => Box::new(htsjdk_cram::reader::LazyReference),
    };
    let header = reader.header.clone();
    let (scratch, mut writer) = scratch_writer(&header, indexed)?;
    let mut pending = None;
    loop {
        match reader.next_container(fasta.as_mut()) {
            Ok(Some(records)) => {
                for record in &records {
                    write_record(&mut writer, record)?;
                }
            }
            Ok(None) => break,
            Err(error) => {
                pending = Some(refusal(error));
                break;
            }
        }
    }
    let bai = finish(writer, indexed)?;
    Ok((scratch, header, bai, pending))
}

/// SAM text read a line at a time, as `read_sam` reads it whole: the header is the leading `@`
/// lines, blank lines are skipped. `None` for text that is not SAM.
fn transcode_sam(path: &Path) -> Result<Option<Transcoded>, Thrown> {
    let file = std::fs::File::open(path).map_err(|error| Thrown::user(format!("{error}")))?;
    let mut reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut line = Vec::new();
    let next = |reader: &mut std::io::BufReader<std::fs::File>, line: &mut Vec<u8>| {
        line.clear();
        reader.read_until(b'\n', line).is_ok()
    };
    let mut text = String::new();
    loop {
        if !next(&mut reader, &mut line) {
            return Ok(None);
        }
        if !line.starts_with(b"@") {
            break;
        }
        let Ok(header_line) = std::str::from_utf8(&line) else {
            return Ok(None);
        };
        text.push_str(header_line.trim_end_matches(['\n', '\r']));
        text.push('\n');
    }
    let header = htsjdk_bam::reader::parse_header_text(&text);
    let contigs: std::collections::HashMap<&str, i32> = header
        .sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.name.as_str(), index as i32))
        .collect();
    let (scratch, mut writer) = scratch_writer(&header, false)?;
    while !line.is_empty() {
        let Ok(record_line) = std::str::from_utf8(&line) else {
            return Ok(None);
        };
        // `str::lines`, which `read_sam` splits by, drops a trailing `\r`.
        let record_line = record_line.trim_end_matches('\n').trim_end_matches('\r');
        if !record_line.is_empty() {
            let Ok(record) =
                htsjdk_bam::text_parse::parse_line(record_line, |name| contigs.get(name).copied())
            else {
                return Ok(None);
            };
            write_record(&mut writer, &record)?;
        }
        if !next(&mut reader, &mut line) {
            return Ok(None);
        }
    }
    finish(writer, false)?;
    Ok(Some((scratch, header, None, None)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The line-at-a-time transcode reads the records and the header `read_sam` reads whole,
    /// blank lines and `\r\n` included.
    #[test]
    fn sam_text_transcodes_to_the_records_read_sam_reads() {
        let text = "@HD\tVN:1.6\tSO:coordinate\r\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:0\tSM:s\n\
                    a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:0\r\n\n\
                    b\t16\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:0\n\
                    c\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII";
        let path =
            std::env::temp_dir().join(format!("gatk-rs-reads-input-{}.sam", std::process::id()));
        std::fs::write(&path, text).unwrap();
        let (scratch, header, bai, _) = transcode_sam(&path).unwrap().expect("SAM text");
        let source = ReadsDataSource::open_unindexed(&scratch.0).unwrap();
        let _ = std::fs::remove_file(&path);
        let (expected_header, expected) =
            htsjdk_bam::sam_file::read_sam(&text.replace('\r', "")).unwrap();
        assert!(bai.is_none());
        assert_eq!(header, expected_header);
        assert_eq!(source.iter_all().unwrap(), expected);
    }

    /// Text that is not SAM is left to the BAM open, whose refusal it was.
    #[test]
    fn text_that_is_not_sam_is_not_transcoded() {
        let path =
            std::env::temp_dir().join(format!("gatk-rs-reads-input-{}.vcf", std::process::id()));
        std::fs::write(&path, "##fileformat=VCFv4.2\n#CHROM\tPOS\n").unwrap();
        let transcoded = transcode_sam(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        assert!(transcoded.is_none());
    }
}
