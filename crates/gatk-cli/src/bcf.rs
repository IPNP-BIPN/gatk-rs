//! BCF as htsjdk 5.0.0 reads and writes it, for the tools here that carry a VCF as text.
//!
//! `VCFFileReader` reads a name ending `.bcf` with `BCF2Codec`, and `VariantContextWriterBuilder`
//! writes one with `BCF2Writer`; both are [`htsjdk_vcf::bcf2`]. The tools here work on VCF text, so
//! a BCF is read as the text `VCFWriter` writes of its records, and a text written to a `.bcf` name
//! is encoded record by record.
//!
//! Text is not enough to write a BCF from a BCF. The BCF2 codec's record has its INFO in the order
//! of the codec's own map, its floats as the 32 bits the file held (VCF text keeps three decimals),
//! and a genotype block `BCF2Writer` copies without decoding it. So each record read here is kept
//! beside the line it was rendered as, and an output line that is that line, over the same
//! samples, is written from the kept record, as Java writes the record the tool passed through.
//! A record read from a VCF is kept the same way, because its INFO reaches the BCF in the order of
//! the codec's map, which follows the file's order where two keys share a bucket, and the text a
//! tool writes has its INFO sorted. A VCF's records are only parsed for this when a BCF is written.
//! A line a tool changed is parsed back and written as a record from `VCFCodec`, which is the
//! `Origin::VcfCodec` the writer replays.
//!
//! shortcut: a record a tool changed is rebuilt from its text, so a BCF read, edited and written
//! as BCF loses its floats' digits past the third and takes `VCFCodec`'s INFO order; upgrade by
//! carrying records rather than text through that tool.

use std::collections::HashMap;
use std::sync::Mutex;

use gatk_tools::main_entry::{Thrown, PORT_FAILURE};
use htsjdk_tribble::index_write::Feature;
use htsjdk_vcf::bcf2::codec::read_bcf;
use htsjdk_vcf::bcf2::writer::{Bcf2Writer, Origin};
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::header::VcfHeader;
use htsjdk_vcf::variant::VariantContext;

use crate::runners::{java_absolute_path, modified_millis, tribble_on_the_fly};

/// `FileExtensions.BCF`, which both the reader and the writer builder match case-sensitively.
pub(crate) fn is_bcf(path: &str) -> bool {
    path.ends_with(".bcf")
}

/// A record read in this run: the samples its genotypes were decoded under, the record, and the
/// codec that read it.
struct Kept {
    samples: Vec<String>,
    record: VariantContext,
    origin: Origin,
}

/// Every record read in this run, by the VCF line it renders as.
static KEPT: Mutex<Option<HashMap<String, Vec<Kept>>>> = Mutex::new(None);

/// The VCFs read in this run whose records are not in [`KEPT`] yet.
static VCFS: Mutex<Vec<String>> = Mutex::new(Vec::new());

/// A VCF a tool reads as text, whose records a BCF written later may be.
pub(crate) fn note_vcf(path: &str) {
    VCFS.lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .push(path.to_string());
}

/// The records of every VCF noted since the last call, each under the line it renders as.
fn keep_noted(kept: &mut HashMap<String, Vec<Kept>>) {
    let paths = std::mem::take(&mut *VCFS.lock().unwrap_or_else(|poisoned| poisoned.into_inner()));
    for path in paths {
        let Ok(bytes) = std::fs::read(&path) else {
            continue;
        };
        let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
            htsjdk_bgzf::read::decompress_all(&bytes).unwrap_or_default()
        } else {
            bytes
        };
        // A file the tool could read is one this parses; one it could not wrote no BCF.
        let Ok(file) = htsjdk_vcf::reader::read_vcf(&String::from_utf8_lossy(&bytes)) else {
            continue;
        };
        let encoder = VcfEncoder::new(&file.header);
        for record in file.records {
            if let Ok(line) = encoder.encode(&record) {
                kept.entry(line).or_default().push(Kept {
                    samples: file.header.samples.clone(),
                    record,
                    origin: Origin::VcfCodec,
                });
            }
        }
    }
}

/// An exception as the library threw it, its class leaked: a handful per run, which ends with it.
fn thrown(class: &str, message: String) -> Thrown {
    crate::picard::thrown(class, message)
}

/// `new VCFFileReader(path)` on a BCF, read whole: the VCF text `VCFWriter` writes of its header
/// and records, or the exception as the reader surfaces it.
pub(crate) fn read_as_text(path: &str) -> Result<String, Thrown> {
    let file = read(path)?;
    let encoder = VcfEncoder::new(&file.header);
    let mut text = file.header.write();
    let mut kept = KEPT.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let kept = kept.get_or_insert_with(HashMap::new);
    for record in file.records {
        let line = encoder
            .encode(&record)
            .map_err(|error| thrown("java.lang.IllegalStateException", format!("{error:?}")))?;
        text.push_str(&line);
        text.push('\n');
        kept.entry(line).or_default().push(Kept {
            samples: file.header.samples.clone(),
            record,
            origin: Origin::Bcf2Codec,
        });
    }
    Ok(text)
}

/// `new VCFFileReader(path)` on a BCF, read whole, or the exception as the reader surfaces it.
pub(crate) fn read(path: &str) -> Result<htsjdk_vcf::bcf2::codec::Bcf2File, Thrown> {
    let bytes = std::fs::read(path)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
    read_bcf(&bytes).map_err(|failure| {
        let error = failure.reader_exception(&format!("file://{}", java_absolute_path(path)));
        thrown(error.class, error.message)
    })
}

/// `BCF2Writer` over a VCF text: the BCF and its `.idx` when `index` gives the dictionary, and its
/// `.md5` when asked.
pub(crate) fn write_text(
    output: &str,
    text: &str,
    index: Option<&[(String, i32)]>,
    md5: bool,
) -> Result<(), Thrown> {
    let (header, records) = records_of(text, true)?;
    write(output, &header, &records, index, md5)
}

/// Where each record of a VCF text starts in the BCF written from it, and where that BCF ends:
/// the positions `IndexingVariantContextWriter` hands an index creator. A record the writer
/// refuses ends the list.
pub(crate) fn positions(text: &str) -> Result<(Vec<i64>, i64), Thrown> {
    let (header, records) = records_of(text, false)?;
    let Ok(mut writer) = Bcf2Writer::new(&header, false) else {
        return Ok((Vec::new(), 0));
    };
    let mut at = writer.header_bytes().len() as i64;
    let mut starts = Vec::with_capacity(records.len());
    for (record, origin) in &records {
        starts.push(at);
        match writer.encode(record, *origin) {
            Ok(bytes) => at += bytes.len() as i64,
            Err(_) => break,
        }
    }
    Ok((starts, at))
}

/// The header and records a VCF text is written from: each line that renders a record read from
/// a BCF over the same samples is that record, which `take` uses up; any other is parsed back.
fn records_of(
    text: &str,
    take: bool,
) -> Result<(VcfHeader, Vec<(VariantContext, Origin)>), Thrown> {
    let file = htsjdk_vcf::reader::read_vcf(text)
        .map_err(|failure| thrown(failure.error.class(), failure.error.message()))?;
    let lines = text
        .split_inclusive('\n')
        .skip_while(|line| line.starts_with('#'))
        .map(|line| line.trim_end_matches('\n'));
    let mut kept = KEPT.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let kept = kept.get_or_insert_with(HashMap::new);
    keep_noted(kept);
    let records = lines
        .zip(file.records)
        .map(|(line, parsed)| {
            let same = kept.get_mut(line).and_then(|candidates| {
                let at = candidates
                    .iter()
                    .position(|candidate| candidate.samples == file.header.samples)?;
                Some(if take {
                    let kept = candidates.remove(at);
                    (kept.record, kept.origin)
                } else {
                    (candidates[at].record.clone(), candidates[at].origin)
                })
            });
            same.unwrap_or((parsed, Origin::VcfCodec))
        })
        .collect();
    Ok((file.header, records))
}

/// `BCF2Writer` over records, each with where it came from, as `IndexingVariantContextWriter`
/// drives it: the indexer gets each record's offset before its bytes, and `close()` the length.
///
/// The builder opens the file before the writer sees a header, and nothing reaches the disk
/// before `close()` through its 128 KiB buffer, so a refusal leaves an empty file.
/// shortcut: an output past 128 KiB that fails keeps what the buffer passed on in Java, and is
/// empty here; upgrade if a row measures one.
pub(crate) fn write(
    output: &str,
    header: &VcfHeader,
    records: &[(VariantContext, Origin)],
    index: Option<&[(String, i32)]>,
    md5: bool,
) -> Result<(), Thrown> {
    let refused = |error: htsjdk_vcf::bcf2::Bcf2Error| {
        let _ = std::fs::write(output, b"");
        thrown(error.class, error.message)
    };
    let mut writer = Bcf2Writer::new(header, false).map_err(refused)?;
    let mut bytes = writer.header_bytes();
    let mut features = Vec::with_capacity(records.len());
    for (record, origin) in records {
        features.push((
            Feature {
                contig: record.contig.clone(),
                start: record.start as i32,
                end: record.stop as i32,
            },
            bytes.len() as i64,
        ));
        bytes.extend(writer.encode(record, *origin).map_err(refused)?);
    }
    std::fs::write(output, &bytes).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    if md5 {
        crate::runners::write_file(
            &format!("{output}.md5"),
            gatk_tools::gather_bam_files::md5_file(&bytes).as_bytes(),
        )?;
    }
    if let Some(dictionary) = index {
        // `Tribble.indexFile`: the output's name with `.idx` appended.
        let size = bytes.len() as i64;
        let idx = tribble_on_the_fly(
            &features,
            size,
            dictionary,
            output,
            size,
            modified_millis(output),
        );
        crate::runners::write_file(&format!("{output}.idx"), &idx)?;
    }
    Ok(())
}

/// What `IndexFactory.createDynamicIndex(path, new BCF2Codec())` indexes: each record's locus at
/// the offset its block starts, after the magic, the version, the header's INT32 length and the
/// header itself, each block being two INT32 lengths and what they count.
pub(crate) fn features(path: &str, bytes: &[u8]) -> Result<Vec<(Feature, i64)>, Thrown> {
    let file = read_bcf(bytes).map_err(|failure| {
        let error = failure.reader_exception(&format!("file://{}", java_absolute_path(path)));
        thrown(error.class, error.message)
    })?;
    // Every length below was read once already by the codec that accepted the file.
    let int = |at: usize| u32::from_le_bytes(bytes[at..at + 4].try_into().expect("four bytes"));
    let mut at = 9 + int(5) as usize;
    let mut found = Vec::with_capacity(file.records.len());
    for record in &file.records {
        found.push((
            Feature {
                contig: record.contig.clone(),
                start: record.start as i32,
                end: record.stop as i32,
            },
            at as i64,
        ));
        at += 8 + int(at) as usize + int(at + 4) as usize;
    }
    Ok(found)
}

