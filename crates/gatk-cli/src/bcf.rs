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
//!
//! A tool that changes a record's INFO hands Java's writer a map it built through copies, puts and
//! removes, each copy presized for what it copies, so the order depends on every step. Such a
//! tool describes its steps to [`rebuilt`], which replays them on the source record's keys and
//! keeps the result under the line the tool writes ([`keep_rebuilt`]). `SelectVariants` does.
//! Any other line a tool changed is parsed back and written as a record from `VCFCodec`, which is
//! the `Origin::VcfCodec` the writer replays.
//!
//! shortcut: a tool that changes INFO without describing its steps writes `VCFCodec`'s order over
//! its text, and a BCF read, edited and written as BCF loses its floats' digits past the third;
//! upgrade by routing that tool's records through [`rebuilt`].

use std::collections::HashMap;
use std::sync::Mutex;

use gatk_engine::java_hash::JavaHashMap;
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};
use htsjdk_tribble::index_write::Feature;
use htsjdk_vcf::bcf2::codec::read_bcf;
use htsjdk_vcf::bcf2::decode::{fully_decode, VCF_CODEC_INFO};
use htsjdk_vcf::bcf2::java_hash::{order_through, Fill};
use htsjdk_vcf::bcf2::writer::{Bcf2Writer, Origin};
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::header::VcfHeader;
use htsjdk_vcf::variant::{Value, VariantContext};

use crate::runners::{java_absolute_path, modified_millis, tribble_on_the_fly};

/// `FileExtensions.BCF`, which both the reader and the writer builder match case-sensitively.
pub(crate) fn is_bcf(path: &str) -> bool {
    path.ends_with(".bcf")
}

/// A record for `BCF2Writer`: the record, what `add` does with it, and the codec that read the
/// record it came from.
pub(crate) type Encoded = (VariantContext, Origin, Origin);

/// A record read in this run: the samples its genotypes were decoded under, the record, and the
/// codec that read it.
struct Kept {
    samples: Vec<String>,
    record: VariantContext,
    origin: Origin,
    /// For a record [`rebuilt`] for its INFO and filters alone, whose genotypes are the line's, the
    /// codec that read the record it was built from.
    rebuilt: Option<Origin>,
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
                    rebuilt: None,
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
            rebuilt: None,
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
    for (record, origin, _) in &records {
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
fn records_of(text: &str, take: bool) -> Result<(VcfHeader, Vec<Encoded>), Thrown> {
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
                let kept = if take {
                    candidates.remove(at)
                } else {
                    let candidate = &candidates[at];
                    Kept {
                        samples: Vec::new(),
                        record: candidate.record.clone(),
                        ..*candidate
                    }
                };
                let Some(read_by) = kept.rebuilt else {
                    return Some((kept.record, kept.origin, kept.origin));
                };
                // The tool's genotypes, as its line gives them, under the INFO its maps built.
                let mut decoded = fully_decode(&parsed, &file.header).ok()?;
                decoded.attributes = kept.record.attributes;
                decoded.filters = kept.record.filters;
                Some((decoded, Origin::Bcf2Codec, read_by))
            });
            same.unwrap_or((parsed, Origin::VcfCodec, Origin::VcfCodec))
        })
        .collect();
    Ok((file.header, records))
}

/// The record a tool read, as the codec that read it built it: the one kept when a BCF was read
/// here, and otherwise `record` itself, which a VCF's text gives with its INFO in file order.
pub(crate) fn source_of(header: &VcfHeader, record: &VariantContext) -> (VariantContext, Origin) {
    let kept = VcfEncoder::new(header)
        .encode(record)
        .ok()
        .and_then(|line| {
            let mut kept = KEPT.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
            kept.as_mut()?
                .get(&line)?
                .iter()
                .find(|candidate| candidate.samples == header.samples)
                .map(|candidate| (candidate.record.clone(), candidate.origin))
        });
    kept.unwrap_or_else(|| (record.clone(), Origin::VcfCodec))
}

/// One thing a tool does to a record's attribute map between reading the record and writing it.
pub(crate) enum Step<'a> {
    /// `new HashMap<>(map)`, or a builder's `attributes(map)`: a new map presized for the old.
    Copy,
    /// `new LinkedHashMap<>(map)`: the old map's order, which a `put` appends to.
    Linked,
    /// `put`, skipped when the written record does not hold the key.
    Put(&'a str),
    /// `remove`, skipped when the written record still holds the key.
    Remove(&'a str),
}

/// `new HashMap<>(source)`: `putMapEntries` presizes an empty table for the source's size, then
/// puts in the source's order.
fn copy_of(order: &[String]) -> JavaHashMap<String, ()> {
    let mut map = if order.is_empty() {
        JavaHashMap::new()
    } else {
        JavaHashMap::with_capacity((order.len() as f32 / 0.75 + 1.0) as usize)
    };
    for key in order {
        map.insert(key.clone(), ());
    }
    map
}

/// The record `BCF2Writer` is handed when a tool builds `written` from `source` through `steps`:
/// `written`'s fields, with its INFO in the order of the map the steps end in and each value the
/// Java object behind it, and with the codec's filter order. `value` gives a key's Java value where
/// the text lost it; otherwise a key a step put carries the written value and any other the
/// source's. `None` where the writer would refuse the record, which the text path then reports.
///
/// The record comes back fully decoded, as `Origin::Bcf2Codec` expects. A record from VCFCodec is
/// decoded by the writer in Java, and `fullyDecode` copies its map twice more.
/// shortcut: a decode that fails on two keys names the one `fully_decode`'s order meets first,
/// which may not be the one Java's map meets first; upgrade if a row measures it.
pub(crate) fn rebuilt(
    header: &VcfHeader,
    source: &(VariantContext, Origin),
    written: &VariantContext,
    steps: &[Step],
    value: impl Fn(&str) -> Option<Value>,
) -> Option<VariantContext> {
    let (source, origin) = source;
    let keys: Vec<String> = source
        .attributes
        .iter()
        .map(|(key, _)| key.clone())
        .collect();
    let mut order = match origin {
        Origin::VcfCodec => order_through(&keys, &VCF_CODEC_INFO)
            .into_iter()
            .map(|at| keys[at].clone())
            .collect(),
        Origin::Bcf2Codec => keys,
    };
    let holds = |key: &str| written.attributes.iter().any(|(name, _)| name == key);
    // The map the steps are on: a hash map, or `None` for a linked one, which is `order` itself.
    let mut map: Option<JavaHashMap<String, ()>> = None;
    let mut linked = false;
    let mut put: Vec<&str> = Vec::new();
    for step in steps {
        if let Some(hashed) = &map {
            if matches!(step, Step::Copy | Step::Linked) {
                order = hashed.keys().cloned().collect();
            }
        }
        match step {
            Step::Copy => {
                map = Some(copy_of(&order));
                linked = false;
            }
            Step::Linked => {
                map = None;
                linked = true;
            }
            Step::Put(key) if holds(key) => {
                match map.as_mut() {
                    Some(hashed) => {
                        hashed.insert(key.to_string(), ());
                    }
                    None if linked => {
                        if !order.iter().any(|name| name == key) {
                            order.push(key.to_string());
                        }
                    }
                    // The source's own map is only an order: a step before a copy is a misuse.
                    None => return None,
                }
                put.push(key);
            }
            Step::Remove(key) if !holds(key) => match map.as_mut() {
                Some(hashed) => {
                    hashed.remove(&key.to_string());
                }
                None if linked => order.retain(|name| name != key),
                None => return None,
            },
            _ => {}
        }
    }
    let order: Vec<String> = map.map_or(order, |map| map.keys().cloned().collect());
    let of = |attributes: &[(String, Value)], key: &str| {
        attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.clone())
    };
    let attributes = order
        .iter()
        .filter_map(|key| {
            let java = value(key).or_else(|| {
                if put.contains(&key.as_str()) {
                    of(&written.attributes, key)
                } else {
                    of(&source.attributes, key).or_else(|| of(&written.attributes, key))
                }
            })?;
            Some((key.clone(), java))
        })
        .collect();
    let filters = source.filters.as_ref().map(|filters| {
        let mut distinct: Vec<String> = Vec::new();
        for filter in filters {
            if !distinct.contains(filter) {
                distinct.push(filter.clone());
            }
        }
        match origin {
            // `new HashSet<>(list)` over the codec's filter list.
            Origin::VcfCodec => {
                let capacity = ((filters.len() as f32 / 0.75) as usize + 1).max(16);
                order_through(&distinct, &[Fill::PutSized(capacity)])
                    .into_iter()
                    .map(|at| distinct[at].clone())
                    .collect()
            }
            Origin::Bcf2Codec => distinct,
        }
    });
    let record = VariantContext {
        attributes,
        filters: filters.clone(),
        ..written.clone()
    };
    let mut decoded = fully_decode(&record, header).ok()?;
    let order = match origin {
        Origin::VcfCodec => {
            // `fullyDecode` builds from the record, a copy, and decodes into a third map.
            let held: Vec<String> = record
                .attributes
                .iter()
                .map(|(key, _)| key.clone())
                .collect();
            let keys: Vec<String> = copy_of(&held).keys().cloned().collect();
            order_through(&keys, &[Fill::PutSized(10), Fill::PutAll])
                .into_iter()
                .map(|at| keys[at].clone())
                .collect()
        }
        Origin::Bcf2Codec => order,
    };
    decoded.attributes = order
        .iter()
        .filter_map(|key| Some((key.clone(), of(&decoded.attributes, key)?)))
        .collect();
    decoded.filters = filters;
    Some(decoded)
}

/// A record [`rebuilt`] for the line `as_written` renders as under `header`, so a BCF written from
/// that text is written from it, ahead of any record read with the same line.
pub(crate) fn keep_rebuilt(
    header: &VcfHeader,
    record: VariantContext,
    read_by: Origin,
    as_written: &VariantContext,
) {
    let Ok(line) = VcfEncoder::new(header).encode(as_written) else {
        return;
    };
    let mut kept = KEPT.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    kept.get_or_insert_with(HashMap::new)
        .entry(line)
        .or_default()
        .insert(
            0,
            Kept {
                samples: header.samples.clone(),
                record,
                origin: Origin::Bcf2Codec,
                rebuilt: Some(read_by),
            },
        );
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
    records: &[Encoded],
    index: Option<&[(String, i32)]>,
    md5: bool,
) -> Result<(), Thrown> {
    let refused = |error: htsjdk_vcf::bcf2::Bcf2Error| {
        let _ = std::fs::write(output, b"");
        thrown(error.class, error.message)
    };
    // A refusal names the record's source: VCFCodec's default name, `Unknown`, which the port's
    // writer prints, or a BCF2Codec record's, which is `null`.
    let named = |read_by: Origin| {
        move |mut error: htsjdk_vcf::bcf2::Bcf2Error| {
            if read_by == Origin::Bcf2Codec {
                error.message = error
                    .message
                    .replace(" from Unknown but ", " from null but ");
            }
            refused(error)
        }
    };
    let mut writer = Bcf2Writer::new(header, false).map_err(refused)?;
    let mut bytes = writer.header_bytes();
    let mut features = Vec::with_capacity(records.len());
    for (record, origin, read_by) in records {
        features.push((
            Feature {
                contig: record.contig.clone(),
                start: record.start as i32,
                end: record.stop as i32,
            },
            bytes.len() as i64,
        ));
        bytes.extend(writer.encode(record, *origin).map_err(named(*read_by))?);
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

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::genotypes_context::GenotypesContext;

    /// With no steps, a VCF record comes back in the order the writer's own `fullyDecode` gives,
    /// which is what the steps then build on.
    #[test]
    fn a_record_with_no_steps_is_ordered_as_the_writer_decodes_it() {
        let text = "##fileformat=VCFv4.2\n\
            ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">\n\
            ##INFO=<ID=AN,Number=1,Type=Integer,Description=\"x\">\n\
            ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n\
            ##INFO=<ID=MQ,Number=1,Type=Float,Description=\"x\">\n\
            ##INFO=<ID=QD,Number=1,Type=Float,Description=\"x\">\n\
            ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"x\">\n\
            ##contig=<ID=chr1,length=100>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
            chr1\t5\t.\tA\tC\t9\t.\tReadPosRankSum=1.5;QD=2.0;MQ=60.0;DP=7;AN=2;AC=1\n";
        let file = htsjdk_vcf::reader::read_vcf(text).unwrap();
        let record = VariantContext {
            genotypes: GenotypesContext::new(Vec::new()),
            ..file.records[0].clone()
        };
        let source = (record.clone(), Origin::VcfCodec);
        let ours = rebuilt(&file.header, &source, &record, &[], |_| None).unwrap();
        let theirs = fully_decode(&record, &file.header).unwrap();
        assert_eq!(ours.attributes, theirs.attributes);
    }
}
