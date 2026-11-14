//! A second reader of every file the goldens carry: noodles against htsjdk-rs (IPNP-BIPN/gatk-rs#1393).
//!
//! Every byte the port writes is checked against the reference's, and those bytes travel in the
//! goldens. What that does not check is the *reading*: the port reads through htsjdk-rs, which this
//! programme also wrote, so a file can equal its golden and still be misread by both of our readers
//! in the same way, or carry an index structure nothing here ever reads back. noodles was written by
//! someone else from the specifications, so each disagreement is either our bug, theirs, or a place
//! where the specification leaves room, and the oracle decides which (docs/a-second-reader.md).
//!
//! What it reads: every base64 BGZF stream (`H4sIBA`), every raw `.bai` (`QkFJ`) and every text VCF
//! (`##fileformat=VCF`) in `crates/*/tests/data/*.txt.gz`. Those include inputs as well as outputs;
//! both are the reference's bytes, and every output among them is a file the port is required to
//! write byte for byte. noodles is a dev-dependency only: it never reaches the byte path.

use std::collections::BTreeMap;
use std::io::{Cursor, Read};
use std::path::{Path, PathBuf};

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::TagValue;
use noodles::sam::alignment::record_buf::data::field::value::Array as NArray;
use noodles::sam::alignment::record_buf::data::field::Value as NTag;
use noodles::{bam, bgzf, csi, sam, tabix, vcf};

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Bgzf,
    Bai,
    Vcf,
}

/// One file found in a golden, named by where it was found.
struct Found {
    at: String,
    bytes: Vec<u8>,
}

#[derive(Default)]
struct Report {
    /// What was read, by kind.
    read: BTreeMap<&'static str, usize>,
    /// Every disagreement, as `kind: where: what`.
    disagreements: Vec<String>,
}

impl Report {
    fn count(&mut self, kind: &'static str) {
        *self.read.entry(kind).or_default() += 1;
    }
    fn differ(&mut self, kind: &str, at: &str, what: String) {
        self.disagreements.push(format!("{kind}: {at}: {what}"));
    }
}

fn goldens() -> Vec<PathBuf> {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut out = Vec::new();
    for entry in std::fs::read_dir(&crates).unwrap() {
        let data = entry.unwrap().path().join("tests").join("data");
        let Ok(files) = std::fs::read_dir(&data) else {
            continue;
        };
        for file in files {
            let path = file.unwrap().path();
            if path.to_string_lossy().ends_with(".txt.gz") {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

fn is_base64(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'+' || byte == b'/' || byte == b'='
}

/// The dumps' `escape`: `\\`, `\t`, `\n`, and `\r` where a dump wrote one.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

/// Base64 runs that start with `start`, decoded, and text VCFs, unescaped.
fn scan(golden: &Path) -> Vec<(Kind, Found)> {
    let name = golden.file_name().unwrap().to_string_lossy().into_owned();
    let text = gatk_corpus::read_golden(golden);
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let at = |what: &str| format!("{name}:{}:{what}", n + 1);
        let bytes = line.as_bytes();
        let mut i = 0;
        while i + 6 <= bytes.len() {
            let start = &bytes[i..i + 6];
            let kind = if start == b"H4sIBA" {
                Some(Kind::Bgzf)
            } else if start == b"QkFJAQ" {
                Some(Kind::Bai)
            } else {
                None
            };
            // A match must start a token, or it is the middle of some other base64 run.
            let starts_token = i == 0 || !is_base64(bytes[i - 1]) || bytes[i - 1] == b'=';
            match kind {
                Some(kind) if starts_token => {
                    let end = (i..bytes.len())
                        .find(|&j| !is_base64(bytes[j]))
                        .unwrap_or(bytes.len());
                    let bytes = gatk_corpus::decode_base64(&line[i..end]);
                    out.push((
                        kind,
                        Found {
                            at: at(&format!("{i}")),
                            bytes,
                        },
                    ));
                    i = end;
                }
                _ => i += 1,
            }
        }
        let mut from = 0;
        while let Some(p) = line[from..].find("##fileformat=VCF") {
            let start = from + p;
            let end = line[start..].find('\t').map_or(line.len(), |e| start + e);
            let bytes = unescape(&line[start..end]).into_bytes();
            out.push((
                Kind::Vcf,
                Found {
                    at: at(&format!("{start}")),
                    bytes,
                },
            ));
            from = end;
        }
    }
    out
}

fn noodles_inflate(bytes: &[u8]) -> std::io::Result<Vec<u8>> {
    let mut out = Vec::new();
    bgzf::io::Reader::new(bytes).read_to_end(&mut out)?;
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// BAM

/// A tag value in the one shape both readers can be put in: htsjdk keeps no integer width and
/// reads `H` as a signed byte array, so the comparison cannot ask for more than that.
#[derive(Debug, PartialEq)]
enum Tag {
    Char(u8),
    Int(i64),
    Float(u32),
    Str(Vec<u8>),
    Ints(Vec<i64>, bool),
    Floats(Vec<u32>),
}

fn ours_tag(value: &TagValue) -> Tag {
    match value {
        TagValue::Char(c) => Tag::Char(*c),
        TagValue::Int(i) => Tag::Int(*i),
        TagValue::Float(f) => Tag::Float(f.to_bits()),
        TagValue::Str(s) => Tag::Str(s.as_bytes().to_vec()),
        TagValue::ByteArray { values, unsigned } => Tag::Ints(
            values
                .iter()
                .map(|&v| if *unsigned { v as u8 as i64 } else { v as i64 })
                .collect(),
            *unsigned,
        ),
        TagValue::ShortArray { values, unsigned } => Tag::Ints(
            values
                .iter()
                .map(|&v| if *unsigned { v as u16 as i64 } else { v as i64 })
                .collect(),
            *unsigned,
        ),
        TagValue::IntArray { values, unsigned } => Tag::Ints(
            values
                .iter()
                .map(|&v| if *unsigned { v as u32 as i64 } else { v as i64 })
                .collect(),
            *unsigned,
        ),
        TagValue::FloatArray(values) => Tag::Floats(values.iter().map(|f| f.to_bits()).collect()),
    }
}

fn noodles_tag(value: &NTag) -> Tag {
    match value {
        NTag::Character(c) => Tag::Char(*c),
        NTag::Int8(v) => Tag::Int(*v as i64),
        NTag::UInt8(v) => Tag::Int(*v as i64),
        NTag::Int16(v) => Tag::Int(*v as i64),
        NTag::UInt16(v) => Tag::Int(*v as i64),
        NTag::Int32(v) => Tag::Int(*v as i64),
        NTag::UInt32(v) => Tag::Int(*v as i64),
        NTag::Float(f) => Tag::Float(f.to_bits()),
        NTag::String(s) => Tag::Str(s.to_vec()),
        // htsjdk decodes `H` into a byte[], which it then reports as a signed `B:c` array.
        NTag::Hex(s) => Tag::Ints(
            s.chunks(2)
                .map(|pair| {
                    let pair = std::str::from_utf8(pair).unwrap();
                    u8::from_str_radix(pair, 16).unwrap() as i8 as i64
                })
                .collect(),
            false,
        ),
        NTag::Array(array) => match array {
            NArray::Int8(v) => Tag::Ints(v.iter().map(|&x| x as i64).collect(), false),
            NArray::UInt8(v) => Tag::Ints(v.iter().map(|&x| x as i64).collect(), true),
            NArray::Int16(v) => Tag::Ints(v.iter().map(|&x| x as i64).collect(), false),
            NArray::UInt16(v) => Tag::Ints(v.iter().map(|&x| x as i64).collect(), true),
            NArray::Int32(v) => Tag::Ints(v.iter().map(|&x| x as i64).collect(), false),
            NArray::UInt32(v) => Tag::Ints(v.iter().map(|&x| x as i64).collect(), true),
            NArray::Float(v) => Tag::Floats(v.iter().map(|f| f.to_bits()).collect()),
        },
    }
}

/// The record as both readers can state it, field by field.
fn ours_fields(r: &BamRecord) -> Vec<(&'static str, String)> {
    let cigar: String = r
        .cigar
        .elements
        .iter()
        .map(|e| format!("{}{}", e.length, e.op.to_char() as char))
        .collect();
    let tags: Vec<String> = r
        .tags
        .iter()
        .map(|(t, v)| format!("{t}:{:?}", ours_tag(v)))
        .collect();
    vec![
        ("name", r.read_name.clone()),
        ("flags", r.flags.to_string()),
        ("reference", r.reference_index.to_string()),
        ("start", r.alignment_start.to_string()),
        ("mapq", r.mapping_quality.to_string()),
        ("cigar", cigar),
        ("mate reference", r.mate_reference_index.to_string()),
        ("mate start", r.mate_alignment_start.to_string()),
        ("tlen", r.inferred_insert_size.to_string()),
        ("bases", String::from_utf8_lossy(&r.read_bases).into_owned()),
        ("qualities", format!("{:?}", r.base_qualities)),
        ("tags", tags.join(" ")),
    ]
}

fn noodles_fields(r: &sam::alignment::RecordBuf) -> Vec<(&'static str, String)> {
    let index = |v: Option<usize>| v.map_or(-1, |i| i as i64).to_string();
    let position = |p: Option<noodles::core::Position>| p.map_or(0, usize::from).to_string();
    let cigar: String = r
        .cigar()
        .as_ref()
        .iter()
        .map(|op| {
            use sam::alignment::record::cigar::op::Kind;
            let c = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            format!("{}{c}", op.len())
        })
        .collect();
    let tags: Vec<String> = r
        .data()
        .iter()
        .map(|(t, v)| {
            let t: [u8; 2] = *t.as_ref();
            format!("{}{}:{:?}", t[0] as char, t[1] as char, noodles_tag(v))
        })
        .collect();
    vec![
        (
            "name",
            r.name()
                .map_or("*".to_string(), |n| String::from_utf8_lossy(n).into_owned()),
        ),
        ("flags", u16::from(r.flags()).to_string()),
        ("reference", index(r.reference_sequence_id())),
        ("start", position(r.alignment_start())),
        (
            "mapq",
            r.mapping_quality().map_or(255, u8::from).to_string(),
        ),
        ("cigar", cigar),
        ("mate reference", index(r.mate_reference_sequence_id())),
        ("mate start", position(r.mate_alignment_start())),
        ("tlen", r.template_length().to_string()),
        (
            "bases",
            String::from_utf8_lossy(r.sequence().as_ref()).into_owned(),
        ),
        ("qualities", format!("{:?}", r.quality_scores().as_ref())),
        ("tags", tags.join(" ")),
    ]
}

/// One BAM, read by both; returns noodles' header and records for the index check.
fn compare_bam(
    found: &Found,
    inflated: &[u8],
    report: &mut Report,
) -> Option<(sam::Header, Vec<BamRecord>)> {
    let ours = match BamReader::new(inflated) {
        Ok(reader) => reader,
        Err(e) => {
            report.differ(
                "bam",
                &found.at,
                format!("htsjdk-rs refuses the header: {e:?}"),
            );
            return None;
        }
    };
    let dictionary = ours.header.binary_sequences.clone();
    let ours: Result<Vec<BamRecord>, _> = ours.collect();

    let mut reader = bam::io::Reader::new(found.bytes.as_slice());
    let header = match reader.read_header() {
        Ok(header) => header,
        Err(e) => {
            report.differ("bam", &found.at, format!("noodles refuses the header: {e}"));
            return None;
        }
    };
    let theirs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| (name.to_string(), usize::from(seq.length())))
        .collect();
    let dictionary_ours: Vec<(String, usize)> = dictionary
        .iter()
        .map(|(n, l)| (n.clone(), *l as usize))
        .collect();
    if theirs != dictionary_ours {
        report.differ(
            "bam",
            &found.at,
            format!("dictionary: ours {dictionary_ours:?}, noodles {theirs:?}"),
        );
    }

    let theirs: Result<Vec<_>, _> = reader.record_bufs(&header).collect();
    let (ours, theirs) = match (ours, theirs) {
        (Ok(o), Ok(t)) => (o, t),
        (o, t) => {
            report.differ(
                "bam",
                &found.at,
                format!(
                    "records: htsjdk-rs {:?}, noodles {:?}",
                    o.map(|v| v.len()),
                    t.map(|v| v.len()).map_err(|e| e.to_string())
                ),
            );
            return None;
        }
    };
    if ours.len() != theirs.len() {
        report.differ(
            "bam",
            &found.at,
            format!("{} records against noodles' {}", ours.len(), theirs.len()),
        );
    }
    for (n, (o, t)) in ours.iter().zip(&theirs).enumerate() {
        for ((field, a), (_, b)) in ours_fields(o).into_iter().zip(noodles_fields(t)) {
            if a != b {
                report.differ(
                    "bam",
                    &found.at,
                    format!("record {n} {field}: htsjdk-rs {a:?}, noodles {b:?}"),
                );
            }
        }
    }
    report.count("bam records");
    *report.read.entry("bam records").or_default() += ours.len() - 1;
    Some((header, ours))
}

// ---------------------------------------------------------------------------------------------
// BAI

/// The BAMs an index may describe: those whose per-reference counts are the index's own.
///
/// Dumps put a `.bai` before its BAM as often as after it, and some carry several of each, so
/// neither neighbour can be trusted. The pseudo-bins' aligned and unaligned counts and the
/// unplaced count narrow it down; several files of one dump can share them, so the caller keeps
/// the candidate the index actually answers for.
fn described<'a>(index: &htsjdk_bam::index::BamIndex, bams: &'a [Bam]) -> Vec<&'a Bam> {
    let wanted: Vec<(i64, i64)> = index
        .references
        .iter()
        .map(|r| {
            r.metadata
                .as_ref()
                .map_or((0, 0), |m| (m.aligned, m.unaligned))
        })
        .collect();
    bams.iter()
        .rev()
        .filter(|(_, _, records)| {
            let mut counts = vec![(0i64, 0i64); wanted.len()];
            let mut unplaced = 0i64;
            for r in records {
                match counts.get_mut(r.reference_index.max(0) as usize) {
                    Some(c) if r.reference_index >= 0 => {
                        if r.read_unmapped() {
                            c.1 += 1
                        } else {
                            c.0 += 1
                        }
                    }
                    _ if r.reference_index < 0 => unplaced += 1,
                    _ => return false,
                }
            }
            counts == wanted && unplaced == index.no_coordinate_records
        })
        .collect()
}

/// A BAM as both readers saw it: where it was, noodles' header, htsjdk-rs's records.
type Bam = (Found, sam::Header, Vec<BamRecord>);

fn compare_bai(found: &Found, bams: &[Bam], report: &mut Report) {
    use csi::binning_index::{BinningIndex, ReferenceSequence as _};

    let ours = htsjdk_bam::index::read_bai(&found.bytes);
    let theirs = bam::bai::io::Reader::new(found.bytes.as_slice()).read_index();
    let (ours, theirs) = match (ours, theirs) {
        (Ok(o), Ok(t)) => (o, t),
        (o, t) => {
            report.differ(
                "bai",
                &found.at,
                format!(
                    "htsjdk-rs {:?}, noodles {:?}",
                    o.map(|_| ()),
                    t.map(|_| ()).map_err(|e| e.to_string())
                ),
            );
            return;
        }
    };
    report.count("bai");

    // The structure, as each reads it.
    let theirs_refs = theirs.reference_sequences();
    if ours.references.len() != theirs_refs.len() {
        report.differ(
            "bai",
            &found.at,
            format!(
                "{} references against noodles' {}",
                ours.references.len(),
                theirs_refs.len()
            ),
        );
        return;
    }
    for (i, (o, t)) in ours.references.iter().zip(theirs_refs).enumerate() {
        let mut ours_bins: Vec<(usize, Vec<(u64, u64)>)> = o
            .bins
            .iter()
            .map(|b| {
                (
                    b.bin_number as usize,
                    b.chunks.iter().map(|c| (c.start, c.end)).collect(),
                )
            })
            .collect();
        ours_bins.sort();
        let mut theirs_bins: Vec<(usize, Vec<(u64, u64)>)> = t
            .bins()
            .iter()
            .map(|(id, bin)| {
                (
                    *id,
                    bin.chunks()
                        .iter()
                        .map(|c| (u64::from(c.start()), u64::from(c.end())))
                        .collect(),
                )
            })
            .collect();
        theirs_bins.sort();
        if ours_bins != theirs_bins {
            report.differ("bai", &found.at, format!("reference {i}: bins differ"));
        }
        let theirs_linear: Vec<u64> = t.index().iter().map(|&v| u64::from(v)).collect();
        if o.linear_index != theirs_linear {
            report.differ(
                "bai",
                &found.at,
                format!(
                    "reference {i}: linear index {:?} against noodles' {theirs_linear:?}",
                    o.linear_index
                ),
            );
        }
        let ours_meta = o.metadata.as_ref().map(|m| {
            (
                m.first_offset,
                m.last_offset,
                m.aligned as u64,
                m.unaligned as u64,
            )
        });
        let theirs_meta = t.metadata().map(|m| {
            (
                u64::from(m.start_position()),
                u64::from(m.end_position()),
                m.mapped_record_count(),
                m.unmapped_record_count(),
            )
        });
        if ours_meta != theirs_meta {
            report.differ(
                "bai",
                &found.at,
                format!("reference {i}: metadata {ours_meta:?} against noodles' {theirs_meta:?}"),
            );
        }
    }
    let theirs_unplaced = theirs.unplaced_unmapped_record_count().unwrap_or(0);
    if ours.no_coordinate_records as u64 != theirs_unplaced {
        report.differ(
            "bai",
            &found.at,
            format!(
                "unplaced {} against noodles' {theirs_unplaced}",
                ours.no_coordinate_records
            ),
        );
    }

    // The index used: every mapped record must come back from a query over its own span.
    let candidates = described(&ours, bams);
    if candidates.is_empty() {
        report.count("bai without its bam");
        return;
    }
    // Pick the candidate cheaply, on the first and last mapped record of each reference, then
    // query every record of the one picked.
    let bam = candidates
        .iter()
        .find(|bam| misses(bam, &theirs, usize::MAX).0.is_empty())
        .unwrap_or(&candidates[0]);
    let mapped = bam
        .2
        .iter()
        .filter(|r| !r.read_unmapped() && r.reference_index >= 0)
        .count();
    // shortcut: about 64 records an index, not all 50,599 of the corpus (81 s, every one returned
    // when measured); query every record again if a reader or the index writer changes.
    let (failures, queried) = misses(bam, &theirs, (mapped / QUERIES_PER_INDEX).max(1));
    *report.read.entry("bai queries").or_default() += queried;
    for failure in failures {
        report.differ("bai query", &found.at, failure);
    }
}

/// How many records of each BAM are queried through its index, besides each reference's ends.
const QUERIES_PER_INDEX: usize = 64;

/// The mapped records of `bam` that noodles does not return from a query over their own span,
/// and how many were queried.
///
/// Queries the first and last mapped record of each reference and every `stride`-th record.
fn misses(bam: &Bam, index: &bam::bai::Index, stride: usize) -> (Vec<String>, usize) {
    let (bam_found, header, records) = bam;
    let mapped: Vec<&BamRecord> = records
        .iter()
        .filter(|r| !r.read_unmapped() && r.reference_index >= 0)
        .collect();
    let picked: Vec<&BamRecord> = mapped
        .iter()
        .enumerate()
        .filter(|(i, r)| {
            let next = mapped.get(i + 1).map(|n| n.reference_index);
            let previous = i.checked_sub(1).map(|p| mapped[p].reference_index);
            i % stride == 0
                || previous != Some(r.reference_index)
                || next != Some(r.reference_index)
        })
        .map(|(_, r)| *r)
        .collect();
    let queried = picked.len();
    let mut reader = bam::io::Reader::new(Cursor::new(bam_found.bytes.as_slice()));
    reader.read_header().unwrap();
    let mut out = Vec::new();
    for r in picked {
        let name = header
            .reference_sequences()
            .get_index(r.reference_index as usize)
            .unwrap()
            .0
            .to_string();
        let start = noodles::core::Position::try_from(r.alignment_start as usize).unwrap();
        let end = r.alignment_end().max(r.alignment_start) as usize;
        let end = noodles::core::Position::try_from(end).unwrap();
        let region = noodles::core::Region::new(name, start..=end);
        let found = reader
            .query(header, index, &region)
            .and_then(|q| q.records().collect::<Result<Vec<_>, _>>())
            .map(|hits| {
                hits.iter().any(|h| {
                    h.name().map(|n| n.to_vec()) == Some(r.read_name.as_bytes().to_vec())
                        && u16::from(h.flags()) == r.flags
                })
            });
        if !matches!(found, Ok(true)) {
            out.push(format!(
                "{} at {region} is not returned ({found:?})",
                r.read_name
            ));
        }
    }
    (out, queried)
}

// ---------------------------------------------------------------------------------------------
// VCF

fn compare_vcf(found: &Found, text: &[u8], report: &mut Report) {
    let Ok(text) = std::str::from_utf8(text) else {
        report.differ("vcf", &found.at, "not UTF-8".to_string());
        return;
    };
    let ours = htsjdk_vcf::reader::read_vcf(text);
    let mut reader = vcf::io::Reader::new(text.as_bytes());
    let header = reader.read_header();
    let (ours, header) = match (ours, header) {
        (Ok(o), Ok(h)) => (o, h),
        (Err(o), Err(_)) => {
            // Both refuse: nothing to compare, and nothing the port writes.
            let _ = o;
            report.count("vcf refused by both");
            return;
        }
        (o, h) => {
            report.differ(
                "vcf",
                &found.at,
                format!(
                    "header: htsjdk-rs {:?}, noodles {:?}",
                    o.map(|_| ()).map_err(|e| e.error.message()),
                    h.map(|_| ()).map_err(|e| e.to_string())
                ),
            );
            return;
        }
    };
    report.count("vcf");
    let theirs: Result<Vec<_>, _> = reader.record_bufs(&header).collect();
    let theirs = match theirs {
        Ok(t) => t,
        Err(e) => {
            report.differ("vcf", &found.at, format!("noodles refuses a record: {e:?}"));
            return;
        }
    };
    if ours.records.len() != theirs.len() {
        report.differ(
            "vcf",
            &found.at,
            format!(
                "{} records against noodles' {}",
                ours.records.len(),
                theirs.len()
            ),
        );
        return;
    }
    *report.read.entry("vcf records").or_default() += theirs.len();
    for (n, (o, t)) in ours.records.iter().zip(&theirs).enumerate() {
        for (field, a, b) in vcf_fields(o, t) {
            if a != b {
                report.differ(
                    "vcf",
                    &found.at,
                    format!("record {n} {field}: htsjdk-rs {a:?}, noodles {b:?}"),
                );
            }
        }
    }
}

fn ours_strings(value: &htsjdk_vcf::variant::Value) -> Vec<String> {
    use htsjdk_vcf::variant::Value;
    match value {
        Value::Missing => vec![".".to_string()],
        Value::Int(i) => vec![i.to_string()],
        Value::Double(d) => vec![d.to_string()],
        Value::Bool(_) => Vec::new(),
        Value::Str(s) => vec![s.clone()],
        Value::List(items) => items.iter().flat_map(ours_strings).collect(),
    }
}

/// noodles' typed values, put back against htsjdk's strings: a number compares by value.
///
/// Both sides are flattened on commas first. htsjdk splits every INFO value on `,` whatever its
/// `Number`, and keeps a FORMAT value it has no slot for as one unsplit string; noodles splits by
/// `Number`. That is two representations of one field, not two readings of it.
fn same_values(ours: &[String], theirs: &[Option<String>], floats: bool) -> bool {
    let ours: Vec<&str> = ours.iter().flat_map(|s| s.split(',')).collect();
    let theirs: Vec<Option<&str>> = theirs
        .iter()
        .flat_map(|s| match s {
            None => vec![None],
            Some(s) => s.split(',').map(Some).collect(),
        })
        .collect();
    ours.len() == theirs.len()
        && ours.iter().zip(&theirs).all(|(a, b)| match b {
            None => *a == ".",
            Some(b) if floats => {
                a.parse::<f32>().map(f32::to_bits) == b.parse::<f32>().map(f32::to_bits)
                    || (a.parse::<f32>().is_ok_and(f32::is_nan) && *b == "NaN")
            }
            Some(b) => a == b,
        })
}

fn info_strings(
    value: Option<&vcf::variant::record_buf::info::field::Value>,
) -> (Vec<Option<String>>, bool) {
    use vcf::variant::record_buf::info::field::value::Array;
    use vcf::variant::record_buf::info::field::Value;
    let s = |v: Option<String>| v;
    match value {
        None => (vec![None], false),
        Some(Value::Flag) => (Vec::new(), false),
        Some(Value::Integer(i)) => (vec![Some(i.to_string())], false),
        Some(Value::Float(f)) => (vec![Some(f.to_string())], true),
        Some(Value::Character(c)) => (vec![Some(c.to_string())], false),
        Some(Value::String(t)) => (vec![Some(t.clone())], false),
        Some(Value::Array(Array::Integer(v))) => (
            v.iter().map(|x| s(x.map(|x| x.to_string()))).collect(),
            false,
        ),
        Some(Value::Array(Array::Float(v))) => (
            v.iter().map(|x| s(x.map(|x| x.to_string()))).collect(),
            true,
        ),
        Some(Value::Array(Array::Character(v))) => (
            v.iter().map(|x| s(x.map(|x| x.to_string()))).collect(),
            false,
        ),
        Some(Value::Array(Array::String(v))) => (v.to_vec(), false),
    }
}

fn sample_strings(
    value: Option<&vcf::variant::record_buf::samples::sample::Value>,
) -> (Vec<Option<String>>, bool) {
    use vcf::variant::record_buf::samples::sample::value::Array;
    use vcf::variant::record_buf::samples::sample::Value;
    match value {
        None => (vec![None], false),
        Some(Value::Integer(i)) => (vec![Some(i.to_string())], false),
        Some(Value::Float(f)) => (vec![Some(f.to_string())], true),
        Some(Value::Character(c)) => (vec![Some(c.to_string())], false),
        Some(Value::String(t)) => (vec![Some(t.clone())], false),
        Some(Value::Genotype(_)) => (Vec::new(), false),
        Some(Value::Array(Array::Integer(v))) => {
            (v.iter().map(|x| x.map(|x| x.to_string())).collect(), false)
        }
        Some(Value::Array(Array::Float(v))) => {
            (v.iter().map(|x| x.map(|x| x.to_string())).collect(), true)
        }
        Some(Value::Array(Array::Character(v))) => {
            (v.iter().map(|x| x.map(|x| x.to_string())).collect(), false)
        }
        Some(Value::Array(Array::String(v))) => (v.to_vec(), false),
    }
}

fn vcf_fields(
    o: &htsjdk_vcf::variant::VariantContext,
    t: &vcf::variant::RecordBuf,
) -> Vec<(String, String, String)> {
    let mut out = Vec::new();
    let mut push = |field: &str, a: String, b: String| out.push((field.to_string(), a, b));
    push(
        "chrom",
        o.contig.clone(),
        t.reference_sequence_name().to_string(),
    );
    push(
        "pos",
        o.start.to_string(),
        t.variant_start().map_or(0, usize::from).to_string(),
    );
    let ids: Vec<&str> = t.ids().as_ref().iter().map(String::as_str).collect();
    push(
        "id",
        o.id.clone(),
        if ids.is_empty() {
            ".".to_string()
        } else {
            ids.join(";")
        },
    );
    let alleles: Vec<String> = o.alleles.iter().map(|a| a.display_string()).collect();
    let mut theirs_alleles = vec![t.reference_bases().to_string()];
    theirs_alleles.extend(t.alternate_bases().as_ref().iter().cloned());
    push("alleles", alleles.join(","), theirs_alleles.join(","));
    let ours_qual = (o.log10_p_error != htsjdk_vcf::variant::NO_LOG10_PERROR)
        .then(|| ((-10.0 * o.log10_p_error) as f32).to_bits());
    push(
        "qual",
        format!("{ours_qual:?}"),
        format!("{:?}", t.quality_score().map(f32::to_bits)),
    );
    let theirs_filters: Vec<String> = t.filters().as_ref().iter().cloned().collect();
    let ours_filters = match &o.filters {
        None => Vec::new(),
        Some(f) if f.is_empty() => vec!["PASS".to_string()],
        Some(f) => f.clone(),
    };
    push("filter", ours_filters.join(";"), theirs_filters.join(";"));

    let ours_info: Vec<&str> = o.attributes.iter().map(|(k, _)| k.as_str()).collect();
    let theirs_info: Vec<&str> = t.info().keys().map(String::as_str).collect();
    push("info keys", ours_info.join(";"), theirs_info.join(";"));
    for (key, value) in &o.attributes {
        let Some(theirs) = t.info().get(key.as_str()) else {
            continue;
        };
        let ours = ours_strings(value);
        let (theirs, floats) = info_strings(theirs);
        if !same_values(&ours, &theirs, floats) {
            push(
                &format!("INFO {key}"),
                format!("{ours:?}"),
                format!("{theirs:?}"),
            );
        }
    }

    let genotypes = o.genotypes.peek();
    let keys: Vec<&str> = t.format().as_ref().iter().map(String::as_str).collect();
    push(
        "samples",
        genotypes.len().to_string(),
        t.samples().values().count().to_string(),
    );
    for (g, sample) in genotypes.iter().zip(t.samples().values()) {
        for (key, value) in keys.iter().zip(sample.values()) {
            let ours: Option<Vec<String>> = match *key {
                "GT" => {
                    let indices: Vec<String> = g
                        .alleles
                        .iter()
                        .map(|a| {
                            if a.is_no_call() {
                                ".".to_string()
                            } else {
                                o.alleles
                                    .iter()
                                    .position(|b| b == a || b.equals_ignore_ref(a))
                                    .map_or("?".to_string(), |i| i.to_string())
                            }
                        })
                        .collect();
                    let theirs = match value {
                        Some(vcf::variant::record_buf::samples::sample::Value::Genotype(gt)) => gt
                            .as_ref()
                            .iter()
                            .map(|a| a.position().map_or(".".to_string(), |p| p.to_string()))
                            .collect::<Vec<_>>(),
                        _ => vec![".".to_string()],
                    };
                    if indices != theirs {
                        push(
                            &format!("{} GT", g.sample_name),
                            indices.join("/"),
                            theirs.join("/"),
                        );
                    }
                    None
                }
                "GQ" => Some(g.gq.map_or(vec![".".into()], |v| vec![v.to_string()])),
                "DP" => Some(g.dp.map_or(vec![".".into()], |v| vec![v.to_string()])),
                "AD" => Some(g.ad.as_ref().map_or(vec![".".into()], |v| {
                    v.iter().map(|x| x.to_string()).collect()
                })),
                "PL" => Some(g.pl.as_ref().map_or(vec![".".into()], |v| {
                    v.iter().map(|x| x.to_string()).collect()
                })),
                "FT" => Some(vec![g.filters.clone().unwrap_or_else(|| ".".into())]),
                other => g
                    .extended
                    .iter()
                    .find(|(k, _)| k == other)
                    .map(|(_, v)| ours_strings(v))
                    .or(Some(vec![".".into()])),
            };
            if let Some(ours) = ours {
                let (theirs, floats) = sample_strings(value.as_ref());
                if !same_values(&ours, &theirs, floats) {
                    push(
                        &format!("{} {key}", g.sample_name),
                        format!("{ours:?}"),
                        format!("{theirs:?}"),
                    );
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------------------------
// Tabix and the BGZF-compressed VCF it indexes

fn check_tabix(found: &Found, vcf_gz: Option<&Found>, report: &mut Report) {
    let index = match tabix::io::Reader::new(found.bytes.as_slice()).read_index() {
        Ok(index) => index,
        Err(e) => {
            report.differ("tbi", &found.at, format!("noodles refuses it: {e}"));
            return;
        }
    };
    report.count("tbi");
    let Some(vcf_gz) = vcf_gz else {
        return;
    };
    let mut reader =
        vcf::io::Reader::new(bgzf::io::Reader::new(Cursor::new(vcf_gz.bytes.as_slice())));
    let Ok(header) = reader.read_header() else {
        return;
    };
    let records: Vec<_> = reader.record_bufs(&header).filter_map(Result::ok).collect();
    for r in &records {
        let Some(start) = r.variant_start() else {
            continue;
        };
        let region = noodles::core::Region::new(r.reference_sequence_name(), start..=start);
        let mut reader =
            vcf::io::Reader::new(bgzf::io::Reader::new(Cursor::new(vcf_gz.bytes.as_slice())));
        reader.read_header().unwrap();
        let hits: Vec<_> = match reader.query(&header, &index, &region) {
            Ok(q) => q.records().filter_map(Result::ok).collect(),
            Err(e) => {
                report.differ("tbi query", &found.at, format!("{region}: {e}"));
                continue;
            }
        };
        report.count("tbi queries");
        let hit = hits.iter().any(|h| {
            h.reference_sequence_name() == r.reference_sequence_name()
                && h.variant_start().and_then(Result::ok) == Some(start)
        });
        if !hit {
            report.differ(
                "tbi query",
                &found.at,
                format!("the record at {region} is not returned"),
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------

/// Every disagreement the gate has found, with the oracle's verdict. A new one fails the test;
/// one that stops occurring fails it too, so the list cannot rot (docs/a-second-reader.md).
const KNOWN: &[(&str, &str, &str)] = &[
    (
        "bgzf: check_terminator_block.txt.gz:23:",
        "block data checksum mismatch",
        "noodles' bug in the sense of the reference: the dump corrupts one block's payload on \
         purpose, noodles checks every CRC32, htsjdk checks none unless asked (BlockGunzipper's \
         checkCrcs is false), and htsjdk-bgzf does as htsjdk does",
    ),
    (
        "vcf: apply_vqsr_site_filtering.txt.gz:",
        "InvalidValue(\"VQSLOD\", InvalidFloat",
        "an input made malformed on purpose (VQSLOD=nonsense); htsjdk decodes INFO lazily and \
         the reference's answer to it is what the golden holds",
    ),
    (
        "vcf: mutect2_gvcf.txt.gz:",
        "InvalidValue(\"AS_MQ\", EmptyValue)",
        "the specification's ambiguity, written by GATK: an allele-specific value per allele with \
         the <NON_REF> one empty (`AS_MQ=57.52,`); the specification has no empty element and \
         noodles refuses it, the port writes the reference's bytes",
    ),
    (
        "vcf: ",
        " FT: htsjdk-rs \"[\\\"\\\"]\", noodles \"[Some(\\\"PASS\\\")]\"",
        "htsjdk-rs's bug: a genotype FT of PASS reads as Some(\"\") where htsjdk's \
         GenotypeBuilder.filters(emptyList) leaves null. is_filtered() agrees with htsjdk, the \
         field does not, and SelectVariants --select-genotype 'FT == \"PASS\"' diverges from \
         the oracle because of it",
    ),
];

#[test]
fn noodles_reads_what_htsjdk_rs_reads() {
    let mut report = Report::default();
    for golden in goldens() {
        let items = scan(&golden);
        let mut bams: Vec<Bam> = Vec::new();
        let mut bais: Vec<Found> = Vec::new();
        // A .tbi follows the BGZF VCF it indexes in every dump that carries both.
        let mut last_vcf_gz: Option<Found> = None;
        for (kind, found) in items {
            match kind {
                Kind::Bai => bais.push(found),
                Kind::Vcf => compare_vcf(&found, &found.bytes, &mut report),
                Kind::Bgzf => {
                    let ours = htsjdk_bgzf::read::decompress_all(&found.bytes);
                    let theirs = noodles_inflate(&found.bytes);
                    let inflated = match (ours, theirs) {
                        (Ok(o), Ok(t)) if o == t => o,
                        // A truncated or corrupt fixture that both refuse is an agreement.
                        (Err(_), Err(_)) => {
                            report.count("bgzf refused by both");
                            continue;
                        }
                        (o, t) => {
                            report.differ(
                                "bgzf",
                                &found.at,
                                format!(
                                    "htsjdk-rs {:?}, noodles {:?}",
                                    o.map(|v| v.len()),
                                    t.map(|v| v.len()).map_err(|e| e.to_string())
                                ),
                            );
                            continue;
                        }
                    };
                    report.count("bgzf");
                    if inflated.starts_with(b"BAM\x01") {
                        if let Some((header, records)) = compare_bam(&found, &inflated, &mut report)
                        {
                            bams.push((found, header, records));
                        }
                    } else if inflated.starts_with(b"TBI\x01") {
                        check_tabix(&found, last_vcf_gz.as_ref(), &mut report);
                    } else if inflated.starts_with(b"CSI\x01") {
                        report.count("csi");
                        if let Err(e) = csi::io::Reader::new(found.bytes.as_slice()).read_index() {
                            report.differ("csi", &found.at, format!("noodles refuses it: {e}"));
                        }
                    } else if inflated.starts_with(b"##fileformat=VCF") {
                        compare_vcf(&found, &inflated, &mut report);
                        last_vcf_gz = Some(found);
                    } else {
                        report.count("bgzf, other content");
                    }
                }
            }
        }
        for found in &bais {
            compare_bai(found, &bams, &mut report);
        }
    }

    let mut known: BTreeMap<usize, usize> = BTreeMap::new();
    let mut unknown = Vec::new();
    for d in &report.disagreements {
        match KNOWN
            .iter()
            .position(|(at, what, _)| d.starts_with(at) && d.contains(what))
        {
            Some(i) => *known.entry(i).or_default() += 1,
            None => unknown.push(d),
        }
    }
    eprintln!("read: {:#?}", report.read);
    for (i, (at, _, verdict)) in KNOWN.iter().enumerate() {
        eprintln!(
            "known, {} times: {at} {verdict}",
            known.get(&i).copied().unwrap_or(0)
        );
    }
    for d in unknown.iter().take(200) {
        eprintln!("{d}");
    }
    assert!(
        unknown.is_empty(),
        "{} disagreements nobody has judged",
        unknown.len()
    );
    let gone: Vec<&str> = (0..KNOWN.len())
        .filter(|i| !known.contains_key(i))
        .map(|i| KNOWN[i].0)
        .collect();
    assert!(
        gone.is_empty(),
        "known disagreements that no longer occur: {gone:?}"
    );
}
