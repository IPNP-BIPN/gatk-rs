//! `MarkDuplicates.doWork`, as `gatk MarkDuplicates` runs it.
//!
//! Which record of a duplicate set keeps `0x400` clear, and what the metrics say, are picard-rs's
//! (`picard-markduplicates`): this is the file around those decisions. The input is read twice,
//! as Picard reads it, one record at a time each pass: the first keeps a small summary per record
//! for the marking, the second writes every record again with its flag. Neither pass holds the
//! reads.
//!
//! What GATK changes is the command line, not the code: the `@PG` and the metrics header carry
//! Barclay's `--NAME value` command line, and `getVersion()` is GATK's.

use std::cell::RefCell;
use std::sync::Arc;

use gatk_barclay::Parser;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_tools::main_entry::{Thrown, PORT_FAILURE, PORT_LIMITATION};
use gatk_tools::read_walker_refusal::{is_block_compressed, BAM_MAGIC};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::{ProgramRecord, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use picard_markduplicates::{
    base_quality_score, duplication_metrics_text, Flags, Options, ReadEndsBuilder, ReadNameRegex,
    Record, ScoringStrategy, Span, TaggingPolicy, Tally,
};

use super::{assert_readable, assert_writable, thrown};
use crate::runners::{argument, arguments, flag, output_compression, scalar};

const DUPLICATE_READ: u16 = 0x400;

fn limitation(what: &str) -> Thrown {
    Thrown::non_user(
        PORT_LIMITATION,
        format!("{what} is a MarkDuplicates feature that this port does not carry yet. This message is the port's own and not GATK's."),
    )
}

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    }
}

/// `Span::of`, on this workspace's cigar: picard-markduplicates pins another htsjdk-rs.
fn span(record: &BamRecord) -> Span {
    let clip = |elements: &mut dyn Iterator<Item = &htsjdk_bam::cigar::CigarElement>| -> i32 {
        let mut clipped = 0;
        for element in elements {
            match element.op {
                Op::S | Op::H => clipped += element.length as i32,
                _ => break,
            }
        }
        clipped
    };
    Span {
        reference_length: record.cigar.reference_length(),
        leading_clip: clip(&mut record.cigar.elements.iter()),
        trailing_clip: clip(&mut record.cigar.elements.iter().rev()),
    }
}

/// `SamFileHeaderMerger.positiveFourDigitBase36Str`, whose digits are upper case.
fn base36(mut n: u32) -> String {
    const DIGITS: &[u8; 36] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let mut out = Vec::new();
    loop {
        out.push(DIGITS[(n % 36) as usize]);
        n /= 36;
        if n == 0 {
            break;
        }
    }
    out.reverse();
    String::from_utf8(out).expect("base36 digits are ASCII")
}

/// `getChainedPgIds`: one `@PG` per program id the records carried, in `HashSet` order (`null`
/// hashes to bucket 0), each id made unique by `PgIdGenerator` and chained to the id it replaces.
fn chained_pg_ids(
    header: &mut SamHeader,
    seen: Vec<Option<String>>,
    record_id: &str,
    version: &str,
    command_line: &str,
    name: Option<&str>,
) -> Vec<(Option<String>, String)> {
    let mut capacity = 16usize;
    while seen.len() > capacity * 3 / 4 {
        capacity *= 2;
    }
    let mut keyed: Vec<(usize, Option<String>)> = seen
        .into_iter()
        .map(|id| {
            let h = id
                .as_deref()
                .map_or(0, gatk_engine::java_hash::string_hash_code) as u32;
            (((h ^ (h >> 16)) as usize) & (capacity - 1), id)
        })
        .collect();
    keyed.sort_by_key(|(bucket, _)| *bucket);

    let mut taken: std::collections::HashSet<String> =
        header.programs.iter().map(|pg| pg.id.clone()).collect();
    let mut counter = taken.len() as u32;
    let mut chained = Vec::new();
    for (_, existing) in keyed {
        let id = if taken.insert(record_id.to_string()) {
            counter += 1;
            record_id.to_string()
        } else {
            loop {
                let candidate = format!("{record_id}.{}", base36(counter));
                counter += 1;
                if taken.insert(candidate.clone()) {
                    break candidate;
                }
            }
        };
        let mut pg = ProgramRecord::new(&id);
        pg.attributes.set("VN", version);
        pg.attributes.set("CL", command_line);
        if let Some(name) = name {
            pg.attributes.set("PN", name);
        }
        if let Some(previous) = &existing {
            pg.attributes.set("PP", previous);
        }
        header.programs.push(pg);
        chained.push((existing, id));
    }
    chained
}

/// `SAMRecordCoordinateComparator.fileOrderCompare` on (reference index, start): an unmapped
/// record (index -1) sorts last.
fn file_order(left: (i32, i32), right: (i32, i32)) -> std::cmp::Ordering {
    let last = |index: i32| if index == -1 { i32::MAX } else { index };
    (last(left.0), left.1).cmp(&(last(right.0), right.1))
}

/// The input, read from the start once per pass: a BAM through `ReadsDataSource`, one record
/// at a time and never inflated whole; anything else as `GatherBamFiles` reads it.
enum Input {
    Bam(ReadsDataSource),
    // shortcut: a SAM or other text input is read whole, as the other Picard runners read it;
    // stream it if a large SAM input ever matters.
    Read(SamHeader, Vec<BamRecord>),
}

impl Input {
    fn open(input: &str) -> Result<Input, Thrown> {
        let path = std::path::Path::new(input);
        let refused = |error: ReadsError| Thrown::user(format!("{error:?}"));
        let bytes = gatk_engine::reads::file_start(path).ok_or_else(|| {
            Thrown::user(gatk_tools::read_walker_refusal::cannot_read(input, false))
        })?;
        let binary = if is_block_compressed(&bytes) {
            ReadsDataSource::check(path)
                .and_then(|()| gatk_engine::reads::inflated_start(path))
                .map(|inflated| inflated.starts_with(&BAM_MAGIC))
                .unwrap_or(false)
        } else {
            bytes.starts_with(&BAM_MAGIC)
        };
        if binary {
            return ReadsDataSource::open_checked(path)
                .map(Input::Bam)
                .map_err(refused);
        }
        let (header, records) = super::gather_bam_files::read_reads(input)?;
        Ok(Input::Read(header, records))
    }

    fn header(&self) -> &SamHeader {
        match self {
            Input::Bam(source) => source.header(),
            Input::Read(header, _) => header,
        }
    }

    fn records(
        &self,
    ) -> Result<Box<dyn Iterator<Item = Result<BamRecord, ReadsError>> + '_>, ReadsError> {
        Ok(match self {
            Input::Bam(source) => Box::new(source.records_all()?),
            Input::Read(_, records) => Box::new(records.iter().cloned().map(Ok)),
        })
    }
}

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let inputs = arguments(parser, "INPUT");
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    let metrics_file = argument(parser, "METRICS_FILE").unwrap_or_default();
    // shortcut: one input, a BAM out, no UMIs or barcodes, no flow mode; each is refused by name,
    // and each is a few lines over picard-rs's bin once a run needs it.
    let [input] = inputs.as_slice() else {
        return Err(limitation("More than one INPUT"));
    };
    if !output.ends_with(".bam") {
        return Err(limitation("An OUTPUT other than a .bam"));
    }
    for unported in [
        "BARCODE_TAG",
        "READ_ONE_BARCODE_TAG",
        "READ_TWO_BARCODE_TAG",
        "MOLECULAR_IDENTIFIER_TAG",
    ] {
        if argument(parser, unported).is_some_and(|value| value != "null") {
            return Err(limitation(unported));
        }
    }
    if flag(parser, "FLOW_MODE") {
        return Err(limitation("FLOW_MODE"));
    }
    assert_readable(input)?;
    assert_writable(&output)?;
    assert_writable(&metrics_file)?;

    let read_name_regex = ReadNameRegex::from_arg(Some(
        // Its default is the "optimized capture" text, so `None` is a `null`: no locations.
        scalar(parser, "READ_NAME_REGEX")
            .as_deref()
            .unwrap_or("null"),
    ))
    .map_err(|message| {
        let (class, message) = message.split_once(": ").unwrap_or(("", &message));
        thrown(class, message.to_string())
    })?;
    let integer = |name: &str, default: i32| {
        scalar(parser, name)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let options = Options {
        scoring: match scalar(parser, "DUPLICATE_SCORING_STRATEGY").as_deref() {
            Some("TOTAL_MAPPED_REFERENCE_LENGTH") => ScoringStrategy::TotalMappedReferenceLength,
            Some("RANDOM") => ScoringStrategy::Random,
            _ => ScoringStrategy::SumOfBaseQualities,
        },
        remove_duplicates: flag(parser, "REMOVE_DUPLICATES"),
        remove_sequencing_duplicates: flag(parser, "REMOVE_SEQUENCING_DUPLICATES"),
        tagging_policy: match scalar(parser, "TAGGING_POLICY").as_deref() {
            Some("OpticalOnly") => TaggingPolicy::OpticalOnly,
            Some("All") => TaggingPolicy::All,
            _ => TaggingPolicy::DontTag,
        },
        clear_dt: flag(parser, "CLEAR_DT"),
        optical_duplicate_pixel_distance: integer("OPTICAL_DUPLICATE_PIXEL_DISTANCE", 100),
        parse_read_names: !read_name_regex.is_off(),
        read_name_regex,
        max_optical_duplicate_set_size: scalar(parser, "MAX_OPTICAL_DUPLICATE_SET_SIZE")
            .and_then(|value| value.parse().ok())
            .unwrap_or(300_000),
        barcode_tag: None,
        read_barcodes: false,
        assume_mate_cigar: false,
    };
    // `ASSUME_SORTED=true` is `ASSUME_SORT_ORDER=coordinate` by the time `openInputs` reads it.
    let assumed = scalar(parser, "ASSUME_SORT_ORDER")
        .filter(|order| order != "null")
        .or_else(|| flag(parser, "ASSUME_SORTED").then(|| "coordinate".to_string()));

    let source = Input::open(input)?;
    let mut header = source.header().clone();
    if let Some(order) = assumed.as_deref() {
        header.set_sort_order(order);
    }
    let header_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    if !matches!(header_order.as_str(), "coordinate" | "queryname") {
        return Err(thrown(
            "picard.PicardException",
            format!(
                "This program requires input that are either coordinate or query sorted \
                 (according to the header, or at least ASSUME_SORT_ORDER and the content.) Found \
                 ASSUME_SORT_ORDER={} and header sortorder={header_order}",
                assumed.as_deref().unwrap_or("null")
            ),
        ));
    }
    let queryname = header_order == "queryname";

    // A duplicate set is cut by library, and `closeEnough` compares read-group indexes. One shared
    // string per read group, so each summary holds a pointer rather than a copy.
    let unknown_library: Arc<str> = Arc::from("Unknown Library");
    let libraries: Vec<Arc<str>> = header
        .read_groups
        .iter()
        .map(|group| match group.attributes.get("LB") {
            Some(library) => Arc::from(library),
            None => unknown_library.clone(),
        })
        .collect();
    let library_of = |record: &BamRecord| -> (Arc<str>, i32) {
        let id = match record.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(s)) => Some(s.as_str()),
            _ => None,
        };
        match id.and_then(|id| header.read_groups.iter().position(|group| group.id == id)) {
            Some(position) => (libraries[position].clone(), position as i32),
            None => (unknown_library.clone(), -1),
        }
    };
    let summary = |record: &BamRecord| -> Record {
        let (library, read_group) = library_of(record);
        Record {
            name: record.read_name.clone(),
            flags: record.flags,
            reference_index: record.reference_index,
            alignment_start: record.alignment_start,
            span: span(record),
            read_length: record.read_bases.len() as u32,
            base_quality_score: base_quality_score(&record.base_qualities),
            mate_reference_index: record.mate_reference_index,
            library,
            read_group,
            barcode: None,
            existing_dt: string_tag(record, b"DT"),
            mate_span: None,
            mate_alignment_start: record.mate_alignment_start,
        }
    };
    let malformed = |error: ReadsError| Thrown::user(format!("{error:?}"));

    // The first pass: a read end per record, and the program ids the records carry.
    let mut seen: Vec<Option<String>> = Vec::new();
    let mut ends = ReadEndsBuilder::new(&options, queryname);
    for record in source.records().map_err(malformed)? {
        let record = record.map_err(malformed)?;
        let id = string_tag(&record, b"PG");
        if !seen.contains(&id) {
            seen.push(id);
        }
        // The first unplaced unmapped read ends the pass, after its `PG` id is collected.
        if record.flags & 0x4 != 0 && record.reference_index == -1 && !queryname {
            break;
        }
        if record.flags & (0x1 | 0x4 | 0x8 | 0x100 | 0x800) == 0x1 && library_of(&record).1 == -1 {
            return Err(thrown(
                "java.lang.NullPointerException",
                "Cannot invoke \"htsjdk.samtools.SAMReadGroupRecord.getReadGroupId()\" because \
                 the return value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null"
                    .to_string(),
            ));
        }
        ends.push(&summary(&record));
    }
    let mut decisions = ends.decide();

    // `createOutHeader`: assuming queryname says the file is queryname GROUPED.
    let mut out_header = header.clone();
    if assumed.as_deref() == Some("queryname") {
        out_header.set_sort_order("unknown");
        out_header.set_group_order("query");
    }
    for comment in arguments(parser, "COMMENT") {
        out_header.comments.push(format!("@CO\t{comment}"));
    }
    let (level, deflater) = output_compression(parser);
    // `COMPRESSION_LEVEL` defaults to `Defaults.COMPRESSION_LEVEL`, which GATK's configuration
    // sets before the tool is built; the declarations were measured on a tool built without it.
    let mut command_line = crate::command_line::expanded("MarkDuplicates", parser);
    if super::compression_level(parser).is_none() {
        command_line = command_line.replacen(
            " --COMPRESSION_LEVEL 5 ",
            &format!(" --COMPRESSION_LEVEL {level} "),
            1,
        );
    }
    // An argument set to `null` reads as `None`. PROGRAM_RECORD_ID and PROGRAM_GROUP_NAME default
    // to the tool's name, so `None` is a `null` that turns them off; the version and the command
    // line default to `null`, which means `getVersion()` and `getCommandLine()`.
    let chained = scalar(parser, "PROGRAM_RECORD_ID").map(|record_id| {
        chained_pg_ids(
            &mut out_header,
            seen,
            &record_id,
            &scalar(parser, "PROGRAM_GROUP_VERSION")
                .unwrap_or_else(|| format!("Version:{}", crate::TOOLKIT_VERSION)),
            &scalar(parser, "PROGRAM_GROUP_COMMAND_LINE").unwrap_or_else(|| command_line.clone()),
            scalar(parser, "PROGRAM_GROUP_NAME").as_deref(),
        )
    });
    let tag_duplicate_set_members = flag(parser, "TAG_DUPLICATE_SET_MEMBERS");
    let add_pg_tag = flag(parser, "ADD_PG_TAG_TO_READS");
    let coordinate_out = out_header.attributes.get("SO") == Some("coordinate");

    let optical_pairs = std::mem::take(&mut decisions.optical_pairs);
    let all_sets = std::mem::take(&mut decisions.all_sets);
    let optical_sets = std::mem::take(&mut decisions.optical_sets);
    let non_optical_sets = std::mem::take(&mut decisions.non_optical_sets);
    let mut flags = Flags::new(decisions, queryname);
    let tally = RefCell::new(Tally::default());
    // A refusal met mid-stream stops the writer through a placeholder error, and is reported in
    // its place.
    let stopped: RefCell<Option<Thrown>> = RefCell::new(None);
    let mut previous: Option<(i32, i32)> = None;
    let target = format!("file://{}", super::absolute(&output));

    let marked = source.records().map_err(malformed)?.filter_map(|record| {
        let record = match record {
            Ok(record) => record,
            Err(error) => return Some(Err(error)),
        };
        let summary = summary(&record);
        let (duplicate, optical, duplicate_set) = flags.next(
            &summary.name,
            summary.secondary_or_supplementary() || summary.unmapped(),
        );
        let (duplicate_type, kept) = tally
            .borrow_mut()
            .add(&summary, &options, duplicate, optical);
        if !kept {
            return None;
        }
        // `SAMFileWriterImpl.assertPresorted`: the writer is presorted, so the record that goes
        // backwards is the one that throws, and `isSorted` prints its key twice.
        if coordinate_out {
            let here = (record.reference_index, record.alignment_start);
            if previous
                .is_some_and(|before| file_order(before, here) == std::cmp::Ordering::Greater)
            {
                let contig = usize::try_from(here.0)
                    .ok()
                    .and_then(|i| out_header.sequences.get(i))
                    .map_or("*", |s| s.name.as_str());
                let key = format!("{contig}:{}", here.1);
                *stopped.borrow_mut() = Some(thrown(
                    "java.lang.IllegalArgumentException",
                    format!(
                        "Alignments added out of order in SAMFileWriterImpl.addAlignment for \
                         {target}. Sort order is coordinate. Offending records are at [{key}] \
                         and [{key}]"
                    ),
                ));
                return Some(Err(ReadsError::Malformed(String::new())));
            }
            previous = Some(here);
        }

        let original_flags = record.flags;
        let mut written = record;
        written.flags &= !DUPLICATE_READ;
        if duplicate {
            written.flags |= DUPLICATE_READ;
        }
        written.tags.remove(Tag::new(b"DT"));
        if let Some(code) = duplicate_type {
            written.tags.insert(Tag::new(b"DT"), TagValue::Str(code));
        }
        if tag_duplicate_set_members && original_flags & (0x4 | 0x100 | 0x800) == 0 {
            if let Some((representative, size)) = duplicate_set {
                written
                    .tags
                    .insert(Tag::new(b"DI"), TagValue::Int(representative as i64));
                written
                    .tags
                    .insert(Tag::new(b"DS"), TagValue::Int(size as i64));
            }
        }
        // `rec.setAttribute(PG, chainedPgIds.get(rec.getStringAttribute(PG)))`.
        if let (true, Some(chained)) = (add_pg_tag, &chained) {
            let old = string_tag(&written, b"PG");
            match chained
                .iter()
                .find(|(id, _)| id.as_deref() == old.as_deref())
            {
                Some((_, id)) => {
                    written
                        .tags
                        .insert(Tag::new(b"PG"), TagValue::Str(id.clone()));
                }
                None => {
                    written.tags.remove(Tag::new(b"PG"));
                }
            }
        }
        Some(Ok(written))
    });

    let create_index = flag(parser, "CREATE_INDEX") && coordinate_out;
    let create_md5 = flag(parser, "CREATE_MD5_FILE");
    let written = gatk_tools::sam_output::write_bam_file(
        std::path::Path::new(&output),
        &out_header,
        marked,
        create_index,
        level,
        deflater,
        create_md5,
    );
    if let Some(refusal) = stopped.take() {
        return Err(refusal);
    }
    let written = written
        .map_err(malformed)?
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    let write = |path: &str, bytes: &[u8]| {
        std::fs::write(path, bytes)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))
    };
    if let Some(bai) = written.index {
        // `IOUtil.basename(output) + ".bai"`: the output's last extension replaced.
        write(
            &std::path::Path::new(&output)
                .with_extension("bai")
                .display()
                .to_string(),
            &bai,
        )?;
    }
    if let Some(md5) = written.md5 {
        write(&format!("{output}.md5"), md5.as_bytes())?;
    }

    let text = duplication_metrics_text(
        &[
            command_line,
            format!("Started on: {}", crate::runners::java_display_now()),
        ],
        &tally.into_inner().finish(Some(&optical_pairs)),
        &header
            .read_groups
            .iter()
            .map(|group| {
                group
                    .attributes
                    .get("LB")
                    .unwrap_or("Unknown Library")
                    .to_string()
            })
            .collect::<Vec<_>>(),
        &all_sets,
        &optical_sets,
        &non_optical_sets,
    );
    write(&metrics_file, text.as_bytes())?;
    Ok(0)
}

#[cfg(test)]
mod tests {
    #[test]
    fn every_argument_converts() {
        let list = gatk_tools::tool_declarations::declarations("MarkDuplicates").unwrap();
        let missing = crate::definitions::missing(list);
        assert!(missing.is_empty(), "{missing:?}");
    }

    #[test]
    fn base36_is_htsjdks_upper_case() {
        assert_eq!(super::base36(0), "0");
        assert_eq!(super::base36(16), "G");
        assert_eq!(super::base36(36), "10");
    }
}
