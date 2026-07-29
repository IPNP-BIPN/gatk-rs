//! Conformance for the ported read filters against GATK 4.6.2.0.
//!
//! The golden is a decision matrix produced by `tools/readfilter-conformance/ReadFilterDump.java`
//! in the pinned container: one row per filter, one character per record, taken by the reference
//! through `SAMRecordToGATKReadAdapter`. The corpus travels in the same file, field by field, so
//! this test judges the records the reference judged rather than a reconstruction of them.
//!
//! Rows whose label carries parameters (`MappingQualityReadFilter(min=30,max=60)`) are instances
//! the reference built; the port rebuilds them from the label, so one list of instantiations
//! drives both sides.
//!
//! What this catches that a unit test does not: `NotProperlyPairedReadFilter` is
//! `isPaired() && !isProperlyPaired()`, not the negation of `ProperlyPairedReadFilter`. The first
//! version of the port used the negation, which keeps every unpaired read. The decision matrix
//! disagreed on five records of nineteen the first time it ran.

use std::io::Read;

use gatk_readfilter::{by_name, with_header, Parameterized, PORTED};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_filters.txt.gz");
    let file = std::fs::File::open(&path).unwrap_or_else(|e| {
        panic!(
            "{}: {e}. Regenerate with tools/conformance/run_suite.py",
            path.display()
        )
    });
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden is not valid gzip");
    text
}

/// The header the reference judged the corpus against.
///
/// It travels in the golden because the resolved filters read the library, sample, platform and
/// contig lengths out of it: a port given a different header would be answering a different
/// question and could agree by accident.
fn header(text: &str) -> SamHeader {
    let mut header = SamHeader::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "sq" => header
                .sequences
                .push(SequenceRecord::new(parts[2], parts[3].parse().unwrap())),
            "rg" => {
                let mut group = ReadGroup::new(parts[1]);
                for field in &parts[2..] {
                    let (key, value) = field.split_once('=').expect("an @RG field is KEY=value");
                    // "null" is how the dump prints an absent attribute; setting it would make the
                    // port match on a string the reference never had.
                    if value != "null" {
                        group.attributes.set(key, value);
                    }
                }
                header.read_groups.push(group);
            }
            _ => {}
        }
    }
    header
}

/// The corpus, in the order the reference judged it.
fn corpus(text: &str) -> Vec<BamRecord> {
    let mut records = Vec::new();
    for line in text.lines() {
        let mut parts = line.splitn(3, '\t');
        if parts.next() != Some("record") {
            continue;
        }
        let index: usize = parts
            .next()
            .expect("a record row has an index")
            .parse()
            .unwrap();
        let fields: Vec<&str> = parts
            .next()
            .expect("a record row has fields")
            .split('|')
            .collect();
        assert_eq!(
            fields.len(),
            11,
            "record {index} has {} fields",
            fields.len()
        );

        // Fields rather than a SAM line, because the corpus contains a record whose flags say
        // mapped while its reference is absent, which is one of the three criteria of
        // GATKRead.isUnmapped and exactly what htsjdk's reader refuses to parse
        // ("RNAME is not specified but flags indicate mapped"). Routing the corpus through SAM
        // text would drop the case the filter most needs.
        let mut record = BamRecord {
            read_name: fields[0].to_string(),
            flags: fields[1].parse().unwrap(),
            reference_index: fields[2].parse().unwrap(),
            alignment_start: fields[3].parse().unwrap(),
            mapping_quality: fields[4].parse().unwrap(),
            mate_reference_index: fields[6].parse().unwrap(),
            mate_alignment_start: fields[7].parse().unwrap(),
            inferred_insert_size: fields[8].parse().unwrap(),
            read_bases: fields[9].as_bytes().to_vec(),
            base_qualities: if fields[10].is_empty() {
                Vec::new()
            } else {
                fields[10].split(',').map(|q| q.parse().unwrap()).collect()
            },
            ..BamRecord::default()
        };
        if fields[5] != "*" {
            record.cigar = htsjdk_bam::text_parse::parse_cigar(fields[5])
                .unwrap_or_else(|e| panic!("record {index} cigar does not parse: {e:?}"));
        }

        assert_eq!(
            records.len(),
            index,
            "records are out of order in the golden"
        );
        records.push(record);
    }

    // Tags travel on their own rows: an OA value ends with a semicolon, so any in-line separator
    // would collide with the data it carries.
    for line in text.lines() {
        let mut parts = line.splitn(4, '\t');
        if parts.next() != Some("tag") {
            continue;
        }
        let index: usize = parts.next().unwrap().parse().unwrap();
        let name = parts.next().expect("a tag row has a name").as_bytes();
        let value = parts.next().expect("a tag row has a value");
        records[index].tags.insert(
            htsjdk_bam::Tag::new(&[name[0], name[1]]),
            htsjdk_bam::tag::TagValue::Str(value.to_string()),
        );
    }
    assert!(!records.is_empty(), "the golden carries no records");
    records
}

#[test]
fn every_filter_matches_the_reference_decision_for_decision() {
    let text = golden();
    let records = corpus(&text);
    let header = header(&text);

    let mut checked = 0;
    let mut compared = 0;
    for line in text.lines() {
        let mut parts = line.splitn(3, '\t');
        if parts.next() != Some("filter") {
            continue;
        }
        let name = parts.next().expect("a filter row has a name");
        let expected = parts.next().expect("a filter row has decisions");
        assert_eq!(
            expected.len(),
            records.len(),
            "{name}: the golden has {} decisions for {} records",
            expected.len(),
            records.len()
        );

        // A filter in the golden that the port does not implement is a failure, not a skip:
        // silently ignoring it is how a suite shrinks without anyone noticing.
        //
        // A parameterised filter carries its parameters in the label, and the port rebuilds the
        // instance from them, so the reference's own instantiation drives the comparison rather
        // than a second list on this side that could drift from it.
        let ours: String = if let Some(filter) = by_name(name) {
            records
                .iter()
                .map(|read| if filter(read) { '1' } else { '0' })
                .collect()
        } else if let Some(filter) = Parameterized::parse(name) {
            records
                .iter()
                .map(|read| if filter.test(read) { '1' } else { '0' })
                .collect()
        } else {
            // The header-dependent family: the label names the filter and its arguments, and the
            // header comes from the golden, so both sides resolve against the same @RG lines.
            let (label, args) = name.split_once('(').expect("a filter label");
            let args = args.strip_suffix(')').expect("a filter label ends with )");
            let values: Vec<String> = args
                .split_once('=')
                .map(|(_, list)| list)
                .unwrap_or("")
                .split('+')
                .filter(|v| !v.is_empty())
                .map(str::to_string)
                .collect();
            records
                .iter()
                .map(|read| {
                    let kept = match label {
                        "HasReadGroupWithHeader" => gatk_readfilter::has_read_group(read),
                        "AlignmentAgreesWithHeaderReadFilter" => {
                            with_header::alignment_agrees_with_header(read, &header)
                        }
                        "WellformedReadFilter" => with_header::wellformed(read, &header),
                        "LibraryReadFilter" => with_header::library(read, &header, &values),
                        "SampleReadFilter" => with_header::sample(read, &header, &values),
                        "PlatformReadFilter" => with_header::platform(read, &header, &values),
                        "PlatformUnitReadFilter" => {
                            with_header::platform_unit(read, &header, &values)
                        }
                        _ => panic!("{name} is in the golden but not ported; add it or remove it"),
                    };
                    if kept {
                        '1'
                    } else {
                        '0'
                    }
                })
                .collect()
        };

        if ours != expected {
            let differing: Vec<usize> = expected
                .chars()
                .zip(ours.chars())
                .enumerate()
                .filter(|(_, (a, b))| a != b)
                .map(|(i, _)| i)
                .collect();
            let names: Vec<&str> = differing
                .iter()
                .map(|&i| records[i].read_name.as_str())
                .collect();
            panic!(
                "{name} differs on {} of {} records: {names:?}\n  reference: {expected}\n  port     : {ours}",
                differing.len(),
                records.len()
            );
        }
        checked += 1;
        compared += records.len();
    }

    assert!(checked > 0, "the golden carries no filter rows");
    println!("{checked} filters, {compared} decisions, all identical");
}

/// Every ported filter must appear in the golden, or it is untested against the reference.
#[test]
fn no_ported_filter_is_missing_from_the_golden() {
    let text = golden();
    let in_golden: Vec<&str> = text
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            (parts.next() == Some("filter")).then(|| parts.next().unwrap())
        })
        .collect();

    let missing: Vec<&&str> = PORTED
        .iter()
        .filter(|name| !in_golden.contains(*name))
        .collect();
    assert!(
        missing.is_empty(),
        "ported but not in the golden, so never compared to the reference: {missing:?}. \
         Add them to ReadFilterDump.filters() and regenerate."
    );
}
