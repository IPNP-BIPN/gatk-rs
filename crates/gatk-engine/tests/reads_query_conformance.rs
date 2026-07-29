//! Conformance for `ReadsDataSource` against GATK 4.6.2.0.
//!
//! The fixture travels in the golden: the BAM and its `.bai` are base64 in the dump, written back
//! to a temporary directory and queried, so the port queries exactly the bytes the reference
//! queried. Each case carries a `spec` row saying what was asked, so this replays the reference's
//! own queries rather than a second copy of the harness's table.
//!
//! What is compared is the whole answer, in order: every returned read's name, start, cigar,
//! flags and assigned position. A count would hide an ordering difference, and the order is the
//! file's, which is what a walker sees.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::ReadsDataSource;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reads_query.txt.gz"),
    )
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

fn parse_interval(text: &str) -> SimpleInterval {
    let (contig, span) = text.split_once(':').expect("contig:start-end");
    let (start, end) = span.split_once('-').expect("start-end");
    SimpleInterval {
        contig: contig.to_string(),
        start: start.parse().expect("start"),
        end: end.parse().expect("end"),
    }
}

/// `ReadsDataSourceDump.describe`: name, start, cigar, flags, assigned contig, assigned start.
///
/// `getStart` is `UNSET_POSITION` for an unmapped read and the harness prints `-1` for it; the
/// *assigned* position is what the record carries, which is how a mate-placed unmapped read shows
/// both that it has no span and that it sits at a coordinate.
fn describe(record: &BamRecord, header: &htsjdk_bam::header::SamHeader) -> String {
    let start = if gatk_engine::read::is_unmapped(record) {
        -1
    } else {
        record.alignment_start
    };
    let contig = if record.reference_index < 0 {
        "*".to_string()
    } else {
        header.sequences[record.reference_index as usize]
            .name
            .clone()
    };
    format!(
        "{}|{}|{}|{}|{}|{}",
        record.read_name,
        start,
        record.cigar.to_text(),
        record.flags,
        contig,
        record.alignment_start,
    )
}

#[test]
fn every_query_returns_the_records_the_reference_returns() {
    let text = golden();

    let dir = std::env::temp_dir().join(format!("gatk-rs-readsconf-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    std::fs::write(&bam, corpus::decode_base64(field(&text, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&text, "bai"))).unwrap();

    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    // The answers, keyed by label, as the reference gave them.
    let mut answers = std::collections::HashMap::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("query\t") {
            let (label, payload) = rest.split_once('\t').unwrap_or((rest, ""));
            answers.insert(label.to_string(), payload.to_string());
        }
    }

    let mut compared = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("spec\t") else {
            continue;
        };
        let (label, spec) = rest.split_once('\t').expect("spec payload");
        let (mode, intervals) = spec.split_once('|').expect("mode|intervals");
        let intervals: Vec<SimpleInterval> = if intervals.is_empty() {
            Vec::new()
        } else {
            intervals.split(',').map(parse_interval).collect()
        };

        let ours = match mode {
            // An unbounded traversal is not a query with no intervals: it is every record in the
            // file, unplaced reads included, which is what a walker with no -L sees.
            "traverse" if intervals.is_empty() => source.iter_all().map(|records| {
                records
                    .iter()
                    .map(|r| describe(r, &header))
                    .collect::<Vec<_>>()
            }),
            // A query and a bounded traversal reach the same filter; the traversal adds the
            // unplaced reads at the end when it was asked for them.
            "query" | "traverse" => source.query(&intervals).map(|records| {
                records
                    .iter()
                    .map(|r| describe(r, &header))
                    .collect::<Vec<_>>()
            }),
            "traverse+unmapped" => source.query(&intervals).and_then(|records| {
                let mut described: Vec<String> =
                    records.iter().map(|r| describe(r, &header)).collect();
                described.extend(
                    source
                        .query_unmapped()?
                        .iter()
                        .map(|r| describe(r, &header)),
                );
                Ok(described)
            }),
            "unmapped" => source.query_unmapped().map(|records| {
                records
                    .iter()
                    .map(|r| describe(r, &header))
                    .collect::<Vec<_>>()
            }),
            other => panic!("unknown spec mode {other}"),
        };

        let expected = &answers[label];
        let ours = match ours {
            // The harness prints a single `E` segment where the reference threw.
            Err(_) => "E".to_string(),
            Ok(records) => records.join("\\n"),
        };
        assert_eq!(&ours, expected, "{label} ({spec})");
        compared += 1;
    }

    std::fs::remove_dir_all(&dir).ok();
    assert!(compared > 0, "the golden carries no queries");
    println!("{compared} queries, all identical to the reference");
}
