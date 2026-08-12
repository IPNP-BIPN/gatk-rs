//! Conformance for `CountBases`, `CountReads` and `FlagStat` against GATK 4.6.2.0, compared as
//! **text**.
//!
//! Golden from `tools/readfilter-conformance/CountingWalkersDump.java`. The three input fixtures and
//! their indexes travel in full, base64, so the port opens the same bytes.
//!
//! # What this suite is for
//!
//! The first three tools whose output is a number rather than a BAM:
//!
//!  * **their default filter is the engine's**, so a malformed read is not counted;
//!  * **`CountBases` counts bases, not span**;
//!  * **`FlagStat`'s percentages are computed in `float` and formatted `#0.00`**, HALF_EVEN;
//!  * **no reads at all is `NaN%`**, a real line and not a failure;
//!  * **`read2` is tested before `read1`**, so a read with both flags counts once, as read2.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::counting_walkers::{count_bases, count_reads, FlagStatus};
use gatk_tools::read_walker;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/counting_walkers.txt.gz"),
    )
}

fn pairs<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .filter_map(|rest| rest.split_once('\t'))
        .collect()
}

/// `count` rows carry `<tool>\t<label>\t<text>`.
fn counts(text: &str) -> Vec<(&str, &str, &str)> {
    pairs(text, "count")
        .into_iter()
        .filter_map(|(tool, rest)| {
            rest.split_once('\t')
                .map(|(label, value)| (tool, label, value))
        })
        .collect()
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn install_fixtures(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for (label, encoded) in pairs(text, "fixture") {
        std::fs::write(
            dir.join(format!("{label}.bam")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture bam");
    }
    for (label, encoded) in pairs(text, "fixtureindex") {
        if encoded == "absent" {
            continue;
        }
        std::fs::write(
            dir.join(format!("{label}.bai")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture index");
    }
}

/// A run's fixture, its interval and whether the default filter was left in place.
fn configuration(label: &str) -> (&str, Vec<&str>, bool) {
    let (fixture, suffix) = match label.split_once('-') {
        Some((fixture, suffix)) => (fixture, suffix),
        None => (label, ""),
    };
    match suffix {
        "" => (fixture, vec![], true),
        "interval" => (fixture, vec!["chr1:1-12"], true),
        "unfiltered" => (fixture, vec![], false),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The reads one labelled run traversed.
fn traverse(dir: &std::path::Path, label: &str) -> (Vec<BamRecord>, SamHeader) {
    let (fixture, intervals, wellformed) = configuration(label);
    let bam = dir.join(format!("{fixture}.bam"));
    let bai = dir.join(format!("{fixture}.bai"));
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    let parsed: Vec<SimpleInterval> = intervals
        .iter()
        .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
        .collect();

    let header_for_filter = header.clone();
    let filter: Box<dyn Fn(&BamRecord) -> bool> = if wellformed {
        Box::new(move |read: &BamRecord| with_header::wellformed(read, &header_for_filter))
    } else {
        Box::new(|_: &BamRecord| true)
    };

    let records = read_walker::traverse(&source, &parsed, filter.as_ref()).expect("the traversal");
    (records, header)
}

/// The contig a reference index names, which is what `read.getContig()` resolves to.
fn contig_of(header: &SamHeader, index: i32) -> Option<&str> {
    usize::try_from(index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.as_str())
}

#[test]
fn every_count_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-counting-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let rows = counts(&text);
    assert!(rows.len() >= 20, "three fixtures and eight runs each");

    let mut compared = 0;
    for (tool, label, expected) in &rows {
        let (records, header) = traverse(&dir, label);
        let ours = match *tool {
            "CountBases" => count_bases(&records).to_string(),
            "CountReads" => count_reads(&records).to_string(),
            "FlagStat" => {
                let mut status = FlagStatus::default();
                for record in &records {
                    let contig = contig_of(&header, record.reference_index);
                    let mate = contig_of(&header, record.mate_reference_index);
                    status.add(record, contig, mate);
                }
                status.to_text()
            }
            other => panic!("no tool {other}"),
        };
        assert_eq!(ours, unescape(expected), "{tool}/{label}");
        compared += 1;
    }
    println!("counting-walkers: {compared} counts compared");
}

/// The engine's default filter is what decides whether a malformed read is counted at all.
#[test]
fn the_default_filter_is_what_the_counts_rest_on() {
    let text = golden();
    let value = |tool: &str, label: &str| -> String {
        counts(&text)
            .into_iter()
            .find(|(t, l, _)| *t == tool && *l == label)
            .unwrap_or_else(|| panic!("no row {tool}/{label}"))
            .2
            .to_string()
    };
    // A file holding nothing but malformed reads counts nothing, and two with the filter off.
    assert_eq!(value("CountReads", "malformed"), "0");
    assert_eq!(value("CountReads", "malformed-unfiltered"), "2");
    // And the empty file's percentage is NaN rather than an error.
    assert!(unescape(&value("FlagStat", "empty")).contains("0 mapped (NaN%)"));
}
