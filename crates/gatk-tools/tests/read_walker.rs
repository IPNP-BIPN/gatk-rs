//! Conformance for the `ReadWalker` traversal against GATK 4.6.2.0.
//!
//! The golden was produced by a real `ReadWalker` subclass run through the real command line, so
//! what is replayed here is the traversal a tool gets. The fixture BAM, its index and the FASTA
//! all travel inside it.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use gatk_readfilter::{not_duplicate, with_header};
use gatk_tools::read_walker;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_walker.txt.gz"),
    )
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
            Some(other) => out.push(other),
            None => {}
        }
    }
    out
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

/// The traversals the harness ran, by the label it wrote them under.
///
/// The table is here rather than parsed out of the golden because a label like `all-nofilter` is
/// a *configuration*, not coordinates: it says which filters the tool was run with, and there is
/// nothing in the row to derive that from.
fn configuration(label: &str) -> (Vec<&'static str>, bool, bool) {
    // (interval strings, default filters enabled, reference available)
    match label {
        "all" => (vec![], true, true),
        "chr1" => (vec!["chr1"], true, true),
        "chr1:1-60" => (vec!["chr1:1-60"], true, true),
        "chr1:100-160" => (vec!["chr1:100-160"], true, true),
        "chr1:1-100+101-200" => (vec!["chr1:1-100", "chr1:101-200"], true, true),
        "chr2" => (vec!["chr2"], true, true),
        "all-nofilter" => (vec![], false, true),
        "all-nodup" => (vec![], true, true),
        "all-noref" => (vec![], true, false),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_traversal_hands_apply_what_the_reference_hands_it() {
    let text = golden();

    let dir = std::env::temp_dir().join(format!("gatk-rs-readwalker-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    let fasta = dir.join("ref.fasta");
    std::fs::write(&bam, corpus::decode_base64(field(&text, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&text, "bai"))).unwrap();
    std::fs::write(&fasta, unescape(field(&text, "fasta"))).unwrap();
    std::fs::write(dir.join("ref.fasta.fai"), unescape(field(&text, "fai"))).unwrap();

    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    // The rows the reference produced, grouped by label and in order.
    let mut expected: Vec<(String, Vec<String>)> = Vec::new();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("apply\t") else {
            continue;
        };
        let mut parts = rest.splitn(3, '\t');
        let label = parts.next().expect("label").to_string();
        let _index = parts.next();
        let payload = parts.next().unwrap_or("").to_string();
        match expected.last_mut() {
            Some((last, rows)) if *last == label => rows.push(payload),
            _ => expected.push((label, vec![payload])),
        }
    }
    assert!(!expected.is_empty(), "the golden carries no apply rows");

    let mut compared = 0;
    for (label, rows) in &expected {
        let (interval_strings, default_filters, has_reference) = configuration(label);
        let intervals: Vec<SimpleInterval> = interval_strings
            .iter()
            .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
            .collect();

        // The default filter of a ReadWalker is WellformedReadFilter alone; `all-nodup` adds
        // NotDuplicateReadFilter on top of it, and `all-nofilter` disables the defaults.
        let header_for_filter = header.clone();
        let filter: Box<dyn Fn(&BamRecord) -> bool> = match (default_filters, label.as_str()) {
            (false, _) => Box::new(|_: &BamRecord| true),
            (true, "all-nodup") => Box::new(move |read: &BamRecord| {
                with_header::wellformed(read, &header_for_filter) && not_duplicate(read)
            }),
            (true, _) => {
                Box::new(move |read: &BamRecord| with_header::wellformed(read, &header_for_filter))
            }
        };

        let mut reference = if has_reference {
            Some(ReferenceFileSource::open(&fasta).expect("the fixture reference opens"))
        } else {
            None
        };
        let applied = read_walker::traverse_with_reference(
            &source,
            reference.as_mut(),
            &intervals,
            filter.as_ref(),
        )
        .expect("the traversal runs");

        let ours: Vec<String> = applied
            .into_iter()
            .map(|mut applied| {
                let window = match applied.context.window() {
                    None => "null".to_string(),
                    Some(window) => format!("{}:{}-{}", window.contig, window.start, window.end),
                };
                let bases = match reference.as_mut() {
                    Some(source) => applied.context.bases(source).expect("bases"),
                    None => Vec::new(),
                };
                format!(
                    "{}|{}|{}|{}|{}|{}",
                    applied.read.read_name,
                    gatk_engine::read_utils::start(&applied.read),
                    applied.read.cigar.to_text(),
                    applied.read.flags,
                    window,
                    String::from_utf8(bases).expect("ASCII bases"),
                )
            })
            .collect();

        assert_eq!(&ours, rows, "{label}");
        compared += ours.len();
    }

    std::fs::remove_dir_all(&dir).ok();
    println!(
        "{compared} apply calls over {} traversals, all identical",
        expected.len()
    );
}
