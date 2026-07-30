//! Conformance for `AssemblyRegion` against the oracle.
//!
//! Goldens from `tools/readfilter-conformance/AssemblyRegionDump.java`.
//!
//! Three things this suite exists to pin, each of which a reasonable port gets wrong:
//!
//! * the padding constructor reports an interval it cannot place on the contig as a **null padded
//!   span**, not as a padding error;
//! * `trim(span, padding)` disagrees with its own javadoc, and the `javadoc-trim` case is the
//!   javadoc's worked example run through the code;
//! * trimming **sorts** the reads it keeps, so the order of a trimmed region is the comparator's
//!   and not the order the reads went in.

use gatk_corpus as corpus;
use gatk_engine::assembly_region::{render_interval, AssemblyRegion, RegionError};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::read_utils;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/assembly_region.txt.gz"),
    )
}

fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
    SimpleInterval::new(contig, start, end).expect("a valid interval")
}

/// The dump's `E:<class>:<message>`.
fn render_error(error: &RegionError) -> String {
    format!("E:{}:{}", error.class(), error.message())
}

/// The dump's `describe`.
fn describe(region: &AssemblyRegion) -> String {
    format!(
        "{}|{}|{}|{}",
        render_interval(region.span()),
        render_interval(region.padded_span()),
        region.is_active(),
        region.size()
    )
}

fn describe_result(result: &Result<AssemblyRegion, RegionError>) -> String {
    match result {
        Ok(region) => describe(region),
        Err(error) => render_error(error),
    }
}

/// The dump's `read` row.
fn render_read(record: &BamRecord) -> String {
    format!(
        "{}|{}|{}|{}",
        record.read_name,
        read_utils::start(record),
        record.cigar.to_text(),
        String::from_utf8_lossy(&record.read_bases)
    )
}

fn row<'a>(text: &'a str, prefix: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(prefix))
        .unwrap_or_else(|| panic!("the golden carries no row {prefix:?}"))
}

fn check_reads(text: &str, label: &str, region: &AssemblyRegion) {
    let expected: Vec<&str> = text
        .lines()
        .filter_map(|line| line.strip_prefix(&format!("read\t{label}\t")))
        .map(|rest| rest.split_once('\t').expect("an index and a read").1)
        .collect();
    let ours: Vec<String> = region.reads().iter().map(render_read).collect();
    assert_eq!(ours, expected, "{label}: the reads the region holds");
}

/// Label, contig, start, end, padding: the `ctor` cases, in the dump's order.
const CTOR: &[(&str, &str, i32, i32, i32)] = &[
    ("pad-inside", "chr1", 500, 600, 100),
    ("pad-zero", "chr1", 500, 600, 0),
    ("pad-off-front", "chr1", 10, 20, 100),
    ("pad-off-back", "chr1", 1950, 1990, 100),
    ("pad-whole-contig", "chr1", 500, 600, 5000),
    ("pad-unknown-contig", "chrX", 500, 600, 10),
    ("pad-one-base", "chr1", 500, 500, 10),
    ("pad-active-past-contig", "chr1", 2500, 2600, 10),
];

/// Label, active span, padded span: the `ctor` cases built from two intervals.
const CTOR_PAIR: &[(&str, &str, i32, i32, &str, i32, i32)] = &[
    ("pair-asymmetric", "chr1", 500, 600, "chr1", 400, 1000),
    ("pair-equal", "chr1", 500, 600, "chr1", 500, 600),
    ("pair-not-contained", "chr1", 500, 600, "chr1", 550, 700),
    ("pair-other-contig", "chr1", 500, 600, "chr2", 400, 700),
];

/// Label, span, padding: the single-span `trim` cases, all applied to the `wide` region.
const TRIM: &[(&str, &str, i32, i32, i32)] = &[
    ("trim-inside", "chr1", 500, 1000, 50),
    ("trim-javadoc-example", "chr1", 150, 225, 25),
    ("trim-zero-padding", "chr1", 500, 1000, 0),
    ("trim-padding-past-region", "chr1", 500, 1000, 5000),
    ("trim-single-base", "chr1", 1000, 1000, 10),
    ("trim-partly-outside", "chr1", 1800, 2500, 10),
    ("trim-disjoint", "chr1", 1, 50, 0),
    ("trim-other-contig", "chr2", 100, 200, 10),
];

const TRIM_PAIR: &[(&str, &str, i32, i32, &str, i32, i32)] = &[
    ("trimpair-inside", "chr1", 500, 1000, "chr1", 400, 1100),
    (
        "trimpair-padded-not-containing",
        "chr1",
        500,
        1000,
        "chr1",
        600,
        1100,
    ),
    (
        "trimpair-padded-beyond-original",
        "chr1",
        500,
        1000,
        "chr1",
        1,
        2000,
    ),
];

/// Rebuild one of the dump's named regions: every trim case starts from an untouched copy, because
/// `trim` clips the reads of the region it is called on.
fn build(label: &str, header: &SamHeader, records: &[BamRecord]) -> AssemblyRegion {
    let (contig, start, end, padding) = match label {
        "wide" => ("chr1", 100, 1900, 100),
        "javadoc" => ("chr1", 100, 200, 50),
        "narrow" => ("chr1", 1000, 1010, 5),
        other => panic!("unknown region {other}"),
    };
    let mut region =
        AssemblyRegion::with_padding(interval(contig, start, end), true, padding, header)
            .expect("the dump builds this region");
    let padded = region.padded_span().clone();
    let mut reads: Vec<BamRecord> = records
        .iter()
        // The three refusals `add` makes, applied up front, exactly as the dump applies them: an
        // all-insertion cigar consumes no reference, so 10I at 1790 ends at 1789 and cannot be an
        // interval at all.
        .filter(|record| {
            !gatk_engine::read::is_unmapped(record)
                && read_utils::start(record) >= 1
                && read_utils::end(record) >= read_utils::start(record)
        })
        .filter(|record| {
            let contig = header
                .sequences
                .get(record.reference_index as usize)
                .map(|sequence| sequence.name.as_str())
                .unwrap_or("*");
            padded.overlaps(contig, read_utils::start(record), read_utils::end(record))
        })
        .cloned()
        .collect();
    reads.sort_by(read_utils::compare_read_coordinate);
    region
        .add_all(reads, header)
        .expect("the sorted reads are accepted");
    region
}

#[test]
fn every_comparison_matches_the_reference() {
    let text = golden();
    let records = corpus::records(&text);
    let mut rows = 0;

    for line in text.lines() {
        let Some(rest) = line.strip_prefix("cmp\t") else {
            continue;
        };
        let mut parts = rest.split('\t');
        let i: usize = parts.next().expect("i").parse().expect("a number");
        let j: usize = parts.next().expect("j").parse().expect("a number");
        let expected = parts.next().expect("a sign");
        let ours = match read_utils::compare_read_coordinate(&records[i], &records[j]) {
            std::cmp::Ordering::Less => "-1",
            std::cmp::Ordering::Equal => "0",
            std::cmp::Ordering::Greater => "1",
        };
        assert_eq!(
            ours, expected,
            "comparing record {i} ({}) with {j} ({})",
            records[i].read_name, records[j].read_name
        );
        rows += 1;
    }

    assert!(rows > 0, "the golden carries no comparisons");
    println!("{rows} comparisons identical");
}

#[test]
fn every_construction_matches_the_reference() {
    let text = golden();
    let header = corpus::header(&text);

    for (label, contig, start, end, padding) in CTOR {
        let result =
            AssemblyRegion::with_padding(interval(contig, *start, *end), true, *padding, &header);
        assert_eq!(
            describe_result(&result),
            row(&text, &format!("ctor\t{label}\t")),
            "{label}"
        );
    }

    for (label, contig, start, end, padded_contig, padded_start, padded_end) in CTOR_PAIR {
        let result = AssemblyRegion::new(
            interval(contig, *start, *end),
            interval(padded_contig, *padded_start, *padded_end),
            false,
        );
        assert_eq!(
            describe_result(&result),
            row(&text, &format!("ctor\t{label}\t")),
            "{label}"
        );
    }
}

#[test]
fn every_trim_matches_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let records = corpus::records(&text);

    for (label, source) in [
        ("wide", "wide"),
        ("javadoc", "javadoc"),
        ("narrow", "narrow"),
    ] {
        let region = build(source, &header, &records);
        assert_eq!(
            describe(&region),
            row(&text, &format!("region\t{label}\t")),
            "{label}: the region before trimming"
        );
        check_reads(&text, label, &region);
    }

    for (label, contig, start, end, padding) in TRIM {
        let region = build("wide", &header, &records);
        let result = region.trim_with_padding(&interval(contig, *start, *end), *padding, &header);
        assert_eq!(
            describe_result(&result),
            row(&text, &format!("trim\t{label}\t")),
            "{label}"
        );
        if let Ok(trimmed) = &result {
            check_reads(&text, label, trimmed);
        }
    }

    for (label, contig, start, end, padded_contig, padded_start, padded_end) in TRIM_PAIR {
        let region = build("wide", &header, &records);
        let result = region.trim(
            &interval(contig, *start, *end),
            &interval(padded_contig, *padded_start, *padded_end),
            &header,
        );
        assert_eq!(
            describe_result(&result),
            row(&text, &format!("trim\t{label}\t")),
            "{label}"
        );
        if let Ok(trimmed) = &result {
            check_reads(&text, label, trimmed);
        }
    }

    // The javadoc's own worked example, run through the code.
    let region = build("javadoc", &header, &records);
    let trimmed = region
        .trim_with_padding(&interval("chr1", 150, 225), 25, &header)
        .expect("the javadoc example trims");
    assert_eq!(
        describe(&trimmed),
        row(&text, "trim\tjavadoc-trim\t"),
        "javadoc-trim"
    );
    check_reads(&text, "javadoc-trim", &trimmed);

    let region = build("narrow", &header, &records);
    let trimmed = region
        .trim_with_padding(&interval("chr1", 1005, 1008), 0, &header)
        .expect("the narrow region trims");
    assert_eq!(
        describe(&trimmed),
        row(&text, "trim\tnarrow-trim\t"),
        "narrow-trim"
    );
    check_reads(&text, "narrow-trim", &trimmed);
}

#[test]
fn every_addition_matches_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let records = corpus::records(&text);

    let outcome = |region: &mut AssemblyRegion, record: &BamRecord| -> String {
        match region.add(record.clone(), &header) {
            Ok(()) => "ok".to_string(),
            Err(error) => render_error(&error),
        }
    };

    // Declaration order, which the corpus deliberately breaks near its end.
    let mut region = AssemblyRegion::with_padding(interval("chr1", 100, 1900), true, 100, &header)
        .expect("a region");
    for (index, record) in records.iter().enumerate() {
        assert_eq!(
            outcome(&mut region, record),
            row(&text, &format!("add\tadd-corpus-order\t{index}\t")),
            "add-corpus-order, record {index} ({})",
            record.read_name
        );
    }
    assert_eq!(describe(&region), row(&text, "region\tadd-corpus-order\t"));
    check_reads(&text, "add-corpus-order", &region);

    // The same records sorted, which is the order that is accepted.
    let mut sorted = records.clone();
    sorted.sort_by(read_utils::compare_read_coordinate);
    let mut region = AssemblyRegion::with_padding(interval("chr1", 100, 1900), true, 100, &header)
        .expect("a region");
    for (index, record) in sorted.iter().enumerate() {
        assert_eq!(
            outcome(&mut region, record),
            row(&text, &format!("add\tadd-sorted\t{index}\t")),
            "add-sorted, record {index} ({})",
            record.read_name
        );
    }
    assert_eq!(describe(&region), row(&text, "region\tadd-sorted\t"));
    check_reads(&text, "add-sorted", &region);

    let named = |name: &str| -> BamRecord {
        records
            .iter()
            .find(|record| record.read_name == name)
            .unwrap_or_else(|| panic!("no corpus record named {name}"))
            .clone()
    };
    for (label, name) in [
        ("add-unmapped", "flag_unmapped"),
        ("add-zero-start", "zero_start"),
        ("add-no-reference", "no_reference"),
        ("add-other-contig", "chr2_mapped"),
    ] {
        let mut region =
            AssemblyRegion::with_padding(interval("chr1", 100, 1900), true, 100, &header)
                .expect("a region");
        assert_eq!(
            outcome(&mut region, &named(name)),
            row(&text, &format!("add\t{label}\t0\t")),
            "{label}"
        );
        assert_eq!(describe(&region), row(&text, &format!("region\t{label}\t")));
    }
}

/// The second read list does not survive a trim.
#[test]
fn trimming_drops_the_hard_clipped_pileup_reads() {
    let text = golden();
    let header = corpus::header(&text);
    let records = corpus::records(&text);

    let mut region = build("wide", &header, &records);
    let reads: Vec<BamRecord> = region.reads().to_vec();
    region
        .add_hard_clipped_pileup_reads(reads, &header)
        .expect("the region's own reads are in order");
    assert_eq!(
        format!(
            "{}\t{}",
            region.size(),
            region.hard_clipped_pileup_reads().len()
        ),
        row(&text, "pileup\tpileup-drop\tbefore\t"),
        "before the trim"
    );

    let trimmed = region
        .trim_with_padding(&interval("chr1", 500, 1000), 50, &header)
        .expect("the region trims");
    assert_eq!(
        format!(
            "{}\t{}",
            trimmed.size(),
            trimmed.hard_clipped_pileup_reads().len()
        ),
        row(&text, "pileup\tpileup-drop\tafter\t"),
        "after the trim"
    );
    assert_eq!(trimmed.hard_clipped_pileup_reads().len(), 0);
}

/// The row that settles the javadoc: the class documents one answer and produces another.
#[test]
fn the_trim_javadoc_describes_a_region_the_code_does_not_produce() {
    let text = golden();
    let produced = row(&text, "trim\tjavadoc-trim\t");
    // The javadoc says "a region from 150-200 with 25 bp of padding", which would be
    // chr1:150-200|chr1:125-225.
    assert!(
        !produced.starts_with("chr1:150-200|chr1:125-225|"),
        "the javadoc's answer and the code's now agree, which they did not: {produced}"
    );
}
