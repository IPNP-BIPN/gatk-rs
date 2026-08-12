//! Conformance for `BaseRecalibrator` against GATK 4.6.2.0, compared **character for character**.
//!
//! Golden from `tools/readfilter-conformance/BaseRecalibratorDump.java`. The output is a
//! `GATKReport` of five tables, so the whole text travels and every run is compared line by line.
//! The input BAM, its index, the reference and both known-sites files travel too.
//!
//! # What this suite is for
//!
//! The tool that closes the BQSR cycle: it writes the table `ApplyBQSR` reads.
//!
//!  * **its default read filters are seven**, six BQSR-specific plus `WellformedReadFilter`, and the
//!    fixture carries one read per filter so the list is visible in what the table does not hold;
//!  * **a BED and a VCF naming the same sites produce the same table**;
//!  * **the two additional covariate tables share one report table**, interleaved by the sort;
//!  * **the rows are ordered by their values**, because the sort is `SORT_BY_COLUMN`.

use gatk_corpus as corpus;
use gatk_engine::base_recalibration_engine::EngineArguments;
use gatk_engine::covariates::RecalibrationArguments;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::base_recalibrator::{self, QUANTIZING_LEVELS};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/base_recalibrator.txt.gz"),
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

fn value<'a>(text: &'a str, kind: &str, label: &str) -> Option<&'a str> {
    pairs(text, kind)
        .into_iter()
        .find(|(name, _)| *name == label)
        .map(|(_, value)| value)
}

/// The harness's escaping, for the text rows.
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

/// The seven default filters: six BQSR-specific, plus `WellformedReadFilter`.
fn default_filter(header: &SamHeader) -> impl Fn(&BamRecord) -> bool + '_ {
    move |read: &BamRecord| {
        gatk_readfilter::mapping_quality_not_zero(read)
            && gatk_readfilter::mapping_quality_available(read)
            && gatk_readfilter::mapped(read)
            && gatk_readfilter::not_secondary_alignment(read)
            && gatk_readfilter::not_duplicate(read)
            && gatk_readfilter::passes_vendor_quality_check(read)
            && with_header::wellformed(read, header)
    }
}

/// The arguments each labelled run was given.
fn arguments(label: &str) -> (EngineArguments, i32) {
    let base = EngineArguments::default();
    match label {
        "bed-sites" | "vcf-sites" => (base, QUANTIZING_LEVELS),
        "indel-tables" => (
            EngineArguments {
                compute_indel_bqsr_tables: true,
                ..base
            },
            QUANTIZING_LEVELS,
        ),
        "baq-enabled" => (
            EngineArguments {
                enable_baq: true,
                ..base
            },
            QUANTIZING_LEVELS,
        ),
        "quantizing-4" => (base, 4),
        "preserve-20" => (
            EngineArguments {
                preserve_qscores_less_than: 20,
                ..base
            },
            QUANTIZING_LEVELS,
        ),
        "context-3" => (
            EngineArguments {
                covariates: RecalibrationArguments {
                    mismatches_context_size: 3,
                    ..RecalibrationArguments::default()
                },
                ..base
            },
            QUANTIZING_LEVELS,
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The known sites both files name: `chr1:10-12`, which the BED writes half-open as `9 12`.
fn known_sites() -> Vec<SimpleInterval> {
    vec![SimpleInterval {
        contig: "chr1".to_string(),
        start: 10,
        end: 12,
    }]
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
        std::fs::write(
            dir.join(format!("{label}.bai")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture index");
    }
}

#[test]
fn every_recalibration_table_is_the_reference_line_for_line() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-baserecal-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let reference = text
        .lines()
        .find_map(|line| line.strip_prefix("reference\t"))
        .expect("the golden carries the reference");

    let tables = pairs(&text, "table");
    assert_eq!(tables.len(), 7, "seven runs");

    let mut compared = 0;
    for (label, expected) in &tables {
        let (engine_arguments, levels) = arguments(label);
        let source = ReadsDataSource::open(&dir.join("input.bam"), &dir.join("input.bai"))
            .expect("the fixture opens");
        let header = source.header().clone();
        let filter = default_filter(&header);

        let ours = base_recalibrator::base_recalibrator(
            &source,
            reference.as_bytes(),
            &known_sites(),
            &engine_arguments,
            levels,
            &filter,
        )
        .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let theirs = unescape(expected);
        let ours_lines: Vec<&str> = ours.lines().collect();
        let their_lines: Vec<&str> = theirs.lines().collect();
        assert_eq!(
            ours_lines.len(),
            their_lines.len(),
            "{label}: line count differs"
        );
        for (n, (ours, theirs)) in ours_lines.iter().zip(&their_lines).enumerate() {
            assert_eq!(ours, theirs, "{label}: line {n}");
        }
        compared += 1;
    }
    println!("base-recalibrator: {compared} tables compared line for line");
}

/// A BED and a VCF naming the same sites must produce the same table.
#[test]
fn the_two_known_sites_formats_agree() {
    let text = golden();
    let bed = value(&text, "table", "bed-sites").expect("the bed run");
    let vcf = value(&text, "table", "vcf-sites").expect("the vcf run");
    assert_eq!(bed, vcf, "the reference's two runs already agree");
}

/// The two additional covariates share one table, and its rows come out interleaved by the sort.
#[test]
fn the_additional_covariates_share_one_interleaved_table() {
    let text = golden();
    let table = unescape(value(&text, "table", "bed-sites").expect("the bed run"));
    assert!(table.contains("#:GATKTable:RecalTable2:"));
    assert!(!table.contains("#:GATKTable:RecalTable3:"));

    // The two covariate names alternate rather than appearing in two runs, because the sort orders
    // by the row's values and the quality score comes before the covariate value.
    let names: Vec<&str> = table
        .lines()
        .filter(|line| line.contains("Context") || line.contains("Cycle"))
        .filter(|line| line.starts_with("unit-rg1"))
        .map(|line| {
            if line.contains("Context") {
                "Context"
            } else {
                "Cycle"
            }
        })
        .collect();
    assert!(names.len() > 4, "the covariate table has rows");
    assert!(
        names.windows(2).any(|pair| pair[0] != pair[1]),
        "the two covariates are interleaved, not written one after the other"
    );
}
