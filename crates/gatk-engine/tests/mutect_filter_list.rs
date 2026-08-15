//! Conformance for every Mutect filter's identity and the header lines around them, against
//! GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/MutectFilterListDump.java`.
//!
//! # What this suite is for
//!
//!  * **the filters disagree on both axes the engine sorts them by**, nine per allele and nine per
//!    site, over three error types;
//!  * **the stats file names only the filters that fired**, so an engine that has filtered nothing
//!    writes metadata and a header row and nothing else, in every mode;
//!  * **the header's list is not the engine's list**: `PASS` and `FAIL` belong to no filter, and
//!    `orientation` and `possible_numt` belong to filters an argument or a mode has to build;
//!  * **and the tool writes over Mutect2's header line under the same key.**
//!
//! Every row is compared and every row is bit-identical. Nothing here is arithmetic.

use gatk_corpus as corpus;
use gatk_engine::filtering_stats::COLUMNS;
use gatk_engine::mutect_filter_list::{
    filter_line, FILTERED_FILTERING_STATUS, FILTERING_STATUS_VCF_KEY, FILTERS, INFO_LINES,
    MUTECT2_FILTERING_STATUS, MUTECT_AS_FILTER_NAMES, MUTECT_FILTER_NAMES,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mutect_filter_list.txt.gz"),
    )
}

/// The `stats <label>-empty` and `filters`/`count` rows one engine produces before it has filtered
/// anything: the writer's header row, and no filter rows at all.
fn engine_rows(label: &str) -> Vec<String> {
    vec![
        // Written twice by the dump, before and after it would have filtered: an engine that
        // filtered nothing writes the same file both times.
        format!("stats\t{label}-empty\t{}", COLUMNS.join("\\t")),
        format!("filters\t{label}\t"),
        format!("count\t{label}\t0"),
        format!("stats\t{label}\t{}", COLUMNS.join("\\t")),
    ]
}

fn ours() -> Vec<String> {
    let mut rows = Vec::new();
    for label in ["default", "mitochondria", "tuned"] {
        rows.extend(engine_rows(label));
    }
    for filter in FILTERS {
        rows.push(format!(
            "filter\t{}\t{},{},{}",
            filter.class,
            filter.filter_name,
            filter.error_type.name(),
            filter.arity.name()
        ));
    }
    rows.push(format!(
        "names\tMUTECT_FILTER_NAMES\t{}",
        MUTECT_FILTER_NAMES.join(",")
    ));
    rows.push(format!(
        "names\tMUTECT_AS_FILTER_NAMES\t{}",
        MUTECT_AS_FILTER_NAMES.join(",")
    ));
    for name in MUTECT_FILTER_NAMES {
        rows.push(format!(
            "filterline\t{name}\t{}",
            filter_line(name).unwrap_or_else(|| panic!("no line for {name}"))
        ));
    }
    for (name, line) in INFO_LINES {
        rows.push(format!("infoline\t{name}\t{line}"));
    }
    rows.push(format!(
        "stats\tfiltering-status-mutect2\t{FILTERING_STATUS_VCF_KEY}={MUTECT2_FILTERING_STATUS}"
    ));
    rows.push(format!(
        "stats\tfiltering-status-key\t{FILTERING_STATUS_VCF_KEY}"
    ));
    rows.push(format!(
        "stats\tfiltering-status-line\t{FILTERING_STATUS_VCF_KEY}={FILTERED_FILTERING_STATUS}"
    ));
    rows
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    // The metadata rows of the three empty stats files are the clustering model's, which the
    // `somatic-clustering-model` suite already pins; this suite compares the rest.
    let expected: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| !line.contains("\t#<METADATA>"))
        .collect();
    assert_eq!(
        text.lines().filter(|line| !line.starts_with('#')).count(),
        137,
        "the golden's row count"
    );

    let mine = ours();
    // Compare as sets of rows keyed by their first two fields, since the dump interleaves the
    // engines' rows and this port groups them.
    for row in &mine {
        assert!(expected.contains(&row.as_str()), "not in the golden: {row}");
    }
    assert_eq!(mine.len(), expected.len(), "every row is accounted for");
}
