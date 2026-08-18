//! `FilterMutectCalls` end to end against GATK 4.6.2.0: the FILTER column and `AS_FilterStatus` of
//! every record, in every run the golden holds.
//!
//! Golden from `tools/readfilter-conformance/FilterMutectCallsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the tool makes four passes**, so the filters applied to the first record come from a model
//!    that has already seen the last: the same record filtered alone and in company answers
//!    differently, and the golden runs that pair;
//!  * **six filters are built only when the run is not mitochondrial**, which is why the same input
//!    filtered twice in the two modes gives two different FILTER columns;
//!  * **the threshold is learned in a pass of its own**, after the parameters have stopped moving;
//!  * **and a threshold strategy other than the default moves every column at once.**

use gatk_corpus as corpus;
use gatk_engine::accumulate_data::AccumulationAllele;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::filtering_engine::{EngineArguments, Record};
use gatk_engine::mutect_filter_list::FilterArguments;
use gatk_engine::somatic_clustering_model::AlternateAllele;
use gatk_engine::threshold_calculator::Strategy;
use gatk_tools::filter_mutect_calls::{
    run as filter_mutect_calls, MissingStatsTable, ToolArguments,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_mutect_calls.txt.gz"),
    )
}

/// One record of the dump's input VCF, in the fields the filters read.
struct Row {
    start: i32,
    tumour: (i32, i32),
    normal: (i32, i32),
    tumour_fraction: f64,
    tumour_log_10_odds: f64,
    normal_artifact_log_10_odds: f64,
    normal_log_10_odds: f64,
    population_af: f64,
    median_base_quality: [i32; 2],
    median_read_position: i32,
    event_count_in_region: i32,
}

/// The eight records of `calls.vcf`, as the dump writes them.
const ROWS: [Row; 8] = [
    Row {
        start: 100,
        tumour: (80, 20),
        normal: (99, 1),
        tumour_fraction: 0.200,
        tumour_log_10_odds: 30.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 30],
        median_read_position: 25,
        event_count_in_region: 1,
    },
    Row {
        start: 200,
        tumour: (78, 22),
        normal: (99, 1),
        tumour_fraction: 0.220,
        tumour_log_10_odds: 40.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 30],
        median_read_position: 25,
        event_count_in_region: 1,
    },
    Row {
        start: 300,
        tumour: (79, 21),
        normal: (99, 1),
        tumour_fraction: 0.210,
        tumour_log_10_odds: 35.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 30],
        median_read_position: 25,
        event_count_in_region: 1,
    },
    Row {
        start: 400,
        tumour: (97, 3),
        normal: (99, 1),
        tumour_fraction: 0.030,
        tumour_log_10_odds: 3.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 30],
        median_read_position: 25,
        event_count_in_region: 1,
    },
    Row {
        start: 500,
        tumour: (80, 20),
        normal: (99, 1),
        tumour_fraction: 0.200,
        tumour_log_10_odds: 30.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 5],
        median_read_position: 25,
        event_count_in_region: 1,
    },
    Row {
        start: 600,
        tumour: (80, 20),
        normal: (99, 1),
        tumour_fraction: 0.200,
        tumour_log_10_odds: 30.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 30],
        median_read_position: 0,
        event_count_in_region: 1,
    },
    Row {
        start: 700,
        tumour: (80, 20),
        normal: (99, 1),
        tumour_fraction: 0.200,
        tumour_log_10_odds: 30.0,
        normal_artifact_log_10_odds: 2.0,
        normal_log_10_odds: 5.0,
        population_af: 6.0,
        median_base_quality: [30, 30],
        median_read_position: 25,
        event_count_in_region: 9,
    },
    Row {
        start: 800,
        tumour: (50, 50),
        normal: (60, 40),
        tumour_fraction: 0.500,
        tumour_log_10_odds: 30.0,
        normal_artifact_log_10_odds: -3.0,
        normal_log_10_odds: -5.0,
        population_af: 0.1,
        median_base_quality: [30, 30],
        median_read_position: 25,
        event_count_in_region: 1,
    },
];

fn record(row: &Row) -> Record {
    Record {
        start: row.start,
        reference_length: 1,
        alternates: vec![AccumulationAllele {
            allele: AlternateAllele {
                length: 1,
                symbolic: false,
            },
            non_ref: false,
        }],
        genotypes: vec![
            GenotypeData {
                tumor: true,
                allele_depths: vec![row.tumour.0, row.tumour.1],
                values: Vec::new(),
            },
            GenotypeData {
                tumor: false,
                allele_depths: vec![row.normal.0, row.normal.1],
                values: Vec::new(),
            },
        ],
        // `AF` is written to three decimals and read back, so it is the rounded value.
        allele_fractions: vec![vec![row.tumour_fraction], vec![0.010]],
        phasing: vec![(None, None), (None, None)],
        tumor_log_10_odds: Some(vec![row.tumour_log_10_odds]),
        normal_artifact_log_10_odds: Some(vec![row.normal_artifact_log_10_odds]),
        normal_log_10_odds: Some(vec![row.normal_log_10_odds]),
        population_af: Some(vec![row.population_af]),
        median_base_quality: Some(row.median_base_quality.to_vec()),
        median_mapping_quality: Some(vec![60, 60]),
        median_fragment_length: Some(vec![300, 300]),
        median_read_position: Some(vec![row.median_read_position]),
        // The input carries no AS_SB_TABLE, no AS_UNIQ_ALT_READ_COUNT, no NCount, no ECNTH, no
        // RPA/RU and no PON, so those filters answer an empty list or a zero.
        unique_alt_read_count: None,
        strand_bias_table: None,
        n_count: None,
        event_count_in_region: Some(row.event_count_in_region),
        event_count_in_haplotype: None,
        repeats_per_allele: None,
        repeat_unit: None,
        in_panel_of_normals: false,
        indel_lengths: None,
    }
}

fn arguments(run: &str) -> ToolArguments {
    let mut arguments = ToolArguments {
        callable_sites: Some(1000000.0),
        ..ToolArguments::default()
    };
    match run {
        "default" | "single-record" => {}
        "mitochondria" => {
            arguments.engine = EngineArguments {
                list: FilterArguments {
                    mitochondria: true,
                    ..FilterArguments::default()
                },
                ..EngineArguments::default()
            }
        }
        // A stats table saying nothing was callable switches the empirical priors off.
        "no-callable-sites" => arguments.callable_sites = None,
        "constant-threshold" => {
            arguments.strategy = Strategy::Constant;
            arguments.initial_posterior_threshold = 0.5;
        }
        other => panic!("no run named {other}"),
    }
    arguments
}

fn records(run: &str) -> Vec<Record> {
    if run == "single-record" {
        vec![record(&ROWS[0])]
    } else {
        ROWS.iter().map(record).collect()
    }
}

/// The FILTER column as htsjdk writes it: `PASS` when empty, the names sorted otherwise.
fn filter_column(names: &[String]) -> String {
    if names.is_empty() {
        return "PASS".to_string();
    }
    let mut sorted = names.to_vec();
    sorted.sort();
    sorted.join(";")
}

#[test]
fn every_filter_column_matches_the_golden() {
    let text = golden();
    #[allow(clippy::type_complexity)]
    let mut expected: std::collections::BTreeMap<String, Vec<(String, String, Vec<String>)>> =
        std::collections::BTreeMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.splitn(3, '\t');
        if fields.next() != Some("vcfline") {
            continue;
        }
        let run = fields.next().expect("a run").to_string();
        let payload = fields.next().expect("a line");
        let columns: Vec<&str> = payload.split("\\t").collect();
        let filter = columns[6].to_string();
        let status = columns[7]
            .split(';')
            .find_map(|entry| entry.strip_prefix("AS_FilterStatus="))
            .expect("an AS_FilterStatus")
            .to_string();
        // The phred-scaled posteriors the tool adds, which are the INFO keys the input did not
        // carry. Everything else in the column came in with the record.
        let mut annotations: Vec<String> = columns[7]
            .split(';')
            .filter(|entry| {
                ["SEQQ=", "CONTQ=", "GERMQ=", "STRANDQ=", "STRQ="]
                    .iter()
                    .any(|key| entry.starts_with(key))
            })
            .map(str::to_string)
            .collect();
        annotations.sort();
        expected
            .entry(run)
            .or_default()
            .push((filter, status, annotations));
    }
    assert_eq!(expected.len(), 5, "the runs that produced records");

    let mut complaints: Vec<String> = Vec::new();
    for (run, rows) in &expected {
        let output = filter_mutect_calls(&records(run), &arguments(run));
        assert_eq!(output.records.len(), rows.len(), "{run}: record count");
        for (index, (filter, status, annotations)) in rows.iter().enumerate() {
            let ours = &output.records[index];
            let ours_filter = filter_column(&ours.filters);
            if ours_filter != *filter || ours.as_filter_status != *status {
                complaints.push(format!(
                    "{run} record {index}: ours {ours_filter}|{}, reference {filter}|{status}",
                    ours.as_filter_status
                ));
            }
            let mut ours_annotations: Vec<String> = ours
                .annotations
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect();
            ours_annotations.sort();
            if ours_annotations != *annotations {
                complaints.push(format!(
                    "{run} record {index}: ours {ours_annotations:?}, reference {annotations:?}"
                ));
            }
        }
    }
    assert!(complaints.is_empty(), "{}", complaints.join("\n"));
}

/// The filtering-stats file each run writes, line by line.
#[test]
fn every_filtering_stats_line_matches_the_golden() {
    let text = golden();
    let mut expected: std::collections::BTreeMap<String, Vec<String>> =
        std::collections::BTreeMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.splitn(3, '\t');
        if fields.next() != Some("filtering") {
            continue;
        }
        let run = fields.next().expect("a run").to_string();
        expected
            .entry(run)
            .or_default()
            .push(fields.next().expect("a line").to_string());
    }

    assert_eq!(
        expected.len(),
        5,
        "the runs that wrote a filtering-stats file"
    );
    assert_eq!(
        expected.values().map(Vec::len).sum::<usize>(),
        151,
        "the golden's filtering rows"
    );

    let mut complaints: Vec<String> = Vec::new();
    for (run, lines) in &expected {
        let output = filter_mutect_calls(&records(run), &arguments(run));
        let ours: Vec<String> = output
            .filtering_stats
            .lines()
            .map(|line| line.replace('\t', "\\t"))
            .collect();
        if ours != *lines {
            for (index, expected_line) in lines.iter().enumerate() {
                match ours.get(index) {
                    Some(ours_line) if ours_line == expected_line => {}
                    Some(ours_line) => complaints.push(format!(
                        "{run} line {index}: ours {ours_line}, reference {expected_line}"
                    )),
                    None => complaints.push(format!(
                        "{run} line {index}: missing, reference {expected_line}"
                    )),
                }
            }
            if ours.len() > lines.len() {
                complaints.push(format!("{run}: {} extra lines", ours.len() - lines.len()));
            }
        }
    }
    assert!(complaints.is_empty(), "{}", complaints.join("\n"));
}

/// The rows that are not outputs: the inputs the dump wrote, and the one refusal.
///
/// Comparing them closes the loop. The eight records above are hard-coded from the dump's source;
/// these assertions are what says the hard-coding is the same input the reference was given.
#[test]
fn the_inputs_and_the_refusal_match_the_golden() {
    let text = golden();
    let mut seen = 0;
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.splitn(3, '\t');
        let kind = fields.next().expect("a kind");
        let label = fields.next().expect("a label");
        let payload = fields.next().expect("a payload");
        match kind {
            "stats" => {
                // The `.stats` table Mutect2 writes: two columns and one row.
                let callable = match label {
                    "calls" | "single" => "1000000.0",
                    "no-callable" => "0.0",
                    other => panic!("no stats table named {other}"),
                };
                assert_eq!(
                    payload,
                    format!("statistic\\tvalue\\ncallable\\t{callable}\\n"),
                    "the {label} stats table"
                );
                seen += 1;
            }
            "input" => {
                // Every record the test encodes must appear in the input, with the depths and the
                // log odds it was encoded with.
                let rows: &[Row] = if label == "single" { &ROWS[..1] } else { &ROWS };
                for row in rows {
                    let start = row.start;
                    assert!(
                        payload.contains(&format!("chr1\\t{start}\\t.\\tA\\tC")),
                        "the {label} input carries a record at {start}"
                    );
                    assert!(
                        payload.contains(&format!("TLOD={:.1}", row.tumour_log_10_odds)),
                        "the {label} input carries TLOD {} ",
                        row.tumour_log_10_odds
                    );
                }
                seen += 1;
            }
            "error" => {
                assert_eq!(label, "missing-stats");
                let refusal = MissingStatsTable {
                    path: "filter-mutect-calls-dump/no-such.stats".to_string(),
                };
                assert_eq!(
                    payload,
                    format!("{}:{}", refusal.class(), refusal.message()),
                    "the refusal"
                );
                seen += 1;
            }
            _ => {}
        }
    }
    assert_eq!(seen, 6, "two inputs, three stats tables and one refusal");
}

/// Nothing in the golden is unaccounted for.
///
/// The three tests above take the `vcfline`, `filtering`, `input`, `stats` and `error` rows, and the
/// `header` rows are taken by `filter_mutect_calls_header`. A new row kind fails here rather than
/// being skipped by every test at once.
#[test]
fn every_row_of_the_golden_is_accounted_for() {
    let text = golden();
    let mut counts: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let kind = line.split('\t').next().expect("a kind");
        *counts
            .entry(match kind {
                "header" | "vcfline" | "filtering" | "input" | "stats" | "error" => kind,
                other => panic!("no test takes the {other} rows"),
            })
            .or_default() += 1;
    }
    assert_eq!(counts["header"], 130);
    assert_eq!(counts["vcfline"], 33);
    assert_eq!(counts["filtering"], 151);
    assert_eq!(counts["input"], 2);
    assert_eq!(counts["stats"], 3);
    assert_eq!(counts["error"], 1);
    assert_eq!(
        counts.values().sum::<usize>(),
        320,
        "the golden's row count"
    );
}
