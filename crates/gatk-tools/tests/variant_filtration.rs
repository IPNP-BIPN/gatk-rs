//! Conformance for `VariantFiltration` against GATK 4.6.2.0, compared as the FILTER column of every
//! record of every run, and as the FT of every sample.
//!
//! Golden from `tools/readfilter-conformance/VariantFiltrationDump.java`.
//!
//! # What this suite is for
//!
//!  * **the writer sorts what the tool ordered**, so the file never shows the insertion order;
//!  * **an empty set is `PASS`, or `.` under `--invalidate-previous-filters`**;
//!  * **the cluster test ignores indels** and a narrow window filters nothing;
//!  * **FT is per record**, present only where a sample of that record was filtered;
//!  * **and the mask means opposite things** with one flag turned over.

use gatk_corpus as corpus;
use gatk_tools::variant_filtration::{
    filter_records, rendered_filters, rendered_genotype_filter, Arguments, GenotypeFields,
    MatchExp, Record,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/variant_filtration.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, scanning once so a real backslash is never read as a tab.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
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

/// One data line of a vcf as a record, with its INFO and its genotypes.
fn parse_record(line: &str) -> Record {
    let field: Vec<&str> = line.split('\t').collect();
    let start: i32 = field[1].parse().expect("a position");
    let reference = field[3];
    let is_snp = reference.len() == 1 && field[4].split(',').all(|alternate| alternate.len() == 1);
    let info = field[7]
        .split(';')
        .filter(|entry| !entry.is_empty() && *entry != ".")
        .filter_map(|entry| entry.split_once('='))
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect();
    let format: Vec<&str> = field[8].split(':').collect();
    let genotypes = field[9..]
        .iter()
        .map(|sample| {
            let values: Vec<&str> = sample.split(':').collect();
            let mut fields = std::collections::HashMap::new();
            let mut filters = Vec::new();
            for (key, value) in format.iter().zip(&values) {
                if *key == "FT" {
                    if *value != "PASS" && *value != "." {
                        filters = value.split(';').map(|name| name.to_string()).collect();
                    }
                } else {
                    fields.insert(key.to_string(), value.to_string());
                }
            }
            GenotypeFields { fields, filters }
        })
        .collect();
    Record {
        contig: field[0].to_string(),
        start,
        stop: start + reference.len() as i32 - 1,
        is_snp,
        filters: match field[6] {
            "." => None,
            "PASS" => Some(Vec::new()),
            names => Some(names.split(';').map(|name| name.to_string()).collect()),
        },
        info,
        genotypes,
    }
}

fn input(text: &str, label: &str) -> Vec<Record> {
    let whole = rows(text, "input")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no input {label}"))[1]
        .to_string();
    unescape(&whole)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(parse_record)
        .collect()
}

/// The output lines of one run, as `POS FILTER sample...`.
fn written(text: &str, run: &str) -> Vec<Vec<String>> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(|field| field.to_string()).collect())
        .collect()
}

/// The arguments and expressions of each run, which the golden does not carry.
fn setup(run: &str) -> (Arguments, Vec<MatchExp>, Vec<MatchExp>, &'static str) {
    let site = |name: &str, text: &str| vec![MatchExp::new(name, text).expect("parses")];
    let base = Arguments::default;
    match run {
        "site-filter" => (base(), site("LowQD", "QD < 2.0"), vec![], "records"),
        "site-filter-inverted" => (
            Arguments {
                invert_filter_expression: true,
                ..base()
            },
            site("LowQD", "QD < 2.0"),
            vec![],
            "records",
        ),
        "two-filters" => (
            base(),
            vec![
                MatchExp::new("LowQD", "QD < 2.0").expect("parses"),
                MatchExp::new("LowDP", "DP < 10").expect("parses"),
            ],
            vec![],
            "records",
        ),
        "missing-values" => (
            Arguments {
                missing_values_evaluate_as_failing: true,
                ..base()
            },
            site("LowQD", "QD < 2.0"),
            vec![],
            "records",
        ),
        "genotype-filter" => (base(), vec![], site("LowGQ", "GQ < 30"), "records"),
        "genotype-filter-nocall" => (
            Arguments {
                set_filtered_genotypes_to_no_call: true,
                ..base()
            },
            vec![],
            site("LowGQ", "GQ < 30"),
            "records",
        ),
        "cluster" => (
            Arguments {
                cluster_size: 3,
                cluster_window: 20,
                ..base()
            },
            vec![],
            vec![],
            "clustered",
        ),
        "cluster-narrow" => (
            Arguments {
                cluster_size: 3,
                cluster_window: 5,
                ..base()
            },
            vec![],
            vec![],
            "clustered",
        ),
        "cluster-disabled" => (
            Arguments {
                cluster_size: 3,
                cluster_window: 0,
                ..base()
            },
            vec![],
            vec![],
            "clustered",
        ),
        "mask" => (
            Arguments {
                mask_name: "InMask".to_string(),
                ..base()
            },
            vec![],
            vec![],
            "records",
        ),
        "mask-inverted" => (
            Arguments {
                mask_name: "NotInMask".to_string(),
                filter_records_not_in_mask: true,
                ..base()
            },
            vec![],
            vec![],
            "records",
        ),
        "no-filters" => (base(), site("Never", "QD > 1000.0"), vec![], "records"),
        "no-filters-invalidated" => (
            Arguments {
                invalidate_previous_filters: true,
                ..base()
            },
            site("Never", "QD > 1000.0"),
            vec![],
            "records",
        ),
        other => panic!("no setup for {other}"),
    }
}

fn runs(text: &str) -> Vec<String> {
    let mut seen: Vec<String> = Vec::new();
    for row in rows(text, "vcfline") {
        if !seen.iter().any(|name| name == row[0]) {
            seen.push(row[0].to_string());
        }
    }
    seen
}

#[test]
fn every_filter_column_is_the_reference_s() {
    let text = golden();
    let all = runs(&text);
    assert!(all.len() >= 13, "every run is in the golden: {}", all.len());

    for run in &all {
        let (arguments, site, genotype, label) = setup(run);
        let records = input(&text, label);
        let mask = if run.starts_with("mask") {
            input(&text, "mask")
                .iter()
                .map(|record| (record.contig.clone(), record.start, record.stop))
                .collect()
        } else {
            Vec::new()
        };

        let ours = filter_records(&records, &site, &genotype, &mask, &arguments);
        let expected = written(&text, run);
        assert_eq!(ours.len(), expected.len(), "count/{run}");
        for (index, filtered) in ours.iter().enumerate() {
            assert_eq!(
                rendered_filters(filtered),
                expected[index][6],
                "filter/{run}/{}",
                expected[index][1]
            );
        }
    }
}

/// The FT column exists only where a sample of that record was filtered.
#[test]
fn the_ft_column_is_per_record() {
    let text = golden();
    let (arguments, site, genotype, label) = setup("genotype-filter");
    let records = input(&text, label);
    let ours = filter_records(&records, &site, &genotype, &[], &arguments);
    let expected = written(&text, "genotype-filter");

    for (index, filtered) in ours.iter().enumerate() {
        let format: Vec<&str> = expected[index][8].split(':').collect();
        let has_ft = format.contains(&"FT");
        assert_eq!(
            rendered_genotype_filter(filtered, 0).is_some(),
            has_ft,
            "FT presence at {}",
            expected[index][1]
        );
        if has_ft {
            let position = format.iter().position(|key| *key == "FT").expect("FT");
            for sample in 0..filtered.genotype_filters.len() {
                let want = expected[index][9 + sample]
                    .split(':')
                    .nth(position)
                    .unwrap();
                assert_eq!(
                    rendered_genotype_filter(filtered, sample).unwrap(),
                    want,
                    "FT/{}/{sample}",
                    expected[index][1]
                );
            }
        }
    }
}

/// The same emptiness, two columns.
#[test]
fn an_empty_set_is_pass_or_a_dot() {
    let text = golden();
    for run in ["no-filters", "no-filters-invalidated"] {
        let (arguments, site, genotype, label) = setup(run);
        let ours = filter_records(&input(&text, label), &site, &genotype, &[], &arguments);
        let expected = written(&text, run);
        for (index, filtered) in ours.iter().enumerate() {
            assert_eq!(rendered_filters(filtered), expected[index][6], "{run}");
        }
    }
    // And the pre-existing filter is wiped by the invalidating run, kept by the other.
    assert_eq!(written(&text, "no-filters")[1][6], "OldFilter");
    assert_eq!(written(&text, "no-filters-invalidated")[1][6], ".");
}

/// One flag turned over, and the mask filters the complement.
#[test]
fn the_mask_means_the_opposite_with_one_flag() {
    let text = golden();
    let inside: Vec<String> = written(&text, "mask")
        .iter()
        .map(|line| line[6].clone())
        .collect();
    let outside: Vec<String> = written(&text, "mask-inverted")
        .iter()
        .map(|line| line[6].clone())
        .collect();
    assert_eq!(inside[1], "InMask;OldFilter");
    assert!(
        outside[1] == "OldFilter",
        "the masked record is the one NOT filtered"
    );
    assert!(outside.iter().filter(|f| f.contains("NotInMask")).count() == 4);
}
