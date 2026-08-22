//! Conformance for `NuMTFilterTool` against GATK 4.6.2.0, compared as the whole output VCF of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/NuMTFilterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the twenty-one cutoffs**, which are Poisson quantiles and the only arithmetic in the tool;
//!  * **the tool doing nothing at its defaults**, and nothing again at zero copies;
//!  * **an ordinary VCF making it throw**, once anything is filtered and not before;
//!  * **a record with no alternate allele being filtered**;
//!  * **the depth compared being the maximum across samples**, not the sum;
//!  * **and an existing allele filter being unioned** while the SITE placeholder is replaced.

use gatk_corpus as corpus;
use gatk_tools::numt_filter::{
    apply, cutoff_for, decode_as_filters, max_alt_depth_cutoff, Arguments, NuMTError, Record,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/numt_filter.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

/// Every `cutoff\t<coverage>\t<copies>\t<value>` row of the golden.
fn cutoffs(text: &str) -> Vec<(f64, f64, i32)> {
    text.lines()
        .filter_map(|line| line.strip_prefix("cutoff\t"))
        .map(|row| {
            let fields: Vec<&str> = row.split('\t').collect();
            (
                fields[0].parse().expect("a coverage"),
                fields[1].parse().expect("a copy count"),
                fields[2].parse().expect("a cutoff"),
            )
        })
        .collect()
}

/// One VCF record read back into what the filter reads from it.
fn parse(line: &str) -> Record {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut alleles = vec![fields[3].to_string()];
    if fields[4] != "." {
        alleles.extend(fields[4].split(',').map(str::to_string));
    }
    let format: Vec<&str> = fields[8].split(':').collect();
    let ad_index = format.iter().position(|key| *key == "AD");
    let allele_depths = fields[9..]
        .iter()
        .map(|genotype| {
            let values: Vec<&str> = genotype.split(':').collect();
            ad_index.and_then(|index| values.get(index)).map(|depths| {
                depths
                    .split(',')
                    .map(|depth| depth.parse().expect("a depth"))
                    .collect()
            })
        })
        .collect();
    let as_filter_status = fields[7]
        .split(';')
        .find_map(|entry| entry.strip_prefix("AS_FilterStatus="))
        .map(str::to_string);
    Record {
        alleles,
        allele_depths,
        filters: if fields[6] == "." {
            Vec::new()
        } else {
            fields[6].split(';').map(str::to_string).collect()
        },
        as_filter_status,
    }
}

/// The record written back out, in the two columns this tool touches.
fn rendered(original: &str, record: &Record) -> String {
    let fields: Vec<&str> = original.split('\t').collect();
    let filters = if record.filters.is_empty() {
        ".".to_string()
    } else {
        record.filters.join(";")
    };
    let info = match &record.as_filter_status {
        Some(status) => format!("AS_FilterStatus={status}"),
        None => ".".to_string(),
    };
    let mut out: Vec<String> = fields.iter().map(|field| field.to_string()).collect();
    out[6] = filters;
    out[7] = info;
    out.join("\t")
}

/// The data lines of a VCF, header dropped.
fn records(vcf: &str) -> Vec<&str> {
    vcf.lines().filter(|line| !line.starts_with('#')).collect()
}

fn run(text: &str, label: &str, arguments: &Arguments) -> Result<Vec<String>, NuMTError> {
    let cutoff = cutoff_for(arguments).expect("a cutoff");
    records(&value(text, "input", label))
        .iter()
        .map(|line| apply(&parse(line), cutoff).map(|record| rendered(line, &record)))
        .collect()
}

fn expected(text: &str, label: &str) -> Vec<String> {
    records(&value(text, "filtered", label))
        .iter()
        .map(|line| line.to_string())
        .collect()
}

#[test]
fn every_cutoff_matches_the_golden() {
    let text = golden();
    let rows = cutoffs(&text);
    assert_eq!(rows.len(), 21, "the golden's cutoffs");
    for (coverage, copies, cutoff) in rows {
        assert_eq!(
            max_alt_depth_cutoff(copies, coverage).expect("a cutoff"),
            cutoff,
            "coverage {coverage}, copies {copies}"
        );
    }
}

#[test]
fn every_filtered_vcf_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, arguments) in [
        ("defaults", Arguments::default()),
        (
            "coverage-30",
            Arguments {
                median_autosomal_coverage: 30.0,
                ..Arguments::default()
            },
        ),
        (
            "copies-0",
            Arguments {
                median_autosomal_coverage: 30.0,
                max_numt_autosomal_copies: 0.0,
            },
        ),
        (
            "coverage-1000",
            Arguments {
                median_autosomal_coverage: 1000.0,
                ..Arguments::default()
            },
        ),
        ("no-status-defaults", Arguments::default()),
        (
            "no-alt",
            Arguments {
                median_autosomal_coverage: 30.0,
                ..Arguments::default()
            },
        ),
    ] {
        assert_eq!(
            run(&text, label, &arguments).expect("a run that is not refused"),
            expected(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "the golden's outputs");
}

/// The coverage defaults to zero, so the cutoff is zero and the output is the input.
#[test]
fn the_tool_does_nothing_at_its_defaults() {
    let text = golden();
    assert_eq!(cutoff_for(&Arguments::default()).expect("a cutoff"), 0);
    assert_eq!(
        run(&text, "defaults", &Arguments::default()).expect("a run"),
        records(&value(&text, "input", "defaults"))
            .iter()
            .map(|line| line.to_string())
            .collect::<Vec<String>>()
    );
    // And zero copies keeps it at zero however deep the coverage.
    assert_eq!(
        cutoff_for(&Arguments {
            median_autosomal_coverage: 30.0,
            max_numt_autosomal_copies: 0.0
        })
        .expect("a cutoff"),
        0
    );
}

/// An absent AS_FilterStatus decodes to an empty list, so the length check fails as soon as
/// anything is filtered, and not before.
#[test]
fn an_ordinary_vcf_throws_once_anything_is_filtered() {
    let text = golden();
    assert!(decode_as_filters(None).is_empty());
    assert_eq!(decode_as_filters(Some("SITE|SITE")).len(), 2);

    let error = run(
        &text,
        "no-status-filtered",
        &Arguments {
            median_autosomal_coverage: 30.0,
            ..Arguments::default()
        },
    )
    .expect_err("a refused run");
    assert_eq!(error, NuMTError::ListsNotTheSameSize);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        unescape(
            text.lines()
                .find_map(|line| line.strip_prefix("error\tno-status-filtered\t"))
                .expect("the golden's refusal")
        )
    );
}

/// A record with no alternate allele has an empty list of decisions, which escapes nothing.
#[test]
fn a_record_with_no_alternate_is_filtered() {
    let text = golden();
    let filtered = run(
        &text,
        "no-alt",
        &Arguments {
            median_autosomal_coverage: 30.0,
            ..Arguments::default()
        },
    )
    .expect("a run");
    assert!(filtered[0].contains("\tpossible_numt\t"));
    // And the attribute is left exactly as it was, since nothing was true.
    assert!(filtered[0].contains("AS_FilterStatus=SITE\t"));
}

/// The maximum across samples, not the sum: two samples of fifty under a cutoff of seventy-nine.
#[test]
fn the_depth_compared_is_the_maximum_not_the_sum() {
    let text = golden();
    let filtered = run(
        &text,
        "coverage-30",
        &Arguments {
            median_autosomal_coverage: 30.0,
            ..Arguments::default()
        },
    )
    .expect("a run");
    let record = filtered
        .iter()
        .find(|line| line.starts_with("chrM\t400\t"))
        .expect("the record");
    assert!(record.contains("\tpossible_numt\t"));
    assert!(record.contains("0/1:100,50\t0/1:100,50"));
    assert_eq!(max_alt_depth_cutoff(4.0, 30.0).expect("a cutoff"), 79);
}

/// An existing allele filter is unioned, a repeat is dropped, and SITE is replaced.
#[test]
fn an_existing_allele_filter_is_unioned() {
    let text = golden();
    let filtered = run(
        &text,
        "coverage-30",
        &Arguments {
            median_autosomal_coverage: 30.0,
            ..Arguments::default()
        },
    )
    .expect("a run");
    let by_position = |position: &str| {
        filtered
            .iter()
            .find(|line| line.starts_with(position))
            .expect("the record")
            .clone()
    };
    assert!(by_position("chrM\t600\t").contains("AS_FilterStatus=weak_evidence,possible_numt"));
    assert!(by_position("chrM\t700\t").contains("AS_FilterStatus=possible_numt\t"));
    assert!(by_position("chrM\t200\t").contains("AS_FilterStatus=possible_numt\t"));
    // A multiallelic record can carry the attribute without carrying the site filter.
    let multiallelic = by_position("chrM\t300\t");
    assert!(multiallelic.contains("AS_FilterStatus=SITE|possible_numt"));
    assert!(multiallelic.contains("\t.\tAS_FilterStatus="));
}
