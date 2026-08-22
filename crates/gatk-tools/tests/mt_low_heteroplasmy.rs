//! Conformance for `MTLowHeteroplasmyFilterTool` against GATK 4.6.2.0, compared as the whole
//! output VCF of every run.
//!
//! Golden from `tools/readfilter-conformance/MTLowHeteroplasmyDump.java`.
//!
//! # What this suite is for
//!
//!  * **`--low-het-threshold` doing nothing**, the field being a compile-time constant;
//!  * **the filter being all or nothing across the file**, three sites passing and four failing;
//!  * **an already-filtered site not counting**, and `PASS` counting;
//!  * **the threshold being strict**;
//!  * **`AF=.` throwing** in the first pass, exactly as a genotype with no `AF` does, while a
//!    missing entry inside a multi-valued array becomes `Double.MAX_VALUE`;
//!  * **and the merge into `AS_FilterStatus` throwing** when there is nothing to merge into.

use gatk_corpus as corpus;
use gatk_tools::mt_low_heteroplasmy::{
    alleles_are_artifacts, is_not_filtered, run, Arguments, LowHetError, Record, LOW_HET_THRESHOLD,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mt_low_heteroplasmy.txt.gz"),
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

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn parse(line: &str) -> Record {
    let fields: Vec<&str> = line.split('\t').collect();
    let alternates = if fields[4] == "." {
        Vec::new()
    } else {
        fields[4].split(',').map(str::to_string).collect()
    };
    let format: Vec<&str> = fields[8].split(':').collect();
    let af_index = format.iter().position(|key| *key == "AF");
    let allele_fractions = fields[9..]
        .iter()
        .map(|genotype| {
            let values: Vec<&str> = genotype.split(':').collect();
            af_index
                .and_then(|index| values.get(index))
                // A bare `.` is an ABSENT attribute, which htsjdk reports as no attribute at all.
                .filter(|text| **text != ".")
                .map(|text| {
                    text.split(',')
                        .map(|entry| {
                            if entry == "." {
                                None
                            } else {
                                Some(entry.parse().expect("a fraction"))
                            }
                        })
                        .collect()
                })
        })
        .collect();
    Record {
        alternates,
        allele_fractions,
        filters: if fields[6] == "." {
            Vec::new()
        } else {
            fields[6].split(';').map(str::to_string).collect()
        },
        as_filter_status: fields[7]
            .split(';')
            .find_map(|entry| entry.strip_prefix("AS_FilterStatus="))
            .map(str::to_string),
    }
}

fn rendered(original: &str, record: &Record) -> String {
    let mut out: Vec<String> = original.split('\t').map(str::to_string).collect();
    out[6] = if record.filters.is_empty() {
        ".".to_string()
    } else {
        record.filters.join(";")
    };
    out[7] = match &record.as_filter_status {
        Some(status) => format!("AS_FilterStatus={status}"),
        None => ".".to_string(),
    };
    out.join("\t")
}

fn records(vcf: &str) -> Vec<&str> {
    vcf.lines().filter(|line| !line.starts_with('#')).collect()
}

fn produced(text: &str, label: &str, arguments: &Arguments) -> Result<Vec<String>, LowHetError> {
    let lines = value(text, "input", label);
    let originals: Vec<String> = records(&lines).iter().map(|l| l.to_string()).collect();
    let parsed: Vec<Record> = originals.iter().map(|line| parse(line)).collect();
    let filtered = run(&parsed, arguments)?;
    Ok(originals
        .iter()
        .zip(&filtered)
        .map(|(line, record)| rendered(line, record))
        .collect())
}

fn expected(text: &str, label: &str) -> Vec<String> {
    records(&value(text, "filtered", label))
        .iter()
        .map(|line| line.to_string())
        .collect()
}

#[test]
fn every_filtered_vcf_matches_the_golden() {
    let text = golden();
    let none = Arguments {
        max_allowed_low_hets: 0,
    };
    let mut compared = 0;
    for (label, arguments) in [
        ("three-low-hets", Arguments::default()),
        ("four-low-hets", Arguments::default()),
        ("one-already-filtered", Arguments::default()),
        ("pass-not-dot", Arguments::default()),
        ("allow-none", none),
        ("threshold-0.6", none),
        ("exactly-the-threshold", Arguments::default()),
        ("af-entry-is-a-dot", none),
        ("multiallelic", none),
    ] {
        assert_eq!(
            produced(&text, label, &arguments).expect("a run that is not refused"),
            expected(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's outputs");
}

/// The threshold argument is inert: a run at 0.6 over 0.05 and 0.5 filters only the 0.05 ones,
/// which is what the default threshold does. The port therefore has no threshold argument at all.
#[test]
fn the_threshold_argument_does_nothing() {
    let text = golden();
    assert_eq!(LOW_HET_THRESHOLD, 0.1);
    let filtered = produced(
        &text,
        "threshold-0.6",
        &Arguments {
            max_allowed_low_hets: 0,
        },
    )
    .expect("a run");
    // The 0.5 record is untouched, though 0.5 is below the 0.6 that was asked for.
    let last = filtered.last().expect("a record");
    assert!(last.contains("0/1:0.5"));
    assert!(last.contains("\t.\tAS_FilterStatus=SITE"));
    assert_eq!(filtered, expected(&text, "threshold-0.6"));
}

/// Three unfiltered low sites pass, four fail, and one record's fate is decided by the others.
#[test]
fn the_filter_is_all_or_nothing_across_the_file() {
    let text = golden();
    let three = produced(&text, "three-low-hets", &Arguments::default()).expect("a run");
    assert!(three.iter().all(|line| !line.contains("mt_many_low_hets")));

    let four = produced(&text, "four-low-hets", &Arguments::default()).expect("a run");
    assert_eq!(
        four.iter()
            .filter(|line| line.contains("mt_many_low_hets"))
            .count(),
        4
    );

    // The fourth site already carries a filter, so only three count and nothing is filtered.
    let guarded = produced(&text, "one-already-filtered", &Arguments::default()).expect("a run");
    assert!(guarded
        .iter()
        .all(|line| !line.contains("mt_many_low_hets")));
}

/// `PASS` is not a filter: htsjdk asks whether the filter set is non-empty, and PASS leaves it
/// empty, so those records count and are filtered like any other.
#[test]
fn pass_counts_exactly_like_a_dot() {
    let text = golden();
    assert!(is_not_filtered(&Record {
        alternates: vec!["C".to_string()],
        allele_fractions: Vec::new(),
        filters: vec!["PASS".to_string()],
        as_filter_status: None,
    }));
    let filtered = produced(&text, "pass-not-dot", &Arguments::default()).expect("a run");
    assert_eq!(
        filtered
            .iter()
            .filter(|line| line.contains("mt_many_low_hets"))
            .count(),
        4
    );
    assert_eq!(filtered, expected(&text, "pass-not-dot"));
}

/// Strictly below: exactly 0.1 is not low, so a file of them is left alone.
#[test]
fn the_threshold_is_strict() {
    let text = golden();
    let filtered = produced(&text, "exactly-the-threshold", &Arguments::default()).expect("a run");
    assert!(filtered
        .iter()
        .all(|line| !line.contains("mt_many_low_hets")));
    assert_eq!(
        alleles_are_artifacts(&Record {
            alternates: vec!["C".to_string()],
            allele_fractions: vec![Some(vec![Some(0.1)])],
            filters: Vec::new(),
            as_filter_status: None,
        }),
        vec![false]
    );
}

#[test]
fn the_three_refusals_match_the_golden() {
    let text = golden();
    for label in ["af-is-a-dot", "genotype-without-af"] {
        let error = produced(
            &text,
            label,
            &Arguments {
                max_allowed_low_hets: if label == "af-is-a-dot" { 0 } else { 3 },
            },
        )
        .expect_err("a refused run");
        assert_eq!(error, LowHetError::NoAlleleFraction, "{label}");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }

    let error = produced(
        &text,
        "no-as-filter-status",
        &Arguments {
            max_allowed_low_hets: 0,
        },
    )
    .expect_err("a refused run");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-as-filter-status")
    );
}

/// A missing entry inside a multi-valued array is Double.MAX_VALUE, so that alternate is never
/// low and the record keeps its site filter off.
#[test]
fn a_missing_entry_inside_the_array_is_never_low() {
    let text = golden();
    let filtered = produced(
        &text,
        "af-entry-is-a-dot",
        &Arguments {
            max_allowed_low_hets: 0,
        },
    )
    .expect("a run");
    assert!(filtered[0].contains("AS_FilterStatus=mt_many_low_hets|SITE"));
    assert!(filtered[0].contains("\t.\tAS_FilterStatus="));
}
