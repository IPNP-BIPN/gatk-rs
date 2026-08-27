//! Conformance for `VCFComparator` against GATK 4.6.2.0, compared as the complaint of every run.
//!
//! Golden from `tools/readfilter-conformance/VCFComparatorDump.java`.
//!
//! # What this suite is for
//!
//!  * **the allele check being guarded by the reversed comparison**, so an allele added to actual
//!    is never checked and one missing from it is;
//!  * **an unmatched variant being guarded on a genotype quality of zero**;
//!  * **the position being wrapped around the message rather than inside it**;
//!  * **different filters and unapplied filters being two different complaints**;
//!  * **`--ignore-attribute` taking one key at a time**;
//!  * **and the expected file being named by its tag rather than by its order.**

use gatk_corpus as corpus;
use gatk_tools::vcf_comparator::{
    check_inputs, compare, has_new_alleles, unmatched, Complaint, InputError, Tolerances, Variant,
    DEFAULT_QUAL_CHANGE_ALLOWED,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_comparator.txt.gz"),
    )
}

/// A run either succeeded or produced one message.
fn result(text: &str, label: &str) -> Result<(), String> {
    if text
        .lines()
        .any(|line| line == format!("ok\t{label}=succeeded"))
    {
        return Ok(());
    }
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries a result for {label}"));
    let (_, message) = row.split_once(':').expect("a class and a message");
    Err(message.to_string())
}

fn class(text: &str, label: &str) -> String {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    row.split_once(':').expect("a class").0.to_string()
}

#[allow(clippy::too_many_arguments)]
fn variant(
    source: &str,
    start: i32,
    id: &str,
    reference: &str,
    alternates: &[&str],
    qual: f64,
    filters: &[&str],
    attributes: &[(&str, &str)],
    genotype_qualities: &[i32],
) -> Variant {
    Variant {
        source: source.to_string(),
        contig: "chr1".to_string(),
        start,
        id: id.to_string(),
        reference: reference.to_string(),
        alternates: alternates.iter().map(|a| a.to_string()).collect(),
        qual: Some(qual),
        filters: filters.iter().map(|f| f.to_string()).collect(),
        attributes: attributes
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        genotype_qualities: genotype_qualities.to_vec(),
    }
}

/// The expected file's first site, which every actual is a change away from.
fn expected_1000() -> Variant {
    variant(
        "expected",
        1000,
        ".",
        "A",
        &["C"],
        100.0,
        &["PASS"],
        &[("AC", "1"), ("DP", "20")],
        &[50, 50],
    )
}

/// The expected file's third site, which carries two alternates.
fn expected_4000() -> Variant {
    variant(
        "expected",
        4000,
        ".",
        "C",
        &["A", "T"],
        400.0,
        &["PASS"],
        &[("AC", "[1, 1]"), ("DP", "50")],
        &[80, 80],
    )
}

fn ignoring(keys: &[&str]) -> Tolerances {
    Tolerances {
        ignored_attributes: keys.iter().map(|k| k.to_string()).collect(),
        ..Tolerances::default()
    }
}

#[test]
fn every_comparison_matches_the_golden() {
    let text = golden();
    let mut compared = 0;

    let check = |label: &str, actual: &Variant, expected: &Variant, tolerances: &Tolerances| {
        let produced = compare(actual, expected, tolerances)
            .map(|complaint| complaint.wrapped(&expected.contig, expected.start));
        match (produced, result(&text, label)) {
            (None, Ok(())) => {}
            (Some(message), Err(golden)) => assert_eq!(message, golden, "{label}"),
            (produced, golden) => panic!("{label}: {produced:?} against {golden:?}"),
        }
    };

    // Identical.
    check(
        "identical",
        &expected_1000(),
        &expected_1000(),
        &Tolerances::default(),
    );
    compared += 1;

    // A different QUAL, alone and under three arguments.
    let qual = Variant {
        source: "actual".to_string(),
        qual: Some(105.0),
        ..expected_1000()
    };
    check(
        "qual-differs",
        &qual,
        &expected_1000(),
        &Tolerances::default(),
    );
    check(
        "qual-tolerated",
        &qual,
        &expected_1000(),
        &Tolerances {
            qual_change_allowed: 10.0,
            ..Tolerances::default()
        },
    );
    check(
        "qual-ignored",
        &qual,
        &expected_1000(),
        &Tolerances {
            ignore_quals: true,
            ..Tolerances::default()
        },
    );
    compared += 3;

    // Filters: a different one, and one applied on only one side.
    let filtered = Variant {
        source: "actual".to_string(),
        filters: vec!["LOW".to_string()],
        ..expected_1000()
    };
    check(
        "filters-differ",
        &filtered,
        &expected_1000(),
        &Tolerances::default(),
    );
    check(
        "filters-ignored",
        &filtered,
        &expected_1000(),
        &Tolerances {
            ignore_filters: true,
            ..Tolerances::default()
        },
    );
    let unfiltered = Variant {
        source: "actual".to_string(),
        filters: Vec::new(),
        ..expected_1000()
    };
    check(
        "filters-unapplied",
        &unfiltered,
        &expected_1000(),
        &Tolerances::default(),
    );
    compared += 3;

    // An extra allele, which changes AC with it.
    let extra_allele = variant(
        "actual",
        1000,
        ".",
        "A",
        &["C", "G"],
        100.0,
        &["PASS"],
        &[("AC", "[1, 0]"), ("DP", "20")],
        &[50, 50],
    );
    check(
        "alleles-differ-hits-ac",
        &extra_allele,
        &expected_1000(),
        &Tolerances::default(),
    );
    check(
        "alleles-differ",
        &extra_allele,
        &expected_1000(),
        &ignoring(&["AC"]),
    );
    check(
        "alleles-allowed",
        &extra_allele,
        &expected_1000(),
        &Tolerances {
            allow_extra_alleles: true,
            ..ignoring(&["AC"])
        },
    );
    compared += 3;

    // An allele MISSING from actual, which is the direction the guard reaches.
    let missing_allele = variant(
        "actual",
        4000,
        ".",
        "C",
        &["A"],
        400.0,
        &["PASS"],
        &[("AC", "1"), ("DP", "50")],
        &[80, 80],
    );
    check(
        "allele-missing",
        &missing_allele,
        &expected_4000(),
        &ignoring(&["AC"]),
    );
    check(
        "allele-missing-allowed",
        &missing_allele,
        &expected_4000(),
        &Tolerances {
            allow_extra_alleles: true,
            ..ignoring(&["AC"])
        },
    );
    compared += 2;

    // An INFO difference under three arguments.
    let info = variant(
        "actual",
        1000,
        ".",
        "A",
        &["C"],
        100.0,
        &["PASS"],
        &[("AC", "1"), ("DP", "25")],
        &[50, 50],
    );
    check(
        "info-differs",
        &info,
        &expected_1000(),
        &Tolerances::default(),
    );
    check(
        "info-ignored-key",
        &info,
        &expected_1000(),
        &ignoring(&["DP"]),
    );
    check(
        "info-ignored-wrong-key",
        &info,
        &expected_1000(),
        &ignoring(&["AC"]),
    );
    check(
        "info-ignored-all",
        &info,
        &expected_1000(),
        &Tolerances {
            ignore_annotations: true,
            ..Tolerances::default()
        },
    );
    compared += 4;

    // A dbSNP id, and positions-only over the extra allele.
    let ids = Variant {
        source: "actual".to_string(),
        id: "rs1".to_string(),
        ..expected_1000()
    };
    check("ids-differ", &ids, &expected_1000(), &Tolerances::default());
    check(
        "positions-only-alleles",
        &extra_allele,
        &expected_1000(),
        &Tolerances {
            positions_only: true,
            ..Tolerances::default()
        },
    );
    compared += 2;

    assert_eq!(compared, 18, "the comparisons the golden carries");
}

/// The guard in front of the allele check has its arguments the other way round.
#[test]
fn the_allele_check_is_guarded_by_the_reversed_comparison() {
    let text = golden();
    let expected = expected_1000();
    let extra_allele = variant(
        "actual",
        1000,
        ".",
        "A",
        &["C", "G"],
        100.0,
        &["PASS"],
        &[("AC", "[1, 0]"), ("DP", "20")],
        &[50, 50],
    );
    // Actual really does carry an allele expected lacks.
    assert!(has_new_alleles(&extra_allele, &expected));
    // And expected carries none that actual lacks, which is what the guard asks.
    assert!(!has_new_alleles(&expected, &extra_allele));
    // So with AC ignored the extra allele passes in silence.
    assert!(compare(&extra_allele, &expected, &ignoring(&["AC"])).is_none());
    assert!(result(&text, "alleles-differ").is_ok());

    // The other direction is refused, and --allow-extra-alleles does NOT silence it.
    let missing = variant(
        "actual",
        4000,
        ".",
        "C",
        &["A"],
        400.0,
        &["PASS"],
        &[("AC", "1"), ("DP", "50")],
        &[80, 80],
    );
    assert!(has_new_alleles(&expected_4000(), &missing));
    assert!(compare(&missing, &expected_4000(), &ignoring(&["AC"])).is_some());
    assert!(result(&text, "allele-missing").is_err());
    assert!(
        result(&text, "allele-missing-allowed").is_err(),
        "the flag is for the direction the guard never reaches"
    );
}

/// The complaint is guarded on a genotype quality of zero.
#[test]
fn an_unmatched_variant_needs_a_genotype_quality_of_zero() {
    let text = golden();
    let confident = variant(
        "actual",
        3000,
        ".",
        "T",
        &["A"],
        300.0,
        &["PASS"],
        &[("AC", "1"), ("DP", "40")],
        &[70, 70],
    );
    assert!(unmatched(&confident, &Tolerances::default()).is_none());
    assert!(result(&text, "extra-variant").is_ok());

    let uncertain = Variant {
        genotype_qualities: vec![0, 0],
        ..confident.clone()
    };
    let complaint = unmatched(&uncertain, &Tolerances::default()).expect("a complaint");
    assert_eq!(
        complaint.wrapped("chr1", 3000),
        result(&text, "extra-variant-gq0").expect_err("a complaint")
    );
    // It carries its own position rather than being wrapped with one.
    assert_eq!(complaint.message(), complaint.wrapped("chr1", 3000));
    assert!(!complaint.message().starts_with("At position"));

    // And --ignore-gq0 silences it.
    assert!(unmatched(
        &uncertain,
        &Tolerances {
            ignore_gq0: true,
            ..Tolerances::default()
        }
    )
    .is_none());
    assert!(result(&text, "extra-variant-gq0-ignored").is_ok());
}

/// A filter that differs and a filter that was never applied are two different complaints.
#[test]
fn filters_have_two_complaints() {
    let text = golden();
    let differs = result(&text, "filters-differ").expect_err("a complaint");
    let unapplied = result(&text, "filters-unapplied").expect_err("a complaint");
    assert_ne!(differs, unapplied);
    assert!(differs.contains("different filters"));
    assert!(unapplied.contains("not applied to both"));
    // The unapplied message carries a DOUBLE space, because its own text begins with one.
    assert!(unapplied.contains("chr1:1000  filters"), "{unapplied}");
}

/// One key at a time, so ignoring the wrong key changes nothing.
#[test]
fn ignore_attribute_takes_one_key() {
    let text = golden();
    assert!(result(&text, "info-ignored-key").is_ok());
    assert_eq!(
        result(&text, "info-ignored-wrong-key"),
        result(&text, "info-differs"),
        "ignoring AC does nothing to a DP difference"
    );
    assert!(result(&text, "info-ignored-all").is_ok());
    // And the complaint names the attribute it found.
    assert!(result(&text, "info-differs")
        .expect_err("a complaint")
        .contains("for DP:"));
}

/// The tag, not the order, and the parser's own mutual exclusion.
#[test]
fn the_expected_file_is_named_by_its_tag() {
    let text = golden();
    assert!(check_inputs(&["expected".to_string(), "actual".to_string()]).is_ok());

    let one = check_inputs(&["expected".to_string()]).expect_err("one input");
    assert_eq!(one, InputError::WrongNumberOfInputs);
    assert_eq!(
        result(&text, "one-input").expect_err("a complaint"),
        one.message()
    );

    let untagged =
        check_inputs(&["first".to_string(), "second".to_string()]).expect_err("no expected");
    assert_eq!(untagged, InputError::NoExpectedInput);
    assert_eq!(
        result(&text, "no-expected").expect_err("a complaint"),
        untagged.message()
    );
    // The COUNT is checked first, so one untagged input reports the count.
    assert_eq!(
        check_inputs(&["first".to_string()]).expect_err("one input"),
        InputError::WrongNumberOfInputs
    );

    // --positions-only is refused beside a tolerance it subsumes, by the parser.
    assert_eq!(
        class(&text, "positions-only-with-quals"),
        "org.broadinstitute.barclay.argparser.CommandLineException"
    );
    let exclusive = InputError::MutuallyExclusive {
        argument: "ignore-quals".to_string(),
        other: "positions-only".to_string(),
    };
    assert!(result(&text, "positions-only-with-quals")
        .expect_err("a complaint")
        .starts_with(&exclusive.message()));
}

/// It is not zero, so a difference of a hair is tolerated and one of five is not.
#[test]
fn the_qual_tolerance_has_a_default() {
    assert_eq!(DEFAULT_QUAL_CHANGE_ALLOWED, 0.001);
    let expected = expected_1000();
    let hair = Variant {
        source: "actual".to_string(),
        qual: Some(100.0005),
        ..expected_1000()
    };
    assert!(compare(&hair, &expected, &Tolerances::default()).is_none());
    let five = Variant {
        source: "actual".to_string(),
        qual: Some(105.0),
        ..expected_1000()
    };
    let complaint = compare(&five, &expected, &Tolerances::default()).expect("a complaint");
    assert_eq!(
        complaint,
        Complaint::QualDiffers {
            difference: 5.0,
            tolerance: 0.001
        }
    );
    // The message names both numbers, and the difference goes through Double.toString.
    assert!(complaint.message().contains("differ by 5.0"));
    assert!(complaint.message().ends_with("more than 0.001"));
}
