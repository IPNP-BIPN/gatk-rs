//! Conformance for `VCFComparator` against GATK 4.6.2.0, every run of the golden replayed through
//! the whole comparison: the merge, the grouping, the trim and the checks.
//!
//! Golden from `tools/readfilter-conformance/VCFComparatorDump.java`, which records each run's
//! ROOT cause, so the messages here are the comparison's own rather than the walker's wrap.
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
//!
//! Four runs give actual an allele its genotypes do not call, so the trim drops it and rebuilds
//! the record's INFO map from the keys an annotation claims: `AC` survives it, which is why the
//! attribute complaint still comes first.

use gatk_corpus as corpus;
use gatk_tools::vcf_comparator::{check_inputs, compare, Failure, Input, InputError, Options};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_comparator.txt.gz"),
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
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some(other) => out.push(other),
            None => {}
        }
    }
    out
}

/// The VCF the dump wrote under `name`.
fn vcf(text: &str, name: &str) -> String {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("vcf\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries vcf/{name}"));
    unescape(row)
}

/// What the reference answered for `label`: success, or the root cause's message.
fn reference(text: &str, label: &str) -> Result<(), String> {
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
    Err(unescape(message))
}

fn input(text: &str, name: &str) -> Input {
    let file = htsjdk_vcf::reader::read_vcf(text).expect("the dump's VCFs parse");
    Input {
        name: name.to_string(),
        samples: file.header.samples.clone(),
        records: file.records,
    }
}

/// The run the dump made: `-V:expected expected.vcf -V:actual <actual>.vcf` and the arguments.
fn port(text: &str, actual: &str, options: &Options) -> Result<(), Failure> {
    let inputs = [
        input(&vcf(text, "expected"), "expected"),
        input(&vcf(text, actual), "actual"),
    ];
    compare(&inputs, &["chr1".to_string()], None, None, options)
        .map(|_| ())
        .map_err(|stopped| stopped.failure)
}

fn with(change: impl FnOnce(&mut Options)) -> Options {
    let mut options = Options::default();
    change(&mut options);
    options
}

/// Every comparison the dump made, with the file it compared and the arguments it gave.
fn runs() -> Vec<(&'static str, &'static str, Options)> {
    let ignore = |key: &str| {
        let key = key.to_string();
        with(move |o| o.ignore_attributes = vec![key])
    };
    vec![
        ("identical", "same", Options::default()),
        ("extra-variant", "extra", Options::default()),
        ("extra-variant-gq0", "extra-gq0", Options::default()),
        (
            "extra-variant-gq0-ignored",
            "extra-gq0",
            with(|o| o.ignore_gq0 = true),
        ),
        ("qual-differs", "qual", Options::default()),
        (
            "qual-tolerated",
            "qual",
            with(|o| o.qual_change_allowed = 10.0),
        ),
        ("qual-ignored", "qual", with(|o| o.ignore_quals = true)),
        ("filters-differ", "filtered", Options::default()),
        (
            "filters-ignored",
            "filtered",
            with(|o| o.ignore_filters = true),
        ),
        ("filters-unapplied", "unfiltered", Options::default()),
        ("alleles-differ-hits-ac", "extra-allele", Options::default()),
        ("alleles-differ", "extra-allele", ignore("AC")),
        (
            "alleles-allowed",
            "extra-allele",
            with(|o| {
                o.ignore_attributes = vec!["AC".to_string()];
                o.allow_extra_alleles = true;
            }),
        ),
        ("allele-missing", "missing-allele", ignore("AC")),
        (
            "allele-missing-allowed",
            "missing-allele",
            with(|o| {
                o.ignore_attributes = vec!["AC".to_string()];
                o.allow_extra_alleles = true;
            }),
        ),
        ("info-differs", "info", Options::default()),
        ("info-ignored-key", "info", ignore("DP")),
        ("info-ignored-wrong-key", "info", ignore("AC")),
        (
            "info-ignored-all",
            "info",
            with(|o| o.ignore_annotations = true),
        ),
        ("ids-differ", "ids", Options::default()),
        (
            "positions-only-alleles",
            "extra-allele",
            with(|o| o.positions_only = true),
        ),
        ("warn-on-errors", "qual", with(|o| o.warn_on_errors = true)),
    ]
}

#[test]
fn every_run_of_the_golden_is_reproduced() {
    let text = golden();
    for (label, actual, options) in runs() {
        let expected = reference(&text, label);
        let answered = port(&text, actual, &options);
        let answered = answered.map_err(|failure| match failure {
            Failure::User(message) => message,
            other => panic!("{label}: {other:?}"),
        });
        assert_eq!(answered, expected, "{label}");
    }
}

#[test]
fn the_expected_file_is_named_by_its_tag_and_counted_first() {
    let text = golden();
    let names = |list: &[&str]| list.iter().map(|name| name.to_string()).collect::<Vec<_>>();
    assert_eq!(
        check_inputs(&names(&["expected"])).map_err(|error| error.message()),
        reference(&text, "one-input")
    );
    assert_eq!(
        check_inputs(&names(&["first", "second"])).map_err(|error| error.message()),
        reference(&text, "no-expected")
    );
    // The count is checked before the tag, so one input tagged anything is the count's refusal.
    assert_eq!(
        check_inputs(&names(&["first"])),
        Err(InputError::WrongNumberOfInputs)
    );
    assert_eq!(check_inputs(&names(&["actual", "expected"])), Ok(()));
}

/// A difference inside the traversal names the record the walker was on, which is the one after
/// the group; a difference in the last group names none.
#[test]
fn a_complaint_is_raised_at_the_next_record_or_after_the_traversal() {
    let text = golden();
    let inputs = [
        input(&vcf(&text, "expected"), "expected"),
        input(&vcf(&text, "qual"), "actual"),
    ];
    let stopped = compare(
        &inputs,
        &["chr1".to_string()],
        None,
        None,
        &Options::default(),
    )
    .expect_err("the QUAL at 1000 differs");
    // Expected is named first, so it supplies the first record of each tied pair; the group at
    // 1000 is compared when expected's record at 2000 arrives.
    assert_eq!(stopped.at, Some((0, 1)));

    let inputs = [
        input(&vcf(&text, "expected"), "expected"),
        input(&vcf(&text, "missing-allele"), "actual"),
    ];
    let stopped = compare(
        &inputs,
        &["chr1".to_string()],
        None,
        None,
        &with(|o| o.ignore_attributes = vec!["AC".to_string()]),
    )
    .expect_err("the site at 4000 lost an allele");
    assert_eq!(stopped.at, None);
}
