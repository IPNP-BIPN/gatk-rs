//! Conformance for the counting annotations, against the oracle.
//!
//! Golden from `tools/annotation-conformance/AnnotationDump.java`, which asks the three real
//! annotation objects through their real `InfoFieldAnnotation` interface.
//!
//! The golden corrected three things in the port, and each one is a rule about Java rather than
//! about annotations:
//!
//! ```text
//! combine  empty            E:...:Raw value for RAW_GT_COUNT has 1 values, expected 3. ...
//! combine  overflowing-sum  .,-2147483648,-2147483648
//! combine  space-before-comma  E:...:malformed RAW_GT_COUNT annotation: 1 , 2, 3
//! ```
//!
//! `"".split(", *")` is a one-element array holding the empty string, so the empty field is
//! reported as **one** value and not zero. The sums are `int` additions with no overflow check, so
//! two maxima wrap to `Integer.MIN_VALUE` and are written out. And only a comma followed by
//! *spaces* splits, so `"1 , 2, 3"` has the right arity and a first field of `"1 "`, failing on the
//! integer rather than on the count. A tab after the comma fails the same way, because a tab is
//! not a space.
//!
//! The exception messages carry no `A USER ERROR has occurred:` banner: that belongs to the command
//! line's printer, not to the exception a caller sees.

use std::io::Read;

use gatk_annotation::chromosome_counts::ChromosomeCounts;
use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_annotation::raw_gt_count::{combine_raw_data, RawGtCount};
use gatk_annotation::sample_list::SampleList;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/annotations.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::from_str(bases, is_ref).expect("an allele")
}

fn gt(sample: &str, alleles: Vec<Allele>) -> Genotype {
    Genotype::new(sample, alleles)
}

fn filtered(mut genotype: Genotype, filter: &str) -> Genotype {
    genotype.filters = Some(filter.to_string());
    genotype
}

fn build(alleles: Vec<Allele>, genotypes: Vec<Genotype>) -> VariantContext {
    let mut vc = VariantContext::new("chr1", 100, alleles);
    vc.stop = 100;
    vc.genotypes = genotypes;
    vc
}

/// The dump's fixtures, in its order.
fn fixtures() -> Vec<(&'static str, VariantContext)> {
    let r = allele("A", true);
    let a1 = allele("C", false);
    let a2 = allele("G", false);
    let n = Allele::no_call();
    vec![
        ("no-genotypes", build(vec![r.clone(), a1.clone()], vec![])),
        ("ref-only-no-genotypes", build(vec![r.clone()], vec![])),
        (
            "one-het",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![r.clone(), a1.clone()])],
            ),
        ),
        (
            "one-hom-var",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![a1.clone(), a1.clone()])],
            ),
        ),
        (
            "all-hom-ref",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    gt("s1", vec![r.clone(), r.clone()]),
                    gt("s2", vec![r.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "two-alts",
            build(
                vec![r.clone(), a1.clone(), a2.clone()],
                vec![
                    gt("s1", vec![r.clone(), a1.clone()]),
                    gt("s2", vec![a1.clone(), a2.clone()]),
                ],
            ),
        ),
        (
            "ref-only-site",
            build(vec![r.clone()], vec![gt("s1", vec![r.clone(), r.clone()])]),
        ),
        (
            "all-no-call",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![n.clone(), n.clone()])],
            ),
        ),
        (
            "half-no-call",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![a1.clone(), n.clone()])],
            ),
        ),
        (
            "mixed-and-hom-ref",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    gt("s1", vec![r.clone(), n.clone()]),
                    gt("s2", vec![a1.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "filtered-het",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    filtered(gt("s1", vec![r.clone(), a1.clone()]), "LowGQ"),
                    gt("s2", vec![r.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "het-two-alts",
            build(
                vec![r.clone(), a1.clone(), a2.clone()],
                vec![
                    gt("s1", vec![a1.clone(), a2]),
                    gt("s2", vec![r.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "name-order",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    gt("b", vec![r.clone(), a1.clone()]),
                    gt("A", vec![r.clone(), a1.clone()]),
                    gt("a", vec![r.clone(), a1.clone()]),
                    gt("B", vec![r.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "three-hets",
            build(
                vec![r.clone(), a1.clone()],
                (0..3)
                    .map(|i| gt(&format!("m{i}"), vec![r.clone(), a1.clone()]))
                    .collect(),
            ),
        ),
    ]
}

/// The dump's `combine` fixtures, in its order.
fn combine_cases() -> Vec<(&'static str, Vec<&'static str>)> {
    vec![
        ("single", vec!["[1, 2, 3]"]),
        ("single-no-brackets", vec!["1,2,3"]),
        ("single-no-spaces", vec!["[1,2,3]"]),
        ("doubled", vec!["[1, 2, 3]", "[1, 2, 3]"]),
        (
            "three-ways",
            vec!["[1, 2, 3]", "[10, 20, 30]", "[100, 200, 300]"],
        ),
        ("zeroes", vec!["[0, 0, 0]"]),
        ("round-trip", vec![".,2,3"]),
        ("leading-space", vec!["  [1, 2, 3]  "]),
        ("space-before-comma", vec!["1 , 2, 3"]),
        ("many-spaces", vec!["1,    2,    3"]),
        ("tab-separated", vec!["1,\t2,\t3"]),
        ("brackets-inside", vec!["[1, 2], [3]"]),
        ("two-values", vec!["[1, 2]"]),
        ("four-values", vec!["[1, 2, 3, 4]"]),
        ("empty", vec![""]),
        ("trailing-comma", vec!["1,2,3,"]),
        ("negative", vec!["[-1, -2, -3]"]),
        ("plus-sign", vec!["[+1, 2, 3]"]),
        ("not-a-number", vec!["[a, 2, 3]"]),
        ("overflow", vec!["[2147483648, 2, 3]"]),
        ("max-int", vec!["[2147483647, 2147483647, 2147483647]"]),
        (
            "overflowing-sum",
            vec!["[2147483647, 2147483647, 2147483647]", "[1, 1, 1]"],
        ),
    ]
}

/// The golden's rows of one kind, keyed on everything before the last tab-separated value.
fn value(text: &str, prefix: &str) -> String {
    let needle = format!("{prefix}\t");
    text.lines()
        .find(|line| line.starts_with(&needle))
        .unwrap_or_else(|| panic!("no row for {prefix:?}"))[needle.len()..]
        .to_string()
}

/// The dump's rendering of one map: `key=value[class]` joined with `;`, a `Double` as raw bits and
/// a list in parentheses.
fn render(entries: &[(String, AnnotationValue)]) -> String {
    entries
        .iter()
        .map(|(key, value)| format!("{key}={}[{}]", render_value(value), value.java_class()))
        .collect::<Vec<_>>()
        .join(";")
}

fn render_value(value: &AnnotationValue) -> String {
    match value {
        AnnotationValue::Double(d) => (d.to_bits() as i64).to_string(),
        AnnotationValue::List(values) => format!(
            "({})",
            values
                .iter()
                .map(render_value)
                .collect::<Vec<_>>()
                .join(",")
        ),
        other => other
            .to_java_string()
            .expect("every non-Double value renders"),
    }
}

#[test]
fn key_names_match_the_reference() {
    let text = golden();
    assert_eq!(
        ChromosomeCounts.key_names().join(","),
        value(&text, "keys\tChromosomeCounts")
    );
    assert_eq!(
        SampleList.key_names().join(","),
        value(&text, "keys\tSampleList")
    );
    assert_eq!(
        RawGtCount.key_names().join(","),
        value(&text, "keys\tRawGtCount")
    );
}

#[test]
fn chromosome_counts_matches_the_reference() {
    let text = golden();
    for (label, vc) in fixtures() {
        let ours = render(&ChromosomeCounts.annotate(None, &vc));
        assert_eq!(
            ours,
            value(&text, &format!("anno\tChromosomeCounts\t{label}")),
            "ChromosomeCounts on {label}"
        );
    }
}

#[test]
fn sample_list_matches_the_reference() {
    let text = golden();
    for (label, vc) in fixtures() {
        let ours = render(&SampleList.annotate(None, &vc));
        assert_eq!(
            ours,
            value(&text, &format!("anno\tSampleList\t{label}")),
            "SampleList on {label}"
        );
    }
}

/// `RawGtCount.annotate` returns null on every fixture, which the dump prints as `null` and the
/// port reports as an empty map. The two are different in Java and indistinguishable to a caller
/// that only writes what it is given, so the test asserts the reference's answer is uniform rather
/// than pretending the port reproduces the null.
#[test]
fn raw_gt_count_annotates_nothing_anywhere() {
    let text = golden();
    for (label, vc) in fixtures() {
        assert_eq!(
            value(&text, &format!("anno\tRawGtCount\t{label}")),
            "null",
            "the reference returned something for {label}"
        );
        assert!(RawGtCount.annotate(None, &vc).is_empty());
    }
}

#[test]
fn every_combine_matches_the_reference() {
    let text = golden();
    for (label, raws) in combine_cases() {
        let owned: Vec<String> = raws.iter().map(|raw| (*raw).to_string()).collect();
        let ours = match combine_raw_data(&owned) {
            Ok(combined) => combined,
            // The dump replaces tabs with spaces in the message it prints, since its own format is
            // tab-separated; the same substitution is applied here rather than to the port.
            Err(error) => format!(
                "E:{}:{}",
                error.class(),
                error.message().replace(['\t', '\n'], " ")
            ),
        };
        assert_eq!(
            ours,
            value(&text, &format!("combine\t{label}")),
            "combineRawData on {label}"
        );
    }
}

/// The rows that a port gets wrong by writing the obvious Rust.
#[test]
fn the_rows_that_rust_gets_wrong_by_default() {
    let text = golden();

    // The sums are unchecked int additions: two maxima wrap and are written out. Rust panics on
    // overflow in a debug build, so this row is the difference between a golden and a crash.
    assert_eq!(
        value(&text, "combine\toverflowing-sum"),
        ".,-2147483648,-2147483648"
    );

    // "".split(", *") is one field, not zero, so the arity error says 1.
    assert!(
        value(&text, "combine\tempty").ends_with("has 1 values, expected 3. Annotation value is ")
    );

    // Only a comma followed by SPACES splits, so this has three fields and a bad first one.
    assert!(value(&text, "combine\tspace-before-comma").contains("malformed"));
    // And a tab is not a space, so the same shape with tabs fails the same way.
    assert!(value(&text, "combine\ttab-separated").contains("malformed"));

    // The combine drops the hom-ref count it just summed, so its own output cannot be fed back in.
    assert_eq!(value(&text, "combine\tdoubled"), ".,4,6");
    assert!(value(&text, "combine\tround-trip").contains("malformed"));

    // AC changes Java class with the alternate count, and the dump reports the class.
    assert!(value(&text, "anno\tChromosomeCounts\tone-het").contains("AC=1[java.lang.Integer]"));
    assert!(
        value(&text, "anno\tChromosomeCounts\ttwo-alts").contains("AC=(2,1)[java.util.ArrayList]")
    );

    // SampleList iterates in name order, not storage order, and drops the hom-ref sample.
    assert_eq!(
        value(&text, "anno\tSampleList\tname-order"),
        "Samples=A,a,b[java.lang.String]"
    );
    // Its value is one String with commas in it, not a list.
    assert!(value(&text, "anno\tSampleList\tthree-hets").ends_with("[java.lang.String]"));

    // And the guard is monomorphism, so a site whose only alt-carrying genotype is filtered says
    // nothing, even though ChromosomeCounts still writes AN and AC for it.
    assert_eq!(value(&text, "anno\tSampleList\tfiltered-het"), "");
    assert!(value(&text, "anno\tChromosomeCounts\tfiltered-het").starts_with("AN=2"));
}
