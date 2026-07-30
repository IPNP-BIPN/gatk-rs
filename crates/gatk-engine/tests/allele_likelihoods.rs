//! Conformance for the read side of `AlleleLikelihoods`, against the oracle.
//!
//! Golden from `tools/genotyper-conformance/AlleleLikelihoodsDump.java`.
//!
//! The likelihoods are compared as **raw bits**: they are read back out of the matrix after a
//! conversion, and a decimal rendering would hide a divergence in the last place.

use std::io::Read;

use gatk_engine::allele_likelihoods::{AlleleLikelihoods, BestAllele};
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_vcf::allele::Allele;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/allele_likelihoods.txt.gz");
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

fn reference() -> Allele {
    allele("A", true)
}

fn alt1() -> Allele {
    allele("C", false)
}

fn alt2() -> Allele {
    allele("G", false)
}

/// The dump's rendering of an allele: the display string, `*` for a reference, `null` for none.
fn show(allele: Option<&Allele>) -> String {
    match allele {
        None => "null".to_string(),
        Some(allele) => format!(
            "{}{}",
            allele.display_string(),
            if allele.is_reference() { "*" } else { "" }
        ),
    }
}

fn bits(value: f64) -> i64 {
    value.to_bits() as i64
}

/// One fixture: the alleles, the per-sample evidence counts in order, and the matrix.
struct Case {
    label: &'static str,
    alleles: Vec<Allele>,
    samples: Vec<(&'static str, usize)>,
    values: Vec<Vec<Vec<f64>>>,
}

fn cases() -> Vec<Case> {
    let (r, a1, a2) = (reference(), alt1(), alt2());
    let case = |label, alleles: Vec<Allele>, samples: Vec<(&'static str, usize)>, values| Case {
        label,
        alleles,
        samples,
        values,
    };
    vec![
        case(
            "one-allele",
            vec![r.clone()],
            vec![("s1", 1)],
            vec![vec![vec![-1.0]]],
        ),
        case(
            "two-alleles",
            vec![r.clone(), a1.clone()],
            vec![("s1", 3)],
            vec![vec![vec![-1.0, -3.0, -0.5], vec![-2.0, -0.5, -0.5]]],
        ),
        case(
            "tie",
            vec![r.clone(), a1.clone(), a2.clone()],
            vec![("s1", 1)],
            vec![vec![vec![-1.0], vec![-1.0], vec![-1.0]]],
        ),
        case(
            "tie-between-alts",
            vec![r.clone(), a1.clone(), a2],
            vec![("s1", 1)],
            vec![vec![vec![-5.0], vec![-1.0], vec![-1.0]]],
        ),
        case(
            "just-informative",
            vec![r.clone(), a1.clone()],
            vec![("s1", 1)],
            vec![vec![vec![-1.0], vec![-1.200_000_1]]],
        ),
        case(
            "just-uninformative",
            vec![r.clone(), a1.clone()],
            vec![("s1", 1)],
            vec![vec![vec![-1.0], vec![-1.199_999_9]]],
        ),
        case(
            "exactly-threshold",
            vec![r.clone(), a1.clone()],
            vec![("s1", 1)],
            vec![vec![vec![-1.0], vec![-1.2]]],
        ),
        case(
            "all-infinite",
            vec![r.clone(), a1.clone()],
            vec![("s1", 1)],
            vec![vec![vec![f64::NEG_INFINITY], vec![f64::NEG_INFINITY]]],
        ),
        case(
            "one-infinite",
            vec![r.clone(), a1.clone()],
            vec![("s1", 1)],
            vec![vec![vec![-1.0], vec![f64::NEG_INFINITY]]],
        ),
        case(
            "with-nan",
            vec![r.clone(), a1.clone()],
            vec![("s1", 1)],
            vec![vec![vec![f64::NAN], vec![-1.0]]],
        ),
        case(
            "two-samples",
            vec![r.clone(), a1.clone()],
            vec![("s1", 2), ("s2", 1)],
            vec![
                vec![vec![-1.0, -2.0], vec![-2.0, -1.0]],
                vec![vec![-1.0], vec![-3.0]],
            ],
        ),
        case(
            "empty-sample",
            vec![r, a1],
            vec![("s1", 0)],
            vec![vec![vec![], vec![]]],
        ),
        case("no-alleles", vec![], vec![("s1", 1)], vec![vec![]]),
    ]
}

impl Case {
    fn build(&self) -> AlleleLikelihoods<String> {
        let names: Vec<String> = self
            .samples
            .iter()
            .map(|(name, _)| (*name).to_string())
            .collect();
        let evidence: Vec<Vec<String>> = self
            .samples
            .iter()
            .map(|(name, count)| (0..*count).map(|i| format!("{name}-r{i}")).collect())
            .collect();
        AlleleLikelihoods::new(
            SampleList::new(&names),
            AlleleList::new(&self.alleles),
            evidence,
            self.values.clone(),
        )
        .expect("a well-formed matrix")
    }
}

/// The dump's `best` row for a matrix, one line per piece of evidence.
fn best_rows(label: &str, likelihoods: &AlleleLikelihoods<String>) -> Vec<String> {
    let mut rows = Vec::new();
    for sample in 0..likelihoods.number_of_samples() {
        let evidence = likelihoods
            .sample_evidence(sample)
            .expect("the sample's evidence")
            .to_vec();
        for (index, best) in likelihoods
            .best_alleles_breaking_ties_for_sample(sample, None)
            .into_iter()
            .enumerate()
        {
            rows.push(render_best(label, &best, &evidence[index]));
        }
    }
    rows
}

fn render_best(label: &str, best: &BestAllele, evidence: &str) -> String {
    format!(
        "best\t{label}\t{}:{}\t{}\t{}\t{}\t{}\t{}\t{}",
        best.sample,
        evidence,
        show(best.allele.as_ref()),
        show(best.second_best_allele.as_ref()),
        bits(best.likelihood),
        bits(best.second_best_likelihood),
        bits(best.confidence),
        best.is_informative(),
    )
}

fn golden_rows(text: &str, prefix: &str) -> Vec<String> {
    text.lines()
        .filter(|line| line.starts_with(prefix))
        .map(String::from)
        .collect()
}

#[test]
fn every_matrix_reports_the_reference_counts() {
    let text = golden();
    for case in cases() {
        let likelihoods = case.build();
        let per_sample: Vec<String> = (0..likelihoods.number_of_samples())
            .map(|s| likelihoods.sample_evidence_count(s).to_string())
            .collect();
        let ours = format!(
            "counts\t{}\t{}\t{}\t{}\t{}",
            case.label,
            likelihoods.number_of_samples(),
            likelihoods.number_of_alleles(),
            likelihoods.evidence_count(),
            per_sample.join(",")
        );
        let expected = golden_rows(&text, &format!("counts\t{}\t", case.label));
        assert_eq!(expected.len(), 1, "one counts row for {}", case.label);
        assert_eq!(ours, expected[0], "counts for {}", case.label);
    }
}

#[test]
fn every_best_allele_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let likelihoods = case.build();
        let ours = best_rows(case.label, &likelihoods);
        let expected = golden_rows(&text, &format!("best\t{}\t", case.label));
        assert_eq!(ours, expected, "the best alleles for {}", case.label);
    }
}

/// The natural-log switch, and the two thresholds that stop agreeing after it.
#[test]
fn the_natural_log_switch_matches_the_reference() {
    let text = golden();
    let row = golden_rows(&text, "natural\tswitch\t");
    assert_eq!(row.len(), 1, "one natural row");
    let fields: Vec<&str> = row[0].split('\t').skip(2).collect();

    let mut likelihoods = AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&[reference(), alt1()]),
        vec![vec!["s1-r0".to_string()]],
        vec![vec![vec![-1.0], vec![-1.1]]],
    )
    .expect("a well-formed matrix");

    let before = likelihoods.best_alleles_breaking_ties(None).remove(0);
    likelihoods
        .switch_to_natural_log()
        .expect("the first switch is allowed");
    let after = likelihoods.best_alleles_breaking_ties(None).remove(0);

    assert_eq!(bits(before.confidence).to_string(), fields[0]);
    assert_eq!(bits(after.confidence).to_string(), fields[1]);
    assert_eq!(before.is_informative().to_string(), fields[2]);
    assert_eq!(after.is_informative().to_string(), fields[3]);
    assert_eq!(likelihoods.is_natural_log().to_string(), fields[4]);

    // Switching twice is refused.
    let twice = golden_rows(&text, "natural\tswitch-twice\t");
    assert_eq!(twice.len(), 1);
    assert!(
        twice[0].contains("E:"),
        "the reference refused the second switch"
    );
    assert!(likelihoods.switch_to_natural_log().is_err());
}

/// The rows a port gets wrong by reading the signature rather than the body.
#[test]
fn the_rows_that_the_signature_hides() {
    let text = golden();
    let best = |label: &str| -> Vec<String> { golden_rows(&text, &format!("best\t{label}\t")) };

    // One allele: the second best is that same allele, with a likelihood of negative infinity,
    // because both indices started at zero and the tail found them equal.
    let one = &best("one-allele")[0];
    assert!(
        one.contains("\tA*\tA*\t"),
        "the second best is the only allele"
    );
    assert!(
        one.contains(&format!("\t{}\t", bits(f64::NEG_INFINITY))),
        "and its likelihood is negative infinity"
    );

    // A three-way tie goes to the earliest allele, which is the reference here.
    assert!(best("tie")[0].contains("\tA*\t"));
    // A tie between the two alternates goes to the earlier alternate.
    assert!(best("tie-between-alts")[0].contains("\tC\t"));

    // Two negative infinities give a confidence of zero rather than a NaN, which is what the
    // equality guard is really for.
    assert!(best("all-infinite")[0].contains(&format!("\t{}\t", bits(0.0))));
}
