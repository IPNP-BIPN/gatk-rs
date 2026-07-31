//! Conformance for `QD`, the genotype summaries and `LikelihoodRankSum`, against the oracle.
//!
//! Golden from `tools/annotation-conformance/SiteStatisticsDump.java`.
//!
//! ```text
//! depth  two-samples          10       not 20: the AD-restricted tally replaced the total
//! qd     just-below-threshold QD=34.90
//! qd     at-threshold         QD=28.73  a random draw, not a function of the data
//! ```
//!
//! The last two rows are the point. At and above 35 the reference replaces `QD` with 30 plus a
//! Gaussian jitter, so the value is not reproducible from the input. The port refuses those rows,
//! and this test asserts that it refuses **exactly** the rows the reference randomised: the raw
//! ratio is recomputed here and compared against the threshold, so the boundary itself is
//! measured even though what lies past it is not.

use std::io::Read;

use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_annotation::site_statistics::{
    qual_by_depth, qual_by_depth_depth, GenotypeSummaries, QualByDepthError,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/site_statistics.txt.gz");
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

/// The `QD` fixtures, by label: the QUAL, the per-sample ADs, a raw approximation, and whether the
/// genotypes are hom-ref.
/// One `QD` fixture: the QUAL as a log10 error, the per-sample ADs, a raw approximation, and
/// whether the genotypes are hom-ref.
struct QdFixture {
    log10: Option<f64>,
    ads: Option<Vec<Vec<i32>>>,
    raw: Option<i32>,
    hom_ref: bool,
}

fn qd_case(label: &str) -> (VariantContext, Option<i32>) {
    let (log10, ads, raw, hom_ref) = {
        let fixture = match label {
            "ordinary" => (Some(-25.0), Some(vec![vec![5, 5]]), None, false),
            "just-below-threshold" => (Some(-34.9), Some(vec![vec![5, 5]]), None, false),
            "at-threshold" => (Some(-35.0), Some(vec![vec![5, 5]]), None, false),
            "above-threshold" => (Some(-100.0), Some(vec![vec![5, 5]]), None, false),
            "no-qual" => (None, Some(vec![vec![5, 5]]), None, false),
            "raw-qual-approx" => (None, Some(vec![vec![5, 5]]), Some(300), false),
            "hom-ref-only" => (Some(-25.0), Some(vec![vec![10, 0]]), None, true),
            "one-alt-read" => (Some(-20.0), Some(vec![vec![9, 1]]), None, false),
            "two-alt-reads" => (Some(-20.0), Some(vec![vec![8, 2]]), None, false),
            "zero-ad" => (Some(-25.0), Some(vec![vec![0, 0]]), None, false),
            "no-ad-with-dp" => (Some(-25.0), None, None, false),
            "two-samples" => (Some(-20.0), Some(vec![vec![9, 1], vec![8, 2]]), None, false),
            other => panic!("{other} has no fixture"),
        };
        let fixture = QdFixture {
            log10: fixture.0,
            ads: fixture.1,
            raw: fixture.2,
            hom_ref: fixture.3,
        };
        (fixture.log10, fixture.ads, fixture.raw, fixture.hom_ref)
    };

    let mut vc = VariantContext::new("chr1", START, vec![allele("A", true), allele("C", false)]);
    vc.stop = START;
    // htsjdk stores "no QUAL" as `NO_LOG10_PERROR`, which is 1.0.
    vc.log10_p_error = log10.unwrap_or(1.0);
    let sample_count = ads.as_ref().map(|a| a.len()).unwrap_or(1);
    for index in 0..sample_count {
        let called = if hom_ref {
            vec![allele("A", true), allele("A", true)]
        } else {
            vec![allele("A", true), allele("C", false)]
        };
        let mut genotype = Genotype::new(&format!("s{index}"), called);
        match &ads {
            Some(ads) => genotype.ad = Some(ads[index].clone()),
            None => genotype.dp = Some(17),
        }
        vc.genotypes.push(genotype);
    }
    (vc, raw)
}

/// The `GenotypeSummaries` fixtures.
fn summaries_case(label: &str) -> VariantContext {
    let (gqs, no_calls): (Vec<i32>, usize) = match label {
        "no-genotypes" => (Vec::new(), 0),
        "one-gq" => (vec![50], 0),
        "two-gqs" => (vec![50, 70], 0),
        "three-gqs" => (vec![10, 50, 99], 0),
        "no-gq" => (Vec::new(), 2),
        "with-no-calls" => (vec![50, 70], 2),
        other => panic!("{other} has no fixture"),
    };
    let mut vc = VariantContext::new("chr1", START, vec![allele("A", true), allele("C", false)]);
    vc.stop = START;
    for (index, gq) in gqs.iter().enumerate() {
        let mut genotype = Genotype::new(
            &format!("s{index}"),
            vec![allele("A", true), allele("C", false)],
        );
        genotype.gq = Some(*gq);
        vc.genotypes.push(genotype);
    }
    for index in 0..no_calls {
        vc.genotypes.push(Genotype::new(
            &format!("n{index}"),
            vec![Allele::no_call(), Allele::no_call()],
        ));
    }
    vc
}

#[test]
fn every_depth_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("depth\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and a depth");
        let (vc, _) = qd_case(label);
        assert_eq!(
            qual_by_depth_depth(&vc, None).to_string(),
            expected,
            "depth on {label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no depth rows");
    println!("{count} depths identical");
}

#[test]
fn qd_matches_below_the_threshold_and_is_refused_at_and_above_it() {
    let text = golden();
    let mut compared = 0;
    let mut refused = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("qd\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and a value");
        let (vc, raw) = qd_case(label);
        match qual_by_depth(&vc, None, raw) {
            Ok(Some(value)) => {
                assert_eq!(format!("QD={value}"), expected, "QD on {label}");
                compared += 1;
            }
            Ok(None) => {
                assert_eq!(expected, "", "QD on {label} should be absent");
                compared += 1;
            }
            Err(QualByDepthError::RandomisedAboveThreshold { raw }) => {
                // The reference wrote *something* here, and it was a random draw. What the suite
                // can check is that the raw ratio really was at or past the threshold, and that
                // the value the reference produced is the jittered one rather than the ratio.
                assert!(raw >= 35.0, "{label} was refused with a raw QD of {raw}");
                assert!(expected.starts_with("QD="), "{label} produced {expected}");
                let written: f64 = expected
                    .trim_start_matches("QD=")
                    .parse()
                    .expect("a number");
                assert_ne!(
                    format!("{written:.2}"),
                    format!("{raw:.2}"),
                    "{label}: the reference wrote the raw ratio after all"
                );
                refused += 1;
            }
        }
    }
    assert!(compared > 0, "the golden carries no QD rows");
    assert_eq!(
        refused, 2,
        "the number of randomised rows changed; the threshold moved or a fixture did"
    );
    println!("{compared} QD values identical, {refused} refused as randomised");
}

#[test]
fn every_summary_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("summaries\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').unwrap_or((rest, ""));
        let vc = summaries_case(label);
        let ours = GenotypeSummaries
            .annotate(None, &vc, None)
            .iter()
            .map(|(key, value)| match value {
                AnnotationValue::Int(number) => {
                    format!("{key}={number}[java.lang.Integer]")
                }
                AnnotationValue::Str(text) => format!("{key}={text}[java.lang.String]"),
                other => panic!("{other:?} has no rendering"),
            })
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(ours, expected, "summaries on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no summary rows");
    println!("{count} genotype summaries identical");
}
