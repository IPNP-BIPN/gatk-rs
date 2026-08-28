//! Conformance for `AlleleFrequencyQC` against GATK 4.6.2.0, compared as the statistic and the
//! p-value it writes, recomputed from the report it read them out of.
//!
//! Golden from `tools/readfilter-conformance/AlleleFrequencyQCDump.java`, which keeps both files:
//! the metrics the tool wrote and the `--debug-file` report the numbers came from. The test takes
//! the report, groups it the way the tool does, and checks that the port arrives at the metrics.
//!
//! # What this suite is for
//!
//!  * **the statistic being a constant variance and not an expected count**;
//!  * **that variance being squared**;
//!  * **the p-value being the upper tail of a chi-squared with the bin count less one**;
//!  * **the bin ladder being fixed rather than the data's**;
//!  * **a comparison site with no call against it contributing its own square**;
//!  * **the rows being cut to `called` before the grouping**;
//!  * **a file of one variant being an ordinary case**;
//!  * **and the sample name coming from a header line rather than a genotype column.**

use gatk_corpus as corpus;
use gatk_tools::allele_frequency_qc::{
    chi_squared_statistic, complains, degrees_of_freedom, metrics, p_value, CALLED,
    DEFAULT_ALLOWED_VARIANCE, DEFAULT_THRESHOLD, METRIC_TYPE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/allele_frequency_qc.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// The metrics row a case wrote, as its two numbers.
fn written(text: &str, label: &str) -> (f64, f64) {
    let table = field(text, "metrics", label).unwrap_or_else(|| panic!("metrics/{label}"));
    let mut lines = table.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let row: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    let column = |name: &str| {
        row[header.iter().position(|h| *h == name).expect(name)]
            .parse::<f64>()
            .expect("a number")
    };
    assert_eq!(
        row[header
            .iter()
            .position(|h| *h == "METRIC_TYPE")
            .expect("a type")],
        METRIC_TYPE
    );
    (column("CHI_SQ_VALUE"), column("METRIC_VALUE"))
}

/// The bins the tool grouped, taken from the report it read: the `called` rows, by allele
/// frequency, each holding one average per eval track in the report's own order.
fn bins(text: &str, label: &str) -> Vec<Vec<f64>> {
    let report = field(text, "report", label).unwrap_or_else(|| panic!("report/{label}"));
    let mut lines = report.lines().filter(|line| !line.trim().is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split_whitespace().collect();
    let index = |name: &str| header.iter().position(|h| *h == name).expect(name);
    let (frequency, filter, average) =
        (index("AlleleFrequency"), index("Filter"), index("avgVarAF"));
    let mut order: Vec<String> = Vec::new();
    let mut grouped: std::collections::HashMap<String, Vec<f64>> = std::collections::HashMap::new();
    for line in lines {
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns[filter] != CALLED {
            continue;
        }
        let key = columns[frequency].to_string();
        if !grouped.contains_key(&key) {
            order.push(key.clone());
        }
        grouped
            .entry(key)
            .or_default()
            .push(columns[average].parse().expect("an average"));
    }
    order
        .into_iter()
        .map(|key| grouped.remove(&key).expect("a bin"))
        .collect()
}

/// The port arrives at both numbers the tool wrote, from the report the tool read.
#[test]
fn the_port_recomputes_what_was_written() {
    let text = golden();
    for (label, variance) in [
        ("het-calls", DEFAULT_ALLOWED_VARIANCE),
        ("one-bin-homozygous", DEFAULT_ALLOWED_VARIANCE),
        ("variance-tenth", 0.1),
        ("variance-hundredth", 0.001),
        ("a-filtered-variant", DEFAULT_ALLOWED_VARIANCE),
        ("a-bin-with-one-entry", DEFAULT_ALLOWED_VARIANCE),
        ("one-variant", DEFAULT_ALLOWED_VARIANCE),
    ] {
        let (their_statistic, their_p) = written(&text, label);
        let (our_statistic, our_p) = metrics(&bins(&text, label), variance);
        assert!(
            (our_statistic - their_statistic).abs() < 1e-6,
            "{label}: {our_statistic} against {their_statistic}"
        );
        // The p-value is written with six fraction digits.
        assert!(
            (our_p - their_p).abs() < 1e-6,
            "{label}: {our_p} against {their_p}"
        );
    }
}

/// The variance is squared, so a tenfold variance divides the statistic by a hundred.
#[test]
fn the_variance_is_squared() {
    let text = golden();
    let (hundredth, _) = written(&text, "one-bin-homozygous");
    let (tenth, _) = written(&text, "variance-tenth");
    let (thousandth, _) = written(&text, "variance-hundredth");
    assert!(
        (hundredth / tenth - 100.0).abs() < 1e-6,
        "{hundredth} against {tenth}"
    );
    assert!((thousandth / hundredth - 100.0).abs() < 1e-6);
    // Which is the port's own division, on one bin of a known difference.
    let one = vec![vec![0.5, 0.0]];
    assert_eq!(chi_squared_statistic(&one, 0.01), 0.25 / 0.0001);
    assert!(
        (chi_squared_statistic(&one, 0.1) * 100.0 - chi_squared_statistic(&one, 0.01)).abs() < 1e-9
    );
    // And a bin of one entry contributes nothing, which is the guard the reference carries.
    assert_eq!(chi_squared_statistic(&[vec![0.5]], 0.01), 0.0);
}

/// The bin ladder is fixed, not the data's, so the degrees of freedom do not move.
#[test]
fn the_bin_ladder_is_fixed() {
    let text = golden();
    let wide = bins(&text, "het-calls").len();
    assert_eq!(wide, 61);
    for label in ["one-variant", "a-filtered-variant", "a-bin-with-one-entry"] {
        assert_eq!(bins(&text, label).len(), wide, "{label}");
    }
    assert_eq!(degrees_of_freedom(wide), 60.0);
    // Every bin holds exactly the two rows, one per eval track, so the shorter-list guard never
    // fires against the reference.
    assert!(bins(&text, "het-calls").iter().all(|bin| bin.len() == 2));
}

/// A comparison site the call set has nothing at contributes its own square, not a degree of
/// freedom.
#[test]
fn an_uncalled_comparison_site_contributes_a_square() {
    let text = golden();
    let (without, _) = written(&text, "one-bin-homozygous");
    let (with, _) = written(&text, "a-bin-with-one-entry");
    assert_eq!(
        bins(&text, "a-bin-with-one-entry").len(),
        bins(&text, "one-bin-homozygous").len()
    );
    // The comparison frequency was 0.001, against a nought on the call set's side.
    assert!((with - without - 0.001_f64.powi(2) / 0.01_f64.powi(2)).abs() < 1e-6);
}

/// The rows are cut to `called` before the grouping, so a filtered variant changes nothing.
#[test]
fn a_filtered_variant_changes_nothing() {
    let text = golden();
    assert_eq!(
        written(&text, "a-filtered-variant"),
        written(&text, "one-bin-homozygous")
    );
    let report = field(&text, "report", "a-filtered-variant").expect("a report");
    // The report does carry the filtered rows: they are dropped by the tool and not by the walk.
    assert!(report.contains("filtered"));
}

/// A file with one variant in it is an ordinary case: nought against a p-value of one.
#[test]
fn one_variant_is_not_a_degenerate_case() {
    let text = golden();
    let (statistic, p) = written(&text, "one-variant");
    assert_eq!(statistic, 0.0);
    assert_eq!(p, 1.0);
    assert_eq!(p_value(0.0, 60.0), 1.0);
    // The p-value is the upper tail, so a large statistic answers nought.
    assert!(p_value(7600.25, 60.0) < 1e-12);
    assert!(p_value(76.0025, 60.0) > 0.07 && p_value(76.0025, 60.0) < 0.08);
    // And the threshold decides only whether the run complains.
    assert!(complains(0.0, DEFAULT_THRESHOLD));
    assert!(!complains(1.0, DEFAULT_THRESHOLD));
    assert!(!complains(0.079611, DEFAULT_THRESHOLD));
    assert_eq!(
        written(&text, "threshold-above-the-pvalue"),
        written(&text, "het-calls")
    );
    assert!(complains(
        written(&text, "threshold-above-the-pvalue").1,
        0.99
    ));
}

/// The sample name comes from a header line, so a VCF without one dies on a null.
#[test]
fn the_sample_name_comes_from_the_header() {
    let text = golden();
    let refusal = field(&text, "error", "no-sample-alias").expect("a refusal");
    assert!(
        refusal.starts_with("java.lang.NullPointerException"),
        "{refusal}"
    );
    assert!(refusal.contains("getOtherHeaderLine"), "{refusal}");
    // The runs that do carry one report it, and it is not a genotype column's name.
    let table = field(&text, "metrics", "het-calls").expect("the metrics");
    assert!(table.contains("NA12878"));
    assert!(!table.contains("\ts1\t"));
}
