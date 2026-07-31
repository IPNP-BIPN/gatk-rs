//! Conformance for `MannWhitneyU`, against the oracle.
//!
//! Golden from `tools/genotyper-conformance/MannWhitneyDump.java`.
//!
//! The rows that carry the claim:
//!
//! ```text
//! mwu  all-tied        FIRST_DOMINATES  p = 0.5 exactly, by two mechanisms cancelling
//! mwu  nine-and-nine   FIRST_DOMINATES  the permutation test
//! mwu  ten-and-nine    FIRST_DOMINATES  the normal approximation, one element later
//! mwu  large-ramp                       where a float rank sum stops being exact
//! ```

use std::io::Read;

use gatk_engine::mann_whitney::{self, TestType};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mann_whitney.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn bits(value: f64) -> String {
    (value.to_bits() as i64).to_string()
}

fn side(label: &str) -> TestType {
    match label {
        "FIRST_DOMINATES" => TestType::FirstDominates,
        "SECOND_DOMINATES" => TestType::SecondDominates,
        "TWO_SIDED" => TestType::TwoSided,
        other => panic!("unknown side {other}"),
    }
}

fn ramp(from: i32, count: usize) -> Vec<f64> {
    (0..count).map(|i| (from + i as i32) as f64).collect()
}

fn constant(value: f64, count: usize) -> Vec<f64> {
    vec![value; count]
}

/// The two series each label was measured on.
fn series(label: &str) -> (Vec<f64>, Vec<f64>) {
    match label {
        "tiny-separated" => (vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]),
        "tiny-overlapping" => (vec![1.0, 3.0, 5.0], vec![2.0, 4.0, 6.0]),
        "tiny-reversed" => (vec![4.0, 5.0, 6.0], vec![1.0, 2.0, 3.0]),
        "one-each" => (vec![1.0], vec![2.0]),
        "one-vs-many" => (vec![1.0], vec![2.0, 3.0, 4.0, 5.0]),
        "empty-first" => (Vec::new(), vec![1.0, 2.0]),
        "empty-second" => (vec![1.0, 2.0], Vec::new()),
        "nine-and-nine" => (ramp(1, 9), ramp(2, 9)),
        "ten-and-nine" => (ramp(1, 10), ramp(2, 9)),
        "nine-and-ten" => (ramp(1, 9), ramp(2, 10)),
        "ten-and-ten" => (ramp(1, 10), ramp(2, 10)),
        "all-tied" => (constant(5.0, 12), constant(5.0, 12)),
        "all-tied-small" => (constant(5.0, 4), constant(5.0, 4)),
        "half-tied" => (vec![1.0, 1.0, 2.0, 2.0, 3.0], vec![1.0, 1.0, 2.0, 2.0, 3.0]),
        "tied-across" => (ramp(1, 12), ramp(1, 12)),
        "one-tie-band" => (vec![1.0, 2.0, 2.0, 2.0, 3.0], vec![4.0, 5.0, 6.0, 7.0, 8.0]),
        "large-ramp" => (ramp(1, 300), ramp(301, 300)),
        "large-interleaved" => (
            (0..600).map(|i| (2 * i) as f64).collect(),
            (0..600).map(|i| (2 * i + 1) as f64).collect(),
        ),
        "large-tied" => (constant(30.0, 200), constant(30.0, 200)),
        "qualities" => (
            vec![
                30.0, 30.0, 31.0, 32.0, 30.0, 29.0, 30.0, 30.0, 30.0, 31.0, 30.0, 30.0,
            ],
            vec![
                30.0, 28.0, 30.0, 27.0, 30.0, 30.0, 26.0, 30.0, 30.0, 30.0, 30.0, 25.0,
            ],
        ),
        "mapping-qualities" => (
            constant(60.0, 15),
            vec![
                60.0, 60.0, 60.0, 59.0, 60.0, 60.0, 60.0, 60.0, 57.0, 60.0, 60.0, 60.0, 60.0, 60.0,
                60.0,
            ],
        ),
        "negative" => (vec![-3.0, -2.0, -1.0], vec![-6.0, -5.0, -4.0]),
        "fractional" => (vec![0.5, 1.5, 2.5], vec![1.0, 2.0, 3.0]),
        other => panic!("{other} has no fixture"),
    }
}

#[test]
fn every_test_result_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("mwu\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let label = fields.next().expect("a label");
        let which = side(fields.next().expect("a side"));
        let expected: Vec<&str> = fields.collect();
        let (series1, series2) = series(label);
        let result = mann_whitney::test(&series1, &series2, which);
        assert_eq!(
            vec![
                bits(result.u),
                bits(result.z),
                bits(result.p),
                bits(result.median_shift)
            ],
            expected,
            "{label} {which:?}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no mwu rows");
    println!("{count} rank-sum results bit-identical");
}

#[test]
fn every_z_score_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("z\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let u = f64::from_bits(fields.next().expect("u").parse::<i64>().expect("bits") as u64);
        let n1: usize = fields.next().expect("n1").parse().expect("a number");
        let n2: usize = fields.next().expect("n2").parse().expect("a number");
        let nties =
            f64::from_bits(fields.next().expect("nties").parse::<i64>().expect("bits") as u64);
        let which = side(fields.next().expect("a side"));
        let expected = fields.next().expect("the z score");
        assert_eq!(
            bits(mann_whitney::calculate_z(u, n1, n2, nties, which)),
            expected,
            "calculateZ({u}, {n1}, {n2}, {nties}, {which:?})"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no z rows");
    println!("{count} Z scores bit-identical");
}

/// Everything tied is exactly one half, by two mechanisms cancelling rather than by one formula.
#[test]
fn everything_tied_lands_exactly_on_a_half() {
    let tied = vec![5.0; 12];
    let result = mann_whitney::test(&tied, &tied, TestType::FirstDominates);
    assert_eq!(result.p, 0.5);
}
