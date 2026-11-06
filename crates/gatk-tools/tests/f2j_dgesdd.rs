//! F2J's `dgesdd('A')` on Gramians from 2 to 128 samples, against the oracle, bit for bit.
//!
//! Golden from `tools/readfilter-conformance/DgesddDump.java`, which calls `org.netlib.lapack.Dgesdd`
//! as Breeze does for Spark's LocalLAPACK `computeSVD`. Every row is a matrix built here from the
//! same integers with the same operations, so the input bits are the reference's; the row holds
//! Spark's rank under `rCond` 1e-9, a hash of the singular values' bits, a hash of U's, and the
//! first and last values. Past 25 samples that is `dbdsdc`'s divide and conquer, past 32 the
//! blocked `dormqr`.

use gatk_corpus as corpus;
use gatk_tools::f2j_gesdd::{dgesdd, reachable};

const PRIME: u64 = 0x100000001b3;

fn golden() -> Option<String> {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/f2j_dgesdd.txt.gz");
    path.exists().then(|| corpus::read_golden(&path))
}

/// The dump's `noise`.
fn noise(sample: usize, interval: usize) -> i64 {
    let i = interval as i64 + 7;
    200 + ((sample as i64 + 1) * i * i * 2654435761 + sample as i64 * 97) % 1009
}

/// The dump's `panel`: `2n + 3` rows of `n` samples.
fn panel(kind: &str, n: usize) -> Vec<Vec<f64>> {
    let families = (n / 5).max(2);
    (0..2 * n + 3)
        .map(|r| {
            let mut row: Vec<f64> = (0..n)
                .map(|s| match kind {
                    "families" => {
                        noise(s % families, r) as f64 / 1009.0 * (1 + s / families) as f64
                    }
                    "flat" if s == 0 => 0.0,
                    "flat" => noise(if s == n - 1 { 1 } else { s }, r) as f64 / 1009.0,
                    _ => noise(s, r) as f64 / 1009.0,
                })
                .collect();
            if kind == "centered" || kind == "families" {
                let mut sum = 0.0;
                for value in &row {
                    sum += value;
                }
                let mean = sum / n as f64;
                for value in &mut row {
                    *value -= mean;
                }
            }
            row
        })
        .collect()
}

/// The dump's `gramian`, column major.
fn gramian(x: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut g = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let mut sum = 0.0;
            for row in x {
                sum += row[i] * row[j];
            }
            g[i + j * n] = sum;
        }
    }
    g
}

fn fold(values: &[f64]) -> u64 {
    values
        .iter()
        .fold(0, |h, value| (h ^ value.to_bits()).wrapping_mul(PRIME))
}

/// What the port writes for one matrix, in the dump's format.
fn row(kind: &str, n: usize) -> String {
    let mut a = gramian(&panel(kind, n), n);
    assert!(reachable(n, &a), "{kind} {n}");
    let (s, u) = dgesdd(n, &mut a).expect("it converges");
    let threshold = 1e-9 * s[0].sqrt();
    let rank = s
        .iter()
        .take_while(|value| value.sqrt() >= threshold)
        .count();
    format!(
        "svd\t{kind}\t{n}\tinfo=0\trank={rank}\ts={:016x}\tu={:016x}\tfirst={:016x}\tlast={:016x}",
        fold(&s),
        fold(&u),
        s[0].to_bits(),
        s[n - 1].to_bits()
    )
}

#[test]
fn every_matrix_has_the_bits_the_reference_has() {
    let Some(text) = golden() else {
        // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
        return;
    };
    let mut rows = 0;
    let mut differ = Vec::new();
    for line in text.lines().filter(|line| line.starts_with("svd\t")) {
        let cells: Vec<&str> = line.split('\t').collect();
        let n: usize = cells[2].parse().expect("a size");
        let produced = row(cells[1], n);
        if produced != line {
            differ.push(format!("reference {line}\nport      {produced}"));
        }
        rows += 1;
    }
    assert_eq!(rows, 300);
    assert!(
        differ.is_empty(),
        "{} of {rows} rows differ:\n{}",
        differ.len(),
        differ.join("\n")
    );
}

/// Past 25 samples the reference raises every diagonal entry under the epsilon to it, so even a
/// Gramian with an exact zero row keeps all its values under `rCond`; up to 25 it does not.
#[test]
fn divide_and_conquer_keeps_every_value_a_zero_sample_would_drop() {
    let rank = |n| {
        let mut a = gramian(&panel("flat", n), n);
        let (s, _) = dgesdd(n, &mut a).expect("it converges");
        let threshold = 1e-9 * s[0].sqrt();
        s.iter()
            .take_while(|value| value.sqrt() >= threshold)
            .count()
    };
    assert!(rank(25) < 25);
    assert_eq!(rank(26), 26);
    assert_eq!(rank(60), 60);
}
