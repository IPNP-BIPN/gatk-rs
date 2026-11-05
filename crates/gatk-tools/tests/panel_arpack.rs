//! Spark's ARPACK path over a panel's Gramian, against GATK 4.6.2.0 at `--spark-master local[1]`.
//!
//! Golden from `tools/readfilter-conformance/PanelArpackDump.java`: the raw bits of the eigenvalues
//! `symmetricEigs` returns, LocalARPACK under a hundred samples and DistARPACK from a hundred, a
//! fingerprint of their vectors, and the rank `computeSVD` takes from them; then the tool itself on
//! three panels, its count and the values above the noise as `PanelRankDump` writes them.
//!
//! The bits are the usual `local[1]` run's: a task result that reaches the driver late sums the
//! Gramian, or one product, in another order (docs/spark-parallel.md).

use gatk_corpus as corpus;
use gatk_tools::create_read_count_panel_of_normals::{
    breeze_gramian_times, number_of_eigensamples, preprocess, spark_gramian, spark_gramian_times,
    spark_rows, standardize, truncated_svd, Arguments, Matrix,
};
use gatk_tools::f2j_arpack::symmetric_eigs;

fn golden() -> Option<String> {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/panel_arpack.txt.gz");
    path.exists().then(|| corpus::read_golden(&path))
}

/// `PanelRankDump.noise`.
fn noise(sample: usize, interval: usize) -> i64 {
    let i = interval as i64 + 7;
    200 + ((sample as i64 + 1) * i * i * 2654435761 + sample as i64 * 97) % 1009
}

/// The dump's `value`: `families` centred profiles, each sample a multiple of its family's.
fn panel(n: usize, intervals: usize, families: usize) -> Matrix {
    let rows: Vec<Vec<f64>> = (0..n)
        .map(|s| {
            (0..intervals)
                .map(|i| {
                    let multiple = 1 + ((s / families) % 3) as i64;
                    ((noise(s % families, i) - 704) * multiple) as f64 / 128.0
                })
                .collect()
        })
        .collect();
    Matrix::new(&rows)
}

fn counts(samples: usize, families: usize, intervals: usize, independent: bool) -> Matrix {
    let rows: Vec<Vec<f64>> = (0..samples)
        .map(|s| {
            (0..intervals)
                .map(|i| {
                    if independent {
                        noise(s, i) as f64
                    } else {
                        (noise(s % families, i) * (1 + (s / families) as i64)) as f64
                    }
                })
                .collect()
        })
        .collect();
    Matrix::new(&rows)
}

fn no_filters() -> Arguments {
    Arguments {
        minimum_interval_median_percentile: 0.0,
        maximum_zeros_in_sample_percentage: 100.0,
        maximum_zeros_in_interval_percentage: 100.0,
        extreme_sample_median_percentile: 0.0,
        extreme_outlier_truncation_percentile: 0.0,
        ..Arguments::default()
    }
}

/// The dump's cases: label, samples, intervals, families.
const CASES: [(&str, usize, usize, usize); 11] = [
    ("local-63-full", 63, 300, 63),
    ("local-70-of-40", 70, 300, 40),
    ("local-81-of-25", 81, 250, 25),
    ("local-99-full-wide", 99, 1000, 99),
    ("local-90-short", 90, 60, 90),
    ("local-70-of-12", 70, 300, 12),
    ("distributed-100-full", 100, 400, 100),
    ("distributed-120-of-60", 120, 400, 60),
    ("distributed-150-of-30", 150, 333, 30),
    ("distributed-110-short", 110, 90, 110),
    ("distributed-120-of-15", 120, 300, 15),
];

fn hex(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

/// What the port answers for a case: the dump's row, or `None` where Spark throws.
fn eigen(n: usize, intervals: usize, families: usize, k: usize) -> Option<(usize, String, String)> {
    let matrix = panel(n, intervals, families);
    let (values, vectors) = if n < 100 {
        let gramian = spark_gramian(&matrix);
        symmetric_eigs(n, k, &mut |x: &[f64], y: &mut [f64]| {
            breeze_gramian_times(&gramian, x, y)
        })?
    } else {
        let rows = spark_rows(&matrix);
        symmetric_eigs(n, k, &mut |x: &[f64], y: &mut [f64]| {
            spark_gramian_times(&rows, x, y)
        })?
    };
    let threshold = 1e-9 * values[0].sqrt();
    let rank = values
        .iter()
        .take(k)
        .take_while(|value| value.sqrt() >= threshold)
        .count();
    let bits: Vec<String> = values.iter().map(|value| hex(*value)).collect();
    let prints: Vec<String> = (0..values.len())
        .map(|c| {
            let column = &vectors[c * n..(c + 1) * n];
            let sum = column.iter().fold(0.0, |sum, value| sum + value);
            format!("{}:{}", hex(column[0]), hex(sum))
        })
        .collect();
    Some((rank, bits.join(","), prints.join(",")))
}

/// Every row of the dump, when the first x86-64 CI run has published it.
#[test]
fn arpack_gives_the_reference_eigenvalues_to_the_bit() {
    let Some(text) = golden() else {
        // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
        return;
    };
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    let find = |kind: &str, label: &str| {
        rows.iter()
            .find(|cells| cells[0] == kind && cells[1] == label)
            .cloned()
    };
    for (label, n, intervals, families) in CASES {
        let port = eigen(n, intervals, families, 20);
        if let Some(error) = find("error", label) {
            assert!(port.is_none(), "{label}: the reference threw {}", error[2]);
            continue;
        }
        let row = find("eigen", label).unwrap_or_else(|| panic!("{label} is in the dump"));
        let vectors = find("vectors", label).expect("a vectors row");
        let (rank, bits, prints) = port.unwrap_or_else(|| panic!("{label}: the port threw"));
        assert_eq!(format!("rank={rank}"), row[5], "{label}");
        assert_eq!(bits, row[6], "{label}: eigenvalues");
        assert_eq!(prints, vectors[2], "{label}: eigenvectors");
    }
    for (label, samples, families, intervals, independent) in [
        ("tool-independent-70", 70, 70, 120, true),
        ("tool-families-70-of-30", 70, 30, 120, false),
        ("tool-independent-120", 120, 120, 150, true),
    ] {
        let mut values = preprocess(
            &counts(samples, families, intervals, independent),
            &no_filters(),
        )
        .values;
        standardize(&mut values).expect("positive medians");
        let k = number_of_eigensamples(20, values.samples);
        let produced = truncated_svd(&values, k).singular_values;
        let row = find("rank", label).unwrap_or_else(|| panic!("{label} is in the dump"));
        assert_eq!(
            row[4],
            format!("eigensamples={}", produced.len()),
            "{label}"
        );
        let expected: Vec<f64> = find("values", label).expect("a values row")[2]
            .split(',')
            .map(|value| value.parse().expect("a value"))
            .collect();
        let signal: Vec<f64> = produced
            .iter()
            .copied()
            .filter(|value| *value >= 1e-6 * produced[0])
            .collect();
        assert_eq!(signal.len(), expected.len(), "{label}");
        for (value, reference) in signal.iter().zip(&expected) {
            assert!(
                (value - reference).abs() <= 1e-9 * reference,
                "{label}: {value} {reference}"
            );
        }
    }
}
