//! The eigensamples a panel of normals writes, against GATK 4.6.2.0 at `--spark-master local[1]`.
//!
//! Golden from `tools/readfilter-conformance/PanelRankDump.java`. The sibling suite,
//! `create-read-count-panel-of-normals`, carries the cap and leaves the count out; this one carries
//! the count and the values above the noise. Not their bits: 2 of 44 fresh `local[1]` runs of the
//! reference summed in another order, which moves the rounding and leaves the count
//! (docs/spark-parallel.md).
//!
//! # What this suite is for
//!
//!  * **`numEigensamples` lower than `min(requested, surviving)`**, on panels whose trailing values
//!    are rounding noise;
//!  * **the values above the noise, to eleven digits**, where the reference takes `dgesdd` and
//!    `dbdsqr`, up to 25 samples;
//!  * **the count alone** past that, where divide and conquer or ARPACK answers and the port does
//!    not take their path, on panels with no value near the threshold.

use gatk_corpus as corpus;
use gatk_tools::create_read_count_panel_of_normals::{
    number_of_eigensamples, preprocess, solver_mode, spark_gramian, standardize, truncated_svd,
    Arguments, Matrix, SolverMode,
};
use std::collections::HashMap;

fn golden() -> Option<String> {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/panel_rank.txt.gz");
    path.exists().then(|| corpus::read_golden(&path))
}

/// The dump's `noise`.
fn noise(sample: usize, interval: usize) -> i64 {
    let i = interval as i64 + 7;
    200 + ((sample as i64 + 1) * i * i * 2654435761 + sample as i64 * 97) % 1009
}

fn independent(samples: usize, intervals: usize) -> Matrix {
    let rows: Vec<Vec<f64>> = (0..samples)
        .map(|s| (0..intervals).map(|i| noise(s, i) as f64).collect())
        .collect();
    Matrix::new(&rows)
}

fn families(samples: usize, families: usize, intervals: usize) -> Matrix {
    let rows: Vec<Vec<f64>> = (0..samples)
        .map(|s| {
            (0..intervals)
                .map(|i| (noise(s % families, i) * (1 + (s / families) as i64)) as f64)
                .collect()
        })
        .collect();
    Matrix::new(&rows)
}

/// The sibling dump's `shape`.
fn shape(sample: usize, i: usize) -> usize {
    let base = 100 + (i * 7) % 23;
    if sample == 7 {
        return if i < 4 { base } else { 0 };
    }
    if sample == 8 {
        return if i < 10 { base * 40 } else { base };
    }
    let weight = 1 + if sample == 6 { 3 } else { sample };
    if i < 20 {
        base * weight
    } else {
        base
    }
}

/// The sibling dump's `sampleCounts`, for the nine samples.
fn nine() -> Matrix {
    let rows: Vec<Vec<f64>> = (0..9)
        .map(|sample| {
            (0..40)
                .map(|i| {
                    let mut value = shape(sample, i) * if sample == 6 { 2 } else { 1 };
                    if i < 2 {
                        value = match sample {
                            7 => 0,
                            6 => 2,
                            _ => 1,
                        };
                    }
                    if i == 20 && sample < 3 {
                        value = 0;
                    }
                    value as f64
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

struct Run {
    label: &'static str,
    counts: Matrix,
    arguments: Arguments,
}

fn run(label: &'static str, counts: Matrix, arguments: Arguments) -> Run {
    Run {
        label,
        counts,
        arguments,
    }
}

fn runs() -> Vec<Run> {
    let requested = |k| Arguments {
        number_of_eigensamples: k,
        ..Arguments::default()
    };
    vec![
        run("nine-defaults", nine(), Arguments::default()),
        run("nine-no-filters", nine(), no_filters()),
        run("nine-k4", nine(), requested(4)),
        run("nine-k5", nine(), requested(5)),
        run("nine-k6", nine(), requested(6)),
        run("nine-k7", nine(), requested(7)),
        run("nine-k8", nine(), requested(8)),
        run("nine-k100", nine(), requested(100)),
        run("nine-k2-lapack", nine(), requested(2)),
        run("nine-k1-arpack", nine(), requested(1)),
        run("families-6-of-3", families(6, 3, 30), no_filters()),
        run("families-12-of-3", families(12, 3, 30), no_filters()),
        run("families-12-of-5", families(12, 5, 60), no_filters()),
        run("families-25-of-4", families(25, 4, 40), no_filters()),
        run("families-8-of-3-wide", families(8, 3, 250), no_filters()),
        run("independent-5-wide", independent(5, 250), no_filters()),
        run("independent-25", independent(25, 40), no_filters()),
        run(
            "independent-30",
            independent(30, 50),
            Arguments {
                number_of_eigensamples: 30,
                ..no_filters()
            },
        ),
        run("independent-70-arpack", independent(70, 90), no_filters()),
    ]
}

/// What the port answers: the request it makes, the panel's samples, the singular values.
fn answer(run: &Run) -> (usize, usize, Vec<f64>) {
    let mut values = preprocess(&run.counts, &run.arguments).values;
    standardize(&mut values).expect("positive medians");
    let surviving = values.samples;
    let k = number_of_eigensamples(run.arguments.number_of_eigensamples, surviving);
    (
        k,
        surviving,
        truncated_svd(&values, k)
            .expect("a decomposition")
            .singular_values,
    )
}

fn exact(n: usize, k: usize) -> bool {
    n <= 25 && solver_mode(n, k) == SolverMode::LocalLapack
}

/// Every row of the dump, when the first x86-64 CI run has published it.
#[test]
fn every_panel_writes_the_eigensamples_the_reference_writes() {
    let Some(text) = golden() else {
        // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
        return;
    };
    let mut eigensamples = HashMap::new();
    let mut signal = HashMap::new();
    for line in text.lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        match cells[0] {
            "rank" => {
                let count = cells[4].strip_prefix("eigensamples=").expect("a count");
                let count: usize = count.parse().expect("a number");
                let surviving = cells[3].strip_prefix("surviving=").expect("a count");
                let surviving: usize = surviving.parse().expect("a number");
                eigensamples.insert(cells[1].to_string(), (count, surviving));
            }
            "values" => {
                let values: Vec<f64> = cells[2]
                    .split(',')
                    .map(|value| value.parse().expect("a value"))
                    .collect();
                signal.insert(cells[1].to_string(), values);
            }
            _ => {}
        }
    }
    let mut lower = 0;
    for run in runs() {
        let (k, surviving, values) = answer(&run);
        let (count, expected_surviving) = eigensamples[run.label];
        assert_eq!(surviving, expected_surviving, "{}", run.label);
        assert_eq!(values.len(), count, "{}", run.label);
        if count < k {
            lower += 1;
        }
        if exact(surviving, k) {
            let expected = &signal[run.label];
            let produced: Vec<f64> = values
                .iter()
                .copied()
                .filter(|value| *value >= 1e-6 * values[0])
                .collect();
            assert_eq!(produced.len(), expected.len(), "{}", run.label);
            for (value, reference) in produced.iter().zip(expected) {
                assert!(
                    (value - reference).abs() <= 1e-9 * reference,
                    "{}: {value} {reference}",
                    run.label
                );
            }
        }
    }
    assert!(lower > 0, "the suite exists for a count under the request");
}

/// `auto`'s three answers, the integer divisions included.
#[test]
fn the_mode_follows_n_and_k_as_auto_does() {
    assert_eq!(solver_mode(8, 2), SolverMode::LocalLapack, "2 < 8 / 3 is 2");
    assert_eq!(solver_mode(9, 2), SolverMode::LocalArpack);
    assert_eq!(solver_mode(8, 8), SolverMode::LocalLapack);
    assert_eq!(solver_mode(62, 20), SolverMode::LocalLapack, "62 / 3 is 20");
    assert_eq!(solver_mode(63, 20), SolverMode::LocalArpack);
    assert_eq!(solver_mode(100, 20), SolverMode::DistributedArpack);
    assert_eq!(solver_mode(100, 51), SolverMode::LocalLapack);
    assert_eq!(solver_mode(100, 50), SolverMode::DistributedArpack);
}

/// The Gramian is grouped by slice, so more rows than slices is a different sum from row after row.
#[test]
fn the_gramian_is_the_product_it_says_it_is() {
    let counts = independent(4, 250);
    let gramian = spark_gramian(&counts);
    for i in 0..4 {
        for j in 0..4 {
            let direct: f64 = (0..250).map(|r| counts.get(i, r) * counts.get(j, r)).sum();
            assert!(
                (gramian[i + j * 4] - direct).abs() < direct * 1e-12,
                "{i} {j}"
            );
            assert_eq!(gramian[i + j * 4], gramian[j + i * 4], "symmetric");
        }
    }
}

/// A panel with as many independent samples as it asks for gets what it asks for.
#[test]
fn a_full_rank_panel_gets_its_request() {
    for (samples, label) in [(5, "under"), (25, "at the limit"), (30, "over")] {
        let run = run("full", independent(samples, 50), no_filters());
        let (k, surviving, values) = answer(&run);
        assert_eq!(k, surviving.min(20), "{label}");
        assert_eq!(values.len(), k, "{label}");
    }
}
