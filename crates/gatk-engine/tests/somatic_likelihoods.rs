//! Conformance for the Dirichlet fixed point against the oracle.
//!
//! Golden from `tools/annotation-conformance/SomaticLikelihoodsDump.java`, as raw bit patterns.
//!
//! # What this suite is for
//!
//! G1.9.1 established that a 1-ulp `exp` is what enters this iteration. The open question — the one
//! #96 cannot be settled without — is what the iteration does with it, and there are two ways it
//! could go wrong that look identical unless they are separated:
//!
//!  * **amplification**, where the difference grows across iterations;
//!  * **a different iteration count**, because convergence is a *threshold* test
//!    (`distance1 / sum < 0.001`) and a difference too small to see in the values can still land on
//!    the other side of it.
//!
//! The second would be a far larger divergence than the first, so the golden carries the count and
//! this test asserts it **exactly**. A row where the counts match and the values differ is
//! amplification; a row where the counts differ is the threshold; and the test says which.
//!
//! # The weights rows have no licence exposure at all
//!
//! `effectiveLogMultinomialWeights` is `digamma` arithmetic — no `exp` anywhere — so those rows
//! must be **bit-identical**, and are asserted as such. They are the control: if they drift, the
//! problem is not the exponential.

use std::io::Read;

use gatk_engine::somatic_likelihoods::{
    allele_fractions_posterior, effective_log_multinomial_weights,
};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/somatic_likelihoods.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn from_bits(hex: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(hex, 16).expect("hex bits"))
}

fn parse_list(text: &str) -> Vec<f64> {
    text.split(',').map(from_bits).collect()
}

fn ulps_apart(ours: f64, theirs: f64) -> i64 {
    if ours.to_bits() == theirs.to_bits() {
        return 0;
    }
    if !ours.is_finite() || !theirs.is_finite() {
        return i64::MAX;
    }
    (ours.to_bits() as i64 - theirs.to_bits() as i64).abs()
}

/// `column(reads, per_allele)` from the dump: the same likelihood column repeated.
fn column(reads: usize, per_allele: &[f64]) -> Vec<Vec<f64>> {
    per_allele.iter().map(|v| vec![*v; reads]).collect()
}

/// `split(first, second)`: reads favouring allele 0, then reads favouring allele 1.
fn split(first: usize, second: usize) -> Vec<Vec<f64>> {
    let total = first + second;
    let mut matrix = vec![vec![0.0; total]; 2];
    for read in 0..total {
        let first_allele = read < first;
        matrix[0][read] = if first_allele { -0.01 } else { -4.0 };
        matrix[1][read] = if first_allele { -4.0 } else { -0.01 };
    }
    matrix
}

fn decay(n: usize) -> Vec<f64> {
    (0..n).map(|i| 1.0 / (1.0 + i as f64)).collect()
}

/// The case each label was run over, transcribed from the dump rather than parsed out of the
/// golden: the inputs are a *configuration*, and reading them back from the results would be
/// deriving the question from the answer.
#[allow(clippy::type_complexity)]
fn case(label: &str) -> (Vec<Vec<f64>>, Vec<f64>, Option<Vec<f64>>) {
    let flat2 = vec![1.0, 1.0];
    let flat3 = vec![1.0, 1.0, 1.0];
    match label {
        "one-read-clean" => (vec![vec![-0.001], vec![-10.0]], flat2, None),
        "one-read-flat" => (vec![vec![-1.0], vec![-1.0]], flat2, None),
        "ten-reads-clean" => (column(10, &[-0.001, -10.0]), flat2, None),
        "ten-reads-split" => (split(5, 5), flat2, None),
        "near-tie" => (split(50, 49), flat2, None),
        "three-alleles" => (
            vec![
                vec![-0.1, -0.2, -5.0],
                vec![-3.0, -0.3, -4.0],
                vec![-6.0, -7.0, -0.05],
            ],
            flat3,
            None,
        ),
        "skewed-prior" => (split(5, 5), vec![0.5, 3.0], None),
        "large-prior" => (split(5, 5), vec![100.0, 100.0], None),
        "weighted-uniform" => (split(5, 5), flat2, Some(vec![1.0; 10])),
        "weighted-decaying" => (split(5, 5), flat2, Some(decay(10))),
        "saturated" => (column(4, &[-1e-9, -700.0]), flat2, None),
        "fifty-reads" => (split(25, 25), flat2, None),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn weights_case(label: &str) -> Vec<f64> {
    match label {
        "weights-flat2" => vec![1.0, 1.0],
        "weights-flat3" => vec![1.0, 1.0, 1.0],
        "weights-skewed" => vec![0.5, 3.0],
        "weights-large" => vec![100.0, 100.0],
        "weights-tiny" => vec![1e-3, 1e-3],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn the_fixed_point_lands_where_the_reference_lands() {
    let text = golden();
    let mut rows = 0usize;
    let mut worst = 0i64;
    let mut worst_row = String::new();
    let mut count_mismatches: Vec<String> = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            // The harness replays the loop to observe the iteration count, and checks its replay
            // against the engine's own answer. A `false` here would mean the golden is measuring
            // something other than the engine.
            "agree" => assert_eq!(
                fields[2], "true",
                "{}: the dump's replay disagreed with the engine it is measuring",
                fields[1]
            ),
            "post" => {
                let label = fields[1];
                let expected_iterations: usize = fields[2].parse().expect("an iteration count");
                let expected = parse_list(fields[3]);
                let (likelihoods, prior, weights) = case(label);

                let ours = allele_fractions_posterior(&likelihoods, &prior, weights.as_deref())
                    .unwrap_or_else(|_| panic!("{label} must converge"));

                // Asserted exactly, and first: a different count is a different amount of work,
                // and every value comparison below would be measuring that instead.
                if ours.iterations != expected_iterations {
                    count_mismatches.push(format!(
                        "{label}: {} iterations against {expected_iterations}",
                        ours.iterations
                    ));
                    continue;
                }

                assert_eq!(ours.values.len(), expected.len(), "{label}: length");
                for (index, (ours, theirs)) in ours.values.iter().zip(&expected).enumerate() {
                    let distance = ulps_apart(*ours, *theirs);
                    assert!(
                        distance <= 1,
                        "{label}[{index}]: {distance} ulp apart — the iteration amplified the \
                         exp difference, which is what G1.9 needed to know"
                    );
                    if distance > worst {
                        worst = distance;
                        worst_row = format!("{label}[{index}]");
                    }
                    rows += 1;
                }
            }
            // digamma only: no exp, so no excuse.
            "weights" => {
                let label = fields[1];
                let expected = parse_list(fields[2]);
                let ours = effective_log_multinomial_weights(&weights_case(label))
                    .unwrap_or_else(|| panic!("{label} must produce weights"));
                for (index, (ours, theirs)) in ours.iter().zip(&expected).enumerate() {
                    assert_eq!(
                        ours.to_bits(),
                        theirs.to_bits(),
                        "{label}[{index}] is digamma arithmetic with no exp in it, so it must be \
                         bit-identical"
                    );
                    rows += 1;
                }
            }
            _ => {}
        }
    }

    assert!(
        count_mismatches.is_empty(),
        "the convergence threshold was crossed differently, which is a larger divergence than any \
         number of ulps: {count_mismatches:?}"
    );
    assert!(rows > 30, "the golden shrank to {rows} values");
    println!(
        "Dirichlet fixed point: {rows} values compared, every iteration count matched, worst \
         divergence {worst} ulp{}",
        if worst == 0 {
            String::new()
        } else {
            format!(" at {worst_row}")
        }
    );
}

/// The loop tests after the step, so it always runs at least once — even for an input already at
/// its fixed point.
#[test]
fn the_loop_body_always_runs_once() {
    let result = allele_fractions_posterior(&[vec![-1.0], vec![-1.0]], &[1.0, 1.0], None)
        .expect("converges");
    assert!(
        result.iterations >= 1,
        "`while (!converged)` with `converged` false on entry is a do-while"
    );
}
