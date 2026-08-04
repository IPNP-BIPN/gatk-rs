//! Conformance for `AllelePseudoDepth` against the oracle.
//!
//! Golden from `tools/annotation-conformance/AllelePseudoDepthDump.java`, comparing the **strings**
//! the reference writes into the genotype rather than the doubles behind them. That is not a
//! weakening of the test, it is the observable: `DD` and `DF` are built with
//! `DecimalFormat.format` and joined with commas, and nothing downstream ever sees the doubles.
//!
//! It is also what makes the suite possible. The chain underneath calls `Math.exp`, which has no
//! exact port (htsjdk-rs #71), so a suite on doubles would be asserting something already measured
//! to be false in general. What it is instead is bounded twice over: G1.9.1 put the exponential
//! within 1 ulp and G1.9.2 measured the fixed point at zero, and the formatter is about twelve
//! orders of magnitude coarser than either.
//!
//! The golden carries the inputs as raw bit patterns alongside the outputs, so this file feeds the
//! port exactly the doubles the reference saw. A decimal round trip through the golden would be a
//! second source of divergence in a suite whose whole subject is the first.
//!
//! # Every case is annotated twice, with the same object
//!
//! `composePriorPseudoCounts` memoises one array per allele count and returns **that array**. On
//! the empty-evidence branch the posteriors are that array, so the closing
//! `posteriors[i] -= prior[i]` zeroes the memo, and the next genotype with the same allele count
//! gets a prior of zeros. The reference's own second answer for `empty-evidence` is `NaN,NaN`,
//! because a prior of zeros makes `normalizeSumToOne` divide by zero.
//!
//! A suite that called the annotation once per case would miss it entirely, and a port that
//! copied the array would look more correct and be less faithful.

use std::io::Read;

use gatk_annotation::allele_pseudo_depth::{AllelePseudoDepth, PseudoDepthError, SampleMatrix};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/allele_pseudo_depth.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// One `case` row: everything needed to reproduce the call.
struct Case {
    label: String,
    is_natural_log: bool,
    prior: f64,
    keep_prior_in_count: bool,
    weight_decay: f64,
    allele_indices: Vec<usize>,
    mapping_qualities: Vec<i32>,
    log_likelihoods: Vec<Vec<f64>>,
    /// A null likelihoods object, which is one of the two silent guards.
    absent: bool,
}

fn parse_bits(field: &str) -> Vec<f64> {
    if field.is_empty() {
        return Vec::new();
    }
    field
        .split(',')
        .map(|hex| f64::from_bits(u64::from_str_radix(hex, 16).expect("hex bits")))
        .collect()
}

fn parse_case(fields: &[&str]) -> Case {
    let absent = fields[2] == "-";
    Case {
        label: fields[1].to_string(),
        is_natural_log: fields[2] == "true",
        prior: fields[3].parse().expect("prior"),
        keep_prior_in_count: fields[4] == "true",
        weight_decay: fields[5].parse().expect("weight decay"),
        allele_indices: if fields[6] == "-" {
            Vec::new()
        } else {
            fields[6]
                .split(',')
                .map(|index| index.parse().expect("allele index"))
                .collect()
        },
        mapping_qualities: if fields[7] == "-" {
            Vec::new()
        } else {
            fields[7]
                .split(',')
                .map(|quality| quality.parse().expect("mapping quality"))
                .collect()
        },
        log_likelihoods: if fields[8] == "-" {
            Vec::new()
        } else {
            // A matrix with no reads still has one empty row per allele, which is how the
            // empty-evidence cases reach the branch they exist for.
            fields[8].split(';').map(parse_bits).collect()
        },
        absent,
    }
}

#[test]
fn the_two_keys_are_the_strings_the_reference_writes() {
    let text = golden();
    let mut cases: Vec<Case> = Vec::new();
    // (label, call) -> the row the reference produced.
    let mut expected: Vec<(String, usize, String)> = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "case" => cases.push(parse_case(&fields)),
            "out" => expected.push((
                fields[1].to_string(),
                fields[2].parse().expect("call"),
                format!("{}\t{}", fields[3], fields[4]),
            )),
            "err" => expected.push((
                fields[1].to_string(),
                fields[2].parse().expect("call"),
                format!("E:{}", fields[3]),
            )),
            _ => {}
        }
    }

    let mut compared = 0usize;
    let mut second_call_differs = 0usize;
    let mut failures = Vec::new();

    for case in &cases {
        let mut annotation =
            AllelePseudoDepth::new(case.prior, case.keep_prior_in_count, case.weight_decay);
        let mut answers = Vec::new();
        for call in 1..=2usize {
            let matrix = SampleMatrix {
                log_likelihoods: &case.log_likelihoods,
                mapping_qualities: &case.mapping_qualities,
                is_natural_log: case.is_natural_log,
            };
            // The same annotation object across both calls, which is what carries the memo.
            let produced = annotation.annotate(
                &case.allele_indices,
                if case.absent { None } else { Some(&matrix) },
            );
            let rendered = match produced {
                // An absent key is dumped as `-`, and is not the same as an empty string.
                Ok(None) => "-\t-".to_string(),
                Ok(Some(depths)) => format!("{}\t{}", depths.depth, depths.fraction),
                Err(PseudoDepthError::EvidenceIndexOutOfBounds { .. }) => {
                    "E:java.lang.IndexOutOfBoundsException".to_string()
                }
                Err(PseudoDepthError::NonFiniteSum) => {
                    "E:java.lang.IllegalArgumentException".to_string()
                }
            };
            let want = expected
                .iter()
                .find(|(label, n, _)| *label == case.label && *n == call)
                .map(|(_, _, row)| row.clone())
                .unwrap_or_else(|| panic!("{} call {call} is missing from the golden", case.label));

            let matches = if let Some(class) = want.strip_prefix("E:") {
                // The message carries an index and a length, which the port reports in its own
                // words; the class is the part the reference's behaviour is.
                rendered.starts_with("E:") && class.starts_with(&rendered[2..])
            } else {
                rendered == want
            };
            if !matches {
                failures.push(format!(
                    "{} call {call}: ours {rendered:?}, reference {want:?}",
                    case.label
                ));
            }
            answers.push(rendered);
            compared += 1;
        }
        if answers[0] != answers[1] {
            second_call_differs += 1;
        }
    }

    assert!(
        failures.is_empty(),
        "{} divergence(s):\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert_eq!(compared, 46, "the golden changed size");
    // If no case ever answers differently the second time, the memo write-through has stopped
    // being exercised and the suite has quietly lost the behaviour it was built for.
    assert!(
        second_call_differs >= 2,
        "only {second_call_differs} case(s) answered differently on the second call; the \
         empty-evidence cases that write through the memo are no longer covered"
    );
    println!(
        "AllelePseudoDepth: {compared} calls compared, {second_call_differs} of them changed \
         answer on the second call through the same annotation object"
    );
}
