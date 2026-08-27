//! Conformance for `CalibrateDragstrModel` against GATK 4.6.2.0, compared as the parameter tables
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/CalibrateDragstrModelDump.java`.
//!
//! Scanning the reference and piling up the reads are not measured or ported. What is compared is
//! the table's shape, the rows the estimation never reached, the GCP column, and the refusals.
//!
//! # What this suite is for
//!
//!  * **the table's shape being the hyper-parameters' and not the data's**;
//!  * **a period with no data keeping the defaults, indistinguishably**;
//!  * **GCP never being estimated**;
//!  * **every row being constant across its repeat lengths**;
//!  * **the sites census, and `--minimum-depth` deciding what is used**;
//!  * **the precomputed error probabilities**;
//!  * **and the four refusals.**

use gatk_corpus as corpus;
use gatk_tools::calibrate_dragstr_model::{
    column_header, estimate_period, flanks, gcp_row, initial_groups, log10_one_minus_pow10,
    log10_prob, min_gp_index, precompute, row, table, value_range, Case, Cases, HyperParameters,
    BLOCKS, LOG10_ONE_HALF,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/calibrate_dragstr_model.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

/// The three blocks of one run's table, each as its rows of numbers.
fn blocks(text: &str, label: &str) -> Vec<Vec<Vec<f64>>> {
    let file = section(text, "out", label);
    let mut blocks: Vec<Vec<Vec<f64>>> = Vec::new();
    for line in file.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if BLOCKS.iter().any(|name| line == format!("{name}:")) {
            blocks.push(Vec::new());
            continue;
        }
        if blocks.is_empty() {
            // The column header, before the first block.
            continue;
        }
        blocks.last_mut().expect("a block").push(
            line.split_whitespace()
                .map(|v| v.parse().expect("a number"))
                .collect(),
        );
    }
    blocks
}

/// The header line of one run's table, which is the repeat lengths.
fn header(text: &str, label: &str) -> String {
    section(text, "out", label)
        .lines()
        .find(|line| !line.starts_with('#') && !line.is_empty())
        .expect("a header")
        .to_string()
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), unescape(message))
}

/// The shape is the hyper-parameters', whatever the reference holds.
#[test]
fn the_tables_shape_is_the_hyper_parameters() {
    let text = golden();
    for (label, periods, repeats) in [
        ("forced", 8, 20),
        ("parallel", 8, 20),
        ("repeat-eight", 8, 8),
        ("min-loci-one", 8, 20),
        ("min-depth-ten", 8, 20),
    ] {
        let blocks = blocks(&text, label);
        assert_eq!(blocks.len(), 3, "{label}");
        for block in &blocks {
            assert_eq!(block.len(), periods, "{label}");
            for row in block {
                assert_eq!(row.len(), repeats, "{label}");
            }
        }
        // And the header names the repeat lengths in that width.
        assert_eq!(header(&text, label), column_header(repeats), "{label}");
    }
    // The STR table itself was composed for three periods, so the eight rows are the argument's
    // and not the file's.
    assert_eq!(HyperParameters::default().max_period, 8);
    assert_eq!(HyperParameters::default().max_repeat_length, 20);
}

/// The file gives no sign of which rows were estimated.
#[test]
fn a_period_with_no_data_keeps_the_defaults() {
    let text = golden();
    let forced = blocks(&text, "forced");
    let (gop, api) = (&forced[0], &forced[2]);
    // The fixture's repeats are periods one to three, and those three rows moved.
    assert_eq!(gop[0][0], 49.50);
    assert_eq!(gop[1][0], 43.25);
    assert_eq!(gop[2][0], 40.75);
    assert_eq!(api[0][0], 40.00);
    assert_eq!(api[1][0], 40.00);
    assert_eq!(api[2][0], 40.00);
    // Periods four to eight kept the defaults, whose API is zero.
    for (index, row) in api.iter().enumerate().skip(3) {
        assert_eq!(row[0], 0.00, "period {}", index + 1);
    }
    // A default row and an estimated one are written the same way, so the file cannot be read to
    // tell them apart: both are a full row of the same number.
    for (gop, api) in gop.iter().zip(api.iter()) {
        assert!(gop.iter().all(|value| *value == gop[0]));
        assert!(api.iter().all(|value| *value == api[0]));
    }
}

/// Every row of it is ten over the period.
#[test]
fn gcp_is_never_estimated() {
    let text = golden();
    for label in ["forced", "parallel", "min-loci-one", "min-depth-ten"] {
        let gcp = &blocks(&text, label)[1];
        for (index, written) in gcp.iter().enumerate() {
            let period = index + 1;
            let produced = gcp_row(period, written.len());
            // The file rounds to two decimals, so the comparison is on the written form.
            assert_eq!(row(&produced), row(written), "{label} period {period}");
        }
        assert_eq!(gcp[0][0], 10.00, "{label}");
        assert_eq!(gcp[1][0], 5.00, "{label}");
        assert_eq!(gcp[7][0], 1.25, "{label}");
    }
}

/// The parallel run and the serial one agree on every number.
#[test]
fn running_in_parallel_changes_nothing() {
    let text = golden();
    assert_eq!(blocks(&text, "forced"), blocks(&text, "parallel"));
    // The two files differ only in the command line their headers repeat.
    let strip = |label: &str| -> Vec<String> {
        section(&text, "out", label)
            .lines()
            .filter(|line| !line.starts_with("# commandLine"))
            .map(str::to_string)
            .collect()
    };
    assert_eq!(strip("forced"), strip("parallel"));
}

/// A site under the minimum depth is skipped rather than used.
#[test]
fn the_minimum_depth_decides_what_is_used() {
    let text = golden();
    let count = |label: &str, status: &str| -> Option<i64> {
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("status\t{label}\t{status}=")))
            .map(|value| value.parse().expect("a number"))
    };
    // The fixture gives every block twelve reads, above the default minimum of ten.
    assert!(count("forced", "used").expect("used sites") > 0);
    assert_eq!(count("forced", "skipped"), None);
    // Raising the minimum to ten changes nothing, it being the default already.
    assert_eq!(count("min-depth-ten", "used"), count("forced", "used"));
    // The rest of the sites are capped out by --down-sample-size, which is a cap per period and
    // repeat length rather than a total.
    assert!(count("forced", "downsampled-out").expect("capped sites") > 0);
    let rows: i64 = text
        .lines()
        .find_map(|line| line.strip_prefix("rows\tforced="))
        .expect("a row count")
        .parse()
        .expect("a number");
    assert_eq!(
        rows,
        count("forced", "used").expect("used")
            + count("forced", "downsampled-out").expect("capped")
    );
}

/// The two tables the estimator computes before it looks at any data.
#[test]
fn the_error_probabilities_are_precomputed_from_the_hyper_parameters() {
    let parameters = HyperParameters::default();
    let precomputed = precompute(&parameters);
    assert_eq!(
        precomputed.log10_p_error.len(),
        parameters.phred_gp_values.len()
    );
    assert_eq!(precomputed.log10_p_error[0].len(), parameters.max_period);
    assert_eq!(
        precomputed.log10_p_error[0][0].len(),
        parameters.max_repeat_length
    );
    // The correct probability falls with the repeat's length in BASES, so period two at three
    // repeats is period one at six.
    let one_at_six = precomputed.log10_p_correct[0][0][5];
    let two_at_three = precomputed.log10_p_correct[0][1][2];
    assert!((one_at_six - two_at_three).abs() < 1e-12);
    // And the two tables always sum to one.
    for i in [0, 10, 40] {
        for period in 0..parameters.max_period {
            for repeats in 0..parameters.max_repeat_length {
                let correct = 10f64.powf(precomputed.log10_p_correct[i][period][repeats]);
                let error = 10f64.powf(precomputed.log10_p_error[i][period][repeats]);
                assert!(
                    (correct + error - 1.0).abs() < 1e-9,
                    "{i} {period} {repeats}"
                );
            }
        }
    }
    // The value ranges include their ends.
    assert_eq!(parameters.phred_gp_values.len(), 41);
    assert_eq!(parameters.phred_api_values.len(), 41);
    assert_eq!(parameters.phred_gop_values.len(), 161);
    assert_eq!(*parameters.phred_gp_values.last().expect("a value"), 50.0);
    assert_eq!(value_range(0.0, 1.0, 3.0), vec![0.0, 1.0, 2.0, 3.0]);
}

/// The floor is read from `--max-repeats`, so changing the shape moves the search's range.
#[test]
fn the_gp_floor_depends_on_the_maximum_repeat_length() {
    let twenty = HyperParameters::default();
    let eight = HyperParameters {
        max_repeat_length: 8,
        ..HyperParameters::default()
    };
    // A longer repeat needs a HIGHER floor: the per-position error has to be smaller for the
    // whole repeat to survive, and a smaller error is a larger Phred value.
    assert!(min_gp_index(&twenty, 1) > min_gp_index(&eight, 1));
    // A longer period makes the repeat longer in bases, so it raises the floor for the same
    // reason.
    assert!(min_gp_index(&twenty, 8) > min_gp_index(&twenty, 1));
    // Every floor is inside the array.
    for period in 1..=twenty.max_period {
        assert!(min_gp_index(&twenty, period) < twenty.phred_gp_values.len());
    }
}

/// `log10(1 - 10^a)`, whose two special arguments are answered before any arithmetic.
#[test]
fn the_one_minus_probability_is_exact_at_its_ends() {
    assert!(log10_one_minus_pow10(1.0).is_nan());
    assert_eq!(log10_one_minus_pow10(0.0), f64::NEG_INFINITY);
    // Half in log10 space: one less a half is a half.
    assert!((log10_one_minus_pow10(LOG10_ONE_HALF) - LOG10_ONE_HALF).abs() < 1e-12);
    // A very small probability leaves almost one behind.
    assert!(log10_one_minus_pow10(-30.0).abs() < 1e-12);
    // And the constant is the logarithm it says it is.
    assert!((10f64.powf(LOG10_ONE_HALF) - 0.5).abs() < 1e-15);
}

/// The three genotype terms, of which the last is only present when every read carried the indel.
#[test]
fn the_likelihood_has_a_homozygous_term_only_when_every_read_agrees() {
    let (error, correct) = (-1.0, -0.05);
    let (hom_ref, het, hom_var) = (-0.2, -0.7, -1.0);
    // Ten reads, none with an indel: the homozygous-variant term is absent.
    let none = log10_prob(10, 0, error, correct, hom_ref, het, hom_var);
    // Ten reads, all with an indel: it is present.
    let all = log10_prob(10, 10, error, correct, hom_ref, het, hom_var);
    assert!(all > log10_prob(10, 9, error, correct, hom_ref, het, hom_var));
    assert!(none.is_finite() && all.is_finite());
    // Removing the term by hand gives the value the nine-indel case is measured against.
    let without = log10_prob(10, 10, error, correct, hom_ref, het, f64::NEG_INFINITY);
    assert!(all > without);
}

/// A period with no data leaves the flanks crossed and the whole range as one group.
#[test]
fn a_period_with_no_data_is_one_group() {
    let parameters = HyperParameters::default();
    let empty = Cases::empty(parameters.max_period, parameters.max_repeat_length);
    let (left, right) = flanks(1, &parameters, &empty);
    assert_eq!(left, parameters.max_repeat_length);
    assert_eq!(right, 1);
    assert!(right < left, "the flanks crossed");
    let groups = initial_groups(&parameters, left, right);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0], 1..=parameters.max_repeat_length);
    // With data everywhere the groups are the flanks and every repeat length between them.
    let mut full = Cases::empty(parameters.max_period, parameters.max_repeat_length);
    for repeats in 1..=parameters.max_repeat_length {
        for _ in 0..parameters.min_loci_count {
            full.add(
                1,
                repeats,
                Case {
                    depth: 12,
                    indels: 0,
                },
            );
        }
    }
    let (left, right) = flanks(1, &parameters, &full);
    assert_eq!(left, 1);
    assert_eq!(right, parameters.max_repeat_length - 1);
    let groups = initial_groups(&parameters, left, right);
    assert_eq!(groups.len(), parameters.max_repeat_length);
    assert_eq!(groups[0], 1..=1);
    assert_eq!(
        *groups.last().expect("a group"),
        parameters.max_repeat_length..=parameters.max_repeat_length
    );
}

/// The search settles on a group, and the estimation is monotone across the row.
#[test]
fn the_estimation_is_monotone_across_the_row() {
    let parameters = HyperParameters {
        max_period: 2,
        max_repeat_length: 6,
        min_loci_count: 4,
        ..HyperParameters::default()
    };
    let precomputed = precompute(&parameters);
    let mut cases = Cases::empty(parameters.max_period, parameters.max_repeat_length);
    // The longer the repeat, the more reads carry an indel, which is the signal the model is for.
    for repeats in 1..=parameters.max_repeat_length {
        for _ in 0..10 {
            cases.add(
                1,
                repeats,
                Case {
                    depth: 12,
                    indels: repeats as i32 - 1,
                },
            );
        }
    }
    let estimated = estimate_period(1, &parameters, &precomputed, &cases);
    assert!(!estimated.is_empty());
    // Both columns fall, or hold, from one group to the next.
    for window in estimated.windows(2) {
        assert!(window[0].1.gp >= window[1].1.gp);
        assert!(window[0].1.api + parameters.api_mono_threshold >= window[1].1.api);
    }
    // The groups cover every repeat length exactly once, in order.
    let mut covered = Vec::new();
    for (range, _) in &estimated {
        for repeats in range.clone() {
            covered.push(repeats);
        }
    }
    assert_eq!(
        covered,
        (1..=parameters.max_repeat_length).collect::<Vec<_>>()
    );
    // And GCP is ten over the period whatever the search did.
    for (_, estimate) in &estimated {
        assert_eq!(estimate.gcp, 10.0);
    }
}

/// The layout of the file, rebuilt from a run's own numbers.
#[test]
fn the_file_is_a_header_and_three_blocks() {
    let text = golden();
    let parameters = HyperParameters::default();
    let written = blocks(&text, "forced");
    let rows: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..parameters.max_period)
        .map(|period| {
            (
                written[0][period].clone(),
                written[1][period].clone(),
                written[2][period].clone(),
            )
        })
        .collect();
    let produced = table(&parameters, &rows);
    // The golden's file is its header comments and then exactly this.
    let body: String = section(&text, "out", "forced")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    assert_eq!(produced, body);
    assert_eq!(BLOCKS, ["GOP", "GCP", "API"]);
}

/// Four refusals, each from a different layer.
#[test]
fn the_four_refusals_are_the_goldens() {
    let text = golden();
    // Two are the parser's, on arguments with a minimum.
    for (label, argument, minimum) in [
        ("down-sample-too-small", "down-sample-size", "512"),
        ("shard-too-small", "shard-size", "100"),
    ] {
        let (class, message) = refusal(&text, label);
        assert_eq!(
            class,
            "org.broadinstitute.barclay.argparser.CommandLineException$OutOfRangeArgumentValue"
        );
        assert!(message.contains(argument), "{message}");
        assert!(message.contains(minimum), "{message}");
    }
    // One is the dictionary check, before any read is looked at.
    let (class, message) = refusal(&text, "wrong-reference");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.GATKException"
    );
    assert!(message.contains("UNEQUAL_COMMON_CONTIGS"), "{message}");
    // And one is not a refusal at all: a --max-period below the STR table's own is an index
    // error, the estimator allocating from the argument and indexing from the file.
    let (class, message) = refusal(&text, "period-two");
    assert_eq!(class, "java.lang.ArrayIndexOutOfBoundsException");
    assert!(
        message.contains("Index 2 out of bounds for length 2"),
        "{message}"
    );
}
