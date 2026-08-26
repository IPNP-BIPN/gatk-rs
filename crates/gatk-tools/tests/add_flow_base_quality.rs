//! Conformance for `AddFlowBaseQuality` against GATK 4.6.2.0, compared as the whole quality string
//! of every read of every run.
//!
//! Golden from `tools/readfilter-conformance/AddFlowBaseQualityDump.java`.
//!
//! # What this suite is for
//!
//!  * **the middle bases of an hmer never being computed**, and coming out at the MAXIMAL quality
//!    rather than the minimal one;
//!  * **both ends of the read being overridden** with the hmer's own key probability;
//!  * **the clamp and the floor**, through `--maximal-quality-score` and `--minimal-error-rate`;
//!  * **`--replace-quality-mode` moving the old qualities to `OQ`** instead of writing `XQ`;
//!  * **and a flow order whose cycle is one throwing out of the enumeration**.

use gatk_corpus as corpus;
use gatk_tools::add_flow_base_quality::{
    add_base_quality, calc_flow_order_length, convert_error_prob_to_phred,
    extract_error_prob_bands, generate_base_error_probability, RawProbs, ReadOutput,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/add_flow_base_quality.txt.gz"),
    )
}

/// The inverse of `ReferenceQueryDump.escape`, which escapes the BACKSLASH first.
///
/// A naive pair of replacements would be wrong here: a quality of 59 is written as `\\`, and these
/// strings are full of them.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn field<'a>(line: &'a str, name: &str) -> &'a str {
    line.split('\t')
        .find_map(|part| part.strip_prefix(&format!("{name}=")))
        .unwrap_or_else(|| panic!("the row carries {name}"))
}

/// One read as it went in: its group, its bases and its original qualities.
struct Read {
    group: String,
    bases: String,
    qualities: String,
}

fn read(text: &str, name: &str) -> Read {
    let prefix = format!("read\t{name}\t");
    let line = text
        .lines()
        .find(|line| line.starts_with(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries read/{name}"));
    Read {
        group: field(line, "group").to_string(),
        bases: field(line, "bases").to_string(),
        qualities: field(line, "qual").to_string(),
    }
}

/// The flow order of a read group, which is the only thing the cycle is computed from.
fn flow_order(text: &str, group: &str) -> String {
    let prefix = format!("group\t{group}\t");
    field(
        text.lines()
            .find(|line| line.starts_with(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries group/{group}")),
        "order",
    )
    .to_string()
}

fn key(text: &str, name: &str) -> Vec<i32> {
    let prefix = format!("flow\t{name}\t");
    field(
        text.lines()
            .find(|line| line.starts_with(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries flow/{name}")),
        "key",
    )
    .split(',')
    .map(|value| value.parse().expect("an hmer length"))
    .collect()
}

/// The three raw probabilities per flow. `none` is the neighbour that does not exist.
fn probs(text: &str, name: &str) -> Vec<RawProbs> {
    let prefix = format!("prob\t{name}\t");
    text.lines()
        .filter(|line| line.starts_with(prefix.as_str()))
        .map(|line| {
            let optional = |part: &str| match field(line, part) {
                "none" => None,
                value => Some(value.parse::<f64>().expect("a probability")),
            };
            RawProbs {
                minus: optional("minus"),
                key: field(line, "key").parse().expect("a probability"),
                plus: optional("plus"),
            }
        })
        .collect()
}

fn attribute(text: &str, kind: &str, label: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}\t{name}=");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(unescape)
}

/// The reads of the four-flow input, in the order the walker saw them.
const READS: &[&str] = &[
    "singles",
    "doubles",
    "long-hmer",
    "leading-zero",
    "trailing-hmer",
    "varied-quality",
];

/// label, minimal error rate, maximal quality score, replace mode.
const RUNS: &[(&str, f64, i32, bool)] = &[
    ("default", 1e-3, 93, false),
    ("replace", 1e-3, 93, true),
    ("max-quality", 1e-3, 10, false),
    ("min-error-rate", 0.1, 93, false),
    ("both", 1e-3, 20, true),
];

fn produced(text: &str, name: &str, run: (&str, f64, i32, bool)) -> ReadOutput {
    let (_, min_error_rate, max_quality_score, replace) = run;
    let read = read(text, name);
    add_base_quality(
        &key(text, name),
        &probs(text, name),
        read.bases.len(),
        &read.qualities,
        &flow_order(text, &read.group),
        min_error_rate,
        max_quality_score,
        replace,
    )
    .expect("a four-flow order never leaves the slice")
}

#[test]
fn every_read_of_every_run_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for run in RUNS {
        for name in READS {
            let output = produced(&text, name, *run);
            assert_eq!(
                output.base_quality_attribute,
                attribute(&text, "xq", run.0, name),
                "XQ of {name} in {}",
                run.0
            );
            assert_eq!(
                output.old_quality_attribute,
                attribute(&text, "oq", run.0, name),
                "OQ of {name} in {}",
                run.0
            );
            assert_eq!(
                Some(output.qualities),
                attribute(&text, "qual", run.0, name),
                "QUAL of {name} in {}",
                run.0
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 30, "six reads over five runs");
}

/// The three middle bases of a five-base homopolymer are never written, so they keep the zero the
/// array was allocated with, and a zero is the MAXIMAL quality rather than the minimal one.
#[test]
fn the_middle_of_an_hmer_comes_out_at_the_maximal_quality() {
    let text = golden();
    let read = read(&text, "long-hmer");
    assert!(read.bases.starts_with("TTTTT"), "a five-base run");

    let bands = extract_error_prob_bands(&probs(&text, "long-hmer"), 1e-3);
    let error_prob = generate_base_error_probability(
        &key(&text, "long-hmer"),
        &bands,
        read.bases.len(),
        calc_flow_order_length(&flow_order(&text, &read.group)),
    )
    .expect("a computed read");
    assert_eq!(
        error_prob[1..4],
        [0.0, 0.0, 0.0],
        "the three bases the cursor jumped"
    );
    assert_eq!(
        convert_error_prob_to_phred(&error_prob, 93)[1..4],
        [93, 93, 93]
    );
    assert_eq!(
        attribute(&text, "xq", "default", "long-hmer").as_deref(),
        Some("!~~~\\YY!")
    );
}

/// The first and the last base of the READ take the hmer's own key probability, and the last one
/// overwrites a value the hmer loop had already written.
#[test]
fn both_ends_of_the_read_are_overridden() {
    let text = golden();
    for name in READS {
        let read = read(&text, name);
        let raw = probs(&text, name);
        let bands = extract_error_prob_bands(&raw, 1e-3);
        let key = key(&text, name);
        let error_prob = generate_base_error_probability(
            &key,
            &bands,
            read.bases.len(),
            calc_flow_order_length(&flow_order(&text, &read.group)),
        )
        .expect("a computed read");

        let first_flow = key.iter().position(|hmer| *hmer != 0).expect("an hmer");
        let last_flow = key.iter().rposition(|hmer| *hmer != 0).expect("an hmer");
        assert_eq!(
            error_prob[0], bands[1][first_flow],
            "the first base of {name}"
        );
        assert_eq!(
            error_prob[read.bases.len() - 1],
            bands[1][last_flow],
            "the last base of {name}"
        );
    }
}

/// The clamp is a minimum against the computed score, so a low maximum flattens every base that is
/// not an end.
#[test]
fn the_maximal_quality_score_clamps_every_base() {
    let text = golden();
    let clamped = attribute(&text, "xq", "max-quality", "long-hmer").expect("an XQ");
    for character in clamped.chars() {
        assert!(
            (character as u8) <= b'!' + 10,
            "no base exceeds the maximum of 10"
        );
    }
    // The middle of the hmer was 93 under the default and is 10 here, so the clamp reaches the
    // bases that were never computed as well as the ones that were.
    assert_eq!(clamped, "!++++++!");
}

/// The floor lifts every band, which lowers the confidence the enumeration can reach.
#[test]
fn the_minimal_error_rate_is_a_floor_under_every_band() {
    let text = golden();
    let raw = probs(&text, "singles");
    let low = extract_error_prob_bands(&raw, 1e-3);
    let high = extract_error_prob_bands(&raw, 0.1);
    for band in 0..3 {
        for flow in 0..raw.len() {
            assert!(high[band][flow] >= low[band][flow]);
            assert!(high[band][flow] >= 0.1);
        }
    }
    // A flow whose key is 0 has no shorter neighbour, and takes the floor rather than a zero.
    let leading = probs(&text, "leading-zero");
    let flow = leading
        .iter()
        .position(|raw| raw.minus.is_none())
        .expect("a flow with no shorter neighbour");
    assert_eq!(extract_error_prob_bands(&leading, 1e-3)[0][flow], 1e-3);
}

/// Replacing writes the computed string into QUAL and the original into OQ; the default writes XQ
/// and leaves QUAL alone.
#[test]
fn replace_quality_mode_moves_the_old_qualities_to_oq() {
    let text = golden();
    for name in READS {
        let original = read(&text, name).qualities;
        let default = produced(&text, name, RUNS[0]);
        let replaced = produced(&text, name, RUNS[1]);

        assert_eq!(default.qualities, original, "{name} keeps its qualities");
        assert!(default.old_quality_attribute.is_none());
        assert_eq!(
            replaced.old_quality_attribute.as_deref(),
            Some(&original[..])
        );
        assert!(replaced.base_quality_attribute.is_none());
        // The same string, written to two different places.
        assert_eq!(
            default.base_quality_attribute.as_deref(),
            Some(&replaced.qualities[..])
        );
    }
}

/// A flow order whose first base repeats immediately gives a cycle of one, a slice of a single
/// flow, and an index one past it.
#[test]
fn a_flow_order_of_cycle_one_leaves_the_slice() {
    let text = golden();
    assert_eq!(calc_flow_order_length("TGCA"), 4);
    assert_eq!(calc_flow_order_length("TTGCA"), 1);
    assert_eq!(flow_order(&text, "rg1"), "TTGCA");

    let read = read(&text, "cycle-one");
    let error = add_base_quality(
        &key(&text, "cycle-one"),
        &probs(&text, "cycle-one"),
        read.bases.len(),
        &read.qualities,
        &flow_order(&text, &read.group),
        1e-3,
        93,
        false,
    )
    .expect_err("a slice of one flow");
    assert_eq!(error.index, 1);
    assert_eq!(error.length, 1);

    let expected = text
        .lines()
        .find(|line| line.starts_with("error\tcycle-one\t"))
        .expect("the golden carries the refusal");
    assert!(
        expected.contains(&format!("{}:{}", error.java_class(), error.message())),
        "{expected}"
    );
    assert!(
        expected.contains("generateSidedHmerBaseErrorProbability"),
        "the frame the reference threw from"
    );
}
