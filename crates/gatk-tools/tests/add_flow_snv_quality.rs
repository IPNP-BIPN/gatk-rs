//! Conformance for `AddFlowSNVQuality` against GATK 4.6.2.0, compared as the whole quality string
//! and every per-base attribute of every read of every run.
//!
//! Golden from `tools/readfilter-conformance/AddFlowSNVQualityDump.java`.
//!
//! # What this suite is for
//!
//!  * **the computed base quality being discarded** and rebuilt from the alternates;
//!  * **the phred conversion rounding**, where `AddFlowBaseQuality` truncates;
//!  * **the four snvq modes** over the same two probabilities;
//!  * **`--max-phred-score` moving the floor as well as the clamp**;
//!  * **`--output-quality-attribute` leaving `QUAL` alone**;
//!  * **and a flow order of cycle one dying in the normalisation**, where the sibling dies in the
//!    enumeration.

use gatk_corpus as corpus;
use gatk_tools::add_flow_base_quality::RawProbs;
use gatk_tools::add_flow_snv_quality::{
    add_base_quality, attr_name_for_non_called_base, convert_error_prob_to_phred, get_snvq,
    SnvError, SnvqMode,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/add_flow_snv_quality.txt.gz"),
    )
}

/// The inverse of `ReferenceQueryDump.escape`, which escapes the BACKSLASH first: a quality of 59
/// is the character `\`, and these strings are full of them.
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

/// `calcFlowOrderLength`: the distance to the second occurrence of the order's first base.
fn cycle(order: &str) -> usize {
    let bytes = order.as_bytes();
    match bytes.iter().skip(1).position(|base| *base == bytes[0]) {
        Some(offset) => offset + 1,
        None => bytes.len(),
    }
}

fn qual(text: &str, label: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("qual\t{label}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries qual/{label}/{name}")),
    )
}

fn attr(text: &str, label: &str, name: &str, tag: &str) -> Option<String> {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("attr\t{label}\t{name}\t{tag}=")))
        .map(unescape)
}

const READS: &[&str] = &[
    "singles",
    "doubles",
    "long-hmer",
    "leading-zero",
    "trailing-hmer",
    "varied-quality",
];

/// label, mode, maximum phred score, output attribute.
const RUNS: &[(&str, SnvqMode, Option<f64>, Option<&str>)] = &[
    ("default", SnvqMode::Geometric, None, None),
    ("legacy", SnvqMode::Legacy, None, None),
    ("optimistic", SnvqMode::Optimistic, None, None),
    ("pessimistic", SnvqMode::Pessimistic, None, None),
    ("max-phred-20", SnvqMode::Geometric, Some(20.0), None),
    ("max-phred-10", SnvqMode::Geometric, Some(10.0), None),
    ("attribute", SnvqMode::Geometric, None, Some("BQ")),
];

fn produced(
    text: &str,
    name: &str,
    run: (&str, SnvqMode, Option<f64>, Option<&str>),
) -> gatk_tools::add_flow_snv_quality::ReadOutput {
    let (_, mode, max_phred, attribute) = run;
    let read = read(text, name);
    let order = flow_order(text, &read.group);
    add_base_quality(
        &key(text, name),
        &probs(text, name),
        read.bases.as_bytes(),
        &read.qualities,
        &order,
        cycle(&order),
        max_phred,
        mode,
        attribute,
    )
    .expect("a four-flow order never leaves the order")
}

#[test]
fn every_quality_and_attribute_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for run in RUNS {
        for name in READS {
            let output = produced(&text, name, *run);
            assert_eq!(
                output.qualities,
                qual(&text, run.0, name),
                "QUAL of {name} in {}",
                run.0
            );
            assert_eq!(
                output.quality_attribute,
                attr(&text, run.0, name, "BQ"),
                "BQ of {name} in {}",
                run.0
            );
            for (tag, value) in &output.snvq_attributes {
                assert_eq!(
                    Some(value.clone()),
                    attr(&text, run.0, name, tag),
                    "{tag} of {name} in {}",
                    run.0
                );
                compared += 1;
            }
        }
    }
    assert_eq!(compared, 7 * 6 * 4, "four attributes per read per run");
}

/// The base quality is one minus the called base's own probability, which is one minus the sum of
/// the alternates: nothing the hmer loop computed survives.
#[test]
fn the_base_quality_is_rebuilt_from_the_alternates() {
    let text = golden();
    for name in READS {
        let output = produced(&text, name, RUNS[0]);
        let read = read(&text, name);
        let order = flow_order(&text, &read.group);
        let cycle = cycle(&order);
        for (offset, called) in read.bases.bytes().enumerate() {
            // The called base's own attribute carries the probability the base quality is one
            // minus, so the two strings are complements at every position.
            let called_tag = attr_name_for_non_called_base(called);
            let called_value = output
                .snvq_attributes
                .iter()
                .find(|(tag, _)| *tag == called_tag)
                .map(|(_, value)| value.as_bytes()[offset])
                .expect("the called base's own row");
            let quality = output.qualities.as_bytes()[offset];
            // A high called-base probability is a low error probability, so the called row is the
            // small number and QUAL is the large one, never the other way round.
            assert!(
                quality >= called_value || quality == b'!',
                "{name} at {offset}: {quality} against {called_value}"
            );
        }
        assert_eq!(output.snvq_attributes.len(), cycle);
    }
}

/// The conversion rounds: 0.5 goes up, where `AddFlowBaseQuality` would have truncated it away.
#[test]
fn the_phred_conversion_rounds() {
    // -10 * log10(p) = 3.5 exactly at this probability, which rounds to 4 and truncates to 3.
    let probability = 10f64.powf(-0.35);
    assert_eq!(convert_error_prob_to_phred(&[probability], 60), vec![4]);
    // And a probability of zero takes the clamp rather than an infinity.
    assert_eq!(convert_error_prob_to_phred(&[0.0], 60), vec![60]);
    assert_eq!(convert_error_prob_to_phred(&[0.0], 20), vec![20]);
}

/// Four formulae over the same two probabilities, and the default is the geometric one.
#[test]
fn the_four_snvq_modes_are_four_formulae() {
    let (slice, p1, p2) = (0.25, 0.4, 0.6);
    assert_eq!(get_snvq(slice, p1, p2, SnvqMode::Legacy), slice);
    assert_eq!(get_snvq(slice, p1, p2, SnvqMode::Optimistic), p1 * p2);
    assert_eq!(
        get_snvq(slice, p1, p2, SnvqMode::Pessimistic),
        1.0 - (1.0 - p1) * (1.0 - p2)
    );
    let geometric = get_snvq(slice, p1, p2, SnvqMode::Geometric);
    assert!(geometric > p1 * p2 && geometric < 1.0 - (1.0 - p1) * (1.0 - p2));

    // And the modes really do move the output.
    let text = golden();
    assert_ne!(
        qual(&text, "default", "singles"),
        qual(&text, "legacy", "singles")
    );
    assert_ne!(
        qual(&text, "optimistic", "singles"),
        qual(&text, "pessimistic", "singles")
    );
}

/// At 10 the floor is 0.1, so three alternates take 0.3 of the mass and every base comes out at 5:
/// the clamp is never reached.
#[test]
fn the_maximum_phred_score_moves_the_floor_as_well() {
    let text = golden();
    let low = qual(&text, "max-phred-10", "singles");
    assert!(
        low.bytes().all(|value| value == b'&'),
        "every base at 5, not at the clamp of 10: {low}"
    );
    assert_eq!(produced(&text, "singles", RUNS[5]).qualities, low);

    // At 20 the floor is 0.01, so three alternates take 0.03 and every base comes out at 15.
    // Neither run reaches its clamp: the floor is what decides both.
    let higher = qual(&text, "max-phred-20", "singles");
    assert!(
        higher.bytes().all(|value| value == b'0'),
        "every base at 15: {higher}"
    );
    assert_eq!(produced(&text, "singles", RUNS[4]).qualities, higher);
    assert_ne!(low.as_bytes()[0] - b'!', 10, "the clamp of the first run");
    assert_ne!(higher.as_bytes()[0] - b'!', 20, "the clamp of the second");
}

/// The attribute run writes the base quality to a tag and leaves QUAL as it was.
#[test]
fn the_output_attribute_leaves_the_quality_alone() {
    let text = golden();
    for name in READS {
        let original = read(&text, name).qualities;
        let output = produced(&text, name, RUNS[6]);
        assert_eq!(output.qualities, original, "{name} keeps its qualities");
        assert_eq!(
            output.quality_attribute,
            Some(qual(&text, "default", name)),
            "{name} writes the default run's quality into the tag"
        );
    }
}

/// A flow order of cycle one gets past the enumeration and dies in the normalisation, where the
/// sibling tool dies in the enumeration.
#[test]
fn a_flow_order_of_cycle_one_dies_in_the_normalisation() {
    let text = golden();
    let read = read(&text, "cycle-one");
    let order = flow_order(&text, &read.group);
    assert_eq!(order, "TTGCA");
    assert_eq!(cycle(&order), 1);

    let error = add_base_quality(
        &key(&text, "cycle-one"),
        &probs(&text, "cycle-one"),
        read.bases.as_bytes(),
        &read.qualities,
        &order,
        cycle(&order),
        None,
        SnvqMode::Geometric,
        None,
    )
    .expect_err("a base outside the cycle");
    assert_eq!(
        error,
        SnvError::CalledBaseNotInOrder {
            index: -1,
            length: 1
        }
    );

    let expected = text
        .lines()
        .find(|line| line.starts_with("error\tcycle-one\t"))
        .expect("the golden carries the refusal");
    assert!(
        expected.contains(&format!("{}:{}", error.java_class(), error.message())),
        "{expected}"
    );
    // The index is -1 here; the sibling's was 1, from a different line.
    assert!(expected.contains("Index -1 out of bounds for length 1"));
}
