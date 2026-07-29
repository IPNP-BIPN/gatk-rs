//! Conformance for `AlignmentStateMachine` against GATK 4.6.2.0.
//!
//! One row per `stepForwardOnGenome` call over 26 cigars, so the machine's corners are compared
//! individually rather than through a summary: which stops exist, what offsets each carries, and
//! which two cigars it refuses outright.
//!
//! The cigar is the whole input. The reads are built from it on both sides, at the same alignment
//! start the harness used, because a fixture BAM would put an encoder between the cigar and the
//! machine and this suite is about the machine.

use gatk_corpus as corpus;
use gatk_engine::alignment_state::{AlignmentStateMachine, MalformedRead};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;

/// The alignment start `AlignmentStateDump.makeRead` uses.
const ALIGNMENT_START: i32 = 101;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/alignment_state.txt.gz"),
    )
}

/// The read the harness built for a cigar: bases as long as the cigar's read length, cycling
/// `ACGT`, all qualities 30.
fn read_for(cigar_text: &str) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar_text).expect("a parsable cigar");
    let length = cigar.read_length() as usize;
    BamRecord {
        read_name: format!("read-{cigar_text}"),
        reference_index: 0,
        alignment_start: ALIGNMENT_START,
        mapping_quality: 60,
        read_bases: (0..length).map(|i| b"ACGT"[i % 4]).collect(),
        base_qualities: vec![30; length],
        cigar,
        ..Default::default()
    }
}

#[test]
fn every_stop_matches_the_reference() {
    let text = golden();

    // The cigars in the order the harness ran them, with their rows and their outcome.
    let mut cigars: Vec<String> = Vec::new();
    let mut steps: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut outcomes: std::collections::HashMap<String, (String, usize)> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("step\t") {
            let mut parts = rest.splitn(3, '\t');
            let cigar = parts.next().expect("a cigar").to_string();
            let _index = parts.next();
            steps
                .entry(cigar)
                .or_default()
                .push(parts.next().unwrap_or("").to_string());
        } else if let Some(rest) = line.strip_prefix("end\t") {
            let mut parts = rest.split('\t');
            let cigar = parts.next().expect("a cigar").to_string();
            let outcome = parts.next().expect("an outcome").to_string();
            let count: usize = parts.next().expect("a count").parse().expect("a number");
            cigars.push(cigar.clone());
            outcomes.insert(cigar, (outcome, count));
        }
    }
    assert!(!cigars.is_empty(), "the golden carries no end rows");

    let mut compared = 0;
    let mut refused = 0;
    for cigar in &cigars {
        let read = read_for(cigar);
        let mut machine = AlignmentStateMachine::new(&read);
        let (expected_outcome, expected_count) = &outcomes[cigar];
        let expected_rows = steps.get(cigar).cloned().unwrap_or_default();

        let mut ours: Vec<String> = Vec::new();
        let mut failure: Option<MalformedRead> = None;
        loop {
            match machine.step_forward_on_genome() {
                Err(error) => {
                    failure = Some(error);
                    break;
                }
                Ok(op) => {
                    ours.push(format!(
                        "{}|{}|{}|{}|{}|{}|{}|{}",
                        // `CigarOperator.toString` prints the SAM character, so the sequence-match
                        // operator appears as `=` and not as the enum constant `EQ`.
                        match op {
                            None => "null".to_string(),
                            Some(Op::Eq) => "=".to_string(),
                            Some(op) => format!("{op:?}").to_uppercase(),
                        },
                        machine.read_offset(),
                        machine.genome_offset(),
                        machine.genome_position(),
                        machine.current_cigar_element_offset(),
                        machine.offset_into_current_cigar_element(),
                        machine.is_left_edge(),
                        machine.is_right_edge(),
                    ));
                    if op.is_none() {
                        break;
                    }
                }
            }
        }

        assert_eq!(ours.len(), *expected_count, "{cigar}: step count");
        assert_eq!(&ours, &expected_rows, "{cigar}");
        match (failure, expected_outcome.as_str()) {
            (None, "ok") => {}
            (Some(_), outcome) if outcome.starts_with("E:") => refused += 1,
            (None, outcome) => panic!("{cigar}: the reference raised {outcome}, the port did not"),
            (Some(error), _) => {
                panic!("{cigar}: the port raised {error:?}, the reference did not")
            }
        }
        compared += ours.len();
    }

    println!(
        "{compared} stops over {} cigars, {refused} refused exactly where the reference refused",
        cigars.len()
    );
}
