//! Conformance for `ReservoirDownsampler` against GATK 4.6.2.0.
//!
//! Two rows per case: which reads the reservoir kept, and where the shared generator was left.
//!
//! The second row is the one that matters. `Utils.getRandomGenerator()` is a single static stream,
//! and the downsampler draws for every read past the target *before* deciding whether the slot is
//! inside the reservoir. So a discarded read still advances the stream, and a port that skipped
//! the draw when the slot loses would keep the same reads and leave the generator somewhere else.
//! The golden shows it directly: `under-target` and `nonrandom` both leave the stream at
//! `1057280359`, its first value, because neither took a draw, while `one-over` leaves it at
//! `-873351126`, its second, because it took exactly one.

use gatk_corpus as corpus;
use gatk_engine::downsampling::{ReservoirDownsampler, SlotSource};
use gatk_engine::java_random::JavaRandom;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reservoir.txt.gz"),
    )
}

/// (read count, target, non-random mode) per labelled case.
fn configuration(label: &str) -> (usize, usize, bool) {
    match label {
        "under-target" => (3, 10, false),
        "at-target" => (10, 10, false),
        "one-over" => (11, 10, false),
        "many-over" => (50, 10, false),
        "very-many-over" => (500, 10, false),
        "target-one" => (20, 1, false),
        "nonrandom" => (50, 10, true),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_reservoir_keeps_what_the_reference_keeps() {
    let text = golden();

    let mut keeps: std::collections::HashMap<String, String> = Default::default();
    let mut stats: std::collections::HashMap<String, String> = Default::default();
    let mut labels: Vec<String> = Vec::new();
    let mut refusal: Option<String> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("keep\t") {
            let (label, names) = rest.split_once('\t').expect("a label and names");
            labels.push(label.to_string());
            keeps.insert(label.to_string(), names.to_string());
        } else if let Some(rest) = line.strip_prefix("stats\t") {
            let (label, rest) = rest.split_once('\t').expect("a label and stats");
            stats.insert(label.to_string(), rest.to_string());
        } else if let Some(rest) = line.strip_prefix("error\t") {
            refusal = rest.split_once('\t').map(|(_, class)| class.to_string());
        }
    }
    assert!(!labels.is_empty(), "the golden carries no keep rows");

    for label in &labels {
        let (count, target, non_random) = configuration(label);
        // The harness resets the shared stream before each case, so the port starts each one at
        // the GATK seed too. Measuring otherwise would measure the order of the cases.
        let mut random = JavaRandom::gatk();
        let names: Vec<String> = (0..count).map(|i| format!("r{i:03}")).collect();

        let mut downsampler: ReservoirDownsampler<String> = ReservoirDownsampler::new(target);
        downsampler.set_non_random_replacement_mode(non_random);
        for name in &names {
            let mut slots = if non_random {
                SlotSource::NonRandom
            } else {
                SlotSource::Random(&mut random)
            };
            downsampler.submit(name, name, &mut slots);
        }
        downsampler.signal_end_of_input();
        let discarded = downsampler.discarded();
        let kept = downsampler.consume_finalized_items();

        let ours = kept
            .iter()
            .map(|name| name.as_str())
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(&ours, &keeps[label], "{label}: the reads kept");

        // The generator's next value, which is what shows the discarded draws happened.
        let expected = &stats[label];
        let ours_stats = format!("{}\t{}\t{}", kept.len(), discarded, random.next_int());
        assert_eq!(
            &ours_stats, expected,
            "{label}: size, discarded, stream position"
        );
    }

    assert_eq!(
        refusal.as_deref(),
        Some("java.lang.IllegalArgumentException"),
        "a target of zero is refused"
    );
    println!(
        "{} reservoirs, their reads and the stream position after each, all identical",
        labels.len()
    );
}
