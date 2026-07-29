//! Conformance for `java.util.Random` against the JDK the oracle pins.
//!
//! Implemented from the published contract rather than from OpenJDK, which is what makes it
//! legitimate: the Javadoc states the algorithm and its constants and requires every method to
//! produce exactly this sequence. The contrast with `java_hash`, whose order is documented as
//! unspecified and had to be measured instead, is recorded in
//! `docs/an-unspecified-order-that-reaches-the-output.md`.
//!
//! The rows that matter most are the last of each seed: one generator with the methods mixed. A
//! port whose methods are individually right but whose *draw counts* are not passes every
//! single-method row and fails that one, because `nextDouble` and `nextLong` each take two draws.

use gatk_corpus as corpus;
use gatk_engine::java_random::JavaRandom;

const COUNT: usize = 24;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/java_random.txt.gz"),
    )
}

#[test]
fn every_sequence_matches_the_reference() {
    let text = golden();

    let mut sequences = 0;
    let mut interleaved = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("seq\t") {
            let mut parts = rest.split('\t');
            let seed: i64 = parts.next().expect("a seed").parse().expect("a number");
            let method = parts.next().expect("a method");
            let expected = parts.next().expect("the values");

            let mut random = JavaRandom::new(seed);
            let ours = match method {
                "nextInt" => (0..COUNT)
                    .map(|_| random.next_int().to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                "nextLong" => (0..COUNT)
                    .map(|_| random.next_long().to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                "nextBoolean" => (0..COUNT)
                    .map(|_| if random.next_boolean() { '1' } else { '0' })
                    .collect::<String>(),
                // Raw bits, so the comparison is of the value and not of `Double.toString`.
                "nextDouble" => (0..COUNT)
                    .map(|_| (random.next_double().to_bits() as i64).to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                "nextFloat" => (0..COUNT)
                    .map(|_| (random.next_float().to_bits() as i32).to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                bounded => {
                    let bound: i32 = bounded
                        .strip_prefix("nextInt(")
                        .and_then(|r| r.strip_suffix(')'))
                        .expect("a bounded method")
                        .parse()
                        .expect("a bound");
                    (0..COUNT)
                        .map(|_| random.next_int_bound(bound).to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                }
            };
            assert_eq!(ours, expected, "seed {seed}, {method}");
            sequences += 1;
        } else if let Some(rest) = line.strip_prefix("interleave\t") {
            let (seed, expected) = rest.split_once('\t').expect("a seed and the values");
            let seed: i64 = seed.parse().expect("a number");
            let mut random = JavaRandom::new(seed);
            let ours = (0..COUNT)
                .map(|i| match i % 5 {
                    0 => format!("i{}", random.next_int()),
                    1 => format!("d{}", random.next_double().to_bits() as i64),
                    2 => format!("b{}", u8::from(random.next_boolean())),
                    3 => format!("l{}", random.next_long()),
                    _ => format!("n{}", random.next_int_bound(37)),
                })
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(ours, expected, "seed {seed}, interleaved");
            interleaved += 1;
        }
    }

    assert!(sequences > 0, "the golden carries no sequence rows");
    println!("{sequences} sequences and {interleaved} interleaved streams, all identical");
}

/// The seed GATK pins, checked on its own so a change to it is a failing test rather than a
/// silently different stream.
#[test]
fn the_gatk_generator_starts_where_the_reference_starts() {
    let text = golden();
    let expected = text
        .lines()
        .find_map(|line| line.strip_prefix("seq\t47382911\tnextInt\t"))
        .expect("the golden carries the GATK seed");
    let mut random = JavaRandom::gatk();
    let ours = (0..COUNT)
        .map(|_| random.next_int().to_string())
        .collect::<Vec<_>>()
        .join(",");
    assert_eq!(ours, expected, "Utils.getRandomGenerator()'s first draws");
}
