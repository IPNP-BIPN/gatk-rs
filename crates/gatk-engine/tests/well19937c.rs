//! Conformance for commons-math3's `Well19937c` against the oracle.
//!
//! This is GATK's *second* static generator. `Utils.getRandomDataGenerator()` is a
//! `RandomDataGenerator` over `new Well19937c(47382911L)`, and `LevelingDownsampler` draws from it
//! where `ReservoirDownsampler` draws from the `java.util.Random` in `java_random.rs`. Two
//! generators, two algorithms, two positions to keep straight.
//!
//! The rows worth naming:
//!
//!  * the `long` and `int` constructor rows for the same numeric seed must **differ**, because the
//!    long one seeds two pool words and the int one seeds a single word, which changes what the
//!    scrambler fills the remaining 622 with;
//!  * the `bytes` rows carry the generator's next value after the call, which is the only place the
//!    unconditional trailing draw is visible: at a length that is a multiple of four the bytes
//!    themselves are identical whether or not that draw is taken;
//!  * the `interleave` rows catch draw counts, since `nextDouble` and `nextLong` each take two.

use gatk_corpus as corpus;
use gatk_engine::well19937c::Well19937c;

const COUNT: usize = 24;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/well19937c.txt.gz"),
    )
}

/// Rebuild the generator the dump used, from the constructor it recorded.
fn make(ctor: &str, seed: &str) -> Well19937c {
    match ctor {
        "long" => Well19937c::from_long(seed.parse().expect("a long seed")),
        "int" => Well19937c::from_int(seed.parse().expect("an int seed")),
        "ints" => Well19937c::from_int_array(
            &seed
                .split(',')
                .map(|word| word.parse().expect("an int"))
                .collect::<Vec<i32>>(),
        ),
        other => panic!("unknown constructor {other}"),
    }
}

#[test]
fn every_sequence_matches_the_reference() {
    let text = golden();

    let mut sequences = 0;
    let mut byte_rows = 0;
    let mut interleaved = 0;

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("seq\t") {
            let mut parts = rest.split('\t');
            let ctor = parts.next().expect("a constructor");
            let seed = parts.next().expect("a seed");
            let method = parts.next().expect("a method");
            let expected = parts.next().expect("the values");

            let mut random = make(ctor, seed);
            let ours = match method {
                "nextInt" => join(&mut random, COUNT, |r| r.next_int().to_string()),
                "nextLong" => join(&mut random, COUNT, |r| r.next_long().to_string()),
                "nextBoolean" => (0..COUNT)
                    .map(|_| if random.next_boolean() { '1' } else { '0' })
                    .collect::<String>(),
                // Raw bits, so this compares the value rather than `Double.toString`.
                "nextDouble" => join(&mut random, COUNT, |r| {
                    (r.next_double().to_bits() as i64).to_string()
                }),
                "nextFloat" => join(&mut random, COUNT, |r| {
                    (r.next_float().to_bits() as i32).to_string()
                }),
                bounded if bounded.starts_with("nextInt(") => {
                    let bound: i32 = inner(bounded, "nextInt(").parse().expect("a bound");
                    join(&mut random, COUNT, |r| r.next_int_bound(bound).to_string())
                }
                bounded if bounded.starts_with("nextLong(") => {
                    let bound: i64 = inner(bounded, "nextLong(").parse().expect("a bound");
                    join(&mut random, COUNT, |r| r.next_long_bound(bound).to_string())
                }
                other => panic!("unknown method {other}"),
            };
            assert_eq!(ours, expected, "{ctor} seed {seed}, {method}");
            sequences += 1;
        } else if let Some(rest) = line.strip_prefix("bytes\t") {
            let mut parts = rest.split('\t');
            let ctor = parts.next().expect("a constructor");
            let seed = parts.next().expect("a seed");
            let length: usize = parts.next().expect("a length").parse().expect("a number");
            let expected_hex = parts.next().expect("the bytes");
            let expected_after = parts.next().expect("the stream position");

            let mut random = make(ctor, seed);
            let mut buffer = vec![0u8; length];
            random.next_bytes(&mut buffer);
            let hex: String = buffer.iter().map(|b| format!("{b:02x}")).collect();
            assert_eq!(hex, expected_hex, "{ctor} seed {seed}, nextBytes({length})");
            // Where the stream ended up. A port that skipped the trailing draw agrees on the line
            // above at every length divisible by four and fails here at all of them.
            assert_eq!(
                random.next_int().to_string(),
                expected_after,
                "{ctor} seed {seed}, stream after nextBytes({length})"
            );
            byte_rows += 1;
        } else if let Some(rest) = line.strip_prefix("interleave\t") {
            let mut parts = rest.split('\t');
            let ctor = parts.next().expect("a constructor");
            let seed = parts.next().expect("a seed");
            let expected = parts.next().expect("the values");

            let mut random = make(ctor, seed);
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
            assert_eq!(ours, expected, "{ctor} seed {seed}, interleaved");
            interleaved += 1;
        }
    }

    assert!(sequences > 0, "the golden carries no sequence rows");
    println!(
        "{sequences} sequences, {byte_rows} nextBytes lengths and {interleaved} interleaved \
         streams, all identical"
    );
}

/// The constructor GATK actually uses, checked on its own so a change to it is a failing test
/// rather than a silently different stream.
#[test]
fn the_gatk_data_generator_starts_where_the_reference_starts() {
    let text = golden();
    let expected = text
        .lines()
        .find_map(|line| line.strip_prefix("seq\tlong\t47382911\tnextInt\t"))
        .expect("the golden carries the GATK seed under the long constructor");
    let mut random = Well19937c::gatk();
    let ours = join(&mut random, COUNT, |r| r.next_int().to_string());
    assert_eq!(
        ours, expected,
        "Utils.getRandomDataGenerator()'s first draws"
    );
}

/// The two constructors must not agree for the same number. If they did, the seeding path would be
/// collapsing the two-word seed into one and the GATK stream would be the wrong one from its first
/// value.
#[test]
fn the_long_and_int_constructors_seed_different_pools() {
    let text = golden();
    let from_long = text
        .lines()
        .find_map(|line| line.strip_prefix("seq\tlong\t47382911\tnextInt\t"))
        .expect("the long-constructor row");
    let from_int = text
        .lines()
        .find_map(|line| line.strip_prefix("seq\tint\t47382911\tnextInt\t"))
        .expect("the int-constructor row");
    assert_ne!(
        from_long, from_int,
        "the oracle says the two constructors seed the same pool, which contradicts \
         AbstractWell's two-word split"
    );
}

fn join(
    random: &mut Well19937c,
    count: usize,
    mut draw: impl FnMut(&mut Well19937c) -> String,
) -> String {
    (0..count)
        .map(|_| draw(random))
        .collect::<Vec<_>>()
        .join(",")
}

fn inner<'a>(method: &'a str, prefix: &str) -> &'a str {
    method
        .strip_prefix(prefix)
        .and_then(|rest| rest.strip_suffix(')'))
        .expect("a bounded method")
}
