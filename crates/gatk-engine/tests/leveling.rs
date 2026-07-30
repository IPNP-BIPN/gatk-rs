//! Conformance for `LevelingDownsampler` and the `nextPermutation` underneath it.
//!
//! Two layers, both measured, because both decide which items survive:
//!
//!  * `nextPermutation(n, k)` shuffles all `n` and keeps `k`, so its cost is `n - 1` draws whatever
//!    `k` is. The golden carries the generator's next value after each call, which is what turns
//!    "the right indices" into "the right indices at the right cost";
//!  * the leveling plan is arithmetic and only the selection is random, so a wrong plan shows as
//!    wrong *sizes* and a wrong selection shows as wrong *names* on the correct sizes.
//!
//! The cases are declared here and in `LevelingDownsamplerDump.java`, and the test refuses to pass
//! if either side carries a label the other does not: a case silently dropped from one side would
//! otherwise read as agreement.

use gatk_corpus as corpus;
use gatk_engine::downsampling::LevelingDownsampler;
use gatk_engine::permutation::{self, PermutationError};
use gatk_engine::well19937c::Well19937c;

/// The same table as the dump's, in the same order.
const CASES: &[(&str, &[usize], i64, usize)] = &[
    ("under-target", &[3, 3, 3], 20, 1),
    ("exactly-target", &[3, 3, 3], 9, 1),
    ("one-over", &[3, 3, 3], 8, 1),
    ("even-cut", &[10, 10, 10], 15, 1),
    ("uneven", &[1, 5, 20], 10, 1),
    ("floor-blocks", &[1, 1, 20], 5, 1),
    ("minimum-blocks", &[4, 4, 4], 3, 3),
    ("one-stack", &[25], 4, 1),
    ("empty-among-others", &[0, 6, 6], 5, 1),
    ("no-stacks", &[], 5, 1),
    ("target-zero-min-one", &[4, 4], 0, 1),
    ("target-zero-min-zero", &[4, 4], 0, 0),
];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/leveling.txt.gz"),
    )
}

/// `Utils.resetRandomGenerator()` for the data generator: the dump calls it before every case, so
/// each case measures the code under test rather than the order the cases ran in.
fn reset() -> Well19937c {
    Well19937c::gatk()
}

#[test]
fn every_permutation_matches_the_reference() {
    let text = golden();
    let mut rows = 0;

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("perm\t") {
            let mut parts = rest.split('\t');
            let n: i32 = parts.next().expect("n").parse().expect("a number");
            let k: i32 = parts.next().expect("k").parse().expect("a number");
            let expected = parts.next().expect("the indices");
            let expected_after = parts.next().expect("the stream position");

            let mut random = reset();
            let indices = permutation::next_permutation(n, k, &mut random)
                .unwrap_or_else(|e| panic!("nextPermutation({n}, {k}) refused: {e}"));
            let ours = indices
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(ours, expected, "nextPermutation({n}, {k})");
            // The cost, not just the answer: n - 1 draws whatever k is.
            assert_eq!(
                random.next_int().to_string(),
                expected_after,
                "stream position after nextPermutation({n}, {k})"
            );
            rows += 1;
        }
    }

    assert!(rows > 0, "the golden carries no permutation rows");
    println!("{rows} permutations, indices and stream positions all identical");
}

#[test]
fn the_permutation_refusals_match_the_reference() {
    let text = golden();
    let mut rows = 0;

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("permerror\t") {
            let mut parts = rest.split('\t');
            let n: i32 = parts.next().expect("n").parse().expect("a number");
            let k: i32 = parts.next().expect("k").parse().expect("a number");
            let expected = parts.next().expect("the class");

            let mut random = reset();
            let ours = match permutation::next_permutation(n, k, &mut random) {
                Ok(_) => "none".to_string(),
                Err(PermutationError::SizeExceedsN { .. }) => {
                    "org.apache.commons.math3.exception.NumberIsTooLargeException".to_string()
                }
                Err(PermutationError::NotStrictlyPositive { .. }) => {
                    "org.apache.commons.math3.exception.NotStrictlyPositiveException".to_string()
                }
            };
            assert_eq!(ours, expected, "nextPermutation({n}, {k})");
            rows += 1;
        }
    }

    assert!(rows > 0, "the golden carries no refusal rows");
    println!("{rows} refusals identical");
}

#[test]
fn every_leveling_matches_the_reference() {
    let text = golden();

    let mut seen: Vec<String> = Vec::new();
    for (label, sizes, target, minimum) in CASES {
        // Both list kinds, which take different removal paths upstream and must keep the same
        // items. This port has one path, so the two rows check the oracle's own agreement as much
        // as ours.
        for kind in ["linked", "array"] {
            let mut random = reset();
            let mut downsampler: LevelingDownsampler<String> =
                LevelingDownsampler::with_minimum(*target, *minimum);
            let mut next = 0;
            for size in *sizes {
                let mut stack = Vec::new();
                for _ in 0..*size {
                    stack.push(format!("s{next:02}"));
                    next += 1;
                }
                downsampler.submit(stack);
            }

            let outcome = downsampler.signal_end_of_input(&mut random);

            if let Err(error) = outcome {
                let expected = row(&text, &format!("levelerror\t{label}\t{kind}\t"))
                    .unwrap_or_else(|| {
                        panic!("{label}/{kind}: the port refused, the golden did not")
                    });
                let ours = match error {
                    PermutationError::SizeExceedsN { .. } => {
                        "org.apache.commons.math3.exception.NumberIsTooLargeException"
                    }
                    PermutationError::NotStrictlyPositive { .. } => {
                        "org.apache.commons.math3.exception.NotStrictlyPositiveException"
                    }
                };
                assert_eq!(ours, expected, "{label}/{kind}");
                seen.push(format!("{label}/{kind}"));
                continue;
            }

            let size = downsampler.size();
            let discarded = downsampler.discarded();
            let groups = downsampler.consume_finalized_items();
            let stacks = groups
                .iter()
                .map(|group| group.join(","))
                .collect::<Vec<_>>()
                .join("|");

            let expected = row(&text, &format!("level\t{label}\t{kind}\t"))
                .unwrap_or_else(|| panic!("{label}/{kind} is missing from the golden"));
            assert_eq!(stacks, expected, "{label}/{kind}: kept items");

            let stats = row(&text, &format!("levelstats\t{label}\t{kind}\t"))
                .unwrap_or_else(|| panic!("{label}/{kind} has no stats row"));
            let ours = format!("{size}\t{discarded}\t{}", random.next_int());
            assert_eq!(
                ours, stats,
                "{label}/{kind}: size, discarded, stream position"
            );
            seen.push(format!("{label}/{kind}"));
        }
    }

    // Neither side may carry a case the other does not.
    let in_golden = text
        .lines()
        .filter(|line| line.starts_with("level\t") || line.starts_with("levelerror\t"))
        .count();
    assert_eq!(
        in_golden,
        seen.len(),
        "the golden holds {in_golden} leveling cases and the test ran {}",
        seen.len()
    );
    println!("{} leveling cases identical", seen.len());
}

fn row(text: &str, prefix: &str) -> Option<String> {
    text.lines()
        .find_map(|line| line.strip_prefix(prefix))
        .map(str::to_string)
}
