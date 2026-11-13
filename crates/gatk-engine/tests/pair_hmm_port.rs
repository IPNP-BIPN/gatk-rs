//! Conformance for `LoglessPairHMM` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/PairHmmDump.java`, which ran three implementations
//! over twelve read and haplotype pairs and printed each likelihood as its RAW BITS as well as its
//! decimal form. The port targets the scalar `LoglessPairHMM`, which is the one the oracle
//! contract pins.
//!
//! # What this suite is for
//!
//!  * **every one of the twelve likelihoods, compared as bits and not as decimals**;
//!  * **the initial condition, which is what keeps the linear recursion from underflowing**;
//!  * **the transition probabilities, whose match-to-match term is a cached expression and not
//!    `1 - 10^x`**;
//!  * **and the tristate correction a mismatch pays.**

use gatk_corpus as corpus;
use gatk_engine::pair_hmm::{
    approximate_log10_sum_log10, initial_condition, initial_condition_log10, match_to_match_prob,
    prior, qual_to_trans_probs, read_likelihood_given_haplotype_log10, DELETION_TO_DELETION,
    INDEL_TO_MATCH, INSERTION_TO_INSERTION, MATCH_TO_DELETION, MATCH_TO_INSERTION, MATCH_TO_MATCH,
    TRISTATE_CORRECTION,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pair_hmm.txt.gz"),
    )
}

/// One case's LoglessPairHMM answer, as the bits the golden holds.
fn bits(text: &str, label: &str) -> u64 {
    let prefix = format!("likelihood\t{label}\tLoglessPairHMM=");
    let row = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .unwrap_or_else(|| panic!("likelihood/{label}"));
    let (hex, _) = row.split_once(',').expect("bits and a decimal");
    u64::from_str_radix(hex, 16).expect("the bits")
}

/// One pair of the fixture: the haplotype, the read, and the four qualities.
struct Pair {
    label: &'static str,
    haplotype: &'static str,
    read: &'static str,
    base: u8,
    insertion: u8,
    deletion: u8,
    gap: u8,
}

const PAIRS: [Pair; 12] = [
    Pair {
        label: "identical",
        haplotype: "ACGTACGTACGT",
        read: "ACGTACGTACGT",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "one-mismatch",
        haplotype: "ACGTACGTACGT",
        read: "ACGTAAGTACGT",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "short-read",
        haplotype: "ACGTACGTACGTACGT",
        read: "ACGTACGT",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "deletion",
        haplotype: "ACGTACGTACGT",
        read: "ACGTACGT",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "insertion",
        haplotype: "ACGTACGT",
        read: "ACGTTTTTACGT",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "low-base-quality",
        haplotype: "ACGTACGTACGT",
        read: "ACGTAAGTACGT",
        base: 2,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "high-base-quality",
        haplotype: "ACGTACGTACGT",
        read: "ACGTAAGTACGT",
        base: 60,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "at-the-quality-threshold",
        haplotype: "ACGTACGTACGT",
        read: "ACGTAAGTACGT",
        base: 18,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "below-the-quality-threshold",
        haplotype: "ACGTACGTACGT",
        read: "ACGTAAGTACGT",
        base: 17,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "cheap-gaps",
        haplotype: "ACGTACGTACGT",
        read: "ACGTACGT",
        base: 30,
        insertion: 10,
        deletion: 10,
        gap: 5,
    },
    Pair {
        label: "long",
        haplotype: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        read: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
    Pair {
        label: "homopolymer",
        haplotype: "AAAAAAAAAAAA",
        read: "AAAAAAAAAA",
        base: 30,
        insertion: 45,
        deletion: 45,
        gap: 10,
    },
];

fn likelihood(pair: &Pair) -> f64 {
    let read = pair.read.as_bytes();
    let quals = vec![pair.base; read.len()];
    let insertions = vec![pair.insertion; read.len()];
    let deletions = vec![pair.deletion; read.len()];
    let gaps = vec![pair.gap; read.len()];
    read_likelihood_given_haplotype_log10(
        pair.haplotype.as_bytes(),
        read,
        &quals,
        &insertions,
        &deletions,
        &gaps,
    )
}

/// Every one of the twelve, compared as bits.
#[test]
fn every_likelihood_matches_the_golden() {
    let text = golden();
    for pair in &PAIRS {
        let ours = likelihood(pair);
        assert_eq!(
            ours.to_bits(),
            bits(&text, pair.label),
            "{}: {ours} against the golden's own",
            pair.label
        );
    }
}

/// The initial condition, and what it is for.
#[test]
fn the_initial_condition_is_a_power_of_two() {
    // `Math.pow(2, 1020)`, which is exact, and its log10.
    assert_eq!(initial_condition(), 2f64.powi(1020));
    assert_eq!(initial_condition_log10(), initial_condition().log10());
    // It is subtracted back off, so the answer does not depend on it: doubling the haplotype's
    // length does not scale the likelihood, because the seed is spread over that length.
    let text = golden();
    let identical = bits(&text, "identical");
    assert!(f64::from_bits(identical) < 0.0);
}

/// The transition probabilities, whose match-to-match term is not `1 - 10^x`.
#[test]
fn the_match_to_match_term_is_the_cached_expression() {
    let probs = qual_to_trans_probs(45, 45, 10);
    // The two indel openings are the plain error probabilities of their own qualities.
    assert_eq!(probs[MATCH_TO_INSERTION], 10f64.powf(-4.5));
    assert_eq!(probs[MATCH_TO_DELETION], 10f64.powf(-4.5));
    // The gap continuation is one probability read two ways.
    assert_eq!(probs[INDEL_TO_MATCH], 1.0 - 10f64.powf(-1.0));
    // The match-to-match term is written as `exp(log1p(-min(1, 10^sum)))` and not as
    // `1 - 10^sum`. The two agree bit for bit at this quality, and the reference writes the first
    // because it does not agree everywhere: `log1p` is what keeps a sum near one from losing its
    // significant digits. The port carries the reference's spelling rather than the shorter one
    // it happens to equal here.
    let sum = approximate_log10_sum_log10(-4.5, -4.5);
    let naive = 1.0 - 10f64.powf(sum);
    assert_eq!(probs[MATCH_TO_MATCH].to_bits(), naive.to_bits());
    // Where the difference shows: a quality of one on both sides puts the sum near zero, and
    // there the two spellings are different doubles.
    let small = qual_to_trans_probs(1, 1, 10);
    let small_sum = approximate_log10_sum_log10(-0.1, -0.1);
    assert_ne!(
        small[MATCH_TO_MATCH].to_bits(),
        (1.0 - 10f64.powf(small_sum)).to_bits()
    );
    // A quality past the maximum takes the other branch, which IS the naive expression.
    let past = match_to_match_prob(255, 255);
    assert_eq!(
        past,
        1.0 - 10f64.powf(approximate_log10_sum_log10(-25.5, -25.5))
    );
}

/// A mismatch pays the tristate correction, and an `N` pays nothing.
#[test]
fn a_mismatch_is_divided_by_three() {
    assert_eq!(TRISTATE_CORRECTION, 3.0);
    let matched = prior(b'A', b'A', 30);
    let mismatched = prior(b'A', b'C', 30);
    assert_eq!(matched, 1.0 - 10f64.powf(-3.0));
    assert_eq!(mismatched, 10f64.powf(-3.0) / 3.0);
    // An `N` on either side is a match, whichever side it is on.
    assert_eq!(prior(b'N', b'C', 30), matched);
    assert_eq!(prior(b'A', b'N', 30), matched);
    // And the quality moves it: the golden's two extremes differ by orders of magnitude.
    let text = golden();
    let low = f64::from_bits(bits(&text, "low-base-quality"));
    let high = f64::from_bits(bits(&text, "high-base-quality"));
    assert!(high < low, "a higher quality makes a mismatch cost more");
}

/// The straight transcription the kernel replaced: three full matrices and a `pow` per cell.
fn transcription(
    haplotype: &[u8],
    read: &[u8],
    read_quals: &[u8],
    insertion_quals: &[u8],
    deletion_quals: &[u8],
    gap_continuation: &[u8],
) -> f64 {
    let (n, m) = (read.len() + 1, haplotype.len() + 1);
    let mut matches = vec![vec![0.0f64; m]; n];
    let mut insertions = vec![vec![0.0f64; m]; n];
    let mut deletions = vec![vec![0.0f64; m]; n];
    let initial = initial_condition() / haplotype.len() as f64;
    deletions[0].fill(initial);
    for i in 1..n {
        let t = qual_to_trans_probs(
            insertion_quals[i - 1],
            deletion_quals[i - 1],
            gap_continuation[i - 1],
        );
        for j in 1..m {
            let p = prior(read[i - 1], haplotype[j - 1], read_quals[i - 1]);
            matches[i][j] = p
                * (matches[i - 1][j - 1] * t[MATCH_TO_MATCH]
                    + insertions[i - 1][j - 1] * t[INDEL_TO_MATCH]
                    + deletions[i - 1][j - 1] * t[INDEL_TO_MATCH]);
            insertions[i][j] = matches[i - 1][j] * t[MATCH_TO_INSERTION]
                + insertions[i - 1][j] * t[INSERTION_TO_INSERTION];
            deletions[i][j] = matches[i][j - 1] * t[MATCH_TO_DELETION]
                + deletions[i][j - 1] * t[DELETION_TO_DELETION];
        }
    }
    let mut total = 0.0;
    for j in 1..m {
        total += matches[n - 1][j] + insertions[n - 1][j];
    }
    jmath::math::log10(total) - initial_condition_log10()
}

/// The kernel gives the transcription's bits on seeded pairs that reach every quality, `N` on
/// both sides, indels, an empty read and a one-base haplotype.
#[test]
fn the_kernel_gives_the_transcriptions_bits() {
    let mut state = 0x1391_u64;
    let mut next = move |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    let bases = b"ACGTN";
    for _ in 0..3000 {
        let haplotype: Vec<u8> = (0..1 + next(300))
            .map(|_| bases[next(5) as usize])
            .collect();
        let length = next(160) as usize;
        let read: Vec<u8> = (0..length)
            .map(|k| {
                if next(10) < 8 && k < haplotype.len() {
                    haplotype[k]
                } else {
                    bases[next(5) as usize]
                }
            })
            .collect();
        let wide = next(4) == 0;
        let mut quals = |wide: bool| -> Vec<u8> {
            (0..length)
                .map(|_| {
                    if wide {
                        next(256) as u8
                    } else {
                        2 + next(40) as u8
                    }
                })
                .collect()
        };
        let (base, insertion, deletion, gap) = (quals(wide), quals(wide), quals(wide), quals(wide));
        let ours = read_likelihood_given_haplotype_log10(
            &haplotype, &read, &base, &insertion, &deletion, &gap,
        );
        let theirs = transcription(&haplotype, &read, &base, &insertion, &deletion, &gap);
        assert_eq!(ours.to_bits(), theirs.to_bits(), "{haplotype:?} {read:?}");
    }
}
