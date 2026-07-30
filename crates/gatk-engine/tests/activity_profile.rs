//! Conformance for the activity profile against the oracle.
//!
//! Goldens from `tools/readfilter-conformance/ActivityProfileDump.java`.
//!
//! The row that states the whole point:
//!
//! ```text
//! region  plain-active-then-inactive     0  chr1:100-119  true
//! region  plain-active-then-inactive     1  chr1:120-139  false
//! region  bandpass-active-then-inactive  0  chr1:100-139  true
//! ```
//!
//! Twenty active positions followed by twenty zeros give two regions without the band pass filter
//! and **one** with it: the Gaussian tails of the active sites are added onto the zeros and lift
//! them over the threshold. The input is identical; the filter changes which bases get assembled.
//!
//! The kernels are compared as raw bits. `Math.exp` is the platform's, and an `exp` that differs in
//! the last ulp moves a cut site, so comparing rounded values would hide exactly the failure this
//! suite exists to catch.

use gatk_corpus as corpus;
use gatk_engine::activity_profile::{make_kernel, ActivityProfile};

/// Label, sigma (`None` for a plain profile), threshold, max region size, and the blocks of
/// probabilities, as `(count, value)`. Same order as the dump.
type Case = (
    &'static str,
    Option<f64>,
    f64,
    usize,
    &'static [(usize, f64)],
);

const CASES: &[Case] = &[
    (
        "plain-active-then-inactive",
        None,
        0.002,
        50,
        &[(20, 0.9), (20, 0.0)],
    ),
    (
        "bandpass-active-then-inactive",
        Some(17.0),
        0.002,
        50,
        &[(20, 0.9), (20, 0.0)],
    ),
    (
        "long-active",
        None,
        0.002,
        10,
        &[(5, 0.9), (1, 0.1), (5, 0.9), (1, 0.05), (30, 0.9)],
    ),
    ("plateau", None, 0.002, 10, &[(4, 0.9), (4, 0.1), (30, 0.9)]),
    ("all-inactive", None, 0.002, 20, &[(60, 0.0)]),
    ("too-short", None, 0.002, 50, &[(5, 0.9)]),
    ("bandpass-zeros", Some(17.0), 0.002, 20, &[(40, 0.0)]),
    ("single-state", None, 0.002, 20, &[(1, 0.9)]),
];

const CONTIG_LENGTH: i32 = 10000;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/activity_profile.txt.gz"),
    )
}

fn bits(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| (value.to_bits() as i64).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[test]
fn every_kernel_matches_the_reference_bit_for_bit() {
    let text = golden();
    let mut rows = 0;

    for line in text.lines() {
        let Some(rest) = line.strip_prefix("kernel\t") else {
            continue;
        };
        let mut parts = rest.split('\t');
        let max_filter_size: i32 = parts.next().expect("a size").parse().expect("a number");
        let sigma: f64 = parts.next().expect("a sigma").parse().expect("a number");
        let filter_size: i32 = parts
            .next()
            .expect("a filter size")
            .parse()
            .expect("a number");
        let band_size: usize = parts
            .next()
            .expect("a band size")
            .parse()
            .expect("a number");
        let expected = parts.next().expect("the kernel");

        let profile = ActivityProfile::band_pass(
            50,
            0.002,
            max_filter_size,
            sigma,
            true,
            "chr1",
            CONTIG_LENGTH,
        );
        assert_eq!(
            profile.filter_size(),
            filter_size,
            "filter size at maxFilterSize {max_filter_size}, sigma {sigma}"
        );
        let kernel = profile.kernel().expect("a band pass profile has a kernel");
        assert_eq!(kernel.len(), band_size, "band size");
        assert_eq!(
            bits(kernel),
            expected,
            "kernel at maxFilterSize {max_filter_size}, sigma {sigma}"
        );
        rows += 1;
    }

    assert!(rows > 0, "the golden carries no kernel rows");
    println!("{rows} kernels identical bit for bit");
}

#[test]
fn every_profile_produces_the_reference_regions() {
    let text = golden();

    for (label, sigma, threshold, max_region_size, blocks) in CASES {
        let mut profile = match sigma {
            None => ActivityProfile::new(50, *threshold, "chr1", CONTIG_LENGTH),
            Some(sigma) => {
                ActivityProfile::band_pass(50, *threshold, 50, *sigma, true, "chr1", CONTIG_LENGTH)
            }
        };

        let mut position = 100;
        let mut index = 0;
        // Popped after every add, as the walker does: the same probabilities popped only at the
        // end give different regions, because a pop consumes the states it used.
        for (count, value) in *blocks {
            for _ in 0..*count {
                profile.add(position, *value);
                position += 1;
                for region in profile.pop_ready_regions(0, 1, *max_region_size, false) {
                    check(&text, label, index, &region);
                    index += 1;
                }
            }
        }

        let expected_probs = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("probs\t{label}\t")))
            .unwrap_or_else(|| panic!("{label}: no probs row"));
        assert_eq!(
            bits(&profile.probabilities()),
            expected_probs,
            "{label}: the probabilities left in the profile"
        );

        for region in profile.pop_ready_regions(0, 1, *max_region_size, true) {
            check(&text, label, index, &region);
            index += 1;
        }

        let summary = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("summary\t{label}\t")))
            .unwrap_or_else(|| panic!("{label}: no summary row"));
        assert_eq!(
            format!(
                "{index}\t{}\t{}",
                profile.size(),
                profile.max_prob_propagation_distance()
            ),
            summary,
            "{label}: regions, states left, propagation distance"
        );
    }

    println!("{} profiles identical", CASES.len());
}

fn check(
    text: &str,
    label: &str,
    index: usize,
    region: &gatk_engine::activity_profile::PoppedRegion,
) {
    let expected = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("region\t{label}\t{index}\t")))
        .unwrap_or_else(|| panic!("{label}: no region row {index}"));
    let ours = format!(
        "{}:{}-{}\t{}",
        region.span.contig, region.span.start, region.span.end, region.is_active
    );
    assert_eq!(ours, expected, "{label}, region {index}");
}

/// The comparison that justifies the suite: identical input, one region or two, decided by the
/// filter alone.
#[test]
fn the_band_pass_filter_changes_which_bases_are_assembled() {
    let text = golden();
    let count = |label: &str| -> usize {
        text.lines()
            .filter(|line| line.starts_with(&format!("region\t{label}\t")))
            .count()
    };
    assert_eq!(count("plain-active-then-inactive"), 2);
    assert_eq!(count("bandpass-active-then-inactive"), 1);

    // And the one region covers the positions that were added as exactly zero.
    let merged = text
        .lines()
        .find(|line| line.starts_with("region\tbandpass-active-then-inactive\t0\t"))
        .expect("the golden carries it");
    assert!(merged.contains("chr1:100-139"), "{merged}");
}

/// A kernel built at the adaptive width is not the full kernel truncated: the normalisation
/// divides by a different sum, so every value differs.
#[test]
fn the_adaptive_kernel_is_rebuilt_rather_than_truncated() {
    let full = make_kernel(50, 1.0);
    let narrow = make_kernel(3, 1.0);
    let middle_of_full = full[50];
    let middle_of_narrow = narrow[3];
    assert_ne!(
        middle_of_full.to_bits(),
        middle_of_narrow.to_bits(),
        "the two kernels agree at the centre, so the normalisation is not being redone"
    );
}
