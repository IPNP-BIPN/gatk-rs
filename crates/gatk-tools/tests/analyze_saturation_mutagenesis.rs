//! Conformance for `AnalyzeSaturationMutagenesis` against GATK 4.6.2.0, compared as the three
//! reports of every configuration.
//!
//! Golden from `tools/readfilter-conformance/AnalyzeSaturationMutagenesisDump.java`.
//!
//! The alignment is not measured or ported. What is compared is the census the tool writes, the
//! arithmetic of its percentages, the trim, the flank test and the codon translation.
//!
//! # What this suite is for
//!
//!  * **the census being a tree whose denominators change from level to level**;
//!  * **one of those denominators being wrong, so a line reads over 100%**;
//!  * **the overlapping line counting reads where its rows count pairs**;
//!  * **the quality trim, and TLEN cutting it away when it is zero**;
//!  * **a variant at a read's edge failing the flank test**;
//!  * **the observation threshold deciding which variants are reported**;
//!  * **the codon translation being indexed in base-four over ACGT**;
//!  * **the two refusals**;
//!  * **and the counters being static, so two runs in one JVM add up.**

use gatk_corpus as corpus;
use gatk_tools::analyze_saturation_mutagenesis::{
    census_line, check_orf, check_translation, codon_value, decimal_format, fragment_trim,
    has_sufficient_flank, is_reported, orf_length, parse_orf, percentage, quality_trim, translate,
    Arguments, Interval, OrfInterval, ReportType, BASE_ORDER, DEFAULT_CODON_TRANSLATION,
    ORF_LENGTH_MESSAGE, TRANSLATION_LENGTH_MESSAGE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/analyze_saturation_mutagenesis.txt.gz"),
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

/// The three footer labels, whose denominator is the base calls rather than the reads.
const FOOTER_LABELS: [&str; 3] = [
    "Total base calls",
    "Base calls evaluated for variants",
    "Base calls unevaluated",
];

/// One census, as its lines split into depth, label, count and percentage.
#[derive(Debug, Clone, PartialEq)]
struct Line {
    depth: usize,
    label: String,
    count: u64,
    percentage: String,
}

fn census(text: &str, label: &str) -> Vec<Line> {
    section(text, "reads", label)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let depth = line.chars().take_while(|c| *c == '>').count();
            let columns: Vec<&str> = line[depth..].split('\t').collect();
            Line {
                depth,
                label: columns[0].trim_end_matches(':').to_string(),
                count: columns[1].parse().expect("a count"),
                percentage: columns[2].trim_end_matches('%').to_string(),
            }
        })
        .collect()
}

fn line<'a>(census: &'a [Line], label: &str) -> &'a Line {
    census
        .iter()
        .find(|line| line.label == label)
        .unwrap_or_else(|| panic!("the census carries {label}"))
}

/// The variant rows of one configuration.
fn variants(text: &str, label: &str) -> Vec<Vec<String>> {
    section(text, "variants", label)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// Every percentage of every census is the count over the denominator its level names.
#[test]
fn every_percentage_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "default",
        "min-obs-one",
        "flank-ten",
        "mapq-thirty",
        "min-length-hundred",
        "unpaired-mode",
        "combine-disjoint",
        "two-intervals",
    ] {
        let census = census(&text, label);
        let total = line(&census, "Total Reads").count;
        let evaluable = line(&census, "Evaluable Reads").count;
        // The three top lines are over every read.
        for name in ["Unmapped Reads", "LowQ Reads", "Evaluable Reads"] {
            let found = census
                .iter()
                .find(|line| line.depth == 1 && line.label == name)
                .unwrap_or_else(|| panic!("{label}: {name}"));
            assert_eq!(
                found.percentage,
                percentage(found.count, total),
                "{label}: {name}"
            );
            assert_eq!(
                census_line(1, name, found.count, total),
                format!(">{name}:\t{}\t{}%", found.count, found.percentage),
                "{label}: {name}"
            );
        }
        // The footer's two `>` lines are over the base calls instead.
        let base_calls = line(&census, "Total base calls").count;
        for name in [
            "Base calls evaluated for variants",
            "Base calls unevaluated",
        ] {
            let found = line(&census, name);
            assert_eq!(
                found.percentage,
                percentage(found.count, base_calls),
                "{label}: {name}"
            );
        }
        // Each category's own line is over the EVALUABLE reads.
        for found in census.iter().filter(|line| line.depth == 2) {
            assert_eq!(
                found.percentage,
                percentage(found.count, evaluable),
                "{label}: {}",
                found.label
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 8, "the configurations the port reproduces");
}

/// The unpaired line counts every read while its denominator counts the evaluable ones.
#[test]
fn one_denominator_lets_a_line_read_over_a_hundred_per_cent() {
    let text = golden();
    let census = census(&text, "unpaired-mode");
    let unpaired = line(&census, "Unpaired reads");
    let evaluable = line(&census, "Evaluable Reads");
    // Seventeen reads over sixteen evaluable ones.
    assert_eq!(unpaired.count, 17);
    assert_eq!(evaluable.count, 16);
    assert_eq!(unpaired.percentage, "106.250");
    assert_eq!(percentage(17, 16), "106.250");
    // The rows beneath it are over the category's own total, which IS seventeen.
    let rows: u64 = census
        .iter()
        .filter(|line| line.depth == 3)
        .map(|line| line.count)
        .sum();
    assert_eq!(rows, 17);
    // And the category's own `LowQ Reads` row is not the top-level one: the top says none.
    assert_eq!(line(&census, "LowQ Reads").count, 0);
    assert_eq!(
        census
            .iter()
            .filter(|line| line.depth == 3 && line.label == "LowQ Reads")
            .map(|line| line.count)
            .sum::<u64>(),
        1
    );
}

/// The line is in reads and the rows beneath it are in pairs.
#[test]
fn the_overlapping_line_counts_reads_and_its_rows_count_pairs() {
    let text = golden();
    let census = census(&text, "default");
    let overlapping = line(&census, "Reads in overlapping pairs evaluated together");
    assert_eq!(overlapping.count, 12);
    // Its rows sum to six, which is half of twelve: they are pairs.
    let start = census
        .iter()
        .position(|line| line.label == "Reads in overlapping pairs evaluated together")
        .expect("the category");
    let rows: u64 = census[start + 1..]
        .iter()
        .take_while(|line| line.depth == 3)
        .map(|line| line.count)
        .sum();
    assert_eq!(rows * 2, overlapping.count);
    // The disjoint category is in reads at both levels, so the two are not written alike.
    let disjoint = line(&census, "Reads in disjoint pairs evaluated separately");
    let start = census
        .iter()
        .position(|line| line.label == "Reads in disjoint pairs evaluated separately")
        .expect("the category");
    let rows: u64 = census[start + 1..]
        .iter()
        .take_while(|line| line.depth == 3)
        .map(|line| line.count)
        .sum();
    assert_eq!(rows, disjoint.count);
}

/// The unpaired line is written only when it has reads; the other two always are.
#[test]
fn two_of_the_three_categories_are_always_written() {
    let text = golden();
    let census = census(&text, "unpaired-mode");
    // In unpaired mode the two pair categories are empty and still written.
    assert_eq!(
        line(&census, "Reads in disjoint pairs evaluated separately").count,
        0
    );
    assert_eq!(
        line(&census, "Reads in overlapping pairs evaluated together").count,
        0
    );
    // And a row whose count is zero is left out of every category.
    assert!(census
        .iter()
        .filter(|line| line.depth == 3)
        .all(|line| line.count != 0));
}

/// The first and last runs of `min_length` high-quality bases.
#[test]
fn the_trim_is_the_first_and_last_high_quality_runs() {
    let arguments = Arguments::default();
    // A read that is high quality throughout is not trimmed at all.
    let good = vec![40u8; 120];
    assert_eq!(
        quality_trim(&good, &arguments),
        Interval { start: 0, end: 120 }
    );
    // A read with a low-quality head is trimmed to after it.
    let mut headed = vec![2u8; 20];
    headed.extend(vec![40u8; 100]);
    assert_eq!(
        quality_trim(&headed, &arguments),
        Interval {
            start: 20,
            end: 120
        }
    );
    // A read with no run long enough at all yields the null interval.
    let noisy: Vec<u8> = (0..120).map(|i| if i % 5 == 0 { 2 } else { 40 }).collect();
    assert_eq!(quality_trim(&noisy, &arguments), Interval::NULL);
    assert_eq!(Interval::NULL.size(), 0);
    // A longer minimum takes more away, which is the run the golden's own configuration made.
    let long = Arguments {
        min_length: 100,
        ..Arguments::default()
    };
    assert_eq!(
        quality_trim(&headed, &long),
        Interval {
            start: 20,
            end: 120
        }
    );
    assert!(quality_trim(&headed, &long).size() >= long.min_length);
}

/// A fragment length of zero cuts the whole trim away.
#[test]
fn a_zero_fragment_length_trims_everything() {
    let arguments = Arguments::default();
    let whole = Interval { start: 0, end: 120 };
    // Not properly paired: the trim survives whatever the fragment length says.
    assert_eq!(
        fragment_trim(whole, 120, false, false, 0, &arguments),
        whole
    );
    // Properly paired with a real fragment: the trim survives too.
    assert_eq!(
        fragment_trim(whole, 120, true, false, 219, &arguments),
        whole
    );
    // Properly paired with no fragment length: nothing is left.
    let cut = fragment_trim(whole, 120, true, false, 0, &arguments);
    assert_eq!(cut.size(), 0);
    // Which is why the golden's reads carry a TLEN: with one, none of them is LOW_QUALITY.
    let text = golden();
    assert_eq!(
        census(&text, "default")
            .iter()
            .find(|line| line.depth == 1 && line.label == "LowQ Reads")
            .expect("the line")
            .count,
        0
    );
}

/// A read whose variant sits at its own edge is counted `Insufficient flank`.
#[test]
fn a_variant_at_a_reads_edge_has_no_flank() {
    let text = golden();
    let arguments = Arguments::default();
    let coverage = Interval { start: 0, end: 120 };
    // At the very first base there is nothing before it.
    assert!(!has_sufficient_flank(0, coverage, &arguments));
    assert!(!has_sufficient_flank(1, coverage, &arguments));
    assert!(has_sufficient_flank(2, coverage, &arguments));
    // And at the very last.
    assert!(!has_sufficient_flank(119, coverage, &arguments));
    assert!(has_sufficient_flank(117, coverage, &arguments));
    // Which is the row the census carries for that read.
    let plain_census = census(&text, "default");
    assert!(plain_census
        .iter()
        .any(|line| line.label == "Insufficient flank" && line.count > 0));
    // A wider flank takes more variants away.
    let wide = Arguments {
        min_flanking_length: 10,
        ..Arguments::default()
    };
    assert!(!has_sufficient_flank(5, coverage, &wide));
    assert!(has_sufficient_flank(10, coverage, &wide));
    let widened: u64 = census(&text, "flank-ten")
        .iter()
        .filter(|line| line.label == "Insufficient flank")
        .map(|line| line.count)
        .sum();
    let plain: u64 = plain_census
        .iter()
        .filter(|line| line.label == "Insufficient flank")
        .map(|line| line.count)
        .sum();
    assert!(widened >= plain);
}

/// The threshold decides how many rows the variant report has.
#[test]
fn the_observation_threshold_decides_which_variants_are_reported() {
    let text = golden();
    let default = variants(&text, "default");
    let one = variants(&text, "min-obs-one");
    assert_eq!(default.len(), 1);
    assert_eq!(one.len(), 2);
    // The row the default keeps is the one seen three times, which is the threshold.
    assert_eq!(default[0][0], "3");
    assert!(is_reported(3, &Arguments::default()));
    assert!(!is_reported(2, &Arguments::default()));
    // The extra row of the lower threshold was seen once.
    let extra = one.iter().find(|row| row[0] == "1").expect("the rare row");
    assert!(is_reported(
        1,
        &Arguments {
            min_variant_observations: 1,
            ..Arguments::default()
        }
    ));
    // Both rows name a base change, a codon change and an amino-acid change.
    for row in [&default[0], extra] {
        assert!(row[4].contains(">"), "{row:?}");
        assert!(row[6].contains(">"), "{row:?}");
        assert!(row.last().expect("a name").len() >= 3, "{row:?}");
    }
    assert_eq!(default[0][4], "61:A>T");
    assert_eq!(default[0][6], "21:ACG>TCG");
    assert_eq!(default[0].last().expect("a name"), "T21S");
}

/// Three bases in base-four over `ACGT`.
#[test]
fn a_codon_is_indexed_in_base_four() {
    assert_eq!(BASE_ORDER, *b"ACGT");
    assert_eq!(codon_value(b"AAA"), Some(0));
    assert_eq!(codon_value(b"AAC"), Some(1));
    assert_eq!(codon_value(b"TTT"), Some(63));
    // `A` is 0, `C` is 1 and `G` is 2, so the codon is 0*16 + 1*4 + 2.
    assert_eq!(codon_value(b"ACG"), Some(6));
    // A base that is not one of the four has no codon at all.
    assert_eq!(codon_value(b"ANG"), None);
    assert_eq!(codon_value(b"AC"), None);
    // The default table is sixty-four codes, and translates the golden's own codons.
    assert_eq!(DEFAULT_CODON_TRANSLATION.chars().count(), 64);
    assert!(check_translation(DEFAULT_CODON_TRANSLATION).is_ok());
    assert_eq!(translate(b"ACG", DEFAULT_CODON_TRANSLATION), Some('T'));
    assert_eq!(translate(b"TCG", DEFAULT_CODON_TRANSLATION), Some('S'));
    // Which is the change the golden's variant row names.
    let text = golden();
    assert_eq!(
        variants(&text, "default")[0].last().expect("a name"),
        "T21S"
    );
    // A table of the wrong length is refused by name.
    assert_eq!(
        check_translation("KNKN").expect_err("too short"),
        TRANSLATION_LENGTH_MESSAGE
    );
}

/// One-based and inclusive, and its total length must divide by three.
#[test]
fn the_orf_is_one_based_and_inclusive() {
    assert_eq!(
        parse_orf("1-300"),
        Some(vec![OrfInterval { start: 1, end: 300 }])
    );
    assert_eq!(orf_length(&[OrfInterval { start: 1, end: 300 }]), 300);
    // Two intervals are spliced, so the total is what has to divide.
    let spliced = parse_orf("1-147,151-300").expect("two intervals");
    assert_eq!(orf_length(&spliced), 147 + 150);
    assert_eq!(orf_length(&spliced) % 3, 0);
    assert!(check_orf(&spliced, 300).is_ok());
    // A splice whose parts do NOT each divide by three still works, the total being what is
    // checked; the golden's own two-interval run is the case where both happen to.
    let uneven = parse_orf("1-146,151-301").expect("two intervals");
    assert_ne!(146 % 3, 0);
    assert_ne!((301 - 151 + 1) % 3, 0);
    assert_eq!(orf_length(&uneven) % 3, 0);
    assert!(check_orf(&uneven, 301).is_ok());
    // A length that does not divide is refused, and so is an end past the reference.
    let odd = parse_orf("1-299").expect("one interval");
    assert_eq!(
        check_orf(&odd, 300).expect_err("not a codon"),
        ORF_LENGTH_MESSAGE
    );
    let past = parse_orf("1-150,154-303").expect("two intervals");
    assert!(check_orf(&past, 300)
        .expect_err("past the end")
        .contains("larger than reference length"));
}

/// Both refusals reach the console as the tool's own words, with no exception class.
#[test]
fn the_two_refusals_are_the_goldens() {
    let text = golden();
    let message = |label: &str| -> String {
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
            .unwrap_or_else(|| panic!("the golden carries error/{label}"))
            .to_string()
    };
    assert_eq!(
        message("orf-not-a-codon"),
        format!("A USER ERROR has occurred: {ORF_LENGTH_MESSAGE}")
    );
    assert_eq!(
        message("short-translation"),
        format!("A USER ERROR has occurred: {TRANSLATION_LENGTH_MESSAGE}")
    );
    // No exception class reaches the caller at all.
    for label in ["orf-not-a-codon", "short-translation"] {
        assert!(!message(label).contains("Exception"), "{label}");
    }
}

/// Two invocations in one JVM add up, the counters being static fields.
#[test]
fn the_counters_are_static_and_never_reset() {
    let text = golden();
    let once = census(&text, "once-in-one-jvm");
    let twice = census(&text, "twice-in-one-jvm");
    assert_eq!(line(&once, "Total Reads").count, 17);
    assert_eq!(line(&twice, "Total Reads").count, 34);
    // The READ counters doubled, so their percentages are unchanged.
    for (once, twice) in once.iter().zip(twice.iter()) {
        assert_eq!(once.label, twice.label);
        if FOOTER_LABELS.contains(&once.label.as_str()) {
            continue;
        }
        assert_eq!(once.count * 2, twice.count, "{}", once.label);
        assert_eq!(once.percentage, twice.percentage, "{}", once.label);
    }
    // NOT EVERY COUNTER ACCUMULATES. The base calls do, being a static long; the coverage does
    // not, being read off the reference that each run rebuilds. So the second run reports twice
    // the base calls over the same coverage, and its percentage comes out half the first's.
    assert_eq!(
        line(&once, "Total base calls").count * 2,
        line(&twice, "Total base calls").count
    );
    assert_eq!(
        line(&once, "Base calls evaluated for variants").count,
        line(&twice, "Base calls evaluated for variants").count
    );
    assert_eq!(
        line(&once, "Base calls evaluated for variants").percentage,
        "66.406"
    );
    assert_eq!(
        line(&twice, "Base calls evaluated for variants").percentage,
        "33.203"
    );
    // And the run that shares no JVM with anything reports the same seventeen.
    assert_eq!(line(&census(&text, "min-obs-one"), "Total Reads").count, 17);
}

/// `DecimalFormat("0.000")`, which rounds half to even.
#[test]
fn the_percentages_are_rounded_half_to_even() {
    assert_eq!(decimal_format(0.0), "0.000");
    assert_eq!(decimal_format(100.0), "100.000");
    assert_eq!(percentage(1, 17), "5.882");
    assert_eq!(percentage(16, 17), "94.118");
    assert_eq!(percentage(17, 16), "106.250");
    assert_eq!(percentage(1, 3), "33.333");
    assert_eq!(percentage(2, 3), "66.667");
    // Every percentage the golden wrote is reproduced.
    let text = golden();
    for label in ["default", "unpaired-mode", "min-length-hundred"] {
        let census = census(&text, label);
        let total = line(&census, "Total Reads").count;
        for found in census
            .iter()
            .filter(|line| line.depth == 1 && !FOOTER_LABELS.contains(&line.label.as_str()))
        {
            assert_eq!(found.percentage, percentage(found.count, total), "{label}");
        }
    }
    // The labels are the tool's own and not the enum's names.
    assert_eq!(ReportType::NoFlank.label(), "Insufficient flank");
    assert_eq!(ReportType::IgnoredMate.label(), "Mate ignored");
    assert_eq!(ReportType::CalledVariant.label(), "Called variants");
}
