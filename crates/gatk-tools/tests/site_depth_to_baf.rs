//! Conformance for `SiteDepthtoBAF` against GATK 4.6.2.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/SiteDepthToBafDump.java`.
//!
//! # What this suite is for
//!
//!  * **the value written not being the fraction measured**: one survivor is replaced by `0.5`,
//!    several are shifted so their median lands there;
//!  * **the chi-squared het test** and the four defaults around it;
//!  * **the whole locus disappearing** when the sample deviation exceeds `--max-std`;
//!  * **`DecimalFormat("#.00")`**, which writes a half as `.50`;
//!  * **and the three refusals**, a reference base that is not [ACGT], a sites file that runs out,
//!    and a sites file that names another locus.

use gatk_corpus as corpus;
use gatk_tools::site_depth_to_baf::{
    format_baf, read_depths, read_sites, run, running_stddev, write, Arguments, BafError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/site_depth_to_baf.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn produced(text: &str, label: &str, arguments: &Arguments) -> String {
    let depths = read_depths(&value(text, "depths", label));
    let sites = read_sites(&value(text, "sites", label));
    write(&run(&depths, &sites, arguments).expect("a run that is not refused"))
}

#[test]
fn every_baf_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, arguments) in [
        ("defaults", Arguments::default()),
        (
            "het-0.9",
            Arguments {
                min_het_probability: 0.9,
                ..Arguments::default()
            },
        ),
        (
            "het-0.05",
            Arguments {
                min_het_probability: 0.05,
                ..Arguments::default()
            },
        ),
        (
            "depth-20",
            Arguments {
                min_total_depth: 20,
                ..Arguments::default()
            },
        ),
        (
            "std-0.5",
            Arguments {
                max_std_dev: 0.5,
                ..Arguments::default()
            },
        ),
        (
            "het-0.01",
            Arguments {
                min_het_probability: 0.01,
                ..Arguments::default()
            },
        ),
        (
            "het-0.01-std-0.5",
            Arguments {
                min_het_probability: 0.01,
                max_std_dev: 0.5,
                ..Arguments::default()
            },
        ),
    ] {
        assert_eq!(
            produced(&text, label, &arguments),
            value(&text, "baf", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 7, "the golden's outputs");
}

/// A locus with one survivor is written as a half whatever it measured: at the defaults the only
/// sample at locus 200 is 11 reference against 9 alternate, which is 0.45, and `.50` is written.
#[test]
fn a_lone_survivor_is_replaced_rather_than_adjusted() {
    let text = golden();
    let defaults = produced(&text, "defaults", &Arguments::default());
    assert!(defaults.contains("chr1\t199\t.50\ts1\n"));

    // The same sample, alone, measures 0.45 before the replacement.
    let depths = read_depths(&value(&text, "depths", "defaults"));
    let alone: Vec<_> = depths
        .iter()
        .filter(|depth| depth.position == 200)
        .cloned()
        .collect();
    assert_eq!(alone.len(), 1);
    assert_eq!(alone[0].counts, [11, 0, 9, 0]);
}

/// Several survivors are all moved by `0.5 - median`, so none of the written values is a measured
/// one: locus 500 measures 0.7, 0.6 and 0.5 and is written `.60`, `.50` and `.40`.
#[test]
fn a_locus_is_shifted_onto_its_own_median() {
    let text = golden();
    let loosened = produced(
        &text,
        "het-0.05",
        &Arguments {
            min_het_probability: 0.05,
            ..Arguments::default()
        },
    );
    assert!(loosened.contains("chr1\t499\t.60\ts1\n"));
    assert!(loosened.contains("chr1\t499\t.50\ts2\n"));
    assert!(loosened.contains("chr1\t499\t.40\ts3\n"));
    // And a locus whose median is already a half moves by nothing at all.
    assert!(loosened.contains("chr1\t99\t.50\ts1\n"));
    assert!(loosened.contains("chr1\t99\t.40\ts2\n"));
    assert!(loosened.contains("chr1\t99\t.60\ts3\n"));
}

/// The deviation drops the locus, not the outlier, and it is the SAMPLE deviation: locus 600
/// measures 0.75, 0.5 and 0.25, whose deviation over `n - 1` is exactly 0.25 and over `n` would be
/// 0.204. Both exceed the default, so the discriminating evidence is the pair of runs at 0.2 and
/// 0.5 rather than the arithmetic alone.
#[test]
fn the_deviation_drops_the_whole_locus() {
    let text = golden();
    let strict = produced(
        &text,
        "het-0.01",
        &Arguments {
            min_het_probability: 0.01,
            ..Arguments::default()
        },
    );
    assert!(!strict.contains("chr1\t599\t"));

    let loose = produced(
        &text,
        "het-0.01-std-0.5",
        &Arguments {
            min_het_probability: 0.01,
            max_std_dev: 0.5,
            ..Arguments::default()
        },
    );
    assert!(loose.contains("chr1\t599\t.75\ts1\n"));
    assert!(loose.contains("chr1\t599\t.50\ts2\n"));
    assert!(loose.contains("chr1\t599\t.25\ts3\n"));

    assert_eq!(running_stddev(&[0.75, 0.5, 0.25]), 0.25);
}

/// `#.00` has no integer digit, so a half is `.50` and not `0.50`, and the rounding is HALF_EVEN.
#[test]
fn the_formatter_drops_the_leading_zero() {
    assert_eq!(format_baf(0.5), ".50");
    assert_eq!(format_baf(0.4), ".40");
    assert_eq!(format_baf(0.75), ".75");
    assert_eq!(format_baf(1.0), "1.00");
    assert_eq!(format_baf(0.0), ".00");
    assert_eq!(format_baf(-0.5), "-.50");
    // Ties to even rather than away from zero.
    assert_eq!(format_baf(0.125), ".12");
    assert_eq!(format_baf(0.135), ".14");
}

#[test]
fn the_three_refusals_match_the_golden() {
    let text = golden();
    for (label, expected) in [
        ("bad-ref", "ref call is not [ACGT] in vcf at chr1:100"),
        (
            "short-sites",
            "baf sites vcf exhausted before site depth data",
        ),
        (
            "wrong-locus",
            "expecting locus chr1:100, but found locus chr1:150 in baf sites vcf",
        ),
    ] {
        let depths = read_depths(&value(&text, "depths", label));
        let sites = read_sites(&value(&text, "sites", label));
        let error: BafError =
            run(&depths, &sites, &Arguments::default()).expect_err("a refused run");
        assert_eq!(error.message(), expected, "{label}");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
}
