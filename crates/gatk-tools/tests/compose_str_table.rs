//! Conformance for `ComposeSTRTableFile` against GATK 4.6.2.0, compared as the whole sites table.
//!
//! Golden from `tools/readfilter-conformance/ComposeSTRTableFileDump.java`.
//!
//! # What this suite is for
//!
//!  * **which period wins**, and the integer division that counts its repeats;
//!  * **the sites overlapping** although no position starts a search twice;
//!  * **the gap an N leaves**, period one being the only period tried at a base that is not ACGT;
//!  * **the mask starting at the contig's index**, and the decimation that follows from it;
//!  * **`--max-repeat` changing the masks** without changing the repeats reported;
//!  * **`--max-period` breaking a long repeat into single bases**;
//!  * **and `-L` bounding where the scan starts and nothing else.**

use gatk_corpus as corpus;
use gatk_tools::compose_str_table::{scan, DecimationTable, Locus, Settings};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/compose_str_table.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn contigs(text: &str) -> Vec<(String, Vec<u8>)> {
    ["chr1", "chr2"]
        .iter()
        .map(|name| {
            (
                name.to_string(),
                text.lines()
                    .find_map(|line| line.strip_prefix(&format!("fixture\t{name}=")))
                    .unwrap_or_else(|| panic!("the golden carries {name}"))
                    .as_bytes()
                    .to_vec(),
            )
        })
        .collect()
}

/// The sites table of one run, parsed back into loci.
fn sites(text: &str, label: &str) -> Vec<Locus> {
    let table = unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("sites\t{label}=")))
            .unwrap_or_else(|| panic!("the golden carries the sites of {label}")),
    );
    table
        .lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let start: i64 = fields[2].parse().expect("a start");
            let end: i64 = fields[3].parse().expect("an end");
            let period: usize = fields[4].parse().expect("a period");
            let mask: i64 = fields[5].parse().expect("a mask");
            // The two length columns are derived, and the suite checks them rather than trusting.
            assert_eq!(fields[7].parse::<i64>().expect("a span"), end - start + 1);
            assert_eq!(
                fields[8].parse::<i64>().expect("a repeat count"),
                (end - start + 1) / period as i64
            );
            assert_eq!(fields[6], format!("{mask:b}"), "the binary column");
            Locus {
                contig_index: fields[0].parse().expect("a contig index"),
                start,
                end,
                period,
                repeats: ((end - start + 1) / period as i64) as usize,
                mask,
            }
        })
        .collect()
}

fn check(
    text: &str,
    label: &str,
    table: DecimationTable,
    max_period: usize,
    max_repeat: usize,
    intervals: &[(String, i64, i64)],
) {
    let result = scan(
        &contigs(text),
        intervals,
        Settings {
            max_period,
            max_repeat,
        },
        &table,
    );
    assert_eq!(result.emitted, sites(text, label), "{label}");
}

#[test]
fn the_default_decimation_table() {
    let text = golden();
    check(
        &text,
        "default",
        DecimationTable::default_table(),
        8,
        20,
        &[],
    );
}

#[test]
fn no_decimation_keeps_the_second_contigs_homopolymer() {
    let text = golden();
    check(&text, "no-decimation", DecimationTable::none(), 8, 20, &[]);
    // The only difference between the two runs is chr2's period-one site, whose mask is 1
    // because the counter starts at the contig's index.
    let kept = sites(&text, "no-decimation");
    let decimated = sites(&text, "default");
    assert_eq!(kept.len(), decimated.len() + 1);
    let lost = kept
        .iter()
        .find(|locus| !decimated.contains(locus))
        .expect("one site went");
    assert_eq!((lost.contig_index, lost.period, lost.mask), (1, 1, 1));
}

#[test]
fn a_shorter_maximum_period_breaks_a_repeat_into_single_bases() {
    let text = golden();
    check(&text, "max-period-two", DecimationTable::none(), 2, 20, &[]);
    // The trinucleotide and four-base regions become one site per base, each with its own mask.
    let rows = sites(&text, "max-period-two");
    assert!(rows.iter().filter(|locus| locus.repeats == 1).count() > 20);
}

#[test]
fn a_lower_maximum_repeat_changes_the_masks_and_not_the_repeats() {
    let text = golden();
    check(
        &text,
        "max-repeat-three",
        DecimationTable::none(),
        8,
        3,
        &[],
    );
    let capped = sites(&text, "max-repeat-three");
    let uncapped = sites(&text, "no-decimation");
    // The same sites, at the same places, with the same repeat counts.
    assert_eq!(
        capped
            .iter()
            .map(|l| (l.start, l.end, l.period, l.repeats))
            .collect::<Vec<_>>(),
        uncapped
            .iter()
            .map(|l| (l.start, l.end, l.period, l.repeats))
            .collect::<Vec<_>>()
    );
    // But not the same masks: the period-one sites of chr1 now share one counter.
    assert_ne!(
        capped.iter().map(|l| l.mask).collect::<Vec<_>>(),
        uncapped.iter().map(|l| l.mask).collect::<Vec<_>>()
    );
}

#[test]
fn an_interval_bounds_where_the_scan_starts_and_nothing_else() {
    let text = golden();
    check(
        &text,
        "interval",
        DecimationTable::none(),
        8,
        20,
        &[("chr1".to_string(), 20, 40)],
    );
    let rows = sites(&text, "interval");
    // A site beginning before the interval and another ending after it.
    assert!(rows.iter().any(|locus| locus.start < 20));
    assert!(rows.iter().any(|locus| locus.end > 40));
}

#[test]
fn the_n_leaves_a_gap_and_the_sites_overlap() {
    let text = golden();
    let rows = sites(&text, "no-decimation");
    // The first two sites of chr1 share the base at 9.
    assert_eq!((rows[0].end, rows[1].start), (9, 9));
    // And nothing is reported at the N, which the fixture puts at 37.
    assert!(!rows
        .iter()
        .any(|locus| locus.contig_index == 0 && locus.start <= 37 && 37 <= locus.end));
}
