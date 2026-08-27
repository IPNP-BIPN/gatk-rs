//! Conformance for `ReblockGVCF` against GATK 4.6.2.0, compared as the records of every run.
//!
//! Golden from `tools/readfilter-conformance/ReblockGVCFDump.java`.
//!
//! The annotations the tool recomputes for itself, `RAW_GT_COUNT` and `RAW_MQandDP`, are in the
//! golden but not reproduced: they come from the annotation engine, which is not ported.
//!
//! # What this suite is for
//!
//!  * **adjacent blocks in one band merging at the lowest quality**;
//!  * **the band edges deciding which blocks merge**;
//!  * **`--drop-low-quals` and `--rgq-threshold` touching different records**;
//!  * **a demoted variant becoming a one-base block that does not merge**;
//!  * **`--floor-blocks` writing the bound and dropping MIN_DP and PL**;
//!  * **and the two annotation arguments not being symmetric.**

use gatk_corpus as corpus;
use gatk_tools::reblock_gvcf::{
    band_floor, band_of, check_annotations_to_keep, demote, is_demoted, reblock, Arguments,
    ReblockError, Record, DEFAULT_GQ_BANDS,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reblock_gvcf.txt.gz"),
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

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), message.to_string())
}

fn field(genotype: &str, format: &str, key: &str) -> Option<i32> {
    let keys: Vec<&str> = format.split(':').collect();
    let values: Vec<&str> = genotype.split(':').collect();
    keys.iter()
        .position(|name| *name == key)
        .and_then(|at| values.get(at))
        .and_then(|text| text.parse().ok())
}

/// The input GVCF, read as the records the tool sees.
fn input(text: &str) -> Vec<Record> {
    section(text, "vcf", "input")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let info: Vec<(String, String)> = columns[7]
                .split(';')
                .filter_map(|part| part.split_once('='))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect();
            let alternates: Vec<String> = columns[4]
                .split(',')
                .filter(|allele| *allele != "<NON_REF>")
                .map(str::to_string)
                .collect();
            let end = info
                .iter()
                .find(|(key, _)| key == "END")
                .map(|(_, value)| value.parse().expect("an end"))
                .unwrap_or_else(|| columns[1].parse().expect("a position"));
            let likelihoods: Vec<i32> = field(columns[9], columns[8], "PL")
                .map(|_| Vec::new())
                .unwrap_or_default();
            let _ = likelihoods;
            let reference_likelihood = columns[9]
                .split(':')
                .nth(
                    columns[8]
                        .split(':')
                        .position(|key| key == "PL")
                        .expect("a PL"),
                )
                .and_then(|text| text.split(',').next())
                .and_then(|text| text.parse().ok())
                .unwrap_or(0);
            Record {
                contig: columns[0].to_string(),
                start: columns[1].parse().expect("a position"),
                end,
                alternates,
                depth: field(columns[9], columns[8], "DP").unwrap_or(0),
                minimum_depth: field(columns[9], columns[8], "MIN_DP")
                    .unwrap_or_else(|| field(columns[9], columns[8], "DP").unwrap_or(0)),
                genotype_quality: field(columns[9], columns[8], "GQ").unwrap_or(0),
                reference_likelihood,
                info: info.into_iter().filter(|(key, _)| key != "END").collect(),
                format: columns[8].split(':').map(str::to_string).collect(),
            }
        })
        .collect()
}

/// One run's records, as position, end and genotype quality.
fn measured(text: &str, label: &str) -> Vec<(i32, i32, i32)> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let start: i32 = columns[1].parse().expect("a position");
            let end = columns[7]
                .split(';')
                .find_map(|part| part.strip_prefix("END="))
                .map(|value| value.parse().expect("an end"))
                .unwrap_or(start);
            (
                start,
                end,
                field(columns[9], columns[8], "GQ").expect("a quality"),
            )
        })
        .collect()
}

fn produced(records: &[Record], arguments: &Arguments) -> Vec<(i32, i32, i32)> {
    reblock(records, arguments)
        .into_iter()
        .map(|record| (record.start, record.end, record.genotype_quality))
        .collect()
}

/// label, arguments.
fn runs() -> Vec<(&'static str, Arguments)> {
    vec![
        ("default", Arguments::default()),
        (
            "one-band",
            Arguments {
                gq_bands: vec![60],
                ..Arguments::default()
            },
        ),
        (
            "many-bands",
            Arguments {
                gq_bands: vec![10, 20, 30, 40, 50, 60],
                ..Arguments::default()
            },
        ),
        (
            "rgq-threshold",
            Arguments {
                rgq_threshold: 10,
                ..Arguments::default()
            },
        ),
        (
            "drop-low-quals",
            Arguments {
                drop_low_quals: true,
                ..Arguments::default()
            },
        ),
        (
            "drop-and-threshold",
            Arguments {
                drop_low_quals: true,
                rgq_threshold: 10,
                ..Arguments::default()
            },
        ),
        ("keep-all-alts", Arguments::default()),
        ("keep-annotation", Arguments::default()),
        ("remove-annotation", Arguments::default()),
    ]
}

#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let records = input(&text);
    let mut compared = 0;
    for (label, arguments) in runs() {
        assert_eq!(
            produced(&records, &arguments),
            measured(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 9, "the runs the port reproduces");
}

/// Three blocks at 25, 35 and 80 become one block at 25 under the default bands.
#[test]
fn adjacent_blocks_in_one_band_merge_at_the_lowest_quality() {
    let text = golden();
    let default = measured(&text, "default");
    assert_eq!(default[0], (1000, 1299, 25), "all three, at the lowest");
    // The three arrived as separate records.
    let records = input(&text);
    assert_eq!(records[0].genotype_quality, 25);
    assert_eq!(records[1].genotype_quality, 35);
    assert_eq!(records[2].genotype_quality, 80);
    // And all three fall in the same default band.
    for record in &records[..3] {
        assert_eq!(band_of(record.genotype_quality, DEFAULT_GQ_BANDS), 1);
    }
    assert_eq!(DEFAULT_GQ_BANDS, &[20, 100]);
}

/// A single bound of 60 puts 0 and 25 together; six bounds keep every block apart.
#[test]
fn the_band_edges_decide_which_blocks_merge() {
    let text = golden();
    let one = measured(&text, "one-band");
    assert_eq!(one[0], (1000, 1199, 25), "25 and 35 together");
    assert_eq!(one[1], (1200, 1299, 80), "80 in the other band");
    // The last two blocks, at 25 and 0, merge under the same single bound.
    assert!(one.contains(&(1401, 1599, 0)));
    assert_eq!(band_of(0, &[60]), band_of(25, &[60]));

    let many = measured(&text, "many-bands");
    assert_eq!(many[0], (1000, 1099, 25));
    assert_eq!(many[1], (1100, 1199, 35));
    assert_eq!(many[2], (1200, 1299, 80));
    assert_ne!(
        band_of(25, &[10, 20, 30, 40, 50, 60]),
        band_of(35, &[10, 20, 30, 40, 50, 60])
    );
}

/// One removes a block, the other demotes a variant, and neither does the other's job.
#[test]
fn the_two_thresholds_touch_different_records() {
    let text = golden();
    let dropped = measured(&text, "drop-low-quals");
    let demoted = measured(&text, "rgq-threshold");
    // Dropping removes the GQ0 block at 1500 and leaves the weak variant at 1400.
    assert!(!dropped.iter().any(|(start, ..)| *start == 1500));
    assert!(dropped
        .iter()
        .any(|(start, end, _)| *start == 1400 && *end == 1400));
    // The threshold leaves the block and demotes the variant.
    assert!(demoted.iter().any(|(start, ..)| *start == 1500));
    assert!(demoted.contains(&(1400, 1400, 0)), "a one-base GQ0 block");
    // The default leaves both alone.
    let default = measured(&text, "default");
    assert!(default.iter().any(|(start, ..)| *start == 1500));
    assert!(default
        .iter()
        .any(|(start, _, quality)| *start == 1400 && *quality == 3));
}

/// It becomes a one-base block of its own and does not merge with the blocks either side.
#[test]
fn a_demoted_variant_does_not_merge() {
    let text = golden();
    let records = input(&text);
    let weak = records
        .iter()
        .find(|record| record.start == 1400)
        .expect("the weak variant");
    assert!(is_demoted(weak, 10));
    assert!(!is_demoted(weak, 5), "the comparison is strictly less than");
    let block = demote(weak);
    assert_eq!(block.start, block.end, "one base");
    assert_eq!(block.genotype_quality, 0);
    assert!(block.is_reference_block());

    // The blocks either side are at 25, in the same default band as a GQ0 block would be...
    let demoted = measured(&text, "rgq-threshold");
    assert!(demoted.contains(&(1301, 1399, 25)));
    assert!(demoted.contains(&(1400, 1400, 0)));
    assert!(demoted.contains(&(1401, 1499, 25)));
    // ...and all three stay three records.
    assert_eq!(
        demoted
            .iter()
            .filter(|(start, ..)| (1301..=1499).contains(start))
            .count(),
        3
    );
}

/// The bound, not the observed quality, and MIN_DP and PL go with it.
#[test]
fn floor_blocks_writes_the_bound_and_drops_two_fields() {
    let text = golden();
    let floored = section(&text, "out", "floor-blocks");
    let first = floored
        .lines()
        .find(|line| line.starts_with("chr1\t1000\t"))
        .expect("the first block");
    let columns: Vec<&str> = first.split('\t').collect();
    assert_eq!(columns[8], "GT:DP:GQ", "MIN_DP and PL are gone");
    assert_eq!(field(columns[9], columns[8], "GQ"), Some(20));
    assert_eq!(band_floor(25, DEFAULT_GQ_BANDS), 20);
    assert_eq!(band_floor(80, DEFAULT_GQ_BANDS), 20);
    assert_eq!(band_floor(0, DEFAULT_GQ_BANDS), 0, "below the first bound");
    // The default writes the observed quality and keeps all five fields.
    let plain = section(&text, "out", "default");
    let plain_first = plain
        .lines()
        .find(|line| line.starts_with("chr1\t1000\t"))
        .expect("the first block");
    let plain_columns: Vec<&str> = plain_first.split('\t').collect();
    assert_eq!(plain_columns[8], "GT:DP:GQ:MIN_DP:PL");
    assert_eq!(field(plain_columns[9], plain_columns[8], "GQ"), Some(25));
}

/// One takes a FORMAT key, the other an INFO one, and asking wrongly is refused.
#[test]
fn the_two_annotation_arguments_are_not_symmetric() {
    let text = golden();
    // The INFO annotation survives when asked for and is gone by default.
    assert!(section(&text, "out", "keep-annotation").contains("EXTRA=note"));
    assert!(!section(&text, "out", "default").contains("EXTRA=note"));
    // The FORMAT one is removed when asked for and present by default.
    assert!(!section(&text, "out", "remove-annotation").contains(":SPARE"));
    assert!(section(&text, "out", "default").contains("GT:AD:DP:GQ:PL:SPARE"));

    let (class, message) = refusal(&text, "keep-format-key");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException"
    );
    let produced = check_annotations_to_keep(
        &["SPARE".to_string()],
        &["END".to_string(), "EXTRA".to_string()],
    )
    .expect_err("a FORMAT key");
    assert_eq!(
        produced,
        ReblockError::NotInHeader {
            key: "SPARE".to_string()
        }
    );
    assert_eq!(produced.message(), message);
    // The INFO key it does declare is accepted.
    assert!(check_annotations_to_keep(
        &["EXTRA".to_string()],
        &["END".to_string(), "EXTRA".to_string()]
    )
    .is_ok());
}

/// A variant already biallelic with <NON_REF> is untouched, so the flag is the control.
#[test]
fn keep_all_alts_changes_nothing_here() {
    let text = golden();
    assert_eq!(
        section(&text, "out", "keep-all-alts"),
        section(&text, "out", "default")
    );
}
