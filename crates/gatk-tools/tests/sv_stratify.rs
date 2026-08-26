//! Conformance for `SVStratify` against GATK 4.6.2.0, compared as the whole set of records written
//! by every run.
//!
//! Golden from `tools/readfilter-conformance/SVStratifyDump.java`.
//!
//! # What this suite is for
//!
//!  * **the strata coming back in `java.util.HashMap` order**, not the configuration's;
//!  * **a record matching two strata being written twice into one file** without `--split-output`;
//!  * **an insertion ignoring both thresholds**;
//!  * **`-1` being the only negative that means null**;
//!  * **and the column-count message printing the same number twice**.

use gatk_corpus as corpus;
use gatk_tools::sv_stratify::{
    apply, check_columns, parse_integer_maybe_null, parse_track_string, split_output_files,
    CallRecord, Engine, Interval, StratifyError, Stratum, SvType, Thresholds, Tracks,
    DEFAULT_STRATUM,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sv_stratify.txt.gz"),
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

fn refusal(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"))
        .to_string()
}

/// One track's intervals. A bed is half-open and zero-based, and `loadIntervals` turns it into the
/// closed one-based interval the engine compares with.
fn track(text: &str, name: &str) -> Vec<Interval> {
    section(text, "track", name)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Interval {
                contig: columns[0].to_string(),
                start: columns[1].parse::<i32>().expect("a bed start") + 1,
                end: columns[2].parse().expect("a bed end"),
            }
        })
        .collect()
}

/// The strata of one configuration table, in the file's order.
fn strata(text: &str, name: &str, tracks: &[String]) -> Result<Vec<Stratum>, StratifyError> {
    let table = section(text, "config", name);
    let mut lines = table.lines().filter(|line| !line.is_empty());
    let columns: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .map(str::to_string)
        .collect();
    check_columns(&columns)?;
    let mut out = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split('\t').collect();
        let get = |column: &str| {
            values[columns
                .iter()
                .position(|name| name == column)
                .expect("a column")]
        };
        out.push(Stratum::new(
            get("NAME"),
            SvType::parse(get("SVTYPE")).expect("an sv type"),
            parse_integer_maybe_null(get("MIN_SIZE")),
            parse_integer_maybe_null(get("MAX_SIZE")),
            parse_track_string(get("TRACKS"), tracks)?,
        )?);
    }
    Ok(out)
}

/// The measured records, read out of the input VCF the golden carries.
fn records(text: &str) -> Vec<CallRecord> {
    section(text, "vcf", "input")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let info: Vec<(&str, &str)> = columns[7]
                .split(';')
                .filter_map(|part| part.split_once('='))
                .collect();
            let field = |key: &str| {
                info.iter()
                    .find(|(name, _)| *name == key)
                    .map(|(_, value)| *value)
            };
            let sv_type = SvType::parse(field("SVTYPE").expect("a type")).expect("a known type");
            let start: i32 = columns[1].parse().expect("a position");
            let end: i32 = field("END").expect("an end").parse().expect("an end");
            // BND and INS carry no meaningful length: the record's own length is absent for them,
            // which is what makes them reach only a stratum with neither bound.
            let length = match sv_type {
                SvType::Bnd | SvType::Ctx | SvType::Ins => None,
                _ => Some(field("SVLEN").expect("a length").parse().expect("a length")),
            };
            CallRecord {
                id: columns[2].to_string(),
                sv_type,
                contig_a: columns[0].to_string(),
                position_a: start,
                contig_b: field("CHR2").unwrap_or(columns[0]).to_string(),
                position_b: field("END2")
                    .map(|value| value.parse().expect("a second position"))
                    .unwrap_or(end),
                length,
            }
        })
        .collect()
}

fn tracks(text: &str) -> Tracks {
    Tracks::new(
        &["RM".to_string(), "SD".to_string()],
        &[track(text, "RM"), track(text, "SD")],
    )
    .expect("two distinct tracks")
}

fn engine(text: &str) -> Engine {
    let names = vec!["RM".to_string(), "SD".to_string()];
    Engine::new(
        strata(text, "main", &names).expect("a valid configuration"),
        tracks(text),
    )
    .expect("no reserved name")
}

fn defaults() -> Thresholds {
    Thresholds {
        overlap_fraction: 0.0,
        num_breakpoint_overlaps: 1,
        num_breakpoint_overlaps_interchrom: 1,
    }
}

/// The `ID` and `STRAT` of every record one run wrote, in order, out of a measured output file.
fn written(text: &str, label: &str, file: &str) -> Vec<(String, String)> {
    let prefix = format!("out\t{label}\t{file}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries out/{label}/{file}")),
    )
    .lines()
    .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
    .map(|line| {
        let columns: Vec<&str> = line.split('\t').collect();
        let stratum = columns[7]
            .split(';')
            .find_map(|part| part.strip_prefix("STRAT="))
            .expect("a stratum")
            .to_string();
        (columns[2].to_string(), stratum)
    })
    .collect()
}

fn produced(
    text: &str,
    thresholds: Thresholds,
    allow_multiple: bool,
    split: bool,
) -> Vec<(String, String)> {
    let engine = engine(text);
    let mut out = Vec::new();
    for record in records(text) {
        for written in
            apply(&engine, &record, thresholds, allow_multiple, split).expect("a stratified record")
        {
            out.push((record.id.clone(), written.stratum));
        }
    }
    out
}

/// Every record of every run that produced one, compared whole.
#[test]
fn every_written_record_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, file, thresholds, allow_multiple) in [
        (
            "overlap-half",
            "overlap-half.vcf",
            Thresholds {
                overlap_fraction: 0.5,
                ..defaults()
            },
            false,
        ),
        (
            "overlap-full",
            "overlap-full.vcf",
            Thresholds {
                overlap_fraction: 1.0,
                ..defaults()
            },
            false,
        ),
        ("multiple-allowed", "multiple-allowed.vcf", defaults(), true),
    ] {
        assert_eq!(
            produced(&text, thresholds, allow_multiple, false),
            written(&text, label, file),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 3, "the runs that wrote records");
}

/// A record matching two strata is written twice into ONE file, because the writer is chosen by
/// the default group unless the output is split.
#[test]
fn two_matches_are_written_twice_into_one_file() {
    let text = golden();
    let rows = written(&text, "multiple-allowed", "multiple-allowed.vcf");
    let small: Vec<&(String, String)> =
        rows.iter().filter(|(id, _)| id == "del_small_rm").collect();
    assert_eq!(small.len(), 2, "the same record, twice");
    assert_eq!(small[0].1, "DEL_small_both");
    assert_eq!(small[1].1, "DEL_small_RM");
    assert_eq!(produced(&text, defaults(), true, false), rows);
}

/// The configuration declares `DEL_small_RM` first and the engine reports `DEL_small_both` first,
/// because the strata are a `java.util.HashMap` keyed by name.
#[test]
fn the_strata_come_back_in_hash_map_order() {
    let text = golden();
    let names = vec!["RM".to_string(), "SD".to_string()];
    let declared: Vec<String> = strata(&text, "main", &names)
        .expect("a configuration")
        .iter()
        .map(|stratum| stratum.name.clone())
        .collect();
    let walked: Vec<String> = engine(&text)
        .strata
        .iter()
        .map(|stratum| stratum.name.clone())
        .collect();
    assert_ne!(declared, walked, "the map does not keep the file's order");
    assert_eq!(
        declared
            .iter()
            .position(|name| name == "DEL_small_RM")
            .expect("declared"),
        0
    );
    assert!(
        walked.iter().position(|name| name == "DEL_small_both")
            < walked.iter().position(|name| name == "DEL_small_RM")
    );

    // The refusal lists them in that same order.
    let error = StratifyError::MultipleMatches {
        id: "del_small_rm".to_string(),
        names: vec!["DEL_small_both".to_string(), "DEL_small_RM".to_string()],
    };
    assert!(refusal(&text, "multiple-refused").ends_with(&error.message()));

    // And so does the order the split-output files are created in, after the default group.
    let files = split_output_files(&engine(&text), "strat");
    assert_eq!(files[0], "strat.default.vcf.gz");
    assert_eq!(files.len(), walked.len() + 1);
}

/// An insertion asks for one end in a track whatever the arguments said, so raising the fraction
/// to 1.0 drops every deletion to the default group and leaves the insertion where it was.
#[test]
fn an_insertion_ignores_both_thresholds() {
    let text = golden();
    let full = written(&text, "overlap-full", "overlap-full.vcf");
    assert!(full.contains(&("ins_rm".to_string(), "INS_RM".to_string())));
    for id in ["del_small_rm", "del_large_rm", "del_both"] {
        assert!(
            full.contains(&(id.to_string(), DEFAULT_STRATUM.to_string())),
            "{id} fell to the default group"
        );
    }
    // And the breakpoint threshold does not reach it either.
    let engine = engine(&text);
    let insertion = records(&text)
        .into_iter()
        .find(|record| record.id == "ins_rm")
        .expect("the insertion");
    for required in [1, 2] {
        let matched = engine
            .matches(
                &insertion,
                Thresholds {
                    num_breakpoint_overlaps: required,
                    ..defaults()
                },
            )
            .expect("a match");
        assert_eq!(matched.len(), 1, "at {required}");
        assert_eq!(matched[0].name, "INS_RM");
    }
}

/// `-1` is null and `-2` is a number, which the constructor then refuses.
#[test]
fn only_minus_one_is_a_null_bound() {
    let text = golden();
    assert_eq!(parse_integer_maybe_null("-1"), None);
    assert_eq!(parse_integer_maybe_null(""), None);
    assert_eq!(parse_integer_maybe_null("NULL"), None);
    assert_eq!(parse_integer_maybe_null("NA"), None);
    assert_eq!(parse_integer_maybe_null("-2"), Some(-2));

    let names = vec!["RM".to_string()];
    let error = strata(&text, "negative-min", &names).expect_err("a negative bound");
    assert_eq!(error, StratifyError::NegativeMin);
    assert!(refusal(&text, "negative-min").ends_with(&error.message()));

    // A null minimum reaches Integer.MIN_VALUE rather than zero, which is what lets a stratum take
    // a one-base record.
    let stratum = Stratum::new("x", SvType::Del, None, Some(5000), vec![]).expect("a stratum");
    assert_eq!(stratum.min_size, i32::MIN);
    assert_eq!(
        Stratum::new("y", SvType::Del, Some(5000), None, vec![])
            .expect("a stratum")
            .max_size,
        i32::MAX
    );
}

/// Every configuration refusal, each matching the message the reference wrote.
#[test]
fn the_configuration_refusals_match_the_golden() {
    let text = golden();
    let names = vec!["RM".to_string(), "SD".to_string()];

    for (label, expected) in [
        ("min-over-max", StratifyError::MinGreaterThanMax),
        (
            "bnd-with-size",
            StratifyError::SizedInterchromosomal {
                name: "BND_sized".to_string(),
            },
        ),
        (
            "unknown-track",
            StratifyError::UnknownTrack {
                name: "XX".to_string(),
            },
        ),
        (
            "missing-column",
            StratifyError::MissingColumn {
                column: "TRACKS".to_string(),
            },
        ),
        ("extra-column", StratifyError::ColumnCount { count: 6 }),
    ] {
        let error = strata(&text, label, &names).expect_err(label);
        assert_eq!(error, expected, "{label}");
        assert!(
            refusal(&text, label).ends_with(&error.message()),
            "{label}: {}",
            error.message()
        );
    }

    // The extra column is reported with the same number on both sides of the message.
    assert_eq!(
        StratifyError::ColumnCount { count: 6 }.message(),
        "Expected 6 columns but found 6"
    );

    // The reserved name is refused when the engine is built, not when the row is parsed.
    let reserved = strata(&text, "reserved-name", &names).expect("a parsable row");
    let error = Engine::new(reserved, tracks(&text)).expect_err("the reserved name");
    assert_eq!(error, StratifyError::ReservedName);
    assert!(refusal(&text, "reserved-name").ends_with(&error.message()));
}

/// The two track-argument refusals and the one threshold refusal.
#[test]
fn the_argument_refusals_match_the_golden() {
    let text = golden();

    let error = Tracks::new(
        &["RM".to_string(), "RM".to_string()],
        &[track(&text, "RM"), track(&text, "SD")],
    )
    .expect_err("a duplicate name");
    assert_eq!(
        error,
        StratifyError::DuplicateTrack {
            name: "RM".to_string()
        }
    );
    assert!(refusal(&text, "duplicate-track").ends_with(&error.message()));

    let error = Tracks::new(&["RM".to_string(), "SD".to_string()], &[track(&text, "RM")])
        .expect_err("a count mismatch");
    assert_eq!(error, StratifyError::TrackCountMismatch);
    assert!(refusal(&text, "track-count-mismatch").ends_with(&error.message()));

    let error = Thresholds {
        overlap_fraction: 0.0,
        num_breakpoint_overlaps: 0,
        num_breakpoint_overlaps_interchrom: 1,
    }
    .check()
    .expect_err("both at zero");
    assert_eq!(error, StratifyError::ThresholdsBothZero);
    assert!(refusal(&text, "both-thresholds-zero").ends_with(&error.message()));
}
