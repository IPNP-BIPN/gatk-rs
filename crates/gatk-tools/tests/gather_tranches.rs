//! Conformance for `GatherTranches` against GATK 4.6.2.0, compared as the whole gathered file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/GatherTranchesDump.java`.
//!
//! # What this suite is for
//!
//!  * **a merged level's Ti/Tv is a ratio of sums**, which `no-known` proves by coming back `NaN`;
//!  * **the shards are pooled by level, not by file**, so `reversed` writes what `two-shards`
//!    writes and `extra-level` merges a level from one shard;
//!  * **the sensitivity match answers with the level before the one that overshot**, and stops at
//!    the first target it cannot advance past: seven requests over four levels write two rows;
//!  * **the requested sensitivities are sorted**, so `unsorted-levels` writes them in order;
//!  * **the output rows carry the requested sensitivity**, not the achieved one, and the filter
//!    name is built from the previous row's target;
//!  * **and the output header says version 5** while the inputs must say version 6.

use gatk_corpus as corpus;
use gatk_engine::tranches::Mode;
use gatk_tools::gather_tranches::{gather, DEFAULT_TRUTH_SENSITIVITY_LEVELS};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gather_tranches.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, label: &str) -> String {
    let prefix = format!("tranches\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {label}")),
    )
}

const HEADER: &str = "# Variant quality score tranches file\n\
                      # Version number 6\n\
                      requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,\
                      model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";

/// One row of a VQSLOD tranche file, written the way the harness wrote it.
#[allow(clippy::too_many_arguments, reason = "the file's own eleven columns")]
fn row(
    requested: f64,
    known: i64,
    novel: i64,
    known_titv: f64,
    novel_titv: f64,
    min_vqslod: f64,
    accessible: i32,
    called: i32,
) -> String {
    let sensitivity = if accessible == 0 {
        0.0
    } else {
        called as f64 / accessible as f64
    };
    format!(
        "{requested:.4},{known},{novel},{known_titv:.4},{novel_titv:.4},{min_vqslod:.4},\
         VQSRTranche,SNP,{accessible},{called},{sensitivity:.4}\n"
    )
}

fn first_shard() -> String {
    format!(
        "{HEADER}{}{}{}{}",
        row(4.0, 100, 20, 2.0, 1.5, 4.0, 1000, 500),
        row(2.0, 200, 50, 2.1, 1.6, 2.0, 1000, 800),
        row(0.0, 300, 90, 2.2, 1.7, 0.0, 1000, 950),
        row(-2.0, 400, 150, 2.3, 1.8, -2.0, 1000, 990),
    )
}

fn second_shard() -> String {
    format!(
        "{HEADER}{}{}{}{}",
        row(4.0, 60, 10, 1.8, 1.4, 4.0, 1000, 450),
        row(2.0, 130, 30, 1.9, 1.5, 2.0, 1000, 780),
        row(0.0, 220, 70, 2.0, 1.6, 0.0, 1000, 940),
        row(-2.0, 330, 120, 2.1, 1.7, -2.0, 1000, 985),
    )
}

fn extra_level() -> String {
    format!(
        "{HEADER}{}{}",
        row(4.0, 10, 2, 2.0, 1.5, 4.0, 1000, 400),
        row(1.0, 90, 25, 2.0, 1.5, 1.0, 1000, 700),
    )
}

fn no_known() -> String {
    format!(
        "{HEADER}{}{}",
        row(4.0, 0, 20, 0.0, 1.5, 4.0, 1000, 500),
        row(2.0, 0, 50, 0.0, 1.6, 2.0, 1000, 800),
    )
}

/// The shards, the requested levels and the mode of each run.
fn run(label: &str) -> (Vec<String>, Vec<f64>, Mode) {
    let default = DEFAULT_TRUTH_SENSITIVITY_LEVELS.to_vec();
    match label {
        "two-shards" => (vec![first_shard(), second_shard()], default, Mode::Snp),
        "reversed" => (vec![second_shard(), first_shard()], default, Mode::Snp),
        "one-shard" => (vec![first_shard()], default, Mode::Snp),
        "extra-level" => (vec![first_shard(), extra_level()], default, Mode::Snp),
        "no-known" => (vec![no_known()], default, Mode::Snp),
        "one-level" => (vec![first_shard(), second_shard()], vec![99.0], Mode::Snp),
        "many-levels" => (
            vec![first_shard(), second_shard()],
            vec![100.0, 99.9, 99.5, 99.0, 98.0, 95.0, 90.0],
            Mode::Snp,
        ),
        "unsorted-levels" => (
            vec![first_shard(), second_shard()],
            vec![90.0, 99.9],
            Mode::Snp,
        ),
        "indel-mode" => (vec![first_shard(), second_shard()], default, Mode::Indel),
        "both-mode" => (vec![first_shard(), second_shard()], default, Mode::Both),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_gathered_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "two-shards",
        "reversed",
        "one-shard",
        "extra-level",
        "no-known",
        "one-level",
        "many-levels",
        "unsorted-levels",
        "indel-mode",
        "both-mode",
    ] {
        let (shards, levels, mode) = run(label);
        let ours = gather(&shards, &levels, mode).expect("a run the files allow");
        assert_eq!(ours, value(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 10, "the golden's gathered files");
}

#[test]
fn a_file_of_the_wrong_version_is_refused_before_a_row_is_read() {
    let text = golden();
    let old = first_shard().replace("Version number 6", "Version number 5");
    let error = gather(&[old], &DEFAULT_TRUTH_SENSITIVITY_LEVELS, Mode::Snp)
        .expect_err("a version the parser refuses");
    assert_eq!(
        format!("error\told-version\t{}:{}", error.class(), error.message()),
        text.lines()
            .find(|line| line.starts_with("error\told-version"))
            .expect("the golden carries the refusal")
    );
}
