//! Conformance for `CheckPileup` against GATK 4.6.2.0, compared as its report and its summary.
//!
//! Golden from `tools/readfilter-conformance/CheckPileupDump.java`. The three input BAMs and their
//! indexes travel in full, base64, and the truth files travel as text.
//!
//! # What this suite is for
//!
//!  * **overlapping pairs are fixed by default**, so the qualities reported are not the reads';
//!  * **a duplicate is dropped by the filters**, which shows up as a size disagreement;
//!  * **the bases are compared case-insensitively and the qualities are not**;
//!  * **the message prints the qualities as raw Phred bytes**, where the report line prints them as
//!    SAM qualities;
//!  * **the file is written before the exception**, so a failing run leaves the line that explains
//!    it;
//!  * **and the summary counts the bases after the filters and the fixing**.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use gatk_engine::sam_pileup::{self, SamPileupFeature};
use gatk_tools::check_pileup::{self, CheckPileupArguments};
use gatk_tools::locus_walker;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/check_pileup.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn reference_bases(text: &str) -> String {
    rows(text, "reference")[0][0].to_string()
}

/// The fixtures and the reference, written where the port can open them.
fn install(text: &str, dir: &std::path::Path) -> std::path::PathBuf {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for row in rows(text, "fixture") {
        std::fs::write(
            dir.join(format!("{}.bam", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture bam");
    }
    for row in rows(text, "fixtureindex") {
        if row[1] == "absent" {
            continue;
        }
        std::fs::write(
            dir.join(format!("{}.bai", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture index");
    }

    let bases = reference_bases(text);
    let fasta = dir.join("reference.fasta");
    std::fs::write(&fasta, format!(">chr1\n{bases}\n")).expect("the reference");
    std::fs::write(
        dir.join("reference.fasta.fai"),
        format!(
            "chr1\t{}\t6\t{}\t{}\n",
            bases.len(),
            bases.len(),
            bases.len() + 1
        ),
    )
    .expect("the reference index");
    fasta
}

/// The truth file of one label, decoded through the pileup codec.
fn truth(text: &str, label: &str) -> Vec<SamPileupFeature> {
    let file = rows(text, "truth")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no truth file {label}"))[1]
        .to_string();
    unescape(&file)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| sam_pileup::decode(line).expect("the truth file parses"))
        .collect()
}

/// The BAM, the truth file and the arguments each labelled run used.
fn configuration(label: &str) -> (&str, &str, CheckPileupArguments) {
    let default = CheckPileupArguments::default();
    let carry_on = CheckPileupArguments {
        continue_after_error: true,
        ..default
    };
    match label {
        "agrees" => ("plain", "agrees", default),
        "wrong-size" => ("plain", "wrong-size", default),
        "wrong-bases" => ("plain", "wrong-bases", default),
        "wrong-quals" => ("plain", "wrong-quals", default),
        "incomplete" => ("plain", "incomplete", default),
        "wrong-bases-continue" => ("plain", "wrong-bases", carry_on),
        "incomplete-continue" => ("plain", "incomplete", carry_on),
        "duplicate" => ("duplicate", "agrees", carry_on),
        "overlapping" => ("overlapping", "agrees", carry_on),
        "overlapping-ignored" => (
            "overlapping",
            "agrees",
            CheckPileupArguments {
                ignore_overlaps: true,
                continue_after_error: true,
            },
        ),
        other => panic!("no run {other}"),
    }
}

/// One labelled run: the report it wrote, the summary it returned, and the first refusal.
fn run(
    text: &str,
    dir: &std::path::Path,
    fasta: &std::path::Path,
    label: &str,
) -> (
    String,
    Option<String>,
    Option<check_pileup::CheckPileupError>,
) {
    let (fixture, truth_label, arguments) = configuration(label);
    let bam = dir.join(format!("{fixture}.bam"));
    let bai = dir.join(format!("{fixture}.bai"));
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();
    let mut reference = ReferenceFileSource::open(fasta).expect("the reference opens");
    let features = truth(text, truth_label);

    let records: Vec<BamRecord> =
        gatk_tools::read_walker::traverse(&source, &[], &|_| true).expect("the traversal");

    // `LocusWalker`'s own filters, then the three samtools ones this tool adds.
    let base_filter = locus_walker::default_filter(&header);
    let filter = |read: &BamRecord| {
        base_filter(read)
            && !gatk_engine::read::is_duplicate(read)
            && !gatk_engine::read::fails_vendor_quality_check(read)
            && !gatk_engine::read::is_secondary_alignment(read)
    };

    let applied = locus_walker::traverse(
        &records,
        &header,
        Some(&mut reference),
        None,
        locus_walker::Options::default(),
        &filter,
    )
    .expect("the locus traversal");

    let bases_of_contig = reference_bases(text);
    let mut report = String::new();
    let mut first_error = None;
    let mut loci = 0i64;
    let mut bases = 0i64;

    for one in &applied {
        let reference_base = bases_of_contig.as_bytes()[(one.context.position - 1) as usize];
        let feature = features
            .iter()
            .find(|f| f.contig == one.context.contig && f.position == one.context.position);
        let (line, error) = check_pileup::apply(&one.context, reference_base, feature, &arguments);
        if let Some(line) = line {
            report.push_str(&line);
        }
        if let Some(error) = error {
            if first_error.is_none() {
                first_error = Some(error);
            }
            if !arguments.continue_after_error {
                // The reference throws here, so the run returns no summary at all: the counters it
                // had reached are never printed anywhere.
                return (report, None, first_error);
            }
        }
        loci += 1;
        bases += one.context.pileup.size() as i64;
    }

    (
        report,
        Some(check_pileup::summary(loci, bases)),
        first_error,
    )
}

#[test]
fn every_report_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-checkpileup-{}", std::process::id()));
    let fasta = install(&text, &dir);

    let reports = rows(&text, "report");
    assert_eq!(reports.len(), 10, "ten runs leave a file");

    for row in &reports {
        let label = row[0];
        let (report, _, _) = run(&text, &dir, &fasta, label);
        let expected = if row.len() > 1 {
            unescape(row[1])
        } else {
            String::new()
        };
        assert_eq!(report, expected, "report/{label}");
    }
}

#[test]
fn every_summary_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-checkpileup-sum-{}", std::process::id()));
    let fasta = install(&text, &dir);

    let summaries = rows(&text, "summary");
    assert_eq!(summaries.len(), 6, "six runs finish");

    for row in &summaries {
        let label = row[0];
        let (_, summary, _) = run(&text, &dir, &fasta, label);
        assert_eq!(
            summary.expect("this run finishes"),
            unescape(row[1]),
            "summary/{label}"
        );
    }
}

#[test]
fn every_refusal_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-checkpileup-err-{}", std::process::id()));
    let fasta = install(&text, &dir);

    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 4, "four runs are refused");

    for row in &errors {
        let label = row[0];
        let (_, _, error) = run(&text, &dir, &fasta, label);
        let error = error.expect("this run is refused");
        assert_eq!(
            format!(
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput:Bad input: {}",
                error.message()
            ),
            row[1],
            "error/{label}"
        );
    }
}

/// The overlap fixing, which is what makes the reported qualities not the reads'.
#[test]
fn an_overlapping_pair_comes_out_with_one_quality_zeroed() {
    let text = golden();
    let fixed = rows(&text, "report")
        .into_iter()
        .find(|row| row[0] == "overlapping")
        .expect("the overlapping run")[1]
        .to_string();
    // `r` is the sum of the two qualities and `!` is zero.
    assert!(unescape(&fixed).contains("chr1 21 C CC r! vs. chr1 21 C CC IJ"));

    // With the fixing off the same run agrees with the file and writes nothing.
    let ignored: Vec<Vec<&str>> = rows(&text, "report")
        .into_iter()
        .filter(|row| row[0] == "overlapping-ignored")
        .collect();
    assert_eq!(ignored.len(), 1);
    assert!(ignored[0].len() == 1 || ignored[0][1].is_empty());
}
