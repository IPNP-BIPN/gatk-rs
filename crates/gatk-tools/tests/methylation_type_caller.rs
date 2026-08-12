//! Conformance for `MethylationTypeCaller` against GATK 4.6.2.0, compared as **VCF text**.
//!
//! Golden from `tools/readfilter-conformance/MethylationTypeCallerDump.java`. The two input BAMs and
//! their indexes travel in full, base64, and the reference is written out from the golden's own row,
//! so the port reads exactly what the reference read.
//!
//! # What this suite is for
//!
//!  * **the strand is chosen by the reference base**, a C from forward reads and a G from reverse
//!    ones;
//!  * **DP is the whole pileup** while the two counts are one strand and two bases;
//!  * **a site with no methylated coverage writes nothing**;
//!  * **the contexts are three bases on their own strand**, truncated near a contig edge;
//!  * **an unmapped read is filtered out and a deletion contributes no base**;
//!  * **and the samples are the read groups' samples, sorted then deduplicated**, with no genotypes
//!    on the records.
//!
//! The run with `--add-output-vcf-command-line` left on is **not** compared: its header carries the
//! run's own date. Its shape is in the golden as a masked row and checked here as a shape.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use gatk_tools::methylation_type_caller::methylation_type_caller;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/methylation_type_caller.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn reference_bases(text: &str) -> String {
    rows(text, "reference")[0][0].to_string()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
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

/// The fixture and the intervals each labelled run used.
fn configuration(label: &str) -> (&str, Option<Vec<SimpleInterval>>) {
    match label {
        "plain" => ("bisulfite", None),
        "interval" => (
            "bisulfite",
            Some(vec![
                SimpleInterval::new("chr1", 18, 24).expect("a valid interval")
            ]),
        ),
        "two-samples" => ("samples", None),
        other => panic!("no comparable run {other}"),
    }
}

fn run(dir: &std::path::Path, fasta: &std::path::Path, label: &str) -> String {
    let (fixture, intervals) = configuration(label);
    let bam = dir.join(format!("{fixture}.bam"));
    let bai = dir.join(format!("{fixture}.bai"));
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let mut reference = ReferenceFileSource::open(fasta).expect("the reference opens");

    methylation_type_caller(
        &source,
        &mut reference,
        intervals.as_deref(),
        // `--add-output-vcf-command-line false`, which is what makes the file comparable.
        Vec::new(),
    )
    .expect("the run finishes")
}

#[test]
fn every_comparable_run_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-methylation-{}", std::process::id()));
    let fasta = install(&text, &dir);

    let mut compared = 0;
    for label in ["plain", "interval", "two-samples"] {
        let expected: Vec<String> = rows(&text, "vcfline")
            .into_iter()
            .filter(|row| row[0] == label)
            .map(|row| unescape(row[1]))
            .collect();
        assert!(!expected.is_empty(), "the golden holds the {label} run");

        let ours: Vec<String> = run(&dir, &fasta, label)
            .lines()
            .map(|line| line.to_string())
            .collect();
        assert_eq!(ours, expected, "{label}");
        compared += ours.len();
    }
    println!("methylation-type-caller: {compared} VCF lines compared");
}

/// The site near the contig end, whose context is shorter than everywhere else.
#[test]
fn a_site_near_the_contig_end_gets_a_shorter_context() {
    let text = golden();
    let lines: Vec<String> = rows(&text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == "plain")
        .map(|row| unescape(row[1]))
        .collect();

    let context_of = |position: &str| -> String {
        lines
            .iter()
            .find(|line| line.starts_with(&format!("chr1\t{position}\t")))
            .unwrap_or_else(|| panic!("no record at {position}"))
            .split('\t')
            .nth(7)
            .expect("the INFO column")
            .split(';')
            .find(|field| field.starts_with("REFERENCE_CONTEXT="))
            .expect("the context")
            .trim_start_matches("REFERENCE_CONTEXT=")
            .to_string()
    };

    // Three bases in the middle of the contig, two at the second to last base, one at the last.
    assert_eq!(context_of("30").len(), 3);
    assert_eq!(context_of("59"), "CC");
    assert_eq!(context_of("60"), "C");
}

/// The command-line header line the default run writes, which carries the run's own date.
#[test]
fn the_command_line_header_carries_a_date_no_golden_can_hold() {
    let text = golden();
    let rows = rows(&text, "commandline");
    assert_eq!(rows.len(), 1, "one run leaves the command line on");
    let line = unescape(rows[0][1]);
    assert!(
        line.starts_with("##GATKCommandLine=<ID=MethylationTypeCaller,"),
        "{line}"
    );
    // The dump masked it, which is the whole reason the comparable runs turn the line off.
    assert!(line.contains("Date=\"MASKED\""), "{line}");
}
