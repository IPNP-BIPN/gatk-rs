//! Conformance for `VariantEval` against GATK 4.6.2.0, every run of the dump replayed through the
//! command line and compared as the whole report the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/VariantEvalDump.java`, which writes one eval set
//! (transitions, transversions, an insertion, a deletion and a multiallelic site) and a comparison
//! set holding three of its positions, and runs the tool twelve times. The golden keeps each report
//! whole and each refusal as its root cause. The thirteenth run, `--list`, ends the process and is
//! a limitation of the port.
//!
//! # What this suite is for
//!
//!  * **the report's formatting, column widths and derived rates, byte for byte**;
//!  * **novelty coming from dbSNP and not from `--comp`**;
//!  * **the standard stratifiers turned off, and one stratifier multiplying the rows**;
//!  * **a select expression adding its own stratum**;
//!  * **and an unknown module or stratifier refused as a command line error.**

use gatk_corpus as corpus;

/// `VariantEvalDump.CONTIG_LENGTH`.
const CONTIG_LENGTH: usize = 199_980;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/variant_eval.txt.gz"),
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

fn scratch() -> std::path::PathBuf {
    let dir = std::env::temp_dir().join("gatk-cli-variant-eval");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// `writeReference`: `ACGT` repeated in lines of sixty, with its index and dictionary.
fn write_reference(dir: &std::path::Path) -> String {
    let mut fasta = String::from(">chr1\n");
    for _ in 0..CONTIG_LENGTH / 60 {
        fasta.push_str(&"ACGT".repeat(15));
        fasta.push('\n');
    }
    std::fs::write(dir.join("reference.fasta"), fasta).expect("the reference");
    std::fs::write(
        dir.join("reference.fasta.fai"),
        format!("chr1\t{CONTIG_LENGTH}\t6\t60\t61\n"),
    )
    .expect("the index");
    std::fs::write(
        dir.join("reference.dict"),
        format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:{CONTIG_LENGTH}\n"),
    )
    .expect("the dictionary");
    dir.join("reference.fasta").to_string_lossy().to_string()
}

#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let dir = scratch();
    let reference = write_reference(&dir);
    let eval = dir.join("eval.vcf");
    std::fs::write(&eval, section(&text, "vcf", "eval")).expect("the eval set");
    let comp = dir.join("comp.vcf");
    std::fs::write(&comp, section(&text, "vcf", "comp")).expect("the comparison set");
    let comp = comp.to_string_lossy().to_string();
    let only =
        |module: &'static str| vec!["--do-not-use-all-standard-modules", "true", "-EV", module];
    let runs: Vec<(&str, Vec<&str>)> = vec![
        ("no-comp", vec![]),
        (
            "dbsnp",
            [vec!["--dbsnp", comp.as_str()], only("CountVariants")].concat(),
        ),
        ("with-comp", vec!["--comp", comp.as_str()]),
        (
            "count-variants",
            [vec!["--comp", comp.as_str()], only("CountVariants")].concat(),
        ),
        (
            "titv",
            [vec!["--comp", comp.as_str()], only("TiTvVariantEvaluator")].concat(),
        ),
        ("indel-length", only("IndelLengthHistogram")),
        ("multiallelic", only("MultiallelicSummary")),
        (
            "no-standard-stratifiers",
            [
                vec!["--comp", comp.as_str()],
                only("CountVariants"),
                vec!["--do-not-use-all-standard-stratifications", "true"],
            ]
            .concat(),
        ),
        (
            "stratify-by-type",
            [only("CountVariants"), vec!["-ST", "VariantType"]].concat(),
        ),
        (
            "select-expression",
            [
                only("CountVariants"),
                vec!["-select", "QUAL > 50", "-select-name", "highqual"],
            ]
            .concat(),
        ),
        ("unknown-module", vec!["-EV", "NoSuchEvaluator"]),
        ("unknown-stratifier", vec!["-ST", "NoSuchStratifier"]),
    ];
    for (case, extra) in runs {
        let output = dir.join(format!("out-{case}.txt"));
        let mut args: Vec<String> = [
            "VariantEval",
            "-O",
            &output.to_string_lossy(),
            "-R",
            &reference,
            "--eval",
            &eval.to_string_lossy(),
        ]
        .iter()
        .map(|arg| arg.to_string())
        .collect();
        args.extend(extra.iter().map(|arg| arg.to_string()));
        let outcome = gatk_cli::run(&args);
        if let Some(row) = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{case}\t")))
        {
            let (_, message) = row.split_once(':').expect("a class and a message");
            assert_eq!(outcome.status, 1, "{case}: {}", outcome.stderr);
            assert!(
                outcome.stderr.contains(message),
                "{case}: {}",
                outcome.stderr
            );
            continue;
        }
        assert_eq!(outcome.status, 0, "{case}: {}", outcome.stderr);
        let written = std::fs::read_to_string(&output).expect("the report");
        assert_eq!(written, section(&text, "out", case), "{case}");
    }
}
