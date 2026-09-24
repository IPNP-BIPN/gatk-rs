//! Conformance for `ReblockGVCF` against GATK 4.6.2.0, every run of the dump replayed through the
//! command line and compared as the records the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/ReblockGVCFDump.java`, which writes one GVCF (three
//! adjacent blocks, a confident variant, a weak one between two blocks, and a GQ0 block) and runs
//! the tool eleven times. The golden keeps each output from its `#CHROM` line on, and the one
//! refusal as its root cause.
//!
//! # What this suite is for
//!
//!  * **adjacent blocks in one band merging at the lowest quality**;
//!  * **the band edges deciding which blocks merge**;
//!  * **`--drop-low-quals` and `--rgq-threshold` touching different records**;
//!  * **a demoted variant becoming a GQ0 block that merges with its neighbours' band**;
//!  * **`--floor-blocks` writing the bound and dropping MIN_DP and PL**;
//!  * **and the two annotation arguments not being symmetric.**

use gatk_corpus as corpus;

/// `ReblockGVCFDump.CONTIG_LENGTH`.
const CONTIG_LENGTH: usize = 199_980;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/reblock_gvcf.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-reblock-gvcf");
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
    let input = dir.join("input.g.vcf");
    std::fs::write(&input, section(&text, "vcf", "input")).expect("the input");
    let bands: Vec<&str> = ["10", "20", "30", "40", "50", "60"]
        .iter()
        .flat_map(|band| ["--gvcf-gq-bands", band])
        .collect();
    let runs: Vec<(&str, Vec<&str>)> = vec![
        ("default", vec![]),
        ("one-band", vec!["--gvcf-gq-bands", "60"]),
        ("many-bands", bands),
        ("rgq-threshold", vec!["--rgq-threshold", "10"]),
        ("drop-low-quals", vec!["--drop-low-quals", "true"]),
        (
            "drop-and-threshold",
            vec!["--drop-low-quals", "true", "--rgq-threshold", "10"],
        ),
        ("keep-all-alts", vec!["--keep-all-alts", "true"]),
        ("floor-blocks", vec!["--floor-blocks", "true"]),
        ("keep-annotation", vec!["--annotations-to-keep", "EXTRA"]),
        (
            "remove-annotation",
            vec!["--format-annotations-to-remove", "SPARE"],
        ),
        ("keep-format-key", vec!["--annotations-to-keep", "SPARE"]),
    ];
    for (case, extra) in runs {
        let output = dir.join(format!("out-{case}.g.vcf"));
        let mut args: Vec<String> = [
            "ReblockGVCF",
            "-V",
            &input.to_string_lossy(),
            "-O",
            &output.to_string_lossy(),
            "-R",
            &reference,
        ]
        .iter()
        .map(|arg| arg.to_string())
        .collect();
        args.extend(extra.iter().map(|arg| arg.to_string()));
        let outcome = gatk_cli::run(&args);
        if case == "keep-format-key" {
            let row = text
                .lines()
                .find_map(|line| line.strip_prefix("error\tkeep-format-key\t"))
                .expect("the golden carries the refusal");
            let (_, message) = row.split_once(':').expect("a class and a message");
            assert_ne!(outcome.status, 0, "{case}");
            assert!(
                outcome.stderr.contains(message),
                "{case}: {}",
                outcome.stderr
            );
            continue;
        }
        assert_eq!(outcome.status, 0, "{case}: {}", outcome.stderr);
        let written = std::fs::read_to_string(&output).expect("the output");
        let body: String = written
            .lines()
            .filter(|line| !line.starts_with("##") && !line.is_empty())
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
    }
}
