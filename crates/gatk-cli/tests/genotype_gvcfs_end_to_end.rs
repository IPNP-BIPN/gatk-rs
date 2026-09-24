//! Conformance for `GenotypeGVCFs` against GATK 4.6.2.0, every run of the dump replayed through the
//! command line and compared as the records the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/GenotypeGVCFsDump.java`, which writes one GVCF (a
//! reference block, a het, a hom-var, a reference-best site, a marginal one and a triallelic one)
//! and runs the tool six times in ONE JVM. The golden keeps each output from its `#CHROM` line on.
//!
//! # Why the runs are replayed in order, in one process
//!
//! Every `QD` in the golden is above 35 before the reference replaces it with a Gaussian draw from
//! `Utils.getRandomGenerator()`, which is seeded once per JVM and never reset between the six runs.
//! The port's generator is likewise one per process, so the runs here go in the dump's order and
//! this file holds the only test that drives the tool: a second test in the same binary could run
//! first and move the stream.
//!
//! # What this suite is for
//!
//!  * **a reference block never written, and `<NON_REF>` removed from every site that is**;
//!  * **a site whose best genotype is the reference dropped, and written by `-all-sites`**;
//!  * **the calling threshold deciding `LowQual`, not emission**;
//!  * **an alternate no sample carries dropped, with the likelihoods re-indexed around it**;
//!  * **`AC`, `AN`, `AF`, `ExcessHet` and `QD` computed from the called genotypes**;
//!  * **and `--sample-ploidy` changing nothing when the likelihoods are diploid.**

use gatk_corpus as corpus;

/// `GenotypeGVCFsDump.CONTIG_LENGTH`.
const CONTIG_LENGTH: usize = 199_980;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/genotype_gvcfs.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-genotype-gvcfs");
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
fn every_run_matches_the_golden_in_the_dumps_order() {
    let text = golden();
    let dir = scratch();
    let reference = write_reference(&dir);
    let input = dir.join("input.g.vcf");
    std::fs::write(&input, section(&text, "vcf", "input")).expect("the input");
    let runs: [(&str, &[&str]); 6] = [
        ("default", &[]),
        ("all-sites", &["--include-non-variant-sites", "true"]),
        (
            "call-threshold-2",
            &["--standard-min-confidence-threshold-for-calling", "2"],
        ),
        (
            "call-threshold-50",
            &["--standard-min-confidence-threshold-for-calling", "50"],
        ),
        (
            "keep-combined",
            &["--keep-combined-raw-annotations", "true"],
        ),
        ("ploidy-one", &["--sample-ploidy", "1"]),
    ];
    for (case, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            "GenotypeGVCFs",
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
