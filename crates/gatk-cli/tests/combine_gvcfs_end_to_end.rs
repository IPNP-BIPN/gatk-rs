//! Conformance for `CombineGVCFs` against GATK 4.6.2.0, every run of the dump replayed through the
//! command line and compared as the records the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/CombineGVCFsDump.java`, which writes three
//! single-sample GVCFs and a reference of repeated `ACGT`, then runs the tool nine times. The
//! golden keeps each output from its `#CHROM` line on, so the header lines above it are not
//! compared here; the covering array compares them.
//!
//! # What this suite is for
//!
//!  * **the records being the union of every input's edges, each sample keeping its quality**;
//!  * **a sample that stops early keeping its column as `./.`**;
//!  * **the likelihoods of the samples with no variant expanded to the merged allele set**;
//!  * **the two band arguments, base-pair resolution winning over the grid**;
//!  * **`--call-genotypes` turning the no-calls into calls**;
//!  * **and the same file twice being refused.**

use gatk_corpus as corpus;

/// `CombineGVCFsDump.CONTIG_LENGTH`.
const CONTIG_LENGTH: usize = 199_980;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/combine_gvcfs.txt.gz"),
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

/// A directory of this test's own, named after the case so two cases cannot collide.
fn scratch(case: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("gatk-cli-combine-gvcfs-{case}"));
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

/// Runs the tool over the named samples with the extra arguments, and returns the output from its
/// `#CHROM` line on, as the dump kept it.
fn combined(text: &str, case: &str, samples: &[&str], extra: &[&str]) -> Result<String, String> {
    let dir = scratch(case);
    let reference = write_reference(&dir);
    let mut args: Vec<String> = vec!["CombineGVCFs".to_string()];
    for sample in samples {
        let path = dir.join(format!("{sample}.g.vcf"));
        std::fs::write(&path, section(text, "vcf", sample)).expect("an input");
        args.push("-V".to_string());
        args.push(path.to_string_lossy().to_string());
    }
    let output = dir.join(format!("out-{case}.g.vcf"));
    args.extend([
        "-O".to_string(),
        output.to_string_lossy().to_string(),
        "-R".to_string(),
        reference,
    ]);
    args.extend(extra.iter().map(|arg| arg.to_string()));
    let outcome = gatk_cli::run(&args);
    if outcome.status != 0 {
        return Err(outcome
            .stderr
            .replace(&dir.to_string_lossy().to_string(), "<dir>"));
    }
    let written = std::fs::read_to_string(&output).expect("the output");
    Ok(written
        .lines()
        .filter(|line| !line.starts_with("##") && !line.is_empty())
        .map(|line| format!("{line}\n"))
        .collect())
}

#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let three = ["s1", "s2", "s3"];
    let runs: [(&str, &[&str], &[&str]); 8] = [
        ("three-samples", &three, &[]),
        ("two-samples", &["s1", "s2"], &[]),
        (
            "base-pair-resolution",
            &three,
            &["--convert-to-base-pair-resolution", "true"],
        ),
        (
            "break-bands-100",
            &three,
            &["--break-bands-at-multiples-of", "100"],
        ),
        (
            "break-bands-50",
            &three,
            &["--break-bands-at-multiples-of", "50"],
        ),
        (
            "both-band-arguments",
            &three,
            &[
                "--convert-to-base-pair-resolution",
                "true",
                "--break-bands-at-multiples-of",
                "100",
            ],
        ),
        ("call-genotypes", &three, &["--call-genotypes", "true"]),
        ("one-sample", &["s1"], &[]),
    ];
    for (case, samples, extra) in runs {
        assert_eq!(
            combined(&text, case, samples, extra).as_deref(),
            Ok(section(&text, "out", case).as_str()),
            "{case}"
        );
    }
}

/// Base-pair resolution given with a grid produces what base-pair resolution produced alone.
#[test]
fn base_pair_resolution_wins_over_the_grid() {
    let text = golden();
    assert_eq!(
        section(&text, "out", "both-band-arguments"),
        section(&text, "out", "base-pair-resolution")
    );
}

/// The same file twice is refused by the feature-input check, before any sample is looked at.
#[test]
fn the_same_file_twice_is_refused() {
    let text = golden();
    let recorded = text
        .lines()
        .find_map(|line| line.strip_prefix("error\tduplicate-sample\t"))
        .expect("the golden carries error/duplicate-sample");
    let (_, message) = recorded.split_once(':').expect("a class and a message");
    let refused =
        combined(&text, "duplicate-sample", &["s1", "s1"], &[]).expect_err("the run is refused");
    assert!(refused.contains(message), "{refused}");
}
