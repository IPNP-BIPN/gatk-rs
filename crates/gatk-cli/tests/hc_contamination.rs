//! Conformance for contamination downsampling in HaplotypeCaller and RampedHaplotypeCaller against
//! GATK 4.6.2.0, every run of the dump replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcContaminationDump.java`, over `HcEndToEndDump`'s
//! two samples: `--contamination-fraction-to-filter` at five fractions, with
//! `--use-filtered-reads-for-annotations` and in GVCF mode, per-sample files with and without a
//! default, the file's refusals, and RampedHaplotypeCaller. Each run starts the random data
//! generator afresh, as a new JVM does. Each output is compared from its ##INFO and ##FORMAT lines
//! to its last record; a refusal by its message, the dump's directory masked as `<dir>`.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_contamination.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-contamination");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// `writeReference`: each contig in lines of sixty, with the index and dictionary.
fn write_reference(dir: &std::path::Path, contigs: &[(&str, String)]) -> String {
    let mut fasta = String::new();
    let mut fai = String::new();
    let mut dict = String::from("@HD\tVN:1.6\n");
    for (name, bases) in contigs {
        fasta.push_str(&format!(">{name}\n"));
        let offset = fasta.len();
        for line in bases.as_bytes().chunks(60) {
            fasta.push_str(std::str::from_utf8(line).expect("ASCII"));
            fasta.push('\n');
        }
        fai.push_str(&format!("{name}\t{}\t{offset}\t60\t61\n", bases.len()));
        dict.push_str(&format!("@SQ\tSN:{name}\tLN:{}\n", bases.len()));
    }
    std::fs::write(dir.join("reference.fasta"), fasta).expect("the reference");
    std::fs::write(dir.join("reference.fasta.fai"), fai).expect("the index");
    std::fs::write(dir.join("reference.dict"), dict).expect("the dictionary");
    dir.join("reference.fasta").to_string_lossy().to_string()
}

/// The BAM and its index, written from the golden's SAM text.
fn write_bam(dir: &std::path::Path, sam: &str) -> String {
    let (header, records) = htsjdk_bam::sam_file::read_sam(sam).expect("the SAM text");
    let mut writer = htsjdk_bam::writer::BamWriter::new(Vec::new(), &header)
        .expect("a writer")
        .with_index();
    for record in &records {
        writer.write(record).expect("a record");
    }
    let (bam, bai) = writer.finish_with_index().expect("the BAM");
    std::fs::write(dir.join("input.bam"), bam).expect("the BAM file");
    std::fs::write(dir.join("input.bai"), bai).expect("the index file");
    dir.join("input.bam").to_string_lossy().to_string()
}

#[test]
fn every_run_matches_the_golden_in_the_dumps_order() {
    let text = golden();
    let dir = scratch();
    let e2e_contigs: Vec<(&str, String)> = ["chr1", "chr2"]
        .iter()
        .map(|&contig| (contig, section(&text, "fasta", &format!("e2e:{contig}"))))
        .collect();
    let sub = dir.join("e2e");
    std::fs::create_dir_all(&sub).expect("a set directory");
    let e2e_reference = write_reference(&sub, &e2e_contigs);
    let e2e_bam = write_bam(&sub, &section(&text, "sam", "e2e"));
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let names = [
        "per-sample",
        "spaced",
        "three-columns",
        "empty-field",
        "unparsable",
        "duplicate",
        "duplicate-other",
        "out-of-range",
        "nan",
    ];
    let files: Vec<String> = names
        .iter()
        .map(|name| {
            let path = dir.join(format!("{name}.txt"));
            std::fs::write(&path, section(&text, "file", name)).expect("a contamination file");
            path.to_string_lossy().to_string()
        })
        .collect();
    let path = |name: &str| files[names.iter().position(|n| *n == name).expect("a file")].as_str();
    let missing = dir.join("absent.txt").to_string_lossy().to_string();
    let fraction = "--contamination-fraction-to-filter";
    let file = "--contamination-fraction-per-sample-file";
    let hc = "HaplotypeCaller";
    let runs: Vec<Run<'_>> = vec![
        ("hc-0.05", hc, e2e, vec![fraction, "0.05"]),
        ("hc-0.1", hc, e2e, vec![fraction, "0.1"]),
        ("hc-0.3", hc, e2e, vec![fraction, "0.3"]),
        ("hc-all", hc, e2e, vec![fraction, "1.0"]),
        ("hc-zero", hc, e2e, vec![fraction, "0.0"]),
        (
            "hc-filtered-annotations",
            hc,
            e2e,
            vec![fraction, "0.2", "--use-filtered-reads-for-annotations"],
        ),
        (
            "hc-gvcf",
            hc,
            e2e,
            vec![
                fraction,
                "0.2",
                "-ERC",
                "GVCF",
                "--sample-name",
                "s1",
                "-L",
                "chr1",
            ],
        ),
        ("hc-per-sample", hc, e2e, vec![file, path("per-sample")]),
        (
            "hc-per-sample-default",
            hc,
            e2e,
            vec![file, path("spaced"), fraction, "0.1"],
        ),
        (
            "hc-duplicate-other",
            hc,
            e2e,
            vec![file, path("duplicate-other")],
        ),
        ("hc-nan", hc, e2e, vec![file, path("nan")]),
        (
            "hc-three-columns",
            hc,
            e2e,
            vec![file, path("three-columns")],
        ),
        ("hc-empty-field", hc, e2e, vec![file, path("empty-field")]),
        ("hc-unparsable", hc, e2e, vec![file, path("unparsable")]),
        ("hc-duplicate", hc, e2e, vec![file, path("duplicate")]),
        ("hc-out-of-range", hc, e2e, vec![file, path("out-of-range")]),
        ("hc-missing", hc, e2e, vec![file, missing.as_str()]),
        (
            "ramped-0.2",
            "RampedHaplotypeCaller",
            e2e,
            vec![fraction, "0.2"],
        ),
        (
            "ramped-gvcf",
            "RampedHaplotypeCaller",
            e2e,
            vec![
                fraction,
                "0.2",
                "-ERC",
                "GVCF",
                "--sample-name",
                "s2",
                "-L",
                "chr2",
            ],
        ),
    ];
    let shown = dir.to_string_lossy().to_string();
    for (case, tool, (bam, reference), extra) in &runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            *tool,
            "-I",
            bam,
            "-O",
            &output.to_string_lossy(),
            "-R",
            reference,
        ]
        .iter()
        .map(|arg| arg.to_string())
        .collect();
        args.extend(extra.iter().map(|arg| arg.to_string()));
        let outcome = gatk_cli::run(&args);
        let refusal = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{case}\t")));
        if let Some(refusal) = refusal {
            let message = unescape(refusal.split_once(':').expect("a class and a message").1)
                .replace("<dir>", &shown);
            assert_ne!(outcome.status, 0, "{case} should be refused");
            assert!(
                outcome.stderr.contains(&message),
                "{case}: {}",
                outcome.stderr
            );
            continue;
        }
        assert_eq!(outcome.status, 0, "{case}: {}", outcome.stderr);
        let written = std::fs::read_to_string(&output).expect("the output");
        let mut in_body = false;
        let body: String = written
            .lines()
            .filter(|line| {
                in_body |= line.starts_with("#CHROM");
                (in_body || line.starts_with("##INFO") || line.starts_with("##FORMAT"))
                    && !line.is_empty()
            })
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
    }
}
