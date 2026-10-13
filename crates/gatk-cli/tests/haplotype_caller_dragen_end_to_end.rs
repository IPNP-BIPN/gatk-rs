//! Conformance for `HaplotypeCaller` in DRAGEN mode against GATK 4.6.2.0, every run of the dump
//! replayed through the command line and compared as the records and INFO lines the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/HcDragenEndToEndDump.java`: one contig with three
//! short tandem repeats, two samples' reads built so that BQD (an alternate only near the reads'
//! ends), FRD (an alternate only on reads of poor mapping quality) and DRAGstr (indels inside the
//! repeats) each move a call, and an XQ tag on every read. Fourteen runs in ONE JVM, in the dump's
//! order, since a `QD` above 35 draws from the one generator per process on both sides.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_dragen_end_to_end.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-haplotype-caller-dragen");
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
    let contigs = vec![("chr1", section(&text, "fasta", "chr1"))];
    let reference = write_reference(&dir, &contigs);
    let bam = write_bam(&dir, &section(&text, "sam", "input"));
    let table = dir.join("dragstr.txt");
    std::fs::write(&table, section(&text, "dragstr", "params")).expect("the DRAGstr table");
    let table = table.to_string_lossy().to_string();
    let runs: Vec<(&str, Vec<&str>)> = vec![
        ("default", vec![]),
        ("dragen-mode", vec!["--dragen-mode"]),
        (
            "dragen-mode-dragstr",
            vec!["--dragen-mode", "--dragstr-params-path", &table],
        ),
        ("apply-bqd", vec!["--apply-bqd"]),
        ("apply-frd", vec!["--apply-frd"]),
        (
            "apply-frd-depth-cap",
            vec![
                "--apply-frd",
                "--max-effective-depth-adjustment-for-frd",
                "4",
            ],
        ),
        ("transform-mq", vec!["--transform-dragen-mapping-quality"]),
        ("dragstr", vec!["--dragstr-params-path", &table]),
        (
            "dragstr-no-priors",
            vec!["--dragstr-params-path", &table, "--dont-use-dragstr-priors"],
        ),
        (
            "dragstr-no-hmm-scores",
            vec![
                "--dragstr-params-path",
                &table,
                "--dont-use-dragstr-pair-hmm-scores",
            ],
        ),
        (
            "dragstr-het-hom-ratio",
            vec![
                "--dragstr-params-path",
                &table,
                "--dragstr-het-hom-ratio",
                "4",
            ],
        ),
        (
            "dragen-mode-gvcf",
            vec!["--dragen-mode", "-ERC", "GVCF", "--sample-name", "s1"],
        ),
        (
            "dragen-mode-explicit-threshold",
            vec![
                "--dragen-mode",
                "--standard-min-confidence-threshold-for-calling",
                "60",
            ],
        ),
        (
            "dragen-mode-and-flow-mode",
            vec!["--dragen-mode", "--flow-mode", "STANDARD"],
        ),
    ];
    for (case, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            "HaplotypeCaller",
            "-I",
            &bam,
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
        let refusal = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{case}\t")));
        if let Some(refusal) = refusal {
            let message = refusal.split_once(':').expect("a class and a message").1;
            assert_ne!(outcome.status, 0, "{case} should be refused");
            assert!(
                outcome.stderr.contains(message),
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
                (in_body || line.starts_with("##INFO")) && !line.is_empty()
            })
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
    }
}
