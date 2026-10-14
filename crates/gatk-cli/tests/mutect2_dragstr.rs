//! Conformance for Mutect2's `--dragstr-params-path` and `--dragstr-het-hom-ratio` against GATK
//! 4.6.2.0, every run of the dump replayed through the command line.
//!
//! Golden from `tools/readfilter-conformance/Mutect2DragstrDump.java`: DRAGstr's default table over
//! reads with an A homopolymer and a CA repeat, with a het-hom ratio, with the PairHMM scores off
//! (a missing table then unread), with a normal, and a missing table refused. Each output is
//! compared from its header lines (the command line aside) on; a refusal by its message.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/mutect2_dragstr.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-mutect2-dragstr");
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
fn write_bam(dir: &std::path::Path, sam: &str, name: &str) -> String {
    let (header, records) = htsjdk_bam::sam_file::read_sam(sam).expect("the SAM text");
    let mut writer = htsjdk_bam::writer::BamWriter::new(Vec::new(), &header)
        .expect("a writer")
        .with_index();
    for record in &records {
        writer.write(record).expect("a record");
    }
    let (bam, bai) = writer.finish_with_index().expect("the BAM");
    std::fs::write(dir.join(format!("{name}.bam")), bam).expect("the BAM file");
    std::fs::write(dir.join(format!("{name}.bai")), bai).expect("the index file");
    dir.join(format!("{name}.bam"))
        .to_string_lossy()
        .to_string()
}

#[test]
fn every_run_matches_the_golden_in_the_dumps_order() {
    let text = golden();
    let dir = scratch();
    let contigs = vec![("chr1", section(&text, "fasta", "chr1"))];
    let reference = write_reference(&dir, &contigs);
    let one = write_bam(&dir, &section(&text, "sam", "input"), "input");
    let both = write_bam(&dir, &section(&text, "sam", "both"), "both");
    let params_path = dir.join("dragstr.txt");
    std::fs::write(&params_path, section(&text, "dragstr", "params")).expect("the table");
    let params = params_path.to_string_lossy().to_string();
    let missing = dir.join("absent.txt").to_string_lossy().to_string();
    let table = "--dragstr-params-path";
    let off = "--dont-use-dragstr-pair-hmm-scores";
    let runs: Vec<(&str, &str, Vec<&str>)> = vec![
        ("plain", &one, vec![]),
        ("dragstr", &one, vec![table, &params]),
        (
            "dragstr-ratio",
            &one,
            vec![table, &params, "--dragstr-het-hom-ratio", "3"],
        ),
        ("dragstr-off", &one, vec![table, &params, off]),
        ("dragstr-off-missing", &one, vec![table, &missing, off]),
        (
            "dragstr-tumor-normal",
            &both,
            vec![table, &params, "--normal-sample", "s2"],
        ),
        ("dragstr-missing", &one, vec![table, &missing]),
    ];
    let shown = dir.to_string_lossy().to_string();
    for (case, bam, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            "Mutect2",
            "-I",
            bam,
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
                (in_body || (line.starts_with("##") && !line.starts_with("##GATKCommandLine")))
                    && !line.is_empty()
            })
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
    }
}
