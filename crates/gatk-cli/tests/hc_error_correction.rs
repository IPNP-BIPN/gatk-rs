//! Conformance for read error correction in the callers against GATK 4.6.2.0, every run of the
//! dump replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcErrorCorrectionDump.java`, over `HcEndToEndDump`'s
//! reads and the same reads with three per cent of their bases miscalled: HaplotypeCaller with the
//! pileup corrector at five thresholds, the k-mer corrector (which corrects nothing but changes the
//! regions' tail trimming) and its refusal of a k-mer length of zero, in GVCF mode and with the
//! linked graph, RampedHaplotypeCaller, and Mutect2 tumour-only and tumour-normal. Each output is
//! compared from its ##INFO and ##FORMAT lines to its last record.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_error_correction.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-error-correction");
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
    let set = |name: &str| -> (String, String) {
        let sub = dir.join(name);
        std::fs::create_dir_all(&sub).expect("a set directory");
        let reference = write_reference(&sub, &e2e_contigs);
        let bam = write_bam(&sub, &section(&text, "sam", name));
        (reference, bam)
    };
    let (e2e_reference, e2e_bam) = set("e2e");
    let (noisy_reference, noisy_bam) = set("noisy");
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let noisy = (noisy_bam.as_str(), noisy_reference.as_str());
    let odds = "--error-correction-log-odds";
    let kmer = "--error-correct-reads";
    let length = "--kmer-length-for-read-error-correction";
    let mut runs: Vec<Run<'_>> = vec![
        ("hc-noisy", "HaplotypeCaller", noisy, vec![]),
        ("hc-kmer-corrector", "HaplotypeCaller", noisy, vec![kmer]),
        (
            "hc-kmer-corrector-tuned",
            "HaplotypeCaller",
            noisy,
            vec![
                kmer,
                length,
                "10",
                "--min-observations-for-kmer-to-be-solid",
                "2",
            ],
        ),
        (
            "hc-kmer-corrector-zero",
            "HaplotypeCaller",
            noisy,
            vec![kmer, length, "0"],
        ),
        (
            "hc-kmer-corrector-zero-pileup",
            "HaplotypeCaller",
            noisy,
            vec![kmer, length, "0", odds, "0.0"],
        ),
    ];
    let thresholds: Vec<(String, &str)> = ["-5.0", "0.0", "2.0", "10.0", "1000.0"]
        .iter()
        .map(|t| (format!("hc-pileup-corrector-{t}"), *t))
        .collect();
    for (label, threshold) in &thresholds {
        runs.push((
            label.as_str(),
            "HaplotypeCaller",
            noisy,
            vec![odds, threshold],
        ));
    }
    runs.extend([
        (
            "hc-pileup-corrector-clean",
            "HaplotypeCaller",
            e2e,
            vec![odds, "1000.0"],
        ),
        (
            "hc-pileup-corrector-gvcf",
            "HaplotypeCaller",
            noisy,
            vec![
                odds,
                "10.0",
                kmer,
                "-ERC",
                "GVCF",
                "--sample-name",
                "s1",
                "-L",
                "chr1",
            ],
        ),
        (
            "hc-pileup-corrector-linked",
            "HaplotypeCaller",
            noisy,
            vec![odds, "10.0", "--linked-de-bruijn-graph"],
        ),
        (
            "ramped-pileup-corrector",
            "RampedHaplotypeCaller",
            noisy,
            vec![odds, "10.0"],
        ),
        (
            "ramped-kmer-corrector",
            "RampedHaplotypeCaller",
            noisy,
            vec![kmer],
        ),
        ("mutect2-noisy", "Mutect2", noisy, vec![]),
        (
            "mutect2-pileup-corrector",
            "Mutect2",
            noisy,
            vec![odds, "10.0"],
        ),
        (
            "mutect2-pileup-corrector-tn",
            "Mutect2",
            noisy,
            vec![odds, "2.0", "--normal-sample", "s2"],
        ),
        ("mutect2-kmer-corrector", "Mutect2", noisy, vec![kmer]),
        (
            "mutect2-kmer-corrector-zero",
            "Mutect2",
            noisy,
            vec![kmer, length, "0"],
        ),
    ]);
    for (case, tool, (bam, reference), extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            tool,
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
                (in_body || line.starts_with("##INFO") || line.starts_with("##FORMAT"))
                    && !line.is_empty()
            })
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
    }
}
