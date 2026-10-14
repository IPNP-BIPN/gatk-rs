//! Conformance for Mutect2's `--f1r2-tar-gz` against GATK 4.6.2.0, every run of the dump replayed
//! through the command line.
//!
//! Golden from `tools/readfilter-conformance/Mutect2F1r2Dump.java`: the collector's defaults and
//! thresholds, a median MQ no read reaches, a normal, two tumours, GVCF mode, an interval at the
//! contig's start, and a name not ending in `.tar.gz`. Each output is compared from its header
//! lines (the command line aside) on, and each archive member by name.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/mutect2_f1r2.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-mutect2-f1r2");
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

/// The archive's file members, by name: a gzip of a ustar stream.
fn members(path: &std::path::Path) -> std::collections::BTreeMap<String, String> {
    use std::io::Read;
    let mut tar = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).expect("the archive"))
        .read_to_end(&mut tar)
        .expect("a gzip");
    let mut out = std::collections::BTreeMap::new();
    let mut offset = 0;
    while offset + 512 <= tar.len() && tar[offset] != 0 {
        let header = &tar[offset..offset + 512];
        let field = |from: usize, len: usize| {
            String::from_utf8_lossy(&header[from..from + len])
                .trim_end_matches('\0')
                .trim()
                .to_string()
        };
        let name = field(0, 100);
        let size = usize::from_str_radix(&field(124, 12), 8).expect("an octal size");
        let kind = header[156];
        offset += 512;
        if kind == b'0' || kind == 0 {
            out.insert(
                name,
                String::from_utf8_lossy(&tar[offset..offset + size]).into_owned(),
            );
        }
        offset += size.div_ceil(512) * 512;
    }
    out
}

#[test]
fn every_run_matches_the_golden_in_the_dumps_order() {
    let text = golden();
    let dir = scratch();
    let contigs = vec![("chr1", section(&text, "fasta", "chr1"))];
    let reference = write_reference(&dir, &contigs);
    let one = write_bam(&dir, &section(&text, "sam", "input"), "input");
    let both = write_bam(&dir, &section(&text, "sam", "both"), "both");
    let runs: Vec<(&str, &str, &str, Vec<&str>)> = vec![
        ("f1r2", &one, "f1r2.tar.gz", vec![]),
        (
            "f1r2-thresholds",
            &one,
            "f1r2-thresholds.tar.gz",
            vec![
                "--f1r2-median-mq",
                "30",
                "--f1r2-min-bq",
                "10",
                "--f1r2-max-depth",
                "20",
            ],
        ),
        (
            "f1r2-strict",
            &one,
            "f1r2-strict.tar.gz",
            vec!["--f1r2-median-mq", "61"],
        ),
        (
            "f1r2-tumor-normal",
            &both,
            "f1r2-tumor-normal.tar.gz",
            vec!["--normal-sample", "s2"],
        ),
        ("f1r2-both-tumors", &both, "f1r2-both-tumors.tar.gz", vec![]),
        (
            "f1r2-gvcf",
            &one,
            "f1r2-gvcf.tar.gz",
            vec!["-ERC", "GVCF", "-L", "chr1:380-460"],
        ),
        (
            "f1r2-interval",
            &one,
            "f1r2-interval.tar.gz",
            vec!["-L", "chr1:1-5"],
        ),
        ("f1r2-bad-name", &one, "f1r2.tgz", vec![]),
    ];
    let shown = dir.to_string_lossy().to_string();
    for (case, bam, archive, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let tar_gz = dir.join(archive);
        let mut args: Vec<String> = [
            "Mutect2",
            "-I",
            bam,
            "-O",
            &output.to_string_lossy(),
            "-R",
            &reference,
            "--f1r2-tar-gz",
            &tar_gz.to_string_lossy(),
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
        let expected: std::collections::BTreeMap<String, String> = text
            .lines()
            .filter_map(|line| line.strip_prefix(&format!("member\t{case}\t")))
            .map(|rest| {
                let (name, content) = rest.split_once('=').expect("a member");
                (name.to_string(), unescape(content))
            })
            .collect();
        assert_eq!(members(&tar_gz), expected, "{case}");
    }
}
