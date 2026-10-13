//! Conformance for `RampedHaplotypeCaller`'s off ramps against GATK 4.6.2.0, every run of the dump
//! replayed through the command line and compared as the VCF and the zip the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/RampsEndToEndDump.java`: `HcEndToEndDump`'s two
//! samples and two contigs, run with each off ramp, then with each on ramp reading back a zip the
//! off ramps wrote, in ONE JVM. The golden carries the contigs, the
//! SAM text, each output from its `#CHROM` line on, each zip's entry names in order, each text
//! entry, and the haplotype BAM and its index as hex; the post-assembler run writes them through
//! the JDK deflater, so the bytes are the same on every machine.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/ramps_end_to_end.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-ramps");
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

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[test]
fn every_ramp_matches_the_golden_in_the_dumps_order() {
    let text = golden();
    let dir = scratch();
    let contigs: Vec<(&str, String)> = ["chr1", "chr2"]
        .iter()
        .map(|name| (*name, section(&text, "fasta", name)))
        .collect();
    let reference = write_reference(&dir, &contigs);
    let bam = write_bam(&dir, &section(&text, "sam", "input"));
    let zip_of = |label: &str| {
        dir.join(format!("{label}.zip"))
            .to_string_lossy()
            .to_string()
    };
    let post_assembler = zip_of("post-assembler-off");
    let pre_filter = zip_of("pre-filter-off");
    let pre_assembler = zip_of("pre-assembler-off");
    let unused = zip_of("unused");
    // (label, whether the run's zip is in the golden, the arguments)
    let runs: Vec<(&str, bool, Vec<&str>)> = vec![
        (
            "post-assembler-off",
            true,
            vec![
                "--off-ramp-type",
                "POST_ASSEMBLER_OFF",
                "--off-ramp-file",
                &post_assembler,
                "--use-jdk-deflater",
                "true",
            ],
        ),
        (
            "pre-assembler-off",
            true,
            vec![
                "--off-ramp-type",
                "PRE_ASSEMBLER_OFF",
                "--off-ramp-file",
                &pre_assembler,
                "-L",
                "chr1:200-600",
                "-L",
                "chr2",
            ],
        ),
        (
            "pre-filter-off",
            true,
            vec![
                "--off-ramp-type",
                "PRE_FILTER_OFF",
                "--off-ramp-file",
                &pre_filter,
                "-L",
                "chr1:1150-1350",
            ],
        ),
        (
            "no-ramp-file",
            false,
            vec!["--off-ramp-type", "PRE_FILTER_OFF"],
        ),
        (
            "post-assembler-on",
            false,
            vec![
                "--on-ramp-type",
                "POST_ASSEMBLER_ON",
                "--on-ramp-file",
                &post_assembler,
            ],
        ),
        (
            "post-filter-on",
            false,
            vec![
                "--on-ramp-type",
                "POST_FILTER_ON",
                "--on-ramp-file",
                &pre_filter,
                "-L",
                "chr1:1150-1350",
            ],
        ),
        (
            "post-filter-on-wide",
            false,
            vec![
                "--on-ramp-type",
                "POST_FILTER_ON",
                "--on-ramp-file",
                &pre_filter,
                "-L",
                "chr1:1000-1500",
            ],
        ),
        (
            "no-on-ramp-file",
            false,
            vec!["--on-ramp-type", "POST_FILTER_ON"],
        ),
        (
            "on-after-off",
            false,
            vec![
                "--on-ramp-type",
                "POST_ASSEMBLER_ON",
                "--on-ramp-file",
                &post_assembler,
                "--off-ramp-type",
                "PRE_FILTER_OFF",
                "--off-ramp-file",
                &unused,
            ],
        ),
    ];
    for (case, zipped, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            "RampedHaplotypeCaller",
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

        // A refusal: the reference's exception, before anything is read.
        let refusal = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{case}\t")));
        if let Some(refusal) = refusal {
            assert_ne!(outcome.status, 0, "{case} is refused");
            let (class, message) = refusal.split_once(':').expect("class:message");
            assert!(
                outcome.stderr.contains(message),
                "{case}: {class}: {}",
                outcome.stderr
            );
            continue;
        }
        assert_eq!(outcome.status, 0, "{case}: {}", outcome.stderr);
        let written = std::fs::read_to_string(&output).expect("the output");
        let body: String = written
            .lines()
            .skip_while(|line| !line.starts_with("#CHROM"))
            .filter(|line| !line.is_empty())
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
        if !zipped {
            continue;
        }

        let zip = dir.join(format!("{case}.zip"));
        let mut archive = zip::ZipArchive::new(std::fs::File::open(&zip).expect("the zip"))
            .expect("a zip archive");
        let names: Vec<String> = (0..archive.len())
            .map(|i| {
                archive
                    .by_index(i)
                    .expect("an entry")
                    .name()
                    .expect("a name")
                    .to_string()
            })
            .collect();
        assert_eq!(
            names.join(","),
            section(&text, "entry", case),
            "{case} entries"
        );
        for (i, name) in names.iter().enumerate() {
            let mut bytes = Vec::new();
            std::io::Read::read_to_end(&mut archive.by_index(i).expect("an entry"), &mut bytes)
                .expect("the entry's bytes");
            let label = format!("{case}/{name}");
            if name.ends_with(".bam") || name.ends_with(".bai") {
                assert_eq!(hex(&bytes), section(&text, "bytes", &label), "{label}");
            } else {
                assert_eq!(
                    String::from_utf8(bytes).expect("a text entry"),
                    section(&text, "content", &label),
                    "{label}"
                );
            }
        }
    }
}
