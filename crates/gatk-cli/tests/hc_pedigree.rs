//! Conformance for the callers' pedigree arguments against GATK 4.6.2.0, every run of the dump
//! replayed through the command line.
//!
//! Golden from `tools/readfilter-conformance/HcPedigreeDump.java`, over `HcDragenEndToEndDump`'s
//! reads with a third sample: --pedigree and --founder-id with PossibleDeNovo,
//! TransmittedSingleton, ExcessHet and InbreedingCoeff, the de novo thresholds, PED files with two
//! families, a comment and a missing parent, a malformed one, and the refusals without a pedigree
//! annotation. Each output is compared from its ##INFO and ##FORMAT lines to its last record.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_pedigree.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-pedigree");
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
    let mut peds = std::collections::HashMap::new();
    for label in ["trio", "mixed", "malformed", "unknown-sex"] {
        let path = dir.join(format!("{label}.ped"));
        std::fs::write(&path, section(&text, "ped", label)).expect("a PED file");
        peds.insert(label, path.to_string_lossy().to_string());
    }
    let (trio, mixed) = (peds["trio"].as_str(), peds["mixed"].as_str());
    let de_novo = ["-A", "PossibleDeNovo", "-A", "TransmittedSingleton"];
    let with = |front: &[&'static str], back: &[&'static str]| -> Vec<&'static str> {
        front.iter().chain(back).copied().collect()
    };
    let leak = |text: &str| -> &'static str { Box::leak(text.to_string().into_boxed_str()) };
    let (trio, mixed) = (leak(trio), leak(mixed));
    let (malformed, unknown_sex) = (leak(&peds["malformed"]), leak(&peds["unknown-sex"]));
    let runs: Vec<(&str, &str, Vec<&str>)> = vec![
        (
            "hc-trio",
            "HaplotypeCaller",
            with(&["--pedigree", trio], &de_novo),
        ),
        (
            "hc-trio-thresholds",
            "HaplotypeCaller",
            with(
                &[
                    "--pedigree",
                    trio,
                    "--denovo-parent-gq-threshold",
                    "99",
                    "--denovo-depth-threshold",
                    "40",
                ],
                &de_novo,
            ),
        ),
        (
            "hc-trio-default-annotations",
            "HaplotypeCaller",
            vec!["-ped", trio],
        ),
        (
            "hc-founder-id",
            "HaplotypeCaller",
            vec!["--founder-id", "s2", "--founder-id", "s3"],
        ),
        (
            "hc-founder-id-one",
            "HaplotypeCaller",
            with(&["--founder-id", "s1"], &de_novo),
        ),
        (
            "hc-founder-id-and-pedigree",
            "HaplotypeCaller",
            vec!["--founder-id", "s1", "--pedigree", mixed],
        ),
        (
            "hc-mixed",
            "HaplotypeCaller",
            with(&["--pedigree", mixed], &de_novo),
        ),
        (
            "hc-gvcf-trio",
            "HaplotypeCaller",
            with(
                &["--pedigree", trio, "-ERC", "GVCF", "--sample-name", "s1"],
                &de_novo,
            ),
        ),
        (
            "hc-malformed",
            "HaplotypeCaller",
            with(&["--pedigree", malformed], &de_novo),
        ),
        (
            "hc-unknown-sex",
            "HaplotypeCaller",
            with(&["--pedigree", unknown_sex], &de_novo),
        ),
        (
            "hc-no-pedigree-annotation",
            "HaplotypeCaller",
            vec!["--pedigree", trio, "--disable-tool-default-annotations"],
        ),
        ("mutect2-pedigree", "Mutect2", vec!["--pedigree", trio]),
        ("mutect2-founder-id", "Mutect2", vec!["--founder-id", "s2"]),
        (
            "mutect2-trio",
            "Mutect2",
            with(&["--pedigree", trio, "-A", "ExcessHet"], &de_novo),
        ),
    ];
    for (case, tool, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            tool,
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
            // `PedReader`'s source is an identity hash, masked in the golden.
            for part in message.split("<hash>") {
                assert!(outcome.stderr.contains(part), "{case}: {}", outcome.stderr);
            }
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
