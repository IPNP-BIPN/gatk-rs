//! Conformance for what a read walker writes, chosen by its output's name, against GATK 4.7.0.0.
//!
//! `GATKTool.createSAMWriter` hands the output's name to htsjdk 5.0.0's
//! `SAMFileWriterFactory.makeWriter`: a `.cram` is htsjdk's CRAM 3.1 with `CRAMBAIIndexer`'s
//! `<output>.bai`, a `.sam` is text with no index, and any other name is a BAM whose index strips
//! only a trailing `.bam`. The golden holds every file each PrintReads run left, and the port's
//! files must be the same bytes, the same names, and no more of them.
//!
//! The runs were made with RELATIVE paths, which the `@PG` command line records, so this test
//! makes them from a directory of its own laid out as the dump's was. A CRAM's file definition
//! also carries the output's absolute URI cut to twenty bytes (bytes 6 to 25), which only the
//! directory decides; those twenty bytes are taken from the golden before the comparison, and
//! nothing else is.

use gatk_corpus as corpus;

const GOLDEN: &str = "tests/data/cram_output.txt.gz";

/// The cases the dump ran, with the arguments it gave each, the `-O` last.
const CASES: &[(&str, &str, &[&str])] = &[
    (
        "cram",
        "out.cram",
        &["-I", "fixtures/reads.bam", "-R", "fixtures/reference.fasta"],
    ),
    (
        "cram-3.0",
        "out.cram",
        &[
            "-I",
            "fixtures/reads.bam",
            "-R",
            "fixtures/reference.fasta",
            "--output-cram-version",
            "3.0",
        ],
    ),
    (
        "cram-3.1",
        "out.cram",
        &[
            "-I",
            "fixtures/reads.bam",
            "-R",
            "fixtures/reference.fasta",
            "--output-cram-version",
            "3.1",
        ],
    ),
    (
        "cram-md5",
        "out.cram",
        &[
            "-I",
            "fixtures/reads.bam",
            "-R",
            "fixtures/reference.fasta",
            "--create-output-bam-md5",
        ],
    ),
    (
        "cram-no-index",
        "out.cram",
        &[
            "-I",
            "fixtures/reads.bam",
            "-R",
            "fixtures/reference.fasta",
            "--create-output-bam-index",
            "false",
        ],
    ),
    (
        "cram-interval",
        "out.cram",
        &[
            "-I",
            "fixtures/reads.bam",
            "-R",
            "fixtures/reference.fasta",
            "-L",
            "chr1:100-400",
        ],
    ),
    (
        "cram-no-reference",
        "out.cram",
        &["-I", "fixtures/reads.bam"],
    ),
    ("sam", "out.sam", &["-I", "fixtures/reads.bam"]),
    (
        "sam-md5",
        "out.sam",
        &["-I", "fixtures/reads.bam", "--create-output-bam-md5"],
    ),
    ("other-extension", "out.tmp", &["-I", "fixtures/reads.bam"]),
    (
        "bad-cram-version",
        "out.bam",
        &["-I", "fixtures/reads.bam", "--output-cram-version", "3.2"],
    ),
    (
        "cram-queryname",
        "out.cram",
        &[
            "-I",
            "inputs/queryname.bam",
            "-R",
            "fixtures/reference.fasta",
        ],
    ),
];

/// The golden's rows of one kind, split on tabs after the kind.
fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| line.strip_prefix(kind)?.strip_prefix('\t'))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

#[test]
fn every_file_is_the_one_the_reference_wrote() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        eprintln!("skipping: {GOLDEN} is frozen from a CI run");
        return;
    }
    let text = corpus::read_golden(&golden);

    // The dump's working directory: the corpus's fixtures, its own queryname input, and one
    // directory per case. This is
    // the only test in this binary, so moving the process into it races nothing.
    let root = std::env::temp_dir().join(format!("gatk-cli-cram-output-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    for row in rows(&text, "input") {
        let path = root.join(row[0]);
        std::fs::create_dir_all(path.parent().expect("a directory")).expect("a scratch directory");
        std::fs::write(&path, corpus::decode_base64(row[1])).expect("a fixture");
    }
    std::env::set_current_dir(&root).expect("the scratch directory");

    let files = rows(&text, "file");
    let errors = rows(&text, "error");
    for (case, output, arguments) in CASES {
        let dir = format!("cram-output/{case}");
        std::fs::create_dir_all(&dir).expect("the case's directory");
        let mut argv = vec!["PrintReads".to_string()];
        argv.extend(arguments.iter().map(|a| a.to_string()));
        argv.extend(["-O".to_string(), format!("{dir}/{output}")]);
        let run = gatk_cli::run(&argv);

        // A refusal is compared on its message; the class is what decides the status.
        match errors.iter().find(|row| row[0] == *case) {
            Some(row) => {
                let message = row[1].split_once(':').expect("class:message").1;
                assert_ne!(run.status, 0, "{case}");
                assert!(run.stderr.contains(message), "{case}: {}", run.stderr);
            }
            None => assert_eq!(run.status, 0, "{case}: {}", run.stderr),
        }

        let expected: Vec<(&str, Vec<u8>)> = files
            .iter()
            .filter(|row| row[0] == *case)
            .map(|row| (row[1], corpus::decode_base64(row[2])))
            .collect();
        let mut produced: Vec<String> = std::fs::read_dir(&dir)
            .expect("the case's directory")
            .map(|entry| {
                entry
                    .expect("an entry")
                    .file_name()
                    .to_string_lossy()
                    .into_owned()
            })
            .collect();
        produced.sort();
        let names: Vec<&str> = expected.iter().map(|(name, _)| *name).collect();
        assert_eq!(produced, names, "{case}: the files left");

        for (name, bytes) in &expected {
            let mut ours = std::fs::read(format!("{dir}/{name}")).expect("a file");
            // A CRAM's digest covers those twenty bytes too: each side's is checked against its
            // own CRAM rather than against the other's.
            if let Some(cram) = name.strip_suffix(".md5").filter(|c| c.ends_with(".cram")) {
                let digest = |bytes: &[u8]| gatk_tools::gather_bam_files::md5_file(bytes);
                let theirs = &expected
                    .iter()
                    .find(|(n, _)| *n == cram)
                    .expect("the CRAM")
                    .1;
                assert_eq!(
                    bytes,
                    digest(theirs).as_bytes(),
                    "{case}: the golden's digest"
                );
                let our_cram = std::fs::read(format!("{dir}/{cram}")).expect("the CRAM");
                assert_eq!(ours, digest(&our_cram).as_bytes(), "{case}: {name}");
                continue;
            }
            if name.ends_with(".cram") && ours.len() >= 26 && bytes.len() >= 26 {
                ours[6..26].copy_from_slice(&bytes[6..26]);
            }
            assert!(
                ours == *bytes,
                "{case}: {name} differs from the reference's"
            );
        }
    }
}
