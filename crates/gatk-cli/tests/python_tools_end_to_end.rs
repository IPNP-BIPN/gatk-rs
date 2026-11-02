//! `TrainVariantAnnotationsModel` and `ScoreVariantAnnotations` past the interpreter.
//!
//! The dump comes from `tools/python-tools/vets_dump.sh`, which runs in
//! `broadinstitute/gatk:4.6.2.0`, because the pinned image has no python. Its golden is frozen from
//! a CI artefact like every other; until then the replay reads `GATK_RS_VETS_DUMP`, or says it has
//! nothing to compare. The rest of the file is what the tools refuse, which needs no golden.
//!
//! The reference stops without `python` or its packages and the port keeps those refusals, so a
//! `python` that imports anything stands in for the environment. `PATH` is process-wide, so this
//! file holds one test.

use sha2::{Digest, Sha256};
use std::path::Path;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/vets_isolation_forest.txt.gz";

fn from_hex(text: &str) -> Vec<u8> {
    text.as_bytes()
        .chunks(2)
        .map(|pair| {
            let digit = |c: u8| (c as char).to_digit(16).expect("a hex digit") as u8;
            digit(pair[0]) * 16 + digit(pair[1])
        })
        .collect()
}

fn to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// A run, or what it printed on stderr when it failed.
fn run(argv: &[&str]) -> Result<(), String> {
    let argv: Vec<String> = argv.iter().map(|s| (*s).to_string()).collect();
    let ran = gatk_cli::run(&argv);
    if ran.status == 0 {
        Ok(())
    } else {
        Err(ran.stderr)
    }
}

/// `dump_vets.py`'s scorer digest, which is blind to the script's name and to the frames.
fn scorer_digest(data: &[u8]) -> String {
    let mut out = data[..2].to_vec();
    let mut at = 2;
    while at < data.len() {
        if data[at] == 0x95 {
            let length =
                u64::from_le_bytes(data[at + 1..at + 9].try_into().expect("eight")) as usize;
            out.extend_from_slice(&data[at + 9..at + 9 + length]);
            at += 9 + length;
        } else {
            out.extend_from_slice(&data[at..]);
            break;
        }
    }
    let needle = b"/tmp/isolation-forest.";
    let start = out
        .windows(needle.len())
        .position(|window| window == needle)
        .expect("a script name");
    let end = start + needle.len() + 23;
    let mut normal = out[..start - 1].to_vec();
    normal.extend_from_slice(b"?NAME");
    normal.extend_from_slice(&out[end..]);
    to_hex(&Sha256::digest(&normal))
}

/// Every pass of a dump, in its directory.
fn replay(dump: &str, root: &Path) {
    let mut failures = Vec::new();
    let mut wanted: Vec<(String, String, String, String)> = Vec::new();
    for line in dump.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "input" => std::fs::write(root.join(fields[1]), from_hex(fields[2])).expect("a file"),
            kind => wanted.push((
                kind.to_string(),
                fields[1].to_string(),
                fields[2].to_string(),
                fields[3].to_string(),
            )),
        }
    }
    let seven = root.join("seven.json");
    std::fs::write(&seven, r#"{"random_state": 1, "n_estimators": 7, "n_jobs": 2}"#).expect("a file");
    for tag in ["default", "seven"] {
        let dir = root.join(tag);
        std::fs::create_dir_all(&dir).expect("a directory");
        let annotations = format!("{tag}.extract.annot.hdf5");
        let prefix = format!("{tag}/train");
        let mut train = vec![
            "TrainVariantAnnotationsModel",
            "--annotations-hdf5",
            &annotations,
            "-O",
            &prefix,
        ];
        if tag == "seven" {
            train.extend(["--hyperparameters-json", "seven.json"]);
        }
        run(&train).expect("the training");
        let out = format!("{tag}/score");
        let mut score = vec![
            "ScoreVariantAnnotations",
            "-V",
            "calls.vcf",
            "-A",
            "QD",
            "-A",
            "MQ",
            "-A",
            "FS",
            "-A",
            "SOR",
            "--resource:train,training=true,calibration=true",
            "resource.vcf",
            "--model-prefix",
            &prefix,
            "-O",
            &out,
            "--add-output-vcf-command-line",
            "false",
            "--create-output-variant-index",
            "false",
            "--do-not-gzip-vcf-output",
        ];
        if tag == "default" {
            score.extend([
                "--snp-calibration-sensitivity-threshold",
                "0.9",
                "--indel-calibration-sensitivity-threshold",
                "0.95",
            ]);
        }
        run(&score).expect("the scoring");
    }
    for (kind, tag, name, value) in wanted {
        let path = root.join(&tag).join(&name);
        let bytes = std::fs::read(&path).expect("a written file");
        let same = match kind.as_str() {
            "file" => bytes == from_hex(&value),
            "scorer" => scorer_digest(&bytes) == value,
            "scores" => {
                let reader =
                    gatk_tools::hdf5_reader::Hdf5Reader::open(bytes).expect("an HDF5 file");
                let (_, scores) = reader.doubles("/data/scores").expect("scores");
                let scores: Vec<u8> = scores.iter().flat_map(|x| x.to_le_bytes()).collect();
                to_hex(&scores) == value
            }
            other => panic!("a row of kind {other}"),
        };
        if !same {
            failures.push(format!("{kind} {tag}/{name}"));
        }
    }
    assert!(failures.is_empty(), "differs: {failures:?}");
}

fn executable(path: &Path, script: &str) {
    use std::os::unix::fs::PermissionsExt;
    std::fs::write(path, script).expect("a script");
    std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).expect("its mode");
}

/// Two columns and a dozen rows the port's own `ExtractVariantAnnotations` writes.
fn annotations_from_a_small_vcf(root: &Path) {
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=QD,Number=1,Type=Float,Description=\"QD\">\n\
         ##INFO=<ID=MQ,Number=1,Type=Float,Description=\"MQ\">\n\
         ##contig=<ID=chr1,length=100000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..40 {
        vcf.push_str(&format!(
            "chr1\t{}\t.\tA\tC\t50\t.\tQD={}.5;MQ={}.0\n",
            100 + i * 10,
            i % 7,
            40 + i % 5
        ));
    }
    std::fs::write(root.join("small.vcf"), &vcf).expect("a file");
    run(&["IndexFeatureFile", "-I", "small.vcf"]).expect("an index");
    run(&[
        "ExtractVariantAnnotations",
        "-V",
        "small.vcf",
        "-A",
        "QD",
        "-A",
        "MQ",
        "--mode",
        "SNP",
        "--resource:train,training=true,calibration=true",
        "small.vcf",
        "-O",
        "small",
        "--add-output-vcf-command-line",
        "false",
        "--create-output-variant-index",
        "false",
        "--do-not-gzip-vcf-output",
    ])
    .expect("an extraction");
}

fn refusals(root: &Path, bin: &Path) {
    let train = |extra: &[&str]| {
        let mut argv = vec![
            "TrainVariantAnnotationsModel",
            "--annotations-hdf5",
            "small.annot.hdf5",
            "--mode",
            "SNP",
            "-O",
            "refused",
        ];
        argv.extend_from_slice(extra);
        run(&argv)
    };
    let limit = |result: Result<(), String>, about: &str| {
        let message = result.expect_err("a refusal");
        assert!(message.contains(about), "{message}");
        assert!(
            message.contains("This message is the port's own and not GATK's."),
            "{message}"
        );
    };
    std::fs::write(root.join("unseeded.json"), "{}").expect("a file");
    limit(
        train(&["--hyperparameters-json", "unseeded.json"]),
        "without a random_state",
    );
    std::fs::write(root.join("sampled.json"), r#"{"random_state":0,"max_samples":10}"#)
        .expect("a file");
    limit(
        train(&["--hyperparameters-json", "sampled.json"]),
        "max_samples",
    );
    limit(
        train(&["--unlabeled-annotations-hdf5", "small.annot.hdf5"]),
        "--unlabeled-annotations-hdf5",
    );
    assert!(
        !root.join("refused.snp.scorer.pkl").exists(),
        "a refusal leaves no model"
    );

    // A model that was not the script's is not read as if it were.
    train(&[]).expect("a model");
    std::fs::write(root.join("refused.snp.scorer.pkl"), b"not a pickle").expect("a file");
    let refused = run(&[
        "ScoreVariantAnnotations",
        "-V",
        "small.vcf",
        "-A",
        "QD",
        "-A",
        "MQ",
        "--mode",
        "SNP",
        "--model-prefix",
        "refused",
        "-O",
        "refused-score",
        "--add-output-vcf-command-line",
        "false",
        "--create-output-variant-index",
        "false",
        "--do-not-gzip-vcf-output",
    ]);
    limit(refused, "refused.snp.scorer.pkl");

    // The reference's own refusal where the environment lacks a package.
    executable(&bin.join("python"), "#!/bin/sh\ncase \"$2\" in *sklearn*) exit 1;; esac\nexit 0\n");
    let message = train(&[]).expect_err("the missing package");
    assert!(
        message.contains("A required Python package (\"sklearn\") could not be imported"),
        "{message}"
    );
}

#[test]
fn past_the_interpreter() {
    let scratch = std::env::temp_dir().join(format!("python-tools-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).expect("a directory");
    let bin = scratch.join("bin");
    std::fs::create_dir_all(&bin).expect("a directory");
    executable(&bin.join("python"), "#!/bin/sh\nexit 0\n");
    let root = scratch.join("work");
    std::fs::create_dir_all(&root).expect("a directory");
    let previous = std::env::current_dir().expect("a directory");
    std::env::set_current_dir(&root).expect("the scratch directory");
    // SAFETY: this file holds one test, so nothing else in the process reads the environment.
    unsafe { std::env::set_var("PATH", &bin) };

    let golden = previous.join(GOLDEN);
    let dump = std::env::var("GATK_RS_VETS_DUMP")
        .ok()
        .map(std::path::PathBuf::from)
        .or_else(|| golden.exists().then_some(golden));
    match dump {
        Some(path) => {
            let bytes = std::fs::read(&path).expect("the dump");
            let text = if bytes.starts_with(&[0x1f, 0x8b]) {
                use std::io::Read;
                let mut text = String::new();
                flate2::read::MultiGzDecoder::new(&bytes[..])
                    .read_to_string(&mut text)
                    .expect("a gzip");
                text
            } else {
                String::from_utf8(bytes).expect("text")
            };
            replay(&text, &root);
        }
        None => eprintln!("nothing to replay: the suite vets-isolation-forest has no golden yet"),
    }

    annotations_from_a_small_vcf(&root);
    refusals(&root, &bin);

    std::env::set_current_dir(&previous).expect("the directory");
    let _ = std::fs::remove_dir_all(&scratch);
}
