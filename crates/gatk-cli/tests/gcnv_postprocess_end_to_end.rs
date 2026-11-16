//! `PostprocessGermlineCNVCalls` past the interpreter, against the frozen reference.
//!
//! The dump comes from `tools/python-tools/gcnv_dump.sh`, which runs the gCNV trio in
//! `broadinstitute/gatk:4.7.0.0` under `gcnv_frozen_env.sh`, then the postprocessing of every
//! sample, and prints the shards' bytes (`input` rows) and every output's sha256 (`file` rows). Its
//! golden is frozen from a CI artefact like every other; until then the replay reads
//! `GATK_RS_GCNV_DUMP`, or says it has nothing to compare (docs/gcnv-postprocess.md).
//!
//! The reference stops without `python` or `gcnvkernel` and the port keeps that refusal, so a
//! `python` that imports anything stands in for the environment. `PATH` is process-wide, so this
//! file holds one test.

use sha2::{Digest, Sha256};
use std::path::Path;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/gcnv_frozen.txt.gz";

fn from_hex(text: &str) -> Vec<u8> {
    text.as_bytes()
        .chunks(2)
        .map(|pair| {
            let digit = |c: u8| (c as char).to_digit(16).expect("a hex digit") as u8;
            digit(pair[0]) * 16 + digit(pair[1])
        })
        .collect()
}

fn executable(path: &Path, script: &str) {
    use std::os::unix::fs::PermissionsExt;
    std::fs::write(path, script).expect("a script");
    std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).expect("its mode");
}

/// `measure_postprocess.sh`'s command line for one sample.
fn postprocess(sample: usize) -> Result<(), String> {
    let index = sample.to_string();
    let out = |name: &str| format!("post-s{sample}/post.{name}");
    let argv: Vec<String> = [
        "PostprocessGermlineCNVCalls",
        "--model-shard-path",
        "gcnv/g-model",
        "--calls-shard-path",
        "gcnv/g-calls",
        "--contig-ploidy-calls",
        "ploidy/p-calls",
        "--sample-index",
        &index,
        "--allosomal-contig",
        "chrX",
        "--autosomal-ref-copy-number",
        "2",
        "--output-genotyped-intervals",
        &out("intervals.vcf"),
        "--output-genotyped-segments",
        &out("segments.vcf"),
        "--output-denoised-copy-ratios",
        &out("denoised.tsv"),
        "--create-output-variant-index",
        "false",
        "--add-output-vcf-command-line",
        "false",
    ]
    .iter()
    .map(|s| (*s).to_string())
    .collect();
    let ran = gatk_cli::run(&argv);
    if ran.status == 0 {
        Ok(())
    } else {
        Err(ran.stderr)
    }
}

fn replay(dump: &str, root: &Path, bin: &Path) {
    let mut wanted = Vec::new();
    for line in dump.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "input" => {
                let path = root.join(fields[1]);
                std::fs::create_dir_all(path.parent().expect("a parent")).expect("a directory");
                std::fs::write(path, from_hex(fields[2])).expect("a file");
            }
            "file" if fields[1].starts_with("post-s") => {
                wanted.push((
                    fields[1].to_string(),
                    fields[2].to_string(),
                    fields[3].to_string(),
                ));
            }
            _ => {}
        }
    }
    assert_eq!(wanted.len(), 18, "three files for each of six samples");
    let mut failures = Vec::new();
    for sample in 0..6 {
        std::fs::create_dir_all(root.join(format!("post-s{sample}"))).expect("a directory");
        postprocess(sample).expect("the postprocessing");
    }
    for (path, size, sha) in wanted {
        let bytes = std::fs::read(root.join(&path)).expect("a written file");
        let digest: String = Sha256::digest(&bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        if bytes.len().to_string() != size || digest != sha {
            failures.push(path);
        }
    }
    assert!(failures.is_empty(), "differs: {failures:?}");

    // The reference's own refusal where the environment lacks the package.
    executable(
        &bin.join("python"),
        "#!/bin/sh\ncase \"$2\" in *gcnvkernel*) exit 1;; esac\nexit 0\n",
    );
    let message = postprocess(0).expect_err("the missing package");
    assert!(
        message.contains("A required Python package (\"gcnvkernel\") could not be imported"),
        "{message}"
    );
}

#[test]
fn past_the_interpreter() {
    let scratch = std::env::temp_dir().join(format!("gcnv-postprocess-{}", std::process::id()));
    let bin = scratch.join("bin");
    let root = scratch.join("work");
    std::fs::create_dir_all(&bin).expect("a directory");
    std::fs::create_dir_all(&root).expect("a directory");
    executable(&bin.join("python"), "#!/bin/sh\nexit 0\n");
    let previous = std::env::current_dir().expect("a directory");
    let golden = previous.join(GOLDEN);
    let dump = std::env::var("GATK_RS_GCNV_DUMP")
        .ok()
        .map(std::path::PathBuf::from)
        .or_else(|| golden.exists().then_some(golden));
    let Some(path) = dump else {
        eprintln!("nothing to replay: the suite gcnv-frozen has no golden yet");
        return;
    };
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
    std::env::set_current_dir(&root).expect("the scratch directory");
    // SAFETY: this file holds one test, so nothing else in the process reads the environment.
    unsafe { std::env::set_var("PATH", &bin) };
    replay(&text, &root, &bin);
    std::env::set_current_dir(&previous).expect("the directory");
    let _ = std::fs::remove_dir_all(&scratch);
}
