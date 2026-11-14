//! `DetermineGermlineContigPloidy` past the interpreter, against the frozen reference: the
//! fixture of `tools/python-tools/make_gcnv_fixture.py` (stdlib only, seeded), the command line of
//! `measure_gcnv.sh`, and every `ploidy/` file of the `gcnv-frozen` golden, byte for byte.
//!
//! The reference stops without `python` or `gcnvkernel` and the port keeps that refusal, so a
//! `python` that imports anything stands in for the environment. `PATH` is process-wide, so this
//! file holds one test.

use sha2::{Digest, Sha256};
use std::path::Path;

const GOLDEN: &str = "../gatk-tools/tests/data/gcnv_frozen.txt.gz";
const FIXTURE: &str = "../../tools/python-tools/make_gcnv_fixture.py";

fn executable(path: &Path, script: &str) {
    use std::os::unix::fs::PermissionsExt;
    std::fs::write(path, script).expect("a script");
    std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).expect("its mode");
}

/// `measure_gcnv.sh`'s command line.
fn ploidy(fixture: &Path) -> Result<(), String> {
    let f = |name: &str| fixture.join(name).display().to_string();
    let mut argv: Vec<String> = vec!["DetermineGermlineContigPloidy".into()];
    for i in 0..6 {
        argv.push("-I".into());
        argv.push(f(&format!("sample{i}.counts.tsv")));
    }
    for arg in [
        "-L",
        &f("intervals.interval_list"),
        "--interval-merging-rule",
        "OVERLAPPING_ONLY",
        "--contig-ploidy-priors",
        &f("priors.tsv"),
        "--output",
        "ploidy",
        "--output-prefix",
        "p",
        "--max-advi-iter-first-epoch",
        "200",
        "--max-advi-iter-subsequent-epochs",
        "50",
        "--num-thermal-advi-iters",
        "100",
        "--max-training-epochs",
        "3",
        "--min-training-epochs",
        "1",
        "--max-calling-iters",
        "3",
    ] {
        argv.push(arg.to_string());
    }
    let ran = gatk_cli::run(&argv);
    if ran.status == 0 {
        Ok(())
    } else {
        Err(ran.stderr)
    }
}

#[test]
fn past_the_interpreter() {
    let previous = std::env::current_dir().expect("a directory");
    let golden = previous.join(GOLDEN);
    if !golden.exists() {
        eprintln!("nothing to replay: the suite gcnv-frozen has no golden yet");
        return;
    }
    let scratch = std::env::temp_dir().join(format!("gcnv-ploidy-{}", std::process::id()));
    let bin = scratch.join("bin");
    let root = scratch.join("work");
    let fixture = scratch.join("fixture");
    for dir in [&bin, &root, &fixture] {
        std::fs::create_dir_all(dir).expect("a directory");
    }
    let made = std::process::Command::new("python3")
        .arg(previous.join(FIXTURE))
        .arg(&fixture)
        .output();
    if !made.as_ref().is_ok_and(|o| o.status.success()) {
        eprintln!("no python3 to make the fixture with: {made:?}");
        return;
    }
    executable(&bin.join("python"), "#!/bin/sh\nexit 0\n");

    let bytes = std::fs::read(&golden).expect("the golden");
    let mut text = String::new();
    {
        use std::io::Read;
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .expect("a gzip");
    }
    let wanted: Vec<(String, String, String)> = text
        .lines()
        .filter(|line| line.starts_with("file\tploidy/"))
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f[1].to_string(), f[2].to_string(), f[3].to_string())
        })
        .collect();
    assert_eq!(wanted.len(), 38, "the ploidy files of the golden");

    std::env::set_current_dir(&root).expect("the scratch directory");
    // SAFETY: this file holds one test, so nothing else in the process reads the environment.
    unsafe { std::env::set_var("PATH", &bin) };
    let ran = ploidy(&fixture);
    let mut failures = Vec::new();
    if ran.is_ok() {
        for (path, size, sha) in &wanted {
            let bytes = std::fs::read(root.join(path)).unwrap_or_default();
            let digest: String = Sha256::digest(&bytes)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            if bytes.len().to_string() != *size || digest != *sha {
                failures.push(path.clone());
            }
        }
    }

    // The reference's own refusal where the environment lacks the package.
    executable(
        &bin.join("python"),
        "#!/bin/sh\ncase \"$2\" in *gcnvkernel*) exit 1;; esac\nexit 0\n",
    );
    let refused = ploidy(&fixture);
    std::env::set_current_dir(&previous).expect("the directory");
    let _ = std::fs::remove_dir_all(&scratch);

    ran.expect("the run");
    assert!(failures.is_empty(), "differs: {failures:?}");
    let message = refused.expect_err("the missing package");
    assert!(
        message.contains("A required Python package (\"gcnvkernel\") could not be imported"),
        "{message}"
    );
}
