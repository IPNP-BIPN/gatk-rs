//! The eight tools that stop at an interpreter the pinned image does not have, replayed against
//! GATK 4.6.2.0 as far as it goes there.
//!
//! Golden from `tools/readfilter-conformance/InterpreterToolsDump.java`, which records each run's
//! command line, what it printed, why it failed and everything it left under its output
//! directory.
//!
//! # What this suite is for
//!
//!  * **PlotDenoisedCopyRatios and PlotModeledSegments**, every input read and checked before
//!    `RScriptExecutor` looks for its executable, the output directory made by then;
//!  * **DetermineGermlineContigPloidy, GermlineCNVCaller and NVScoreVariants**, whose `onStartup`
//!    imports a package, so that only the parser refuses ahead of the missing python;
//!  * **PostprocessGermlineCNVCalls**, where `GATKTool.onStartup` comes before that import, and
//!    `-L` reads the shards' interval lists;
//!  * **TrainVariantAnnotationsModel**, the matrices each backend reads and refuses before it
//!    needs python;
//!  * **and ScoreVariantAnnotations**, a walker: the startup, the header, the backend, and under
//!    `PYTHON_SCRIPT` a whole pass and the annotations file, or a VCF when nothing was extracted.
//!
//! # How the runs are replayed
//!
//! Every path in the golden is relative, so the test changes into a scratch directory and replays
//! each `argv` row as it is. The HDF5 files the rows read are the port's own: the `setup-` rows
//! make them through `ExtractVariantAnnotations`, and their contents are compared as text since
//! two runs of the library write different bytes. `PATH` is emptied, because a machine with a
//! python is not the machine the reference was measured on; the last check puts one on it.

use gatk_corpus as corpus;
use gatk_tools::hdf5_reader::Hdf5Reader;
use gatk_tools::hdf5_writer::Hdf5File;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/interpreter_tools.txt.gz";

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// Every path under `dir` relative to it, directories marked, in the dump's order.
fn under(dir: &std::path::Path) -> Vec<(String, Option<Vec<u8>>)> {
    let mut found = Vec::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(next) = pending.pop() {
        let Ok(entries) = std::fs::read_dir(&next) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            let name = path
                .strip_prefix(dir)
                .expect("under the directory")
                .to_string_lossy()
                .to_string();
            if path.is_dir() {
                found.push((name, None));
                pending.push(path);
            } else {
                found.push((name, Some(std::fs::read(&path).expect("a written file"))));
            }
        }
    }
    found.sort();
    found
}

/// `Double.toString`, as `Arrays.toString` prints an element.
fn double(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

/// The dump's `annotations`: names, the matrix, and three labels, a missing label being a dash.
fn annotations(bytes: &[u8]) -> String {
    let reader = Hdf5Reader::open(bytes.to_vec()).expect("an HDF5 file");
    let names = reader.strings("/annotations/names").expect("the names");
    let columns = names.len();
    let (_, values) = reader
        .doubles("/annotations/chunk_0")
        .expect("the first chunk");
    let rows: Vec<String> = values
        .chunks(columns.max(1))
        .map(|row| {
            format!(
                "[{}]",
                row.iter()
                    .map(|v| double(*v))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        })
        .collect();
    let mut text = format!(
        "names=[{}];annotations=[{}]",
        names.join(", "),
        rows.join(", ")
    );
    for label in ["snp", "training", "calibration"] {
        let values = match reader.doubles(&format!("/labels/{label}")) {
            Some((_, values)) => format!(
                "[{}]",
                values
                    .iter()
                    .map(|v| (*v == 1.0).to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            None => "-".to_string(),
        };
        text.push_str(&format!(";{label}={values}"));
    }
    text
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn scores_file(values: &[f64]) -> Vec<u8> {
    let mut file = Hdf5File::new();
    file.double_array("/data/scores", values);
    file.to_bytes()
}

/// The first failing message of a golden row, or none.
fn refusal_message(error: &str) -> &str {
    // The dump prints a null message as `null`; `Throwable.toString` prints none.
    match error.split_once(':') {
        Some((_, "null")) | None => "",
        Some((_, message)) => message,
    }
}

fn replay(text: &str, root: &std::path::Path) {
    std::fs::create_dir_all(root.join("interpreter-dump")).expect("a scratch directory");
    std::env::set_current_dir(root).expect("the scratch directory");
    let cwd = std::env::current_dir()
        .expect("a cwd")
        .display()
        .to_string();

    let mut labels = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        match fields[..] {
            ["folder", name] => {
                std::fs::create_dir_all(format!("interpreter-dump/{name}")).expect("a folder");
            }
            ["input", name, encoded] => std::fs::write(
                format!("interpreter-dump/{name}"),
                corpus::decode_base64(encoded),
            )
            .expect("a fixture"),
            ["argv", label, _] => labels.push(label.to_string()),
            _ => {}
        }
    }
    for name in [
        "calls.vcf",
        "resource.vcf",
        "resource-snps.vcf",
        "resource-every-snp.vcf",
    ] {
        let indexed = gatk_cli::run(&[
            "IndexFeatureFile".to_string(),
            "-I".to_string(),
            format!("interpreter-dump/{name}"),
        ]);
        assert_eq!(indexed.status, 0, "{}", indexed.stderr);
    }
    std::fs::write(
        "interpreter-dump/calibrated.snp.calibrationScores.hdf5",
        scores_file(&[1.0, 2.0, 3.0]),
    )
    .expect("the scores");
    std::fs::write(
        "interpreter-dump/infinite.snp.calibrationScores.hdf5",
        scores_file(&[1.0, f64::INFINITY]),
    )
    .expect("the scores");
    assert!(!labels.is_empty(), "the golden holds no runs");

    let row = |kind: &str, label: &str| -> Option<String> {
        let prefix = format!("{kind}\t{label}\t");
        text.lines()
            .find(|line| line.starts_with(&prefix))
            .map(|line| unescape(&line[prefix.len()..]))
    };
    let mut failures = Vec::new();
    for label in &labels {
        let argv: Vec<String> = row("argv", label)
            .expect("a command line")
            .split('\t')
            .map(str::to_string)
            .collect();
        let out_dir = root.join("interpreter-dump/out").join(label);
        std::fs::create_dir_all(&out_dir).expect("an output directory");
        let ran = gatk_cli::run(&argv);
        let stdout = ran.stdout.replace(&cwd, "<cwd>");
        let expected = row("stdout", label).expect("a stdout row");
        if stdout != expected {
            failures.push(format!(
                "{label}: stdout {stdout:?}, reference {expected:?}"
            ));
        }
        let stderr = ran.stderr.replace(&cwd, "<cwd>");
        match row("error", label) {
            Some(error) => {
                let class = error
                    .split_once(':')
                    .map_or(error.as_str(), |(class, _)| class);
                let status = if class.starts_with("org.broadinstitute.barclay") {
                    1
                } else if class.contains("hellbender.exceptions.UserException") {
                    2
                } else {
                    3
                };
                if ran.status != status
                    || !stderr.contains(refusal_message(&error))
                    || (status == 3 && !stderr.contains(class))
                {
                    failures.push(format!(
                        "{label}: status {}, reference refused with {error:?}; stderr {stderr:?}",
                        ran.status
                    ));
                }
            }
            None if ran.status != 0 => {
                failures.push(format!("{label}: status {}: {stderr}", ran.status));
            }
            None => {}
        }

        // What the run left behind, refused or not.
        let mut reference: Vec<(String, String)> = Vec::new();
        for (kind, width) in [("dir", 3), ("index", 3), ("file", 4), ("hdf5", 4)] {
            let prefix = format!("{kind}\t{label}\t");
            for line in text.lines().filter_map(|line| line.strip_prefix(&prefix)) {
                let (name, value) = if width == 3 {
                    (line, String::new())
                } else {
                    let (name, value) = line.split_once('\t').expect("a file row");
                    (name, value.to_string())
                };
                let value = if kind == "file" {
                    hex(&corpus::decode_base64(&value))
                } else {
                    value
                };
                reference.push((name.to_string(), format!("{kind}:{value}")));
            }
        }
        reference.sort();
        let produced: Vec<(String, String)> = under(&out_dir)
            .into_iter()
            .map(|(name, bytes)| {
                let value = match bytes {
                    None => "dir:".to_string(),
                    // A Tribble index holds the absolute path and a modification time.
                    Some(_) if name.ends_with(".idx") => "index:".to_string(),
                    Some(bytes) if name.ends_with(".hdf5") => {
                        format!("hdf5:{}", annotations(&bytes))
                    }
                    Some(bytes) => format!("file:{}", hex(&bytes)),
                };
                (name, value)
            })
            .collect();
        let names = |files: &[(String, String)]| -> Vec<String> {
            files.iter().map(|(name, _)| name.clone()).collect()
        };
        if names(&produced) != names(&reference) {
            failures.push(format!(
                "{label}: files {:?}, reference {:?}",
                names(&produced),
                names(&reference)
            ));
            continue;
        }
        for ((name, ours), (_, theirs)) in produced.iter().zip(&reference) {
            if ours != theirs {
                failures.push(format!("{label}/{name}: {ours}, reference {theirs}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// With the interpreter on the path the reference would run it, which the port does not; the
/// message says so and is the port's own. The refusals ahead of it stay the reference's.
fn with_interpreters(root: &std::path::Path) {
    let bin = root.join("bin");
    std::fs::create_dir_all(&bin).expect("a directory");
    std::fs::write(bin.join("python"), "").expect("a python");
    std::fs::write(bin.join("Rscript"), "").expect("an Rscript");
    // SAFETY: this file holds one test, so nothing else in the process reads the environment.
    unsafe { std::env::set_var("PATH", &bin) };
    let to = |list: &[&str]| -> Vec<String> { list.iter().map(|s| (*s).to_string()).collect() };
    for argv in [
        to(&[
            "PlotDenoisedCopyRatios",
            "--standardized-copy-ratios",
            "interpreter-dump/standardized.tsv",
            "--denoised-copy-ratios",
            "interpreter-dump/denoised.tsv",
            "--sequence-dictionary",
            "interpreter-dump/ref.dict",
            "--output-prefix",
            "sample",
            "--output",
            "interpreter-dump/out/with-r",
        ]),
        to(&[
            "NVScoreVariants",
            "-V",
            "nothing.vcf",
            "-R",
            "nothing.fasta",
            "-O",
            "out.vcf",
        ]),
    ] {
        let ran = gatk_cli::run(&argv);
        assert_ne!(ran.status, 0);
        assert!(
            ran.stderr
                .contains("This message is the port's own and not GATK's."),
            "{}",
            ran.stderr
        );
    }
}

#[test]
fn every_run_replays_as_the_reference_ran_it() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
        return;
    }
    let text = corpus::read_golden(&golden);
    let root = std::env::temp_dir().join("gatk-cli-interpreter-tools");
    let _ = std::fs::remove_dir_all(&root);
    // SAFETY: this file holds one test, so nothing else in the process reads the environment.
    unsafe { std::env::set_var("PATH", "") };
    replay(&text, &root);
    with_interpreters(&root);
}
