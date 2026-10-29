//! The eleven Spark tools this port runs, replayed end to end against GATK 4.6.2.0 with
//! `--spark-master local[1]`.
//!
//! Golden from `tools/readfilter-conformance/SparkWalkersDump.java`, which records each run's
//! command line, what it printed, every file it wrote and, where it failed, why.
//!
//! # What this suite is for
//!
//!  * **the counts printed rather than returned**, and `-O` receiving the same text;
//!  * **Disq's guessed first record**: with no `.sbi` the reads start where `BamRecordGuesser`
//!    first accepts ten chained records, and a mapped read with no cigar is never one, so the
//!    shared fixture's chr1 reads are not read at all and `CountReadsSpark` counts 3 of 11;
//!  * **PrintReadsSpark's BAM, `.bai` and `.sbi`**, and its sharded directory;
//!  * **BQSR on one partition**: `BaseRecalibratorSpark`'s table over the reads Disq's guess
//!    leaves, `ApplyBQSRSpark` applying it through `PrintReadsSpark`'s writer, and
//!    `BQSRPipelineSpark` chaining the two without writing the table;
//!  * **SortSamSpark**, that writer asked to sort to the header: a `.bai` for a coordinate sort
//!    and only the `.sbi` for a queryname one;
//!  * **PrintVariantsSpark's VCF**, as text with no index, or block-compressed with a `.tbi` whose
//!    last chunk ends where the terminator starts;
//!  * **and PileupSpark's text-file directory**, Hadoop's checksums included.
//!
//! # How the runs are replayed
//!
//! Every path in the golden is relative, so `PrintReadsSpark`'s `@PG` command line is the same
//! here as in the dump: the test changes into a scratch directory and replays each `argv` row
//! as it is. That is why this file holds one test, since the working directory is the process's.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/spark_walkers.txt.gz"),
    )
}

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

/// Every file under `dir`, by path relative to it, in the dump's order.
fn files_under(dir: &std::path::Path) -> Vec<(String, Vec<u8>)> {
    let mut found = Vec::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(next) = pending.pop() {
        let Ok(entries) = std::fs::read_dir(&next) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                pending.push(path);
            } else {
                let name = path
                    .strip_prefix(dir)
                    .expect("under the directory")
                    .to_string_lossy()
                    .to_string();
                found.push((name, std::fs::read(&path).expect("a written file")));
            }
        }
    }
    found.sort();
    found
}

#[test]
fn every_run_replays_as_the_reference_ran_it() {
    let text = golden();
    let root = std::env::temp_dir().join("gatk-cli-spark-walkers");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(root.join("sparkwalkers-dump")).expect("a scratch directory");
    std::env::set_current_dir(&root).expect("the scratch directory");
    let cwd = std::env::current_dir()
        .expect("a cwd")
        .display()
        .to_string();

    let mut labels = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        match fields[..] {
            ["input", name, encoded] => std::fs::write(
                format!("sparkwalkers-dump/{name}"),
                corpus::decode_base64(encoded),
            )
            .expect("a fixture"),
            ["argv", label, _] => labels.push(label.to_string()),
            _ => {}
        }
    }
    // The dump indexed this copy with IndexFeatureFile; the port's index is measured byte for
    // byte by `index-feature-file`, so it is made here the same way.
    std::fs::copy(
        "sparkwalkers-dump/plain.vcf",
        "sparkwalkers-dump/indexed.vcf",
    )
    .expect("a copy");
    let indexed = gatk_cli::run(&[
        "IndexFeatureFile".to_string(),
        "-I".to_string(),
        "sparkwalkers-dump/indexed.vcf".to_string(),
    ]);
    assert_eq!(indexed.status, 0, "{}", indexed.stderr);
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
        let out_dir = root.join("sparkwalkers-dump/out").join(label);
        if argv.iter().any(|arg| arg.contains("/out/")) {
            std::fs::create_dir_all(&out_dir).expect("an output directory");
        }
        let ran = gatk_cli::run(&argv);
        let stdout = ran.stdout.replace(&cwd, "<cwd>");
        let expected = row("stdout", label).expect("a stdout row");
        if stdout != expected {
            failures.push(format!(
                "{label}: stdout {stdout:?}, reference {expected:?}"
            ));
        }
        match row("error", label) {
            Some(error) => {
                let message = error.split_once(':').map_or("", |(_, message)| message);
                let stderr = ran.stderr.replace(&cwd, "<cwd>");
                if ran.status == 0 || !stderr.contains(message) {
                    failures.push(format!(
                        "{label}: status {}, reference refused with {error:?}; stderr {stderr:?}",
                        ran.status
                    ));
                }
            }
            None if ran.status != 0 => {
                failures.push(format!("{label}: status {}: {}", ran.status, ran.stderr));
            }
            None => {}
        }
        let prefix = format!("file\t{label}\t");
        let reference: Vec<(String, Vec<u8>)> = text
            .lines()
            .filter_map(|line| line.strip_prefix(&prefix))
            .map(|rest| {
                let (name, encoded) = rest.split_once('\t').expect("a file row");
                (name.to_string(), corpus::decode_base64(encoded))
            })
            .collect();
        let produced = files_under(&out_dir);
        let names = |files: &[(String, Vec<u8>)]| -> Vec<String> {
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
                failures.push(format!(
                    "{label}/{name}: {} bytes, reference {} bytes",
                    ours.len(),
                    theirs.len()
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
