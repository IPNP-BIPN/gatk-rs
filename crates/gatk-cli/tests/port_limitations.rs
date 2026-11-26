//! Conformance for what the port used to refuse as its own limitation, against GATK 4.7.0.0.
//!
//! Golden from `tools/argument-conformance/PortLimitationDump.java`: each case is one command line
//! the reference ran from a directory of its own, with the refusal it threw and every file it
//! left. The golden carries the command line, so this test replays exactly what the dump ran, and
//! a case added to the dump needs nothing here.

use gatk_corpus as corpus;
use gatk_tools::main_entry::{self, Failure};

const GOLDEN: &str = "tests/data/port_limitations.txt.gz";

/// The golden's rows of one kind, split on tabs after the kind, each field unescaped.
fn rows(text: &str, kind: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter_map(|line| line.strip_prefix(kind)?.strip_prefix('\t'))
        .map(|rest| rest.split('\t').map(unescape).collect())
        .collect()
}

fn unescape(field: &str) -> String {
    let mut out = String::with_capacity(field.len());
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

#[test]
fn every_case_is_the_reference_s_answer() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        eprintln!("skipping: {GOLDEN} is frozen from a CI run");
        return;
    }
    let text = corpus::read_golden(&golden);

    // The dump's working directory, its fixtures included. This is the only test in this binary,
    // so moving the process into it races nothing.
    let root =
        std::env::temp_dir().join(format!("gatk-cli-port-limitations-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).expect("a scratch directory");
    for row in rows(&text, "input") {
        let path = root.join(&row[0]);
        std::fs::create_dir_all(path.parent().expect("a directory")).expect("a scratch directory");
        std::fs::write(&path, corpus::decode_base64(&row[1])).expect("a fixture");
    }
    for row in rows(&text, "directory") {
        std::fs::create_dir_all(root.join(&row[0])).expect("a directory input");
    }
    std::env::set_current_dir(&root).expect("the scratch directory");

    // A message naming a file by its absolute path names the dump's directory; ours is this one.
    let theirs = rows(&text, "cwd").remove(0).remove(0);
    let ours = std::env::current_dir()
        .expect("the scratch directory")
        .display()
        .to_string();

    let files = rows(&text, "file");
    let dirs = rows(&text, "dir");
    let errors = rows(&text, "error");
    let kinds = rows(&text, "kind");
    for row in rows(&text, "argv") {
        let (case, argv) = row.split_first().expect("a case");
        let dir = format!("port-limitations/{case}");
        std::fs::create_dir_all(&dir).expect("the case's directory");
        let run = gatk_cli::run(argv);

        match errors.iter().find(|error| error[0] == *case) {
            Some(error) => {
                let (class, message) = error[1].split_once(':').expect("class:message");
                let message = message.replace(&theirs, &ours);
                let kind = &kinds.iter().find(|kind| kind[0] == *case).expect("a kind")[1];
                // `mainEntry`'s status for each handler, and what `handleNonUserException` prints.
                let (failure, printed) = match kind.as_str() {
                    "commandline" => (Failure::CommandLine, message),
                    "user" => (Failure::User, message),
                    _ => (Failure::Other, format!("{class}: {message}")),
                };
                assert_eq!(
                    run.status,
                    main_entry::exit_status(failure),
                    "{case}: {}",
                    run.stderr
                );
                assert!(run.stderr.contains(&printed), "{case}: {}", run.stderr);
            }
            None => assert_eq!(run.status, 0, "{case}: {}", run.stderr),
        }

        // Everything the run left, directories included, against what the reference left.
        let mut expected: Vec<String> = files
            .iter()
            .chain(dirs.iter())
            .filter(|entry| entry[0] == *case)
            .map(|entry| entry[1].clone())
            .collect();
        expected.sort();
        let mut produced = Vec::new();
        let mut pending = vec![std::path::PathBuf::from(&dir)];
        while let Some(at) = pending.pop() {
            for entry in std::fs::read_dir(&at).expect("a directory") {
                let path = entry.expect("an entry").path();
                if path.is_dir() {
                    pending.push(path.clone());
                }
                let relative = path.strip_prefix(&dir).expect("under the case");
                produced.push(relative.display().to_string());
            }
        }
        produced.sort();
        assert_eq!(produced, expected, "{case}: what the run left");
        for file in files.iter().filter(|file| file[0] == *case) {
            let bytes = std::fs::read(format!("{dir}/{}", file[1])).expect("a file");
            assert!(
                bytes == corpus::decode_base64(&file[2]),
                "{case}: {} differs from the reference's",
                file[1]
            );
        }
    }
}

/// AlleleFrequencyQC's refusal for fewer than two bins, which the runner writes out by hand.
#[test]
fn the_chi_squared_refusal_is_commons_math_s() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        return;
    }
    let probes = rows(&corpus::read_golden(&golden), "probe");
    let class = "org.apache.commons.math3.exception.NotStrictlyPositiveException";
    for (case, shape) in [("chi-squared0", "0"), ("chi-squared-1", "-0.5")] {
        let probe = probes
            .iter()
            .find(|probe| probe[0] == case)
            .expect("a probe");
        assert_eq!(probe[1], format!("{class}:shape ({shape})"), "{case}");
    }
}
