//! Conformance for what the port used to refuse as its own limitation, against GATK 4.7.0.0.
//!
//! Golden from `tools/argument-conformance/PortLimitationDump.java`: each case is one command line
//! the reference ran from a directory of its own, with the refusal it threw and every file it
//! left. The golden carries the command line, so this test replays exactly what the dump ran, and
//! a case added to the dump needs nothing here.

use gatk_corpus as corpus;

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
    std::env::set_current_dir(&root).expect("the scratch directory");

    let files = rows(&text, "file");
    let errors = rows(&text, "error");
    for row in rows(&text, "argv") {
        let (case, argv) = row.split_first().expect("a case");
        let dir = format!("port-limitations/{case}");
        std::fs::create_dir_all(&dir).expect("the case's directory");
        let run = gatk_cli::run(argv);

        match errors.iter().find(|error| error[0] == *case) {
            Some(error) => {
                let message = error[1].split_once(':').expect("class:message").1;
                assert_ne!(run.status, 0, "{case}");
                assert!(run.stderr.contains(message), "{case}: {}", run.stderr);
            }
            None => assert_eq!(run.status, 0, "{case}: {}", run.stderr),
        }

        let expected: Vec<(&str, Vec<u8>)> = files
            .iter()
            .filter(|file| file[0] == *case)
            .map(|file| (file[1].as_str(), corpus::decode_base64(&file[2])))
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
            let ours = std::fs::read(format!("{dir}/{name}")).expect("a file");
            assert!(
                ours == *bytes,
                "{case}: {name} differs from the reference's"
            );
        }
    }
}
