//! `GenotypeGVCFs`, `GnarlyGenotyper` and `SelectVariants` over a `gendb://` workspace, replayed
//! against GATK 4.7.0.0.
//!
//! Golden from `tools/readfilter-conformance/GenomicsDBQueryDump.java`, which imports the shared
//! corpus's `cg_a`, `cg_b` and `cg_c` into a workspace and records each caller's output VCF or why
//! it failed. The replay imports the same inputs with the port's own `GenomicsDBImport` and runs
//! the same command lines. The import's `Date`, which the workspace's header carries into every
//! output, is masked on both sides, and so is the run's directory.
//!
//! All of it is skipped until a CI run has produced the golden. The refusals that are the port's
//! own are checked beside the replay, since the reference has nothing to say about them.

use gatk_corpus as corpus;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/genomicsdb_query.txt.gz";

const PORT_OWN: &str = "This message is the port's own and not GATK's.";

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

/// The dump's `masked`: the run's directory and every `Date="..."`.
fn masked(text: &str, cwd: &str) -> String {
    let text = text.replace(&format!("{cwd}/"), "<cwd>/");
    let mut out = String::with_capacity(text.len());
    let mut rest = text.as_str();
    while let Some(at) = rest.find("Date=\"") {
        let (before, after) = rest.split_at(at + "Date=\"".len());
        out.push_str(before);
        out.push_str("<masked>");
        let end = after.find('"').expect("a closed Date");
        rest = &after[end..];
    }
    out.push_str(rest);
    out
}

fn replay(text: &str, root: &std::path::Path) {
    std::fs::create_dir_all(root.join("genomicsdb-query-dump/out")).expect("a scratch directory");
    std::env::set_current_dir(root).expect("the scratch directory");
    let cwd = std::env::current_dir()
        .expect("a cwd")
        .display()
        .to_string();

    let mut labels = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        match fields[..] {
            ["input", name, encoded] => std::fs::write(
                format!("genomicsdb-query-dump/{name}"),
                corpus::decode_base64(encoded),
            )
            .expect("a fixture"),
            ["argv", label, _] => labels.push(label.to_string()),
            _ => {}
        }
    }
    assert!(!labels.is_empty(), "the golden holds no runs");
    for name in ["cg_a", "cg_b", "cg_c"] {
        // The importer reads the whole input and asks only that its index be there.
        std::fs::write(format!("genomicsdb-query-dump/{name}.g.vcf.idx"), b"").expect("an index");
    }

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
        let ran = gatk_cli::run(&argv);
        let stderr = masked(&ran.stderr, &cwd);
        if let Some(error) = row("error", label) {
            let (class, message) = error.split_once(':').unwrap_or((error.as_str(), ""));
            let status = if class.starts_with("org.broadinstitute.barclay") {
                1
            } else if class.contains("hellbender.exceptions.UserException") {
                2
            } else {
                3
            };
            if ran.status != status
                || (message != "null" && !stderr.contains(message))
                || (status == 3 && !stderr.contains(class))
            {
                failures.push(format!(
                    "{label}: status {}, reference refused with {error:?}; stderr {stderr:?}",
                    ran.status
                ));
            }
            continue;
        }
        if ran.status != 0 {
            failures.push(format!("{label}: status {}: {stderr}", ran.status));
            continue;
        }
        let Some(expected) = row("out", label) else {
            continue;
        };
        let output = format!("genomicsdb-query-dump/out/{label}.vcf");
        let produced = masked(&std::fs::read_to_string(&output).unwrap_or_default(), &cwd);
        if produced != expected {
            let at = produced
                .lines()
                .zip(expected.lines())
                .position(|(ours, theirs)| ours != theirs)
                .unwrap_or(produced.lines().count().min(expected.lines().count()));
            failures.push(format!(
                "{label}: line {} is {:?}, reference {:?}",
                at + 1,
                produced.lines().nth(at),
                expected.lines().nth(at)
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// What this port refuses says so and is not a user error of GATK's.
fn refusals() {
    let dump = "genomicsdb-query-dump";
    let to = |list: &[&str]| -> Vec<String> { list.iter().map(|s| (*s).to_string()).collect() };
    // Two import batches write two fragments to each array.
    let imported = gatk_cli::run(&to(&[
        "GenomicsDBImport",
        "-V",
        &format!("{dump}/cg_a.g.vcf"),
        "-V",
        &format!("{dump}/cg_b.g.vcf"),
        "-L",
        "chr1",
        "--batch-size",
        "1",
        "--genomicsdb-workspace-path",
        &format!("{dump}/ws-batches"),
    ]));
    assert_eq!(imported.status, 0, "{}", imported.stderr);
    let reference = format!("{dump}/cgv_ref.fasta");
    let select = |variant: &str, extra: &[&str]| -> Vec<String> {
        let mut argv = to(&[
            "SelectVariants",
            "-R",
            &reference,
            "-V",
            variant,
            "-O",
            &format!("{dump}/out/refused.vcf"),
        ]);
        argv.extend(to(extra));
        argv
    };
    let workspace = format!("gendb://{dump}/ws");
    for (label, argv) in [
        (
            "bcf-codec",
            select(&workspace, &["--genomicsdb-use-bcf-codec"]),
        ),
        ("cloud", select("gendb.gs://bucket/ws", &[])),
        (
            "two-fragments",
            select(&format!("gendb://{dump}/ws-batches"), &[]),
        ),
    ] {
        let ran = gatk_cli::run(&argv);
        assert_ne!(ran.status, 0, "{label}");
        assert!(ran.stderr.contains(PORT_OWN), "{label}: {}", ran.stderr);
    }
}

#[test]
fn the_callers_read_the_workspace_as_the_reference_does() {
    let root = std::env::temp_dir().join("gatk-cli-genomicsdb-query");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).expect("a scratch directory");
    let root = root.canonicalize().expect("a scratch directory");
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        eprintln!("skipping: {GOLDEN} is frozen from a CI run");
        return;
    }
    let text = corpus::read_golden(&golden);
    replay(&text, &root);
    refusals();
}
