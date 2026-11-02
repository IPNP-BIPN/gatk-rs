//! `GenomicsDBImport` replayed against GATK 4.6.2.0: its refusals, and every byte of the workspaces
//! it writes.
//!
//! Golden from `tools/readfilter-conformance/GenomicsDBImportDump.java`, which records each run's
//! command line, what it printed or why it failed, and everything in the workspace. Three things
//! differ from one run of the reference to the next and are masked on both sides: a fragment's
//! name, a batch's row-bounds file's name, and the `Date` in `vcfheader.vcf`. The array's absolute
//! path in `__array_schema.tdb` moves with the directory and is reduced to the array's own name.
//!
//! All of it is skipped until a CI run has produced the golden, and the suite that freezes the
//! tool's argument declarations with it: the parser cannot read this tool's arguments before that.
//! The refusals that are the port's own are checked beside the replay, since the reference has
//! nothing to say about them.

use gatk_corpus as corpus;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/genomicsdb_import.txt.gz";

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

fn is_dir_named(path: &std::path::Path, prefix: &str) -> bool {
    path.is_dir()
        && path
            .file_name()
            .is_some_and(|name| name.to_string_lossy().starts_with(prefix))
}

/// The dump's `stableBytes`.
fn stable_bytes(path: &std::path::Path) -> Vec<u8> {
    let raw = std::fs::read(path).expect("a written file");
    let name = path.file_name().expect("a name").to_string_lossy();
    if name == "vcfheader.vcf" {
        let text = String::from_utf8(raw).expect("a text header");
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
        return out.into_bytes();
    }
    if name == "__array_schema.tdb" {
        let length = i32::from_le_bytes(raw[8..12].try_into().expect("four bytes")) as usize;
        let array = String::from_utf8(raw[12..12 + length].to_vec()).expect("a path");
        let mut out = raw[..8].to_vec();
        out.extend_from_slice(array.rsplit('/').next().expect("a name").as_bytes());
        out.extend_from_slice(&raw[12 + length..]);
        return out;
    }
    raw
}

/// The dump's `stableName` walk: each path under the workspace, fragments and row bounds numbered.
fn workspace_rows(workspace: &std::path::Path) -> Vec<(String, Option<Vec<u8>>)> {
    fn walk(dir: &std::path::Path, prefix: &str, found: &mut Vec<(String, Option<Vec<u8>>)>) {
        let entries: Vec<std::path::PathBuf> = std::fs::read_dir(dir)
            .expect("a directory")
            .flatten()
            .map(|entry| entry.path())
            .collect();
        let mut fragments: Vec<&std::path::PathBuf> = entries
            .iter()
            .filter(|path| is_dir_named(path, "__"))
            .collect();
        fragments.sort_by_key(|path| {
            let name = path.file_name().expect("a name").to_string_lossy();
            name.rsplit('_')
                .next()
                .and_then(|stamp| stamp.parse::<u128>().ok())
                .expect("a fragment's timestamp")
        });
        let mut metas: Vec<&std::path::PathBuf> = entries
            .iter()
            .filter(|path| {
                !path.is_dir()
                    && path
                        .file_name()
                        .is_some_and(|name| name.to_string_lossy().starts_with("genomicsdb_meta_"))
            })
            .collect();
        metas.sort_by_key(|path| std::fs::read(path).expect("a row-bounds file"));
        for path in &entries {
            let name = path
                .file_name()
                .expect("a name")
                .to_string_lossy()
                .to_string();
            let stable = if let Some(index) = fragments.iter().position(|one| *one == path) {
                format!("fragment-{index}")
            } else if let Some(index) = metas.iter().position(|one| *one == path) {
                format!("genomicsdb_meta_{index}.json")
            } else {
                name
            };
            let relative = if prefix.is_empty() {
                stable
            } else {
                format!("{prefix}/{stable}")
            };
            if path.is_dir() {
                found.push((relative.clone(), None));
                walk(path, &relative, found);
            } else {
                found.push((relative, Some(stable_bytes(path))));
            }
        }
    }
    let mut found = Vec::new();
    if workspace.is_dir() {
        walk(workspace, "", &mut found);
    }
    found.sort();
    found
}

fn replay(text: &str, root: &std::path::Path) {
    std::fs::create_dir_all(root.join("genomicsdb-dump")).expect("a scratch directory");
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
                format!("genomicsdb-dump/{name}"),
                corpus::decode_base64(encoded),
            )
            .expect("a fixture"),
            ["argv", label, _] => labels.push(label.to_string()),
            _ => {}
        }
    }
    assert!(!labels.is_empty(), "the golden holds no runs");
    for name in [
        "a", "b", "c", "d", "big", "eq", "again", "two", "mnp", "extra", "nopass", "qual", "values",
    ] {
        // The runner reads the whole input and asks only that its index be there.
        std::fs::write(format!("genomicsdb-dump/{name}.g.vcf.idx"), b"").expect("an index");
    }

    std::fs::write("genomicsdb-dump/bgz.vcf.gz.tbi", b"").expect("an index");

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
        let out_dir = root.join("genomicsdb-dump/out").join(label);
        std::fs::create_dir_all(&out_dir).expect("an output directory");
        let ran = gatk_cli::run(&argv);
        let stderr = ran.stderr.replace(&format!("{cwd}/"), "<cwd>/");
        match row("error", label) {
            Some(error) => {
                let (class, message) = error.split_once(':').unwrap_or((error.as_str(), ""));
                let status = if class.starts_with("org.broadinstitute.barclay") {
                    1
                } else if class.contains("hellbender.exceptions.UserException")
                    || class.contains("UnableToCreateGenomicsDBWorkspace")
                {
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
            None if ran.status != 0 => {
                failures.push(format!("{label}: status {}: {stderr}", ran.status));
                continue;
            }
            None => {}
        }
        let returned = row("returned", label).expect("a returned row");
        let expected = format!(
            "{}Tool returned:\n{returned}\n",
            row("stdout", label).expect("a stdout row")
        );
        let stdout = ran.stdout.replace(&format!("{cwd}/"), "<cwd>/");
        if stdout != expected {
            failures.push(format!(
                "{label}: stdout {stdout:?}, reference {expected:?}"
            ));
        }

        let mut reference: Vec<(String, Option<Vec<u8>>)> = Vec::new();
        let dir_prefix = format!("dir\t{label}\t");
        let file_prefix = format!("file\t{label}\t");
        for line in text.lines() {
            if let Some(name) = line.strip_prefix(&dir_prefix) {
                reference.push((name.to_string(), None));
            } else if let Some(rest) = line.strip_prefix(&file_prefix) {
                let (name, encoded) = rest.split_once('\t').expect("a file row");
                reference.push((name.to_string(), Some(corpus::decode_base64(encoded))));
            }
        }
        reference.sort();
        let produced = workspace_rows(&out_dir.join("ws"));
        let names = |rows: &[(String, Option<Vec<u8>>)]| -> Vec<String> {
            rows.iter().map(|(name, _)| name.clone()).collect()
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
                failures.push(format!("{label}/{name} differs from the reference"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

const HEADER: &str = "##fileformat=VCFv4.2\n\
    ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
    ##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
    ##contig=<ID=chr1,length=1000>\n\
    ##contig=<ID=chr2,length=500>\n";

fn vcf(sample: &str, records: &str) -> String {
    format!("{HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n{records}")
}

/// What this port refuses says so, is not a user error of GATK's, and leaves no workspace.
fn refusals(root: &std::path::Path) {
    let dir = root.join("refusals");
    std::fs::create_dir_all(&dir).expect("a directory");
    std::env::set_current_dir(&dir).expect("the directory");
    let record = "chr1\t5\t.\tA\tC,<NON_REF>\t.\t.\t.\tGT\t0/1\n";
    for (name, text) in [
        ("a.g.vcf", vcf("A", record)),
        ("b.g.vcf", vcf("B", record)),
        (
            "undeclared.g.vcf",
            vcf(
                "U",
                "chr1\t5\t.\tA\tC,<NON_REF>\t.\t.\tNOT_IN_HEADER=1\tGT\t0/1\n",
            ),
        ),
    ] {
        std::fs::write(name, text).expect("a fixture");
        std::fs::write(format!("{name}.idx"), b"").expect("an index");
    }
    let to = |list: &[&str]| -> Vec<String> { list.iter().map(|s| (*s).to_string()).collect() };
    for (label, rest) in [
        (
            "update",
            to(&[
                "--genomicsdb-update-workspace-path",
                "ws-update",
                "-V",
                "a.g.vcf",
            ]),
        ),
        (
            "merge-contigs",
            to(&[
                "--genomicsdb-workspace-path",
                "ws-merge",
                "-V",
                "a.g.vcf",
                "-L",
                "chr1",
                "-L",
                "chr2",
                "--merge-contigs-into-num-partitions",
                "1",
            ]),
        ),
        (
            "consolidate-batches",
            to(&[
                "--genomicsdb-workspace-path",
                "ws-consolidate",
                "-V",
                "a.g.vcf",
                "-V",
                "b.g.vcf",
                "-L",
                "chr1",
                "--consolidate",
                "--batch-size",
                "1",
            ]),
        ),
        (
            "undeclared-field",
            to(&[
                "--genomicsdb-workspace-path",
                "ws-undeclared",
                "-V",
                "undeclared.g.vcf",
                "-L",
                "chr1",
            ]),
        ),
    ] {
        let mut argv = vec!["GenomicsDBImport".to_string()];
        argv.extend(rest.clone());
        let ran = gatk_cli::run(&argv);
        assert_ne!(ran.status, 0, "{label}");
        assert!(ran.stderr.contains(PORT_OWN), "{label}: {}", ran.stderr);
        for argument in rest.iter().filter(|argument| argument.starts_with("ws-")) {
            assert!(
                !std::path::Path::new(argument).exists(),
                "{label} left {argument} behind"
            );
        }
    }
}

#[test]
fn the_workspace_is_the_references_and_the_refusals_say_whose_they_are() {
    let root = std::env::temp_dir().join("gatk-cli-genomicsdb-import");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).expect("a scratch directory");
    let root = root.canonicalize().expect("a scratch directory");
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        eprintln!("skipping: {GOLDEN} and the tool's declarations are frozen from a CI run");
        return;
    }
    refusals(&root);
    let text = corpus::read_golden(&golden);
    replay(&text, &root);
}
