//! What the SV suites share: a golden's runs, replayed from a scratch directory.

use gatk_corpus as corpus;

/// Above this the golden holds a file's length and CRC.
const WHOLE_FILE_LIMIT: usize = 16384;

fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for byte in bytes {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xEDB8_8320 & (crc & 1).wrapping_neg());
        }
    }
    !crc
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

/// Every path under `dir` relative to it, directories marked, without Hadoop's checksums.
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
            } else if !name.ends_with(".crc") {
                found.push((name, Some(std::fs::read(&path).expect("a written file"))));
            }
        }
    }
    found.sort();
    found
}

/// `produced N contigs in S secs.`, with its duration gone.
fn without_the_duration(bytes: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(bytes);
    let mut out = String::new();
    for line in text.split_inclusive('\n') {
        match (line.rfind(" contigs in "), line.rfind(" secs.")) {
            (Some(from), Some(to)) if from < to => {
                out.push_str(&line[..from]);
                out.push_str(" contigs in <S>");
                out.push_str(&line[to..]);
            }
            _ => out.push_str(line),
        }
    }
    out.into_bytes()
}

/// The date of a VCF's `##GATKCommandLine`.
fn without_the_date(bytes: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(bytes);
    let mut out = String::new();
    for line in text.split_inclusive('\n') {
        match line
            .find(",Date=\"")
            .filter(|_| line.starts_with("##GATKCommandLine"))
        {
            Some(at) => {
                out.push_str(&line[..at]);
                out.push_str(",Date=\"<D>\">");
                out.push_str(if line.ends_with('\n') { "\n" } else { "" });
            }
            None => out.push_str(line),
        }
    }
    out.into_bytes()
}

pub fn replay(text: &str, root: &std::path::Path) {
    std::fs::create_dir_all(root.join("sv-evidence-dump")).expect("a scratch directory");
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
                format!("sv-evidence-dump/{name}"),
                corpus::decode_base64(encoded),
            )
            .expect("a fixture"),
            ["argv", label, _] => labels.push(label.to_string()),
            _ => {}
        }
    }
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
        let out_dir = root.join("sv-evidence-dump/out").join(label);
        std::fs::create_dir_all(&out_dir).expect("an output directory");
        // The dump makes this one ahead of every run but the one that shows why.
        for pair in argv.windows(2) {
            if pair[0] == "--fastq-dir" && label != "find-fastq-no-dir" {
                std::fs::create_dir_all(&pair[1]).expect("a FASTQ directory");
            }
        }
        let ran = gatk_cli::run(&argv);
        let stdout = ran.stdout.replace(&cwd, "<cwd>");
        let expected = row("stdout", label).unwrap_or_default();
        if stdout != expected {
            failures.push(format!(
                "{label}: stdout {stdout:?}, reference {expected:?}"
            ));
        }
        let stderr = ran.stderr.replace(&cwd, "<cwd>");
        match row("error", label) {
            Some(error) => {
                // The dump prints a null message as `null`; `Throwable.toString` prints none.
                let message = match error.split_once(':') {
                    Some((_, "null")) | None => "",
                    Some((_, message)) => message,
                };
                if ran.status == 0 || !stderr.contains(message) {
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
        let mut reference: Vec<(String, Option<Vec<u8>>)> = Vec::new();
        let dir_prefix = format!("dir\t{label}\t");
        let file_prefix = format!("file\t{label}\t");
        let crc_prefix = format!("crc\t{label}\t");
        for line in text.lines() {
            if let Some(name) = line.strip_prefix(&dir_prefix) {
                reference.push((name.to_string(), None));
            } else if let Some(rest) = line.strip_prefix(&file_prefix) {
                let (name, encoded) = rest.split_once('\t').unwrap_or((rest, ""));
                reference.push((name.to_string(), Some(corpus::decode_base64(encoded))));
            } else if let Some(rest) = line.strip_prefix(&crc_prefix) {
                let (name, summary) = rest.split_once('\t').expect("a crc row");
                reference.push((name.to_string(), Some(summary.as_bytes().to_vec())));
            }
        }
        reference.sort();
        let produced: Vec<(String, Option<Vec<u8>>)> = under(&out_dir)
            .into_iter()
            .map(|(name, bytes)| {
                let bytes = if name == "intervals.txt" {
                    bytes.map(|bytes| without_the_duration(&bytes))
                } else if name.ends_with(".vcf") {
                    bytes.map(|bytes| without_the_date(&bytes))
                } else {
                    bytes
                };
                let bytes = bytes.map(|bytes| {
                    if bytes.len() > WHOLE_FILE_LIMIT {
                        format!("{}\t{}", bytes.len(), crc32(&bytes)).into_bytes()
                    } else {
                        bytes
                    }
                });
                (name, bytes)
            })
            .collect();
        let names = |files: &[(String, Option<Vec<u8>>)]| -> Vec<String> {
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
                failures.push(format!("{label}: {name} differs"));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} runs differ:\n{}",
        failures.len(),
        labels.len(),
        failures.join("\n")
    );
}
