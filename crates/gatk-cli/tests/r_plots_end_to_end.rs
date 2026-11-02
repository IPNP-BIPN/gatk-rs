//! The three tools that hand their work to an R script, replayed against GATK 4.6.2.0 with R on
//! both sides.
//!
//! Reference from `tools/readfilter-conformance/RPlotsDump.java`, in the Broad's own image of the
//! release, since the pinned oracle has no `Rscript`. The port runs the rewrite under
//! `crates/gatk-cli/src/r/`.
//!
//! # What this suite is for
//!
//! It does not compare a pixel: the rewrite is on ggplot2, and cairo's rasterisation and the font
//! metrics differ besides. What it holds is everything around the image:
//!
//!  * **the side files, byte for byte**, `PlotDenoisedCopyRatios`' four MADs, including the `NA` a
//!    one-point file writes;
//!  * **the status of every run and the first line of its refusal**, down to the bare
//!    `quit(status=1)` of a script run on its own, and the plot that fails on a limit that is not
//!    finite with its MADs already written;
//!  * **the files each run leaves**, in its output and in the working directory, where the
//!    rewrite adds exactly what `EXPORTS` names and only on a run that succeeds;
//!  * **each PNG's pixel size and resolution**;
//!  * **and the numbers plotted**: the exported `.data.rds`, read back through `RPlotsReadBack.R`,
//!    against what the reference's own scripts handed their graphics calls, recorded through
//!    `RPlotsTrace.R`. That pins the log2-to-linear transform, the midpoint rounding, the offsets,
//!    each panel's limit, the running index with its shared point and its run off the table, and
//!    the inverse logit, the `called` selection and the inner join of the allele-frequency script.
//!
//! # What it needs
//!
//! An `Rscript` with `optparse`, `data.table`, `ggplot2`, `cowplot` and `dplyr`, which the `test`
//! job does not install, so a machine without one prints that and passes; `GATK_RS_REQUIRE_R` makes
//! that a failure, and the `oracle-with-r` job sets it. The reference is the frozen golden, or
//! `GATK_RS_R_PLOTS_DUMP`, the path of a dump just regenerated, which is how that job compares
//! before anything is frozen.

use gatk_corpus as corpus;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/r_plots.txt.gz";

/// A package nothing could ever install.
const NO_SUCH_PACKAGE: &str = "gatkRsNoSuchPackage";

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

fn manifest_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn r_source(name: &str) -> String {
    manifest_path(&format!("src/r/{name}")).display().to_string()
}

fn read_back_script() -> String {
    manifest_path("../../tools/readfilter-conformance/RPlotsReadBack.R")
        .display()
        .to_string()
}

/// Which script a run went through.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Denoised,
    Modeled,
    Qc,
}

impl Kind {
    fn of_tool(tool: &str) -> Kind {
        match tool {
            "PlotDenoisedCopyRatios" => Kind::Denoised,
            "PlotModeledSegments" => Kind::Modeled,
            _ => Kind::Qc,
        }
    }

    fn of_script(name: &str) -> Kind {
        match name {
            "denoised" => Kind::Denoised,
            "modeled" => Kind::Modeled,
            _ => Kind::Qc,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Kind::Denoised => "denoised",
            Kind::Modeled => "modeled",
            Kind::Qc => "qc",
        }
    }

    /// `Rscript -e` over the port's own files, in the order the binary sources them.
    fn expression(self) -> String {
        let mut names = vec!["RPlotBootstrap.R"];
        match self {
            Kind::Denoised => names.extend(["CNVPlottingLibrary.R", "PlotDenoisedCopyRatios.R"]),
            Kind::Modeled => names.extend(["CNVPlottingLibrary.R", "PlotModeledSegments.R"]),
            Kind::Qc => names.push("plotAlleleFrequencyQC.R"),
        }
        names
            .iter()
            .map(|name| format!("source('{}');", r_source(name)))
            .collect()
    }

    /// What the rewrite writes that the reference does not, each beside the file it derives from.
    fn exports(self, prefix: &str) -> Vec<String> {
        match self {
            Kind::Denoised => vec![
                format!("{prefix}sample.denoised.data.rds"),
                format!("{prefix}sample.denoised.ggplot.rds"),
            ],
            Kind::Modeled => vec![
                format!("{prefix}sample.modeled.data.rds"),
                format!("{prefix}sample.modeled.ggplot.rds"),
            ],
            Kind::Qc => vec![
                format!("{prefix}qc.data.rds"),
                format!("{prefix}qc.af.ggplot.rds"),
                format!("{prefix}qc.af_differences.ggplot.rds"),
                format!("{prefix}qc.af_combined.png"),
            ],
        }
    }
}

/// Every file under `dir`, relative to it.
fn files_under(dir: &Path) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
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
                found.insert(
                    path.strip_prefix(dir)
                        .expect("under the directory")
                        .to_string_lossy()
                        .to_string(),
                );
            }
        }
    }
    found
}

/// What a run added to the working directory, as `cwd/` rows, gone afterwards.
fn take_working_directory() -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    for entry in std::fs::read_dir(".").expect("the scratch directory").flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if name != "r-plots-dump" {
            std::fs::remove_file(entry.path()).expect("a file in the scratch directory");
            found.insert(format!("cwd/{name}"));
        }
    }
    found
}

/// IHDR's width and height and pHYs's pixels per metre, as `RPlotsDump.pngShape` prints them.
fn png_shape(path: &Path) -> String {
    let bytes = std::fs::read(path).expect("a PNG");
    let number = |at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().expect("four bytes"));
    let mut per_metre = "-".to_string();
    let mut at = 8;
    while at + 8 <= bytes.len() {
        let kind = &bytes[at + 4..at + 8];
        if kind == b"pHYs" {
            per_metre = number(at + 8).to_string();
        }
        if kind == b"IDAT" {
            break;
        }
        at += 12 + number(at) as usize;
    }
    format!("{}x{}\t{per_metre}", number(16), number(20))
}

/// The rows of one kind and label, each with the fields after the label.
fn rows<'a>(text: &'a str, kind: &str, label: &str) -> Vec<&'a str> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .collect()
}

/// The rows the port's export reads back as.
fn read_back(kind: Kind, frame: &Path, root: &Path) -> Result<String, String> {
    let into = root.join("r-plots-dump/read-back.txt");
    let ran = Command::new("Rscript")
        .arg(read_back_script())
        .args([kind.name(), &frame.display().to_string()])
        .arg(&into)
        .output()
        .expect("Rscript");
    if !ran.status.success() {
        return Err(String::from_utf8_lossy(&ran.stderr).to_string());
    }
    std::fs::read_to_string(&into).map_err(|error| error.to_string())
}

/// The export is a ggplot and comes back as one.
fn reads_as_ggplot(figure: &Path) -> bool {
    Command::new("Rscript")
        .args([
            "-e",
            "x = readRDS(commandArgs(TRUE)[1]); if (!inherits(x, 'ggplot')) quit(status=1)",
        ])
        .arg(figure)
        .output()
        .expect("Rscript")
        .status
        .success()
}

/// One replayed run, as far as comparing its outputs goes.
struct Run<'a> {
    label: &'a str,
    kind: Kind,
    out_dir: &'a Path,
    prefix: &'a str,
    succeeded: bool,
}

/// What a run produced, against the reference's list and what each file holds there.
fn compare_outputs(text: &str, run: &Run, root: &Path, failures: &mut Vec<String>) {
    let Run {
        label,
        kind,
        out_dir,
        prefix,
        succeeded,
    } = *run;
    let mut produced = files_under(out_dir);
    produced.extend(take_working_directory());
    let mut expected: BTreeSet<String> = rows(text, "produced", label)
        .iter()
        .map(|row| (*row).to_string())
        .collect();
    let exports = if succeeded { kind.exports(prefix) } else { Vec::new() };
    expected.extend(exports.iter().cloned());
    if produced != expected {
        failures.push(format!(
            "{label}: files {produced:?}, expected {expected:?}"
        ));
        return;
    }
    for row in rows(text, "side", label) {
        let (name, value) = row.split_once('\t').expect("a side row");
        let ours = std::fs::read_to_string(out_dir.join(name)).expect("a side file");
        if ours != unescape(value) {
            failures.push(format!("{label}/{name}: {ours:?}, reference {value:?}"));
        }
    }
    for row in rows(text, "png", label) {
        let (name, shape) = row.split_once('\t').expect("a png row");
        let ours = png_shape(&out_dir.join(name));
        if ours != shape {
            failures.push(format!("{label}/{name}: {ours:?}, reference {shape:?}"));
        }
    }
    for row in rows(text, "metrics", label) {
        let name = format!("{prefix}qc.metrics");
        let table: Vec<String> = std::fs::read_to_string(out_dir.join(&name))
            .expect("the metrics")
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect();
        if table.join("\n") != unescape(row) {
            failures.push(format!("{label}: metrics {table:?}, reference {row:?}"));
        }
    }
    if !succeeded {
        return;
    }
    // The numbers, read back out of the export.
    let frame = out_dir.join(&exports[0]);
    match read_back(kind, &frame, root) {
        Ok(ours) => {
            let ours: Vec<&str> = ours.lines().collect();
            let reference: Vec<String> =
                rows(text, "trace", label).iter().map(|row| unescape(row)).collect();
            if ours != reference {
                failures.push(format!(
                    "{label}: the export reads back as\n{}\nand the reference plotted\n{}",
                    ours.join("\n"),
                    reference.join("\n")
                ));
            }
        }
        Err(error) => failures.push(format!("{label}: reading the export back: {error}")),
    }
    for export in exports.iter().filter(|name| name.ends_with(".ggplot.rds")) {
        if !reads_as_ggplot(&out_dir.join(export)) {
            failures.push(format!("{label}/{export}: not a ggplot"));
        }
    }
}

fn replay(text: &str, root: &Path, failures: &mut Vec<String>) {
    std::fs::create_dir_all(root.join("r-plots-dump/out")).expect("a scratch directory");
    std::env::set_current_dir(root).expect("the scratch directory");
    let cwd = std::env::current_dir()
        .expect("a cwd")
        .display()
        .to_string();

    for line in text.lines() {
        if let ["input", name, content] = line.splitn(3, '\t').collect::<Vec<&str>>()[..] {
            std::fs::write(format!("r-plots-dump/{name}"), unescape(content)).expect("a fixture");
        }
    }
    for name in ["eval.vcf", "comp.vcf"] {
        let indexed = gatk_cli::run(&[
            "IndexFeatureFile".to_string(),
            "-I".to_string(),
            format!("r-plots-dump/{name}"),
        ]);
        assert_eq!(indexed.status, 0, "{}", indexed.stderr);
    }

    let labels: Vec<String> = text
        .lines()
        .filter_map(|line| line.strip_prefix("argv\t"))
        .map(|rest| rest.split('\t').next().expect("a label").to_string())
        .collect();
    assert!(!labels.is_empty(), "the reference holds no runs");
    for label in &labels {
        let argv: Vec<String> = unescape(rows(text, "argv", label)[0])
            .split('\t')
            .map(str::to_string)
            .collect();
        let out_dir = root.join("r-plots-dump/out").join(label);
        std::fs::create_dir_all(&out_dir).expect("an output directory");
        let ran = gatk_cli::run(&argv);
        let stdout = ran.stdout.replace(&cwd, "<cwd>");
        let expected = unescape(rows(text, "stdout", label)[0]);
        if stdout != expected {
            failures.push(format!(
                "{label}: stdout {stdout:?}, reference {expected:?}"
            ));
        }
        let stderr = ran.stderr.replace(&cwd, "<cwd>");
        let refusal = rows(text, "error", label).first().map(|row| unescape(row));
        let status = match &refusal {
            None => 0,
            Some(error) => {
                let class = error.split_once(':').map_or(error.as_str(), |(class, _)| class);
                if class.starts_with("org.broadinstitute.barclay") {
                    1
                } else if class.contains("hellbender.exceptions.UserException") {
                    2
                } else {
                    3
                }
            }
        };
        let message = refusal
            .as_deref()
            .and_then(|error| error.split_once(':'))
            .map_or("", |(_, message)| message);
        if ran.status != status || !stderr.contains(message) {
            failures.push(format!(
                "{label}: status {}, reference {status} refusing with {refusal:?}; stderr {stderr}",
                ran.status
            ));
        }
        let run = Run {
            label,
            kind: Kind::of_tool(&argv[0]),
            out_dir: &out_dir,
            prefix: if argv[0] == "AlleleFrequencyQC" {
                ""
            } else {
                "plots/"
            },
            succeeded: status == 0,
        };
        compare_outputs(text, &run, root, failures);
    }

    // The scripts on their own, on arguments Java's own checks would have refused first.
    let scripts: Vec<&str> = text
        .lines()
        .filter_map(|line| line.strip_prefix("script\t"))
        .collect();
    assert!(!scripts.is_empty(), "the reference holds no script runs");
    for script in scripts {
        let fields: Vec<&str> = script.splitn(3, '\t').collect();
        let [label, kind, arguments] = fields[..] else {
            panic!("a script row");
        };
        let kind = Kind::of_script(kind);
        let out_dir = root.join("r-plots-dump/out").join(label);
        std::fs::create_dir_all(&out_dir).expect("an output directory");
        let ran = Command::new("Rscript")
            .arg("-e")
            .arg(kind.expression())
            .args(unescape(arguments).split('\t'))
            .output()
            .expect("Rscript");
        let status = ran.status.code().expect("a status").to_string();
        let expected = rows(text, "script-status", label)[0];
        if status != expected {
            failures.push(format!(
                "{label}: status {status}, reference {expected}: {}",
                String::from_utf8_lossy(&ran.stderr)
            ));
            continue;
        }
        let run = Run {
            label,
            kind,
            out_dir: &out_dir,
            prefix: "",
            succeeded: status == "0",
        };
        compare_outputs(text, &run, root, failures);
    }
}

/// A package that cannot be had is named on stderr and the run ends with status 1, since the
/// status is the only channel the Java side reads.
fn bootstrap_refusal(root: &Path, failures: &mut Vec<String>) {
    let library = root.join("r-plots-dump/created-on-demand/library");
    let ran = Command::new("Rscript")
        .arg("-e")
        .arg(format!(
            "source('{}');EnsurePackages('{NO_SUCH_PACKAGE}')",
            r_source("RPlotBootstrap.R")
        ))
        .env("GATK_R_REPOSITORY", "file:///nothing/here")
        .env("GATK_R_LIBRARY", &library)
        .output()
        .expect("Rscript");
    let stderr = String::from_utf8_lossy(&ran.stderr);
    if ran.status.code() != Some(1) || !stderr.contains(NO_SUCH_PACKAGE) {
        failures.push(format!(
            "the bootstrap: status {:?}, stderr {stderr}",
            ran.status.code()
        ));
    }
    if !library.is_dir() {
        failures.push("the bootstrap: the user library was not created".to_string());
    }
}

#[test]
fn the_plots_around_their_images_are_the_references() {
    let golden = manifest_path(GOLDEN);
    let text = match std::env::var("GATK_RS_R_PLOTS_DUMP") {
        Ok(dump) => std::fs::read_to_string(dump).expect("the dump"),
        Err(_) if golden.exists() => corpus::read_golden(&golden),
        Err(_) => {
            // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
            eprintln!("skipping: {GOLDEN} is not committed yet");
            return;
        }
    };
    let has_rscript = Command::new("Rscript")
        .arg("--version")
        .output()
        .is_ok_and(|ran| ran.status.success());
    if !has_rscript {
        assert!(
            std::env::var_os("GATK_RS_REQUIRE_R").is_none(),
            "GATK_RS_REQUIRE_R is set and there is no Rscript"
        );
        eprintln!("skipping: no Rscript on this machine");
        return;
    }

    let root = std::env::temp_dir().join(format!("gatk-rs-r-plots-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).expect("a scratch directory");
    let root = root.canonicalize().expect("a scratch directory");

    let mut failures = Vec::new();
    replay(&text, &root, &mut failures);
    bootstrap_refusal(&root, &mut failures);
    std::env::set_current_dir(manifest_path(".")).expect("the crate");
    let _ = std::fs::remove_dir_all(&root);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
