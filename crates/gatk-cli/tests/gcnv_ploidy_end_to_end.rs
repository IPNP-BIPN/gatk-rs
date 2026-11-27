//! `DetermineGermlineContigPloidy` against the frozen reference: the fixture of
//! `tools/python-tools/make_gcnv_fixture.py` (stdlib only, seeded), the command line of
//! `measure_gcnv.sh`, and every `ploidy/` file of the `gcnv-frozen` golden, byte for byte.
//!
//! The port hands gcnvkernel GATK's own script, as the reference does. Here a `python` stands in
//! for it and keeps what the Java side handed over: the argument list and the two temporary tables.
//! The native port of the script, `gatk_tools::gcnv_ploidy`, then runs on exactly those, so the
//! check is the CLI's Java and the library together. `PATH` is process-wide, so this file holds one
//! test.

use gatk_tools::gcnv_ploidy::{Cohort, Settings};
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

/// The script's arguments as python reads them: `--name=value`.
fn settings(argv: &str) -> Settings {
    let get = |name: &str| -> String {
        argv.lines()
            .find_map(|l| l.strip_prefix(&format!("--{name}=")))
            .unwrap_or_else(|| panic!("--{name}"))
            .to_string()
    };
    let f = |name: &str| get(name).parse::<f64>().expect("a float");
    let i = |name: &str| get(name).parse::<i64>().expect("an integer");
    Settings {
        mapping_error_rate: f("mapping_error_rate"),
        psi_s_scale: f("psi_s_scale"),
        mean_bias_sd: f("mean_bias_sd"),
        psi_j_scale: f("psi_j_scale"),
        learning_rate: f("learning_rate"),
        adamax_beta1: f("adamax_beta1"),
        adamax_beta2: f("adamax_beta2"),
        log_emission_samples_per_round: i("log_emission_samples_per_round"),
        log_emission_sampling_rounds: i("log_emission_sampling_rounds"),
        log_emission_sampling_median_rel_error: f("log_emission_sampling_median_rel_error"),
        max_advi_iter_first_epoch: i("max_advi_iter_first_epoch"),
        max_advi_iter_subsequent_epochs: i("max_advi_iter_subsequent_epochs"),
        min_training_epochs: i("min_training_epochs"),
        max_training_epochs: i("max_training_epochs"),
        initial_temperature: f("initial_temperature"),
        num_thermal_advi_iters: i("num_thermal_advi_iters"),
        convergence_snr_averaging_window: i("convergence_snr_averaging_window"),
        convergence_snr_trigger_threshold: f("convergence_snr_trigger_threshold"),
        convergence_snr_countdown_window: i("convergence_snr_countdown_window"),
        max_calling_iters: i("max_calling_iters"),
        caller_update_convergence_threshold: f("caller_update_convergence_threshold"),
        caller_external_admixing_rate: f("caller_external_admixing_rate"),
        disable_annealing: get("disable_annealing") == "true",
        disable_caller: get("disable_caller") == "true",
        disable_sampler: get("disable_sampler") == "true",
    }
}

/// The samples-by-coverage-per-contig table and the intervals, as the script reads them.
fn cohort(coverage: &str, intervals: &str, priors_text: String) -> Cohort {
    let mut lines = coverage.lines().filter(|l| !l.starts_with('@'));
    let contigs: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .skip(1)
        .map(String::from)
        .collect();
    let samples = lines
        .map(|l| {
            let mut f = l.split('\t');
            let name = f.next().expect("a name").to_string();
            (name, f.map(|c| c.parse().expect("a count")).collect())
        })
        .collect();
    let intervals_per_contig = contigs
        .iter()
        .map(|c| {
            intervals
                .lines()
                .filter(|l| !l.starts_with('@') && l.split('\t').next() == Some(c.as_str()))
                .count() as u64
        })
        .collect();
    Cohort {
        contigs,
        intervals_per_contig,
        samples,
        priors_text,
    }
}

fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
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
    let kept = scratch.join("kept");
    for dir in [&bin, &root, &fixture, &kept] {
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
    // The import check passes; the script's run keeps its arguments and the two tables.
    executable(
        &bin.join("python"),
        &format!(
            "#!/bin/sh\n[ \"$1\" = -c ] && exit 0\nkept={}\nfor a in \"$@\"; do echo \"$a\"; done > $kept/argv.txt\n\
             for a in \"$@\"; do case \"$a\" in\n\
             --sample_coverage_metadata=*) /bin/cp \"${{a#*=}}\" $kept/coverage.tsv;;\n\
             --interval_list=*) /bin/cp \"${{a#*=}}\" $kept/intervals.tsv;;\n\
             --contig_ploidy_prior_table=*) /bin/cp \"${{a#*=}}\" $kept/priors.tsv;;\n\
             esac; done\n",
            kept.display()
        ),
    );

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

    // The reference's own refusal where the environment lacks the package.
    executable(
        &bin.join("python"),
        "#!/bin/sh\ncase \"$2\" in *gcnvkernel*) exit 1;; esac\nexit 0\n",
    );
    let refused = ploidy(&fixture);
    std::env::set_current_dir(&previous).expect("the directory");
    ran.expect("the run");

    let read = |name: &str| std::fs::read_to_string(kept.join(name)).expect("a kept file");
    let argv = read("argv.txt");
    let files = gatk_tools::gcnv_ploidy::run(
        &cohort(
            &read("coverage.tsv"),
            &read("intervals.tsv"),
            read("priors.tsv"),
        ),
        &settings(&argv),
    )
    .expect("the native model");
    // The script copies these two of its inputs into the model.
    let mut got: Vec<(String, String)> = vec![
        (
            "ploidy/p-model/interval_list.tsv".into(),
            read("intervals.tsv"),
        ),
        (
            "ploidy/p-model/contig_ploidy_prior.tsv".into(),
            read("priors.tsv"),
        ),
    ];
    got.extend(
        files
            .model
            .into_iter()
            .map(|(n, t)| (format!("ploidy/p-model/{n}"), t)),
    );
    got.extend(
        files
            .calls
            .into_iter()
            .map(|(n, t)| (format!("ploidy/p-calls/{n}"), t)),
    );
    let _ = std::fs::remove_dir_all(&scratch);

    let mut failures = Vec::new();
    for (path, size, sha) in &wanted {
        let ours = got
            .iter()
            .find(|(name, _)| name == path)
            .map(|(_, t)| t.as_bytes());
        if ours.is_none_or(|b| b.len().to_string() != *size || sha256(b) != *sha) {
            failures.push(path.clone());
        }
    }
    assert!(failures.is_empty(), "differs: {failures:?}");
    assert!(
        argv.lines().any(|l| l == "--psi_s_scale=1.000000e-04"),
        "{argv}"
    );
    let message = refused.expect_err("the missing package");
    assert!(
        message.contains("A required Python package (\"gcnvkernel\") could not be imported"),
        "{message}"
    );
}
