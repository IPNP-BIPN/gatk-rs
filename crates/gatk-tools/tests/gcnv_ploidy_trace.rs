//! The ploidy model against a trace of the reference's compiled graphs, call by call: every ADVI
//! step's loss, `mu` and `rho`, every log-emission draw, every caller update.
//!
//! The trace comes from `tools/python-tools/ploidy_trace.sh` (`trace/` under its output). It is a
//! debugging aid, not a golden: without `GCNV_PLOIDY_TRACE` this test says so and passes; the
//! golden check is `crates/gatk-cli/tests/gcnv_ploidy_end_to_end.rs`.

use gatk_tools::gcnv_ploidy::{run_watched, Cohort, Event, Settings};
use std::path::Path;

fn find(dir: &Path, prefix: &str) -> String {
    let entry = std::fs::read_dir(dir)
        .expect("the trace directory")
        .map(|e| e.expect("an entry").path())
        .find(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with(prefix))
        })
        .unwrap_or_else(|| panic!("no {prefix}* in the trace"));
    std::fs::read_to_string(entry).expect("a file")
}

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

fn cohort(dir: &Path) -> Cohort {
    let coverage = find(dir, "in.samples-by-coverage");
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
    let intervals = find(dir, "in.intervals");
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
        priors_text: find(dir, "in.priors"),
    }
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

#[test]
fn the_graphs_are_replayed_call_by_call() {
    let Some(dir) = std::env::var_os("GCNV_PLOIDY_TRACE") else {
        eprintln!("GCNV_PLOIDY_TRACE is not set: nothing to replay");
        return;
    };
    let dir = Path::new(&dir);
    let trace = std::fs::read_to_string(dir.join("trace.tsv")).expect("trace.tsv");
    let mut wanted: Vec<(String, Vec<u64>)> = Vec::new();
    // Each step's loss, then its updates: `rho` before `mu`; the port reports loss, mu, rho.
    let mut rho: Vec<u64> = Vec::new();
    for line in trace.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let values = || -> Vec<u64> {
            fields[2]
                .split(' ')
                .filter(|v| !v.is_empty())
                .map(|v| u64::from_str_radix(v, 16).expect("hex bits"))
                .collect()
        };
        match fields[0] {
            "f2.out0" => wanted.push(("step".into(), values())),
            "f2.upd.2.rho" => rho = values(),
            "f2.upd.4.mu" => {
                let step = &mut wanted.last_mut().expect("a step").1;
                step.extend(values());
                step.append(&mut rho);
            }
            "f3.out0" => wanted.push(("draw".into(), values())),
            "f4.out0" => wanted.push(("call".into(), values())),
            _ => {}
        }
    }
    let argv = std::fs::read_to_string(dir.join("argv.txt")).expect("argv.txt");
    let mut got: Vec<(String, Vec<u64>)> = Vec::new();
    run_watched(&cohort(dir), &settings(&argv), &mut |event| match event {
        Event::Step { loss, mu, rho } => {
            let mut v = vec![loss.to_bits()];
            v.extend(bits(mu));
            v.extend(bits(rho));
            got.push(("step".into(), v));
        }
        Event::Draw(d) => got.push(("draw".into(), bits(d))),
        Event::Call(c) => got.push(("call".into(), bits(c))),
    })
    .expect("the run");
    for (i, (w, g)) in wanted.iter().zip(&got).enumerate() {
        if w != g {
            let first = w.1.iter().zip(&g.1).position(|(a, b)| a != b);
            panic!(
                "event {i} ({}) differs at value {first:?}: want {:?}, got {:?}",
                w.0,
                first.map(|p| f64::from_bits(w.1[p])),
                first.map(|p| f64::from_bits(g.1[p]))
            );
        }
    }
    assert_eq!(wanted.len(), got.len(), "as many events");
}
