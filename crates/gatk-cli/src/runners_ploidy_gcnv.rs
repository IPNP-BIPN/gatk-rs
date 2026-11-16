//! `DetermineGermlineContigPloidy` past its import check, cohort mode: `doWork`'s Java (the
//! intervals, the counts summed per contig, the arguments handed to the script) and then
//! `cohort_determine_ploidy_and_depth.py`, ported in [`gatk_tools::gcnv_ploidy`] against the frozen
//! reference of `docs/gcnv-ploidy.md`.
//!
//! Case mode (`--model`) runs another script, which is not ported, and says so in the port's words.

use super::*;
use gatk_tools::gcnv_ploidy::{self, Cohort, Settings};

/// `String.format("%e", value)` read back by Python's `float`, which is what the script gets.
fn through_percent_e(value: f64) -> f64 {
    gatk_engine::pair_hmm_likelihood_engine::java_scientific(value)
        .parse()
        .unwrap_or(value)
}

fn double(parser: &Parser, name: &str, default: f64) -> f64 {
    through_percent_e(
        scalar(parser, name)
            .and_then(|text| text.parse().ok())
            .unwrap_or(default),
    )
}

fn integer(parser: &Parser, name: &str, default: i64) -> i64 {
    scalar(parser, name)
        .and_then(|text| text.parse().ok())
        .unwrap_or(default)
}

/// A read-counts TSV: its `@` header, then `CONTIG START END COUNT`.
struct Counts {
    header: String,
    sample: String,
    rows: Vec<(String, i64, i64, i64)>,
}

fn read_counts(path: &str) -> Result<Counts, Thrown> {
    if path.ends_with(".hdf5") || path.ends_with(".h5") {
        return Err(port_limit("Reading an HDF5 read-counts file"));
    }
    let text =
        std::fs::read_to_string(path).map_err(|error| Thrown::user(format!("{path}: {error}")))?;
    let header: String = text
        .lines()
        .take_while(|line| line.starts_with('@'))
        .map(|line| format!("{line}\n"))
        .collect();
    let sample = header
        .lines()
        .filter(|line| line.starts_with("@RG"))
        .flat_map(|line| line.split('\t'))
        .find_map(|field| field.strip_prefix("SM:"))
        .unwrap_or_default()
        .to_string();
    let rows = text
        .lines()
        .filter(|line| !line.starts_with('@') && !line.trim().is_empty())
        .skip(1)
        .filter_map(|line| {
            let c: Vec<&str> = line.split('\t').collect();
            Some((
                c.first()?.to_string(),
                c.get(1)?.parse().ok()?,
                c.get(2)?.parse().ok()?,
                c.get(3)?.parse().ok()?,
            ))
        })
        .collect();
    Ok(Counts {
        header,
        sample,
        rows,
    })
}

/// `DetermineGermlineContigPloidy.doWork`, cohort mode.
pub fn run(parser: &Parser) -> Outcome {
    let inputs = arguments(parser, "input");
    let mut seen = std::collections::HashSet::new();
    if !inputs.iter().all(|path| seen.insert(path.clone())) {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "List of input read-count files cannot contain duplicates.",
        ));
    }
    let output = require(parser, "output");
    let prefix = require(parser, "output-prefix");
    if argument(parser, "model").is_some() {
        return Err(port_limit(
            "Case mode (--model, case_determine_ploidy_and_depth.py)",
        ));
    }
    if inputs.len() < 2 {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "At least two samples must be provided in COHORT mode.",
        ));
    }
    let Some(priors) = argument(parser, "contig-ploidy-priors") else {
        return Err(Thrown::user(
            "Bad input: Contig-ploidy priors must be provided in cohort mode.",
        ));
    };
    validate_copy_number_intervals(parser)?;

    // `resolveIntervals` against the first counts file's dictionary.
    let first = read_counts(&inputs[0])?;
    let dictionary_text: String = std::iter::once("@HD\tVN:1.6\n".to_string())
        .chain(
            first
                .header
                .lines()
                .filter(|line| line.starts_with("@SQ"))
                .map(|line| format!("{line}\n")),
        )
        .collect();
    let header = SamHeader {
        sequences: tables::read_interval_collection(&inputs[0])?,
        ..SamHeader::default()
    };
    let intervals: Vec<(String, i64, i64)> = match interval_arguments(parser, &header)? {
        Some(parameters) => parameters
            .intervals
            .iter()
            .map(|i| (i.contig.clone(), i64::from(i.start), i64::from(i.end)))
            .collect(),
        None => first
            .rows
            .iter()
            .map(|(c, s, e, _)| (c.clone(), *s, *e))
            .collect(),
    };
    let wanted: std::collections::HashSet<&(String, i64, i64)> = intervals.iter().collect();
    let mut contigs: Vec<String> = Vec::new();
    for (contig, _, _) in &intervals {
        if !contigs.contains(contig) {
            contigs.push(contig.clone());
        }
    }

    // `writeSamplesByCoveragePerContig`: each file subset to the intervals, summed per contig.
    let mut samples = Vec::new();
    for path in &inputs {
        let counts = read_counts(path)?;
        let mut per_contig = vec![0i64; contigs.len()];
        let mut kept = 0;
        for (contig, start, end, count) in &counts.rows {
            if wanted.contains(&(contig.clone(), *start, *end)) {
                kept += 1;
                let j = contigs.iter().position(|c| c == contig).expect("a contig");
                per_contig[j] += count;
            }
        }
        if kept != wanted.len() {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!(
                    "Intervals for read-count file {path} do not contain all specified intervals."
                ),
            ));
        }
        samples.push((counts.sample, per_contig));
    }

    let priors_text = std::fs::read_to_string(&priors)
        .map_err(|error| Thrown::user(format!("{priors}: {error}")))?;
    let interval_text: String = dictionary_text.clone()
        + "CONTIG\tSTART\tEND\n"
        + &intervals
            .iter()
            .map(|(c, s, e)| format!("{c}\t{s}\t{e}\n"))
            .collect::<String>();

    let settings = Settings {
        mapping_error_rate: double(parser, "mapping-error-rate", 0.3),
        psi_s_scale: double(parser, "sample-psi-scale", 0.0001),
        mean_bias_sd: double(parser, "mean-bias-standard-deviation", 1.0),
        psi_j_scale: double(parser, "global-psi-scale", 0.001),
        learning_rate: double(parser, "learning-rate", 0.05),
        adamax_beta1: double(parser, "adamax-beta-1", 0.9),
        adamax_beta2: double(parser, "adamax-beta-2", 0.999),
        log_emission_samples_per_round: integer(parser, "log-emission-samples-per-round", 2000),
        log_emission_sampling_rounds: integer(parser, "log-emission-sampling-rounds", 100),
        log_emission_sampling_median_rel_error: double(
            parser,
            "log-emission-sampling-median-rel-error",
            0.0005,
        ),
        max_advi_iter_first_epoch: integer(parser, "max-advi-iter-first-epoch", 1000),
        max_advi_iter_subsequent_epochs: integer(parser, "max-advi-iter-subsequent-epochs", 1000),
        min_training_epochs: integer(parser, "min-training-epochs", 20),
        max_training_epochs: integer(parser, "max-training-epochs", 100),
        initial_temperature: double(parser, "initial-temperature", 2.0),
        num_thermal_advi_iters: integer(parser, "num-thermal-advi-iters", 5000),
        convergence_snr_averaging_window: integer(parser, "convergence-snr-averaging-window", 5000),
        convergence_snr_trigger_threshold: double(parser, "convergence-snr-trigger-threshold", 0.1),
        convergence_snr_countdown_window: integer(parser, "convergence-snr-countdown-window", 10),
        max_calling_iters: integer(parser, "max-calling-iters", 1),
        caller_update_convergence_threshold: double(
            parser,
            "caller-update-convergence-threshold",
            0.001,
        ),
        caller_external_admixing_rate: double(parser, "caller-external-admixing-rate", 0.75),
        disable_annealing: scalar(parser, "disable-annealing").as_deref() == Some("true"),
        disable_caller: scalar(parser, "disable-caller").as_deref() == Some("true"),
        disable_sampler: scalar(parser, "disable-sampler").as_deref() == Some("true"),
    };
    let cohort = Cohort {
        contigs,
        intervals_per_contig: {
            let mut t = Vec::new();
            for (contig, _, _) in &intervals {
                let j = cohort_index(&t, contig);
                if j == t.len() {
                    t.push((contig.clone(), 0u64));
                }
                t[j].1 += 1;
            }
            t.into_iter().map(|(_, n)| n).collect()
        },
        samples,
        priors_text: priors_text.clone(),
    };
    let files = gcnv_ploidy::run(&cohort, &settings).map_err(|why| {
        Thrown::non_user(
            PORT_LIMITATION,
            format!("{why}. This message is the port's own and not GATK's."),
        )
    })?;

    let model = format!("{output}/{prefix}-model");
    let calls = format!("{output}/{prefix}-calls");
    let write = |path: String, text: &str| -> Result<(), Thrown> {
        if let Some(parent) = std::path::Path::new(&path).parent() {
            std::fs::create_dir_all(parent)
                .map_err(|error| Thrown::user(format!("{}: {error}", parent.display())))?;
        }
        std::fs::write(&path, text).map_err(|error| Thrown::user(format!("{path}: {error}")))
    };
    write(format!("{model}/interval_list.tsv"), &interval_text)?;
    write(format!("{model}/contig_ploidy_prior.tsv"), &priors_text)?;
    for (name, text) in &files.model {
        write(format!("{model}/{name}"), text)?;
    }
    for (name, text) in &files.calls {
        write(format!("{calls}/{name}"), text)?;
    }
    Ok(None)
}

fn cohort_index(t: &[(String, u64)], contig: &str) -> usize {
    t.iter().position(|(c, _)| c == contig).unwrap_or(t.len())
}
