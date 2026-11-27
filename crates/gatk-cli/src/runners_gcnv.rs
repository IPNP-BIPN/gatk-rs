//! `DetermineGermlineContigPloidy.doWork` and `GermlineCNVCaller.doWork`: the Java around
//! gcnvkernel. The arguments are checked in the reference's order, the intervals resolved, the
//! read counts subset and written to temporary files, and GATK's own script is run on them with the
//! argument list the Java composes. Everything the tools write is the script's.

use super::*;

const COHORT_PLOIDY: PythonResource = PythonResource {
    name: "cohort_determine_ploidy_and_depth.py",
    text: include_str!("python/cohort_determine_ploidy_and_depth.py"),
};
const CASE_PLOIDY: PythonResource = PythonResource {
    name: "case_determine_ploidy_and_depth.py",
    text: include_str!("python/case_determine_ploidy_and_depth.py"),
};
const COHORT_CALLING: PythonResource = PythonResource {
    name: "cohort_denoising_calling.py",
    text: include_str!("python/cohort_denoising_calling.py"),
};
const CASE_CALLING: PythonResource = PythonResource {
    name: "case_denoising_calling.py",
    text: include_str!("python/case_denoising_calling.py"),
};

/// `GermlineCNVCaller.STARTING_SEED` and gcnvkernel's exit status for an inference that diverged.
const STARTING_SEED: i64 = 1984;
const DIVERGED_INFERENCE_EXIT_CODE: i32 = 239;

/// `String.format("%e", value)`.
fn percent_e(value: f64) -> String {
    if !value.is_finite() {
        return gatk_engine::pair_hmm_likelihood_engine::java_scientific(value);
    }
    gatk_engine::java_format::format_scientific_upper(value, 6).replace('E', "e")
}

fn double(parser: &Parser, name: &str, default: f64) -> f64 {
    double_or(parser, name, default)
}

fn integer(parser: &Parser, name: &str, default: i64) -> i64 {
    scalar(parser, name)
        .and_then(|text| text.parse().ok())
        .unwrap_or(default)
}

/// A boolean argument whose default may be true.
fn boolean(parser: &Parser, name: &str, default: bool) -> bool {
    scalar(parser, name).map_or(default, |text| text == "true")
}

/// `ParamUtils.isPositive` and `isPositiveOrZero`.
fn positive(value: f64, or_zero: bool, message: String) -> Result<(), Thrown> {
    if value > 0.0 || (or_zero && value == 0.0) {
        return Ok(());
    }
    Err(illegal(message))
}

/// `HybridADVIArgumentCollection`'s arguments, in `generatePythonArguments`' order: python's name,
/// the command line's, and whether it goes through `%e` rather than `%d`.
const ADVI: [(&str, &str, bool); 19] = [
    ("learning_rate", "learning-rate", true),
    ("adamax_beta1", "adamax-beta-1", true),
    ("adamax_beta2", "adamax-beta-2", true),
    (
        "log_emission_samples_per_round",
        "log-emission-samples-per-round",
        false,
    ),
    (
        "log_emission_sampling_rounds",
        "log-emission-sampling-rounds",
        false,
    ),
    (
        "log_emission_sampling_median_rel_error",
        "log-emission-sampling-median-rel-error",
        true,
    ),
    (
        "max_advi_iter_first_epoch",
        "max-advi-iter-first-epoch",
        false,
    ),
    (
        "max_advi_iter_subsequent_epochs",
        "max-advi-iter-subsequent-epochs",
        false,
    ),
    ("min_training_epochs", "min-training-epochs", false),
    ("max_training_epochs", "max-training-epochs", false),
    ("initial_temperature", "initial-temperature", true),
    ("num_thermal_advi_iters", "num-thermal-advi-iters", false),
    (
        "convergence_snr_averaging_window",
        "convergence-snr-averaging-window",
        false,
    ),
    (
        "convergence_snr_trigger_threshold",
        "convergence-snr-trigger-threshold",
        true,
    ),
    (
        "convergence_snr_countdown_window",
        "convergence-snr-countdown-window",
        false,
    ),
    ("max_calling_iters", "max-calling-iters", false),
    (
        "caller_update_convergence_threshold",
        "caller-update-convergence-threshold",
        true,
    ),
    (
        "caller_internal_admixing_rate",
        "caller-internal-admixing-rate",
        true,
    ),
    (
        "caller_external_admixing_rate",
        "caller-external-admixing-rate",
        true,
    ),
];

/// `GermlineContigPloidyHybridADVIArgumentCollection`'s defaults.
const PLOIDY_ADVI: [f64; 19] = [
    0.05, 0.9, 0.999, 2000.0, 100.0, 0.0005, 1000.0, 1000.0, 20.0, 100.0, 2.0, 5000.0, 5000.0, 0.1,
    10.0, 1.0, 0.001, 0.75, 0.75,
];

/// `GermlineCNVHybridADVIArgumentCollection`'s defaults.
const CALLER_ADVI: [f64; 19] = [
    0.01, 0.9, 0.99, 50.0, 10.0, 0.005, 5000.0, 200.0, 10.0, 50.0, 1.5, 2500.0, 500.0, 0.1, 10.0,
    10.0, 0.001, 0.75, 1.0,
];

/// `HybridADVIArgumentCollection.generatePythonArguments`. Its `validate` checks nothing.
fn advi_arguments(parser: &Parser, defaults: &[f64; 19]) -> Vec<String> {
    let mut arguments: Vec<String> = ADVI
        .iter()
        .zip(defaults)
        .map(|((python, name, scientific), &default)| {
            if *scientific {
                format!("--{python}={}", percent_e(double(parser, name, default)))
            } else {
                format!("--{python}={}", integer(parser, name, default as i64))
            }
        })
        .collect();
    for (python, name) in [
        ("disable_caller", "disable-caller"),
        ("disable_sampler", "disable-sampler"),
        ("disable_annealing", "disable-annealing"),
    ] {
        arguments.push(format!("--{python}={}", boolean(parser, name, false)));
    }
    arguments
}

/// `IOUtils.assertFileIsReadable(IOUtils.getPath(path))`.
fn assert_file_is_readable(path: &str) -> Result<(), Thrown> {
    let refuse = |why: &str| {
        Thrown::user(format!(
            "Couldn't read file file://{}. Error was: {why}",
            java_absolute_path(path)
        ))
    };
    let Ok(metadata) = std::fs::metadata(path) else {
        return Err(refuse("It doesn't exist."));
    };
    if !metadata.is_file() {
        return Err(refuse("It isn't a regular file"));
    }
    if std::fs::File::open(path).is_err() {
        return Err(refuse("It is not readable, check the file permissions"));
    }
    Ok(())
}

/// `CopyNumberArgumentValidationUtils.validateInputs(File...)`: a file must be readable, and
/// anything else is let through.
fn validate_input(path: Option<&str>) -> Result<(), Thrown> {
    match path {
        Some(path) if std::path::Path::new(path).is_file() => tables::can_read_file(path),
        _ => Ok(()),
    }
}

/// `validateArguments`' common tail: no duplicates, the inputs, the prefix and the output.
fn validate_paths(
    parser: &Parser,
    duplicates: &str,
    files: &[Option<String>],
) -> Result<(Vec<String>, String, String), Thrown> {
    let inputs = arguments(parser, "input");
    let mut seen = std::collections::HashSet::new();
    if !inputs.iter().all(|path| seen.insert(path.clone())) {
        return Err(illegal(duplicates));
    }
    for path in &inputs {
        assert_file_is_readable(path)?;
    }
    for path in files {
        validate_input(path.as_deref())?;
    }
    let prefix = require(parser, "output-prefix");
    non_empty(&prefix)?;
    let output = require(parser, "output");
    prepare_output_directory(&output)?;
    Ok((inputs, prefix, output))
}

/// A `SimpleCountCollection` TSV.
struct ReadCounts {
    sample: String,
    dictionary: Vec<SequenceRecord>,
    records: Vec<(Span, String)>,
}

/// `SimpleCountCollection.read`, for the TSV the collection writes.
fn read_counts(path: &str) -> Result<ReadCounts, Thrown> {
    if path.ends_with(".hdf5") || path.ends_with(".h5") {
        return Err(port_limit("Reading an HDF5 read-counts file"));
    }
    let (dictionary, spans) = tables::read_interval_spans(path)?;
    let (sample, columns, rows) = tables::read_sample_rows(path)?;
    let count = columns.iter().position(|c| c == "COUNT").ok_or_else(|| {
        Thrown::user(
            "Bad input: Bad header in file.  Not all mandatory columns are present.  Missing: COUNT",
        )
    })?;
    let records = spans
        .into_iter()
        .zip(rows)
        .map(|(span, row)| (span, row[count].clone()))
        .collect();
    Ok(ReadCounts {
        sample,
        dictionary,
        records,
    })
}

/// A `SimpleIntervalCollection`'s text: the dictionary, the columns, the intervals.
fn interval_text(dictionary: &[SequenceRecord], intervals: &[Span]) -> String {
    let mut text = SamHeader {
        sequences: dictionary.to_vec(),
        ..SamHeader::default()
    }
    .encode();
    text.push_str("CONTIG\tSTART\tEND\n");
    for span in intervals {
        text.push_str(&format!("{}\t{}\t{}\n", span.contig, span.start, span.end));
    }
    text
}

/// `CopyNumberArgumentValidationUtils.resolveIntervals`: `-L` against the first counts file's
/// dictionary, or that file's own intervals.
fn resolve_intervals(
    parser: &Parser,
    first: &str,
) -> Result<(Vec<SequenceRecord>, Vec<Span>), Thrown> {
    assert_file_is_readable(first)?;
    let specified = intervals_specified(parser);
    if specified {
        validate_copy_number_intervals(parser)?;
    }
    let counts = read_counts(first)?;
    let header = SamHeader {
        sequences: counts.dictionary.clone(),
        ..SamHeader::default()
    };
    let intervals = match interval_arguments(parser, &header)? {
        Some(parameters) => parameters
            .intervals
            .iter()
            .map(|i| Span {
                contig: i.contig.clone(),
                start: i.start,
                end: i.end,
            })
            .collect(),
        None => counts.records.into_iter().map(|(span, _)| span).collect(),
    };
    Ok((counts.dictionary, intervals))
}

fn intervals_specified(parser: &Parser) -> bool {
    !arguments(parser, "intervals").is_empty() || !arguments(parser, "exclude-intervals").is_empty()
}

fn key(span: &Span) -> (String, i32, i32) {
    (span.contig.clone(), span.start, span.end)
}

/// `streamOfSubsettedAndValidatedReadCounts`: each file's counts at the intervals, in its order.
fn subset_read_counts(paths: &[String], intervals: &[Span]) -> Result<Vec<ReadCounts>, Thrown> {
    let wanted: std::collections::HashSet<(String, i32, i32)> = intervals.iter().map(key).collect();
    paths
        .iter()
        .map(|path| {
            let mut counts = read_counts(path)?;
            counts
                .records
                .retain(|(span, _)| wanted.contains(&key(span)));
            if counts.records.len() != wanted.len() {
                return Err(illegal(format!(
                    "Intervals for read-count file {path} do not contain all specified intervals."
                )));
            }
            Ok(counts)
        })
        .collect()
}

/// Where the tool's outputs go: `addTrailingSlashIfNecessary(outputDir.getAbsolutePath())`.
fn output_path(output: &str, prefix: &str, suffix: &str) -> String {
    let directory = java_absolute_path(output);
    let slash = if directory.ends_with('/') { "" } else { "/" };
    canonical_path(&format!("{directory}{slash}{prefix}{suffix}"))
}

/// `DetermineGermlineContigPloidy.doWork`.
pub(super) fn ploidy(parser: &Parser) -> Outcome {
    // `validateArguments`.
    let mapping_error_rate = double(parser, "mapping-error-rate", 0.3);
    let sample_psi_scale = double(parser, "sample-psi-scale", 0.0001);
    let mean_bias_sd = double(parser, "mean-bias-standard-deviation", 1.0);
    let global_psi_scale = double(parser, "global-psi-scale", 0.001);
    for (value, what, name) in [
        (
            mean_bias_sd,
            "Prior standard deviation of the contig-level mean coverage bias",
            "mean-bias-standard-deviation",
        ),
        (
            mapping_error_rate,
            "Typical mapping error rate",
            "mapping-error-rate",
        ),
        (
            global_psi_scale,
            "Prior scale of contig coverage unexplained variance",
            "global-psi-scale",
        ),
        (
            sample_psi_scale,
            "Prior scale of the sample-specific correction to the coverage unexplained variance",
            "sample-psi-scale",
        ),
    ] {
        positive(value, false, format!("{what} ({name}) must be positive."))?;
    }
    let priors = argument(parser, "contig-ploidy-priors");
    let model = argument(parser, "model");
    let (inputs, prefix, output) = validate_paths(
        parser,
        "List of input read-count files cannot contain duplicates.",
        &[priors.clone(), model.clone()],
    )?;

    // `setModeAndResolveIntervals`.
    let (dictionary, intervals) = match &model {
        Some(model) => {
            if priors.is_some() {
                return Err(Thrown::user(
                    "Bad input: Invalid combination of inputs: Running in case mode, but \
                     contig-ploidy priors were provided.",
                ));
            }
            if intervals_specified(parser) {
                return Err(Thrown::user(
                    "Bad input: Invalid combination of inputs: Running in CASE mode, but \
                     intervals were provided.",
                ));
            }
            let file = format!("{model}/interval_list.tsv");
            validate_input(Some(&file))?;
            tables::read_interval_spans(&file)?
        }
        None => {
            if inputs.len() < 2 {
                return Err(illegal(
                    "At least two samples must be provided in COHORT mode.",
                ));
            }
            if priors.is_none() {
                return Err(Thrown::user(
                    "Bad input: Contig-ploidy priors must be provided in cohort mode.",
                ));
            }
            resolve_intervals(parser, &inputs[0])?
        }
    };

    let tmp_dir = java_tmp_dir(parser);
    let intervals_file = TempFile::new(&tmp_dir, "intervals", ".tsv")?;
    write_file(
        &intervals_file.0,
        interval_text(&dictionary, &intervals).as_bytes(),
    )?;

    // `writeSamplesByCoveragePerContig`: each file subset, its counts summed per contig.
    let mut contigs: Vec<&str> = Vec::new();
    for span in &intervals {
        if !contigs.contains(&span.contig.as_str()) {
            contigs.push(&span.contig);
        }
    }
    let mut coverage = interval_text(&dictionary, &[]).replace("CONTIG\tSTART\tEND\n", "");
    coverage.push_str(&format!("SAMPLE_NAME\t{}\n", contigs.join("\t")));
    for counts in subset_read_counts(&inputs, &intervals)? {
        let mut sums = vec![0i32; contigs.len()];
        for (span, count) in &counts.records {
            let at = contigs
                .iter()
                .position(|c| *c == span.contig)
                .expect("a contig");
            let count: i32 = count.parse().map_err(|_| {
                Thrown::user(format!(
                    "Bad input: Expected an int value but found {count}"
                ))
            })?;
            // `Collectors.summingInt`.
            sums[at] = sums[at].wrapping_add(count);
        }
        let sums: Vec<String> = sums.iter().map(i32::to_string).collect();
        coverage.push_str(&format!("{}\t{}\n", counts.sample, sums.join("\t")));
    }
    let coverage_file = TempFile::new(&tmp_dir, "samples-by-coverage-per-contig", ".tsv")?;
    write_file(&coverage_file.0, coverage.as_bytes())?;

    // `executeDeterminePloidyAndDepthPythonScript`.
    let mut args = vec![
        format!(
            "--sample_coverage_metadata={}",
            canonical_path(&coverage_file.0)
        ),
        format!(
            "--output_calls_path={}",
            output_path(&output, &prefix, "-calls")
        ),
        format!("--mapping_error_rate={}", percent_e(mapping_error_rate)),
        format!("--psi_s_scale={}", percent_e(sample_psi_scale)),
    ];
    if model.is_none() {
        args.push(format!("--mean_bias_sd={}", percent_e(mean_bias_sd)));
        args.push(format!("--psi_j_scale={}", percent_e(global_psi_scale)));
    }
    args.extend(advi_arguments(parser, &PLOIDY_ADVI));
    let script = match &model {
        None => {
            args.push(format!(
                "--interval_list={}",
                canonical_path(&intervals_file.0)
            ));
            args.push(format!(
                "--contig_ploidy_prior_table={}",
                canonical_path(priors.as_deref().unwrap_or_default())
            ));
            args.push(format!(
                "--output_model_path={}",
                output_path(&output, &prefix, "-model")
            ));
            &COHORT_PLOIDY
        }
        Some(model) => {
            args.push(format!("--input_model_path={}", canonical_path(model)));
            &CASE_PLOIDY
        }
    };
    let script = python_resource(&tmp_dir, script)?;
    run_python_script(&script.0, &args)?;
    Ok(None)
}

/// `GermlineCallingArgumentCollection` and `GermlineDenoisingModelArgumentCollection`'s
/// `validate`: the arguments a model fixes refused in CASE mode, then each bound.
fn validate_caller_collections(parser: &Parser, case: bool) -> Result<(), Thrown> {
    let set = |name: &str| {
        parser
            .definitions()
            .iter()
            .any(|definition| definition.long_name() == name && definition.has_been_set())
    };
    let hidden = |names: &[&str]| -> Result<(), Thrown> {
        for name in names.iter().filter(|_| case) {
            if set(name) {
                return Err(illegal(format!(
                    "Argument '--{name}' cannot be set in the CASE mode."
                )));
            }
        }
        Ok(())
    };
    hidden(&["p-active", "class-coherence-length"])?;
    for (name, default, what) in [
        (
            "cnv-coherence-length",
            10000.0,
            "Coherence length of CNV events",
        ),
        (
            "class-coherence-length",
            10000.0,
            "Coherence length of CNV class domains",
        ),
        ("max-copy-number", 5.0, "Highest allowed copy-number"),
    ] {
        positive(
            double(parser, name, default),
            false,
            format!("{what} ({name}) must be positive."),
        )?;
    }
    hidden(&[
        "max-bias-factors",
        "interval-psi-scale",
        "log-mean-bias-standard-deviation",
        "init-ard-rel-unexplained-variance",
        "enable-bias-factors",
        "num-gc-bins",
        "gc-curve-standard-deviation",
    ])?;
    for (name, default, what) in [
        ("max-bias-factors", 6.0, "Maximum number of bias factors"),
        ("mapping-error-rate", 0.01, "Mapping error rate"),
        ("interval-psi-scale", 0.01, "Typical scale of interval-specific unexplained variance"),
        ("sample-psi-scale", 0.01, "Typical scale of sample-specific correction to the unexplained variance"),
        ("depth-correction-tau", 10000.0, "Precision of read depth pinning to its global value"),
        ("log-mean-bias-standard-deviation", 0.1, "Standard deviation of log mean bias"),
        ("init-ard-rel-unexplained-variance", 0.1, "Initial value of ARD prior precision relative to the scale of interval-specific unexplained variance"),
        ("gc-curve-standard-deviation", 1.0, "Prior standard deviation of the GC curve from flat"),
    ] {
        positive(double(parser, name, default), false, format!("{what} ({name}) must be positive."))?;
    }
    positive(
        double(parser, "active-class-padding-hybrid-mode", 50000.0),
        true,
        "CNV-active interval padding in HYBRID copy-number posterior expectation mode \
         (active-class-padding-hybrid-mode) must be non-negative."
            .to_string(),
    )
}

/// `validateAnnotatedIntervalsSubset` and the collection's `write`: the annotated intervals at the
/// resolved ones, in the annotated file's order, under the counts' dictionary.
fn annotated_subset(
    path: &str,
    dictionary: &[SequenceRecord],
    intervals: &[Span],
) -> Result<String, Thrown> {
    tables::can_read_file(path)?;
    let spans = tables::read_interval_spans(path)?.1;
    let text = std::fs::read_to_string(path).map_err(|error| Thrown::user(error.to_string()))?;
    let mut lines = text.lines().filter(|line| !line.starts_with('@'));
    let columns: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    // `CopyNumberAnnotations.ANNOTATIONS`, in that order, each read as a double.
    let keys: Vec<(&str, usize)> = ["GC_CONTENT", "MAPPABILITY", "SEGMENTAL_DUPLICATION_CONTENT"]
        .into_iter()
        .filter_map(|name| columns.iter().position(|c| *c == name).map(|at| (name, at)))
        .collect();
    let wanted: std::collections::HashSet<(String, i32, i32)> = intervals.iter().map(key).collect();
    let mut out = interval_text(dictionary, &[]).replace("CONTIG\tSTART\tEND\n", "");
    let names: Vec<&str> = keys.iter().map(|(name, _)| *name).collect();
    out.push_str(&format!("CONTIG\tSTART\tEND\t{}\n", names.join("\t")).replace("\t\n", "\n"));
    let mut kept = 0;
    for (span, line) in spans.iter().zip(lines) {
        if !wanted.contains(&key(span)) {
            continue;
        }
        kept += 1;
        let values: Vec<&str> = line.split('\t').collect();
        let mut row = format!("{}\t{}\t{}", span.contig, span.start, span.end);
        for (_, at) in &keys {
            let value = gatk_tools::variant_recalibrator_model::java_parse_double(values[*at])
                .ok_or_else(|| {
                    Thrown::user(format!(
                        "Bad input: Expected a double value but found {}",
                        values[*at]
                    ))
                })?;
            row.push('\t');
            row.push_str(&gatk_engine::java_format::format_decimals(value, 6));
        }
        out.push_str(&row);
        out.push('\n');
    }
    if kept != wanted.len() {
        return Err(illegal(
            "Annotated intervals do not contain all specified intervals.",
        ));
    }
    Ok(out)
}

/// `GermlineCNVCaller.doWork`.
pub(super) fn caller(parser: &Parser) -> Outcome {
    let case = scalar(parser, "run-mode").as_deref() == Some("CASE");
    // `validateArguments`.
    validate_caller_collections(parser, case)?;
    let ploidy_calls = argument(parser, "contig-ploidy-calls");
    let model = argument(parser, "model");
    let annotated = argument(parser, "annotated-intervals");
    let (inputs, prefix, output) = validate_paths(
        parser,
        "List of input read-count file paths cannot contain duplicates.",
        &[ploidy_calls.clone(), model.clone(), annotated.clone()],
    )?;

    // `resolveIntervals`.
    let tmp_dir = java_tmp_dir(parser);
    // The temporary intervals file lives as long as the run: `deleteOnExit`.
    let (intervals, modeling_intervals, _intervals_file) = match &model {
        Some(model) => {
            let file = format!("{model}/interval_list.tsv");
            validate_input(Some(&file))?;
            (tables::read_interval_spans(&file)?.1, file, None)
        }
        None => {
            let (dictionary, intervals) = resolve_intervals(parser, &inputs[0])?;
            let file = TempFile::new(&tmp_dir, "intervals", ".tsv")?;
            let text = match &annotated {
                Some(path) => annotated_subset(path, &dictionary, &intervals)?,
                None => interval_text(&dictionary, &intervals),
            };
            write_file(&file.0, text.as_bytes())?;
            (intervals, file.0.clone(), Some(file))
        }
    };
    if !case {
        if inputs.len() < 2 {
            return Err(illegal(
                "At least two samples must be provided in COHORT mode.",
            ));
        }
    } else {
        if model.is_none() {
            return Err(illegal(
                "An input denoising-model directory must be provided in CASE mode.",
            ));
        }
        if intervals_specified(parser) {
            return Err(Thrown::user(
                "Bad input: Invalid combination of inputs: Running in CASE mode, but intervals \
                 were provided.",
            ));
        }
        if annotated.is_some() {
            return Err(Thrown::user(
                "Bad input: Invalid combination of inputs: Running in CASE mode,but annotated \
                 intervals were provided.",
            ));
        }
    }

    // `writeIntervalSubsetReadCountFiles`.
    let mut count_files = Vec::new();
    for counts in subset_read_counts(&inputs, &intervals)? {
        let file = TempFile::new(&tmp_dir, &format!("{}.rc", counts.sample), ".tsv")?;
        let mut text = tables::header(&counts.sample, Some(&counts.dictionary));
        text.push_str("CONTIG\tSTART\tEND\tCOUNT\n");
        for (span, count) in &counts.records {
            text.push_str(&format!(
                "{}\t{}\t{}\t{count}\n",
                span.contig, span.start, span.end
            ));
        }
        write_file(&file.0, text.as_bytes())?;
        count_files.push(file);
    }

    // `composePythonArguments`.
    let compose = |seed: i64| -> Vec<String> {
        let mut args = vec![
            format!(
                "--ploidy_calls_path={}",
                canonical_path(ploidy_calls.as_deref().unwrap_or_default())
            ),
            format!(
                "--output_calls_path={}",
                output_path(&output, &prefix, "-calls")
            ),
            format!(
                "--output_tracking_path={}",
                output_path(&output, &prefix, "-tracking")
            ),
        ];
        if let Some(model) = &model {
            args.push(format!("--input_model_path={}", canonical_path(model)));
        }
        args.push(format!("--random_seed={seed}"));
        if !case {
            args.push(format!(
                "--modeling_interval_list={}",
                canonical_path(&modeling_intervals)
            ));
            args.push(format!(
                "--output_model_path={}",
                output_path(&output, &prefix, "-model")
            ));
            let gc = if annotated.is_some() { "True" } else { "False" };
            args.push(format!("--enable_explicit_gc_bias_modeling={gc}"));
        }
        args.push("--read_count_tsv_files".to_string());
        args.extend(count_files.iter().map(|file| canonical_path(&file.0)));
        // `GermlineDenoisingModelArgumentCollection.generatePythonArguments`.
        let mode = match scalar(parser, "copy-number-posterior-expectation-mode").as_deref() {
            Some("MAP") => "map",
            Some("EXACT") => "exact",
            _ => "hybrid",
        };
        args.extend([
            format!(
                "--psi_s_scale={}",
                percent_e(double(parser, "sample-psi-scale", 0.01))
            ),
            format!(
                "--mapping_error_rate={}",
                percent_e(double(parser, "mapping-error-rate", 0.01))
            ),
            format!(
                "--depth_correction_tau={}",
                percent_e(double(parser, "depth-correction-tau", 10000.0))
            ),
            format!("--q_c_expectation_mode={mode}"),
            format!(
                "--num_samples_copy_ratio_approx={}",
                integer(parser, "num-samples-copy-ratio-approx", 200)
            ),
        ]);
        if !case {
            let bias = if boolean(parser, "enable-bias-factors", true) {
                "True"
            } else {
                "False"
            };
            args.extend([
                format!(
                    "--max_bias_factors={}",
                    integer(parser, "max-bias-factors", 6)
                ),
                format!(
                    "--psi_t_scale={}",
                    percent_e(double(parser, "interval-psi-scale", 0.01))
                ),
                format!(
                    "--log_mean_bias_std={}",
                    percent_e(double(parser, "log-mean-bias-standard-deviation", 0.1))
                ),
                format!(
                    "--init_ard_rel_unexplained_variance={}",
                    percent_e(double(parser, "init-ard-rel-unexplained-variance", 0.1))
                ),
                format!("--num_gc_bins={}", integer(parser, "num-gc-bins", 20)),
                format!(
                    "--gc_curve_sd={}",
                    percent_e(double(parser, "gc-curve-standard-deviation", 1.0))
                ),
                format!(
                    "--active_class_padding_hybrid_mode={}",
                    integer(parser, "active-class-padding-hybrid-mode", 50000)
                ),
                format!("--enable_bias_factors={bias}"),
                "--disable_bias_factors_in_active_class=False".to_string(),
            ]);
        }
        // `GermlineCallingArgumentCollection.generatePythonArguments`.
        args.extend([
            format!("--p_alt={}", percent_e(double(parser, "p-alt", 5e-4))),
            format!(
                "--cnv_coherence_length={}",
                percent_e(double(parser, "cnv-coherence-length", 10000.0))
            ),
            format!(
                "--max_copy_number={}",
                integer(parser, "max-copy-number", 5)
            ),
        ]);
        if !case {
            args.extend([
                format!(
                    "--p_active={}",
                    gatk_engine::java_format::format_decimals(double(parser, "p-active", 0.1), 6)
                ),
                format!(
                    "--class_coherence_length={}",
                    gatk_engine::java_format::format_decimals(
                        double(parser, "class-coherence-length", 10000.0),
                        6
                    )
                ),
            ]);
        }
        args.extend(advi_arguments(parser, &CALLER_ADVI));
        args
    };

    let resource = if case { &CASE_CALLING } else { &COHORT_CALLING };
    let run = |seed: i64| -> Result<(std::process::Output, String), Thrown> {
        let script = python_resource(&tmp_dir, resource)?;
        python_output(&script.0, &compose(seed))
    };
    let (mut output, mut command_line) = run(STARTING_SEED)?;
    // Restarted once with another seed where the inference diverged.
    if output.status.code() == Some(DIVERGED_INFERENCE_EXIT_CODE) {
        let seed = gatk_engine::java_random::JavaRandom::new(STARTING_SEED).next_int();
        (output, command_line) = run(i64::from(seed))?;
    }
    script_failure(PYTHON_SCRIPT_EXCEPTION, &command_line, &output)?;
    Ok(None)
}
