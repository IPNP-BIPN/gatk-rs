//! The tools that stop at an interpreter: two R plots, the gCNV tools, `NVScoreVariants` and the
//! scalable VQSR's two.
//!
//! The pinned image has neither `Rscript` nor `python`, so what can be measured is everything the
//! reference does on the way there and the user error it ends on. That much is ported, in the
//! reference's order. On a machine that has the interpreter the reference would go on to run a
//! script this port does not carry, and the port says so in its own words rather than writing
//! anything that could pass for the answer.
//!
//! * **the plots** check their arguments, read and cross-check every input, and only then does
//!   `RScriptExecutor.exec` look for its executable;
//! * **`DetermineGermlineContigPloidy`, `GermlineCNVCaller`, `NVScoreVariants`** import a package
//!   in `onStartup`, before `doWork`, so only the parser can refuse first;
//! * **`PostprocessGermlineCNVCalls`** runs `GATKTool.onStartup` ahead of that import, and `-L`
//!   there resolves against the dictionary of the shards' interval lists;
//! * **`TrainVariantAnnotationsModel`** reads the annotation names only to compare them with the
//!   unlabeled file's, and `JAVA_BGMM` and `PYTHON_SCRIPT` read every matrix before they fail;
//! * **`ScoreVariantAnnotations`** is `ExtractVariantAnnotations`' walker with a backend after the
//!   header, a pass and an annotations file before its first script, and a VCF when the pass
//!   extracted nothing and no script is needed.
use super::extract_variant_annotations_runner::{
    labelled_walker_start, matrix, metadata, reduced, validate_output_files, LabelledStart,
};
use super::*;
use gatk_tools::extract_variant_annotations as extract;
use gatk_tools::hdf5_reader::Hdf5Reader;
use gatk_tools::model_segments_models::{illegal, Span};
use gatk_tools::model_segments_tool as tables;
use htsjdk_bam::header::SequenceRecord;

const NOT_IMPLEMENTED: &str = "org.apache.commons.lang3.NotImplementedException";
const BGMM_NOT_IMPLEMENTED: &str = "BGMM module will be implemented in separate PR.";
const GATK_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.GATKException";

/// `RuntimeUtils.which`: a file of that name in any directory of the path.
fn on_the_path(executable: &str) -> bool {
    std::env::var_os("PATH").is_some_and(|path| {
        std::env::split_paths(&path).any(|directory| directory.join(executable).exists())
    })
}

/// `ScriptExecutor.executableMissing`, or the port's limit where the executable exists.
fn stop_at(executable: &str, what_is_not_ported: &str) -> Thrown {
    if on_the_path(executable) {
        return Thrown::non_user(
            PORT_LIMITATION,
            format!("{what_is_not_ported} is not ported. This message is the port's own and not GATK's."),
        );
    }
    Thrown::user(format!(
        "Unable to execute {executable} command: Please add the {executable} directory to your \
         environment ${{PATH}}"
    ))
}

/// `PythonScriptExecutor.checkPythonEnvironmentForPackage`, whose constructor is what fails.
fn import_stop(package: &str) -> Thrown {
    stop_at("python", &format!("A Python environment with {package}"))
}

/// `IOUtils.canReadFile`, which `Utils.nonNull` guards.
fn can_read_file(path: Option<&str>) -> Result<(), Thrown> {
    let Some(path) = path else {
        return Err(illegal("Unexpected null file reference."));
    };
    tables::can_read_file(path)
}

/// `CopyNumberArgumentValidationUtils.validateAndPrepareOutputDirectories`, one level only.
fn prepare_output_directory(path: &str) -> Result<(), Thrown> {
    if std::path::Path::new(path).exists() {
        return Ok(());
    }
    std::fs::create_dir(path).map_err(|_| {
        Thrown::user(format!(
            "Couldn't write file {} because : The output directory does not exist and could not \
             be created.",
            java_absolute_path(path)
        ))
    })
}

/// `Utils.nonEmpty(String)`.
fn non_empty(text: &str) -> Result<(), Thrown> {
    if text.is_empty() {
        return Err(illegal(
            "The string is empty: string must not be null or empty",
        ));
    }
    Ok(())
}

/// `ReferenceUtils.loadFastaDictionary(File)`.
fn load_fasta_dictionary(path: &str) -> Result<Vec<SequenceRecord>, Thrown> {
    let absolute = java_absolute_path(path);
    let text = std::fs::read_to_string(path)
        .map_err(|_| Thrown::non_user("java.lang.RuntimeException", absolute.clone()))?;
    let header = htsjdk_bam::reader::parse_header_text(&text);
    if header.sequences.is_empty() {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile",
            message: Some(format!(
                "Unknown file is malformed: Could not read sequence dictionary from given fasta \
                 file {absolute}"
            )),
        });
    }
    Ok(header.sequences)
}

/// `PlottingUtils.getContigLengthMap`.
fn contig_length_map(
    dictionary: &[SequenceRecord],
    minimum: i32,
) -> Result<Vec<(String, i32)>, Thrown> {
    if dictionary
        .iter()
        .any(|sequence| sequence.name.contains("CONTIG_DELIMITER"))
    {
        return Err(illegal("Contig names cannot contain \"CONTIG_DELIMITER\"."));
    }
    let map: Vec<(String, i32)> = dictionary
        .iter()
        .filter(|sequence| sequence.length >= minimum)
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();
    if map.is_empty() {
        return Err(illegal(
            "There must be at least one contig above the threshold length in the sequence \
             dictionary.",
        ));
    }
    Ok(map)
}

/// `PlottingUtils.validateContigs`: a position past the plotted contig's end. A contig the
/// dictionary leaves out is only the reference's warning.
fn validate_contigs(map: &[(String, i32)], spans: &[Span], file: &str) -> Result<(), Thrown> {
    let beyond = spans.iter().any(|span| {
        map.iter()
            .any(|(contig, length)| *contig == span.contig && span.end > *length)
    });
    if beyond {
        return Err(illegal(format!(
            "Position present in the file {file} exceeds contig length in the sequence dictionary."
        )));
    }
    Ok(())
}

fn minimum_contig_length(parser: &Parser) -> i32 {
    number_or(parser, "minimum-contig-length", 1_000_000)
}

fn require(parser: &Parser, name: &str) -> String {
    argument(parser, name).expect("a required argument")
}

/// `PlotDenoisedCopyRatios.doWork`.
pub fn plot_denoised_copy_ratios(parser: &Parser) -> Outcome {
    let standardized = require(parser, "standardized-copy-ratios");
    let denoised = require(parser, "denoised-copy-ratios");
    let dictionary = require(parser, "sequence-dictionary");
    // `validateArguments`: a missing input passes `validateInputs` and is refused when it is read.
    non_empty(&argument(parser, "output-prefix").unwrap_or_default())?;
    prepare_output_directory(&require(parser, "output"))?;

    let standardized_ratios = tables::read_copy_ratios(&standardized)?;
    let denoised_ratios = tables::read_copy_ratios(&denoised)?;
    if standardized_ratios.spans() != denoised_ratios.spans() {
        return Err(illegal("Intervals in input files must be identical."));
    }
    if standardized_ratios.metadata != denoised_ratios.metadata {
        return Err(illegal("Metadata do not match."));
    }

    let to_plot = load_fasta_dictionary(&dictionary)?;
    let map = contig_length_map(&to_plot, minimum_contig_length(parser))?;
    validate_contigs(&map, &standardized_ratios.spans(), &standardized)?;
    validate_contigs(&map, &denoised_ratios.spans(), &denoised)?;
    Err(stop_at("Rscript", "PlotDenoisedCopyRatios.R"))
}

/// `PlotModeledSegments.validateNumPointsPerContig`, one kind of point: every contig the segments
/// name, with the points the other file has there.
fn validate_points(
    segments: &[(String, i32)],
    points: &[String],
    message: &str,
) -> Result<(), Thrown> {
    let mut contigs: Vec<&String> = Vec::new();
    for (contig, _) in segments {
        if !contigs.contains(&contig) {
            contigs.push(contig);
        }
    }
    let consistent = contigs.iter().all(|contig| {
        let claimed: i32 = segments
            .iter()
            .filter(|(other, _)| other == *contig)
            .map(|(_, count)| *count)
            .sum();
        let found = points.iter().filter(|other| other == contig).count() as i32;
        claimed == found
    });
    if consistent {
        Ok(())
    } else {
        Err(illegal(message))
    }
}

/// `PlotModeledSegments.doWork`.
pub fn plot_modeled_segments(parser: &Parser) -> Outcome {
    let denoised = argument(parser, "denoised-copy-ratios");
    let allelic = argument(parser, "allelic-counts");
    let segments = require(parser, "segments");
    let dictionary = require(parser, "sequence-dictionary");
    if denoised.is_none() && allelic.is_none() {
        return Err(illegal(
            "Must provide at least a denoised-copy-ratios file or an allelic-counts file.",
        ));
    }
    non_empty(&argument(parser, "output-prefix").unwrap_or_default())?;
    prepare_output_directory(&require(parser, "output"))?;

    let ratios = denoised
        .as_deref()
        .map(tables::read_copy_ratios)
        .transpose()?;
    let counts = allelic
        .as_deref()
        .map(tables::read_allelic_counts)
        .transpose()?;
    let modeled = tables::read_modeled_segments(&segments)?;

    // `getValidatedMetadata`.
    if ratios
        .as_ref()
        .is_some_and(|r| r.metadata != modeled.metadata)
        || counts
            .as_ref()
            .is_some_and(|c| c.metadata != modeled.metadata)
    {
        return Err(illegal("Metadata do not match."));
    }

    if let Some(ratios) = &ratios {
        let claimed: Vec<(String, i32)> = modeled
            .records
            .iter()
            .map(|record| (record.span.contig.clone(), record.num_points_copy_ratio))
            .collect();
        let points: Vec<String> = ratios.spans().into_iter().map(|span| span.contig).collect();
        validate_points(
            &claimed,
            &points,
            "Number of denoised-copy-ratio points in input modeled-segments file is inconsistent \
             with that in input denoised-copy-ratios file.",
        )?;
    }
    if let Some(counts) = &counts {
        let claimed: Vec<(String, i32)> = modeled
            .records
            .iter()
            .map(|record| {
                (
                    record.span.contig.clone(),
                    record.num_points_allele_fraction,
                )
            })
            .collect();
        let points: Vec<String> = counts.sites().into_iter().map(|site| site.0).collect();
        validate_points(
            &claimed,
            &points,
            "Number of allelic-count points in input modeled-segments file is inconsistent with \
             that in input heterozygous allelic-counts file.",
        )?;
    }

    let to_plot = load_fasta_dictionary(&dictionary)?;
    let map = contig_length_map(&to_plot, minimum_contig_length(parser))?;
    if let (Some(ratios), Some(path)) = (&ratios, &denoised) {
        validate_contigs(&map, &ratios.spans(), path)?;
    }
    if let (Some(counts), Some(path)) = (&counts, &allelic) {
        let spans: Vec<Span> = counts
            .sites()
            .into_iter()
            .map(|(contig, position)| Span {
                contig,
                start: position,
                end: position,
            })
            .collect();
        validate_contigs(&map, &spans, path)?;
    }
    let spans: Vec<Span> = modeled.records.iter().map(|r| r.span.clone()).collect();
    validate_contigs(&map, &spans, &segments)?;
    Err(stop_at("Rscript", "PlotModeledSegments.R"))
}

/// `DetermineGermlineContigPloidy.onStartup`.
pub fn determine_germline_contig_ploidy(_parser: &Parser) -> Outcome {
    Err(import_stop("gcnvkernel"))
}

/// `GermlineCNVCaller.onStartup`.
pub fn germline_cnv_caller(_parser: &Parser) -> Outcome {
    Err(import_stop("gcnvkernel"))
}

/// `NVScoreVariants.onStartup`.
pub fn nv_score_variants(_parser: &Parser) -> Outcome {
    Err(import_stop("scorevariants"))
}

/// `new SimpleIntervalCollection(<shard>/interval_list.tsv)`, for its dictionary.
fn shard_dictionary(shard: &str) -> Result<Vec<SequenceRecord>, Thrown> {
    tables::read_interval_collection(&format!("{shard}/interval_list.tsv"))
}

/// `PostprocessGermlineCNVCalls.getBestAvailableSequenceDictionary`.
fn shards_dictionary(parser: &Parser) -> Result<SamHeader, Thrown> {
    let calls: Vec<Vec<SequenceRecord>> = arguments(parser, "calls-shard-path")
        .iter()
        .map(|shard| shard_dictionary(shard))
        .collect::<Result<_, _>>()?;
    let models: Vec<Vec<SequenceRecord>> = arguments(parser, "model-shard-path")
        .iter()
        .map(|shard| shard_dictionary(shard))
        .collect::<Result<_, _>>()?;
    let first = calls.first().cloned().unwrap_or_default();
    if calls.iter().any(|dictionary| *dictionary != first) {
        return Err(illegal(
            "The SAM sequence dictionary is not the same for all of the call shards.",
        ));
    }
    if models.iter().any(|dictionary| *dictionary != first) {
        return Err(illegal(
            "The SAM sequence dictionary is either not the same for all of the model shards, or \
             is different from the SAM sequence dictionary of calls shards.",
        ));
    }
    Ok(SamHeader {
        sequences: first,
        ..SamHeader::default()
    })
}

/// `PostprocessGermlineCNVCalls.onStartup`: `GATKTool`'s, and then the import.
pub fn postprocess_germline_cnv_calls(parser: &Parser) -> Outcome {
    // `loadMasterSequenceDictionary`, `initializeReference`, `initializeReads`.
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    let inputs = arguments(parser, "input");
    let reads: Vec<Vec<SequenceRecord>> = inputs
        .iter()
        .map(|path| reads_dictionary(parser, path, inputs.len()))
        .collect::<Result<_, _>>()?;

    // `initializeIntervals`, where this tool's own override answers before any other source.
    if !arguments(parser, "intervals").is_empty()
        || !arguments(parser, "exclude-intervals").is_empty()
    {
        let best = match &master {
            Some(master) => master.clone(),
            None => shards_dictionary(parser)?,
        };
        interval_arguments(parser, &best)?;
    }

    // `validateSequenceDictionaries`.
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            for one in &reads {
                validate_against_master(master, "reads", one)?;
            }
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
        }
        if let Some(reference) = &reference {
            for one in &reads {
                gatk_tools::sequence_dictionary::validate(
                    "reference",
                    &reference.sequences,
                    "reads",
                    one,
                    false,
                    false,
                )
                .map_err(|refusal| Thrown {
                    failure: Failure::User,
                    exception: refusal.java_class(),
                    message: Some(refusal.message()),
                })?;
            }
        }
        // The feature dictionaries are compared with the best available one, asked for even
        // where this tool has no features, which is what reads the shards on every run.
        if master.is_none() {
            shards_dictionary(parser)?;
        }
    }

    Err(import_stop("gcnvkernel"))
}

/// `new HDF5File(file, READ_ONLY)`, reported as each reader of the file reports it.
fn open_hdf5(path: &str, reading: &str) -> Result<Hdf5Reader, Thrown> {
    let bytes = std::fs::read(path).unwrap_or_default();
    let absolute = java_absolute_path(path);
    Hdf5Reader::open(bytes).ok_or_else(|| {
        hdf5_refusal(
            path,
            reading,
            &format!("exception when opening '{absolute}' with READ_ONLY mode: Not an HDF5 file"),
        )
    })
}

/// A `HDF5LibException` as `readLabel` and the others wrap it.
fn hdf5_refusal(path: &str, reading: &str, what: &str) -> Thrown {
    Thrown::non_user(
        GATK_EXCEPTION,
        format!(
            "Exception encountered during reading of {reading} from {}: \
             org.broadinstitute.hdf5.HDF5LibException: {what}",
            java_absolute_path(path)
        ),
    )
}

/// `HDF5File`'s refusal of a dataset that is not there, which names the path as it was given.
fn hdf5_no_dataset(path: &str, reading: &str, dataset: &str) -> Thrown {
    hdf5_refusal(
        path,
        reading,
        &format!(
            "exception when reading from data-set '{dataset}' in file '{path}': Object not found"
        ),
    )
}

/// `LabeledVariantAnnotationsData.readAnnotationNames`.
fn read_annotation_names(path: &str) -> Result<Vec<String>, Thrown> {
    let reader = open_hdf5(path, "annotation names")?;
    reader
        .strings("/annotations/names")
        .ok_or_else(|| hdf5_no_dataset(path, "annotation names", "/annotations/names"))
}

/// `LabeledVariantAnnotationsData.readAnnotations`: `HDF5Utils.readChunkedDoubleMatrix`.
fn read_annotation_rows(path: &str) -> Result<Vec<Vec<f64>>, Thrown> {
    let reader = open_hdf5(path, "annotations")?;
    let refuse = || hdf5_no_dataset(path, "annotations", "/annotations");
    let scalar = |name: &str| -> Result<usize, Thrown> {
        let (_, values) = reader
            .doubles(&format!("/annotations/{name}"))
            .ok_or_else(refuse)?;
        values.first().map(|v| *v as usize).ok_or_else(refuse)
    };
    let (chunks, columns) = (scalar("num_chunks")?, scalar("num_columns")?);
    let mut rows = Vec::new();
    for chunk in 0..chunks {
        let (_, values) = reader
            .doubles(&format!("/annotations/chunk_{chunk}"))
            .ok_or_else(refuse)?;
        rows.extend(values.chunks(columns.max(1)).map(<[f64]>::to_vec));
    }
    Ok(rows)
}

/// `LabeledVariantAnnotationsData.readLabel`.
fn read_label(path: &str, label: &str) -> Result<Vec<bool>, Thrown> {
    let reading = format!("label {label}");
    let reader = open_hdf5(path, &reading)?;
    let (_, values) = reader
        .doubles(&format!("/labels/{label}"))
        .ok_or_else(|| hdf5_no_dataset(path, &reading, &format!("/labels/{label}")))?;
    Ok(values.iter().map(|value| *value == 1.0).collect())
}

/// `TrainVariantAnnotationsModel.readAndValidateAnnotations`, for the one refusal it can make.
fn validate_annotations(
    names: &[String],
    rows: &[&Vec<f64>],
    prefix_with_tag: &str,
) -> Result<(), Thrown> {
    if names.is_empty() {
        return Err(illegal("Number of annotation names must be positive."));
    }
    if rows.is_empty() {
        return Err(illegal("Number of data points must be positive."));
    }
    if names.len() != rows[0].len() {
        return Err(illegal(
            "Number of annotation names must match the number of features in the annotation data.",
        ));
    }
    let missing: Vec<&str> = (0..names.len())
        .filter(|column| rows.iter().all(|row| row[*column].is_nan()))
        .map(|column| names[column].as_str())
        .collect();
    if !missing.is_empty() {
        return Err(bad_input(format!(
            "All values of the following annotations are missing in the training data for the \
             {prefix_with_tag} model: [{}]. Consider repeating the extraction step without \
             specifying these annotations. ",
            missing.join(", ")
        )));
    }
    Ok(())
}

/// `LabeledVariantAnnotationsData.subsetAnnotationsToTemporaryFile`'s three checks and its subset.
fn subset<'a>(
    names: &[String],
    rows: &'a [Vec<f64>],
    keep: &[bool],
) -> Result<Vec<&'a Vec<f64>>, Thrown> {
    if names.is_empty() {
        return Err(illegal("Number of annotation names must be positive."));
    }
    if rows.is_empty() {
        return Err(illegal(
            "Number of annotation data points must be positive.",
        ));
    }
    if names.len() != rows[0].len() {
        return Err(illegal(
            "Number of annotation names must match number of features in annotation data.",
        ));
    }
    Ok(rows
        .iter()
        .zip(keep)
        .filter(|(_, keep)| **keep)
        .map(|(row, _)| row)
        .collect())
}

/// The `--mode` and `--model-backend` constant names, as the parser leaves them.
fn modes_of(parser: &Parser) -> Vec<String> {
    let listed = arguments(parser, "mode");
    if !listed.is_empty() {
        return listed;
    }
    match scalar(parser, "mode") {
        Some(value) => value
            .trim_matches(|c| c == '[' || c == ']')
            .split(", ")
            .map(str::to_string)
            .collect(),
        None => vec!["SNP".to_string(), "INDEL".to_string()],
    }
}

fn backend_of(parser: &Parser) -> String {
    scalar(parser, "model-backend").unwrap_or_else(|| "PYTHON_IFOREST".to_string())
}

/// `TrainVariantAnnotationsModel.doWork`.
pub fn train_variant_annotations_model(parser: &Parser) -> Outcome {
    let annotations = require(parser, "annotations-hdf5");
    let prefix = require(parser, "output");
    let unlabeled = argument(parser, "unlabeled-annotations-hdf5");
    let script = argument(parser, "python-script");
    let hyperparameters = argument(parser, "hyperparameters-json");
    let backend = backend_of(parser);

    // `validateArgumentsAndSetModes`.
    can_read_file(Some(&annotations))?;
    if let Some(unlabeled) = &unlabeled {
        can_read_file(Some(unlabeled))?;
        if read_annotation_names(&annotations)? != read_annotation_names(unlabeled)? {
            return Err(illegal(
                "Annotation names must be identical for positive and unlabeled annotations.",
            ));
        }
    }
    let script_given = |name: &str| -> Result<(), Thrown> {
        if script.is_some() {
            return Err(illegal(format!(
                "Python script should not be provided when using {name} backend."
            )));
        }
        Ok(())
    };
    match backend.as_str() {
        "JAVA_BGMM" => {
            script_given("JAVA_BGMM")?;
            can_read_file(hyperparameters.as_deref())?;
        }
        "PYTHON_SCRIPT" => {
            if hyperparameters.is_none() {
                return Err(illegal(
                    "Hyperparameters JSON must be provided when using PYTHON_SCRIPT backend.",
                ));
            }
            can_read_file(script.as_deref())?;
            can_read_file(hyperparameters.as_deref())?;
        }
        _ => {
            script_given("PYTHON_IFOREST")?;
            if hyperparameters.is_some() {
                can_read_file(hyperparameters.as_deref())?;
            }
            return Err(import_stop("argparse"));
        }
    }

    // `doModelingWorkForVariantType`, SNP before INDEL whatever `--mode` says. The first model
    // needs what is missing, so the first type with training sites is where this ends.
    let modes = modes_of(parser);
    for (kind, is_snp_model) in [("SNP", true), ("INDEL", false)] {
        if !modes.iter().any(|mode| mode == kind) {
            continue;
        }
        let names = read_annotation_names(&annotations)?;
        let rows = read_annotation_rows(&annotations)?;
        let is_training = read_label(&annotations, extract::TRAINING_LABEL)?;
        let _is_calibration = read_label(&annotations, extract::CALIBRATION_LABEL)?;
        let is_snp = read_label(&annotations, extract::SNP_LABEL)?;
        let keep: Vec<bool> = is_training
            .iter()
            .zip(&is_snp)
            .map(|(training, snp)| *training && (*snp == is_snp_model))
            .collect();
        if !keep.iter().any(|keep| *keep) {
            return Err(bad_input(format!(
                "Attempted to train {kind} model, but no suitable training sites were found in \
                 the provided annotations."
            )));
        }
        let training = subset(&names, &rows, &keep)?;
        let tagged = format!("{prefix}.{}", kind.to_lowercase());
        let mut unlabeled_rows = None;
        if let Some(unlabeled) = &unlabeled {
            let all = read_annotation_rows(unlabeled)?;
            let unlabeled_snp = read_label(unlabeled, extract::SNP_LABEL)?;
            let keep: Vec<bool> = unlabeled_snp
                .iter()
                .map(|snp| *snp == is_snp_model)
                .collect();
            if !keep.iter().any(|keep| *keep) {
                return Err(bad_input(format!(
                    "Attempted to train {kind} model, but no suitable unlabeled sites were found \
                     in the provided annotations."
                )));
            }
            let kept = subset(&names, &all, &keep)?
                .into_iter()
                .cloned()
                .collect::<Vec<_>>();
            unlabeled_rows = Some(kept);
        }
        // `trainAndSerializeModel`.
        validate_annotations(&names, &training, &tagged)?;
        if let Some(rows) = &unlabeled_rows {
            validate_annotations(&names, &rows.iter().collect::<Vec<_>>(), &tagged)?;
        }
        if backend == "JAVA_BGMM" {
            return Err(Thrown::non_user(NOT_IMPLEMENTED, BGMM_NOT_IMPLEMENTED));
        }
        return Err(stop_at("python", "Running --python-script"));
    }
    Ok(None)
}

/// `ScoreVariantAnnotations.readCalibrationScoresAndCreateConverter`, whether there is one.
fn calibration_exists(prefix: &str, tag: &str) -> Result<bool, Thrown> {
    let path = format!("{prefix}.{tag}.calibrationScores.hdf5");
    if std::fs::File::open(&path).is_err() {
        return Ok(false);
    }
    let reader = open_hdf5(&path, "scores")?;
    let (_, scores) = reader
        .doubles("/data/scores")
        .ok_or_else(|| hdf5_no_dataset(&path, "scores", "/data/scores"))?;
    if !scores.iter().all(|score| score.is_finite()) {
        return Err(illegal("Calibration scores must all be finite."));
    }
    Ok(true)
}

/// `ScoreVariantAnnotations.constructVCFHeader`, over the input's lines.
fn scored_header(
    parser: &Parser,
    input: &htsjdk_vcf::header::VcfHeader,
    sorted_labels: &[String],
    keys: &ScoreKeys,
) -> htsjdk_vcf::header::VcfHeader {
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    let info = |id: &str, number: Cardinality, line_type: LineType, description: String| {
        HeaderLine::Compound {
            key: "INFO".to_string(),
            id: id.to_string(),
            number,
            line_type,
            description,
            extra: Vec::new(),
        }
    };
    let mut lines = input.lines.clone();
    lines.push(info(
        &keys.score,
        Cardinality::Fixed(1),
        LineType::Float,
        "Score according to the model applied by ScoreVariantAnnotations".to_string(),
    ));
    lines.push(info(
        &keys.calibration,
        Cardinality::Fixed(1),
        LineType::Float,
        format!(
            "Calibration sensitivity corresponding to the value of {}",
            keys.score
        ),
    ));
    lines.push(HeaderLine::Filter {
        id: keys.filter.clone(),
        description: "Low score (corresponding to high calibration sensitivity)".to_string(),
    });
    if let Some(snp) = &keys.snp {
        lines.push(info(
            snp,
            Cardinality::Fixed(0),
            LineType::Flag,
            "This site was considered a SNP during filtering".to_string(),
        ));
    }
    for label in sorted_labels {
        lines.push(info(
            label,
            Cardinality::Fixed(0),
            LineType::Flag,
            format!("This site was labeled as {label} according to resources"),
        ));
    }
    lines.extend(default_tool_vcf_header_lines(
        parser,
        "ScoreVariantAnnotations",
    ));
    htsjdk_vcf::header::VcfHeader {
        lines,
        samples: input.samples.clone(),
    }
}

/// The four names the score's INFO keys and filter go under. A `null` value is the reference's
/// `null`, which only `--snp-key` may be.
struct ScoreKeys {
    snp: Option<String>,
    score: String,
    calibration: String,
    filter: String,
}

impl ScoreKeys {
    /// The names as `new VCFInfoHeaderLine` and `new VCFFilterHeaderLine` take them, in
    /// `constructVCFHeader`'s order, each refusing a `null`.
    fn of(parser: &Parser) -> Result<ScoreKeys, Thrown> {
        let score = argument(parser, "score-key").ok_or_else(|| {
            illegal(
                "Invalid VCFCompoundHeaderLine: key=INFO name=null type=Float desc=Score according \
                 to the model applied by ScoreVariantAnnotations lineType=INFO",
            )
        })?;
        let calibration = argument(parser, "calibration-sensitivity-key").ok_or_else(|| {
            illegal(format!(
                "Invalid VCFCompoundHeaderLine: key=INFO name=null type=Float desc=Calibration \
                 sensitivity corresponding to the value of {score} lineType=INFO"
            ))
        })?;
        let filter = argument(parser, "low-score-filter-name")
            .ok_or_else(|| illegal("Invalid VCFSimpleHeaderLine: key=FILTER name=null"))?;
        Ok(ScoreKeys {
            snp: argument(parser, "snp-key"),
            score,
            calibration,
            filter,
        })
    }
}

/// `ScoreVariantAnnotations`, up to its first script.
pub fn score_variant_annotations(parser: &Parser) -> Outcome {
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "ScoreVariantAnnotations")?;
    let prefix = require(parser, "output");
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;
    let LabelledStart {
        annotations_path,
        vcf_path,
        resources,
        arguments: walker_arguments,
        allele_specific,
        names,
        sorted_labels,
    } = labelled_walker_start(parser, &prefix, &file.header.lines)?;

    // `createVCFWriter` opens the file, and `writeHeader` leaves its bytes in the writer's buffer.
    write_file(&vcf_path, b"")?;
    let keys = ScoreKeys::of(parser)?;
    let header = scored_header(parser, &file.header, &sorted_labels, &keys);

    // `afterOnTraversalStart`, whose other three `Utils.nonNull` the header has answered.
    if argument(parser, "double-format").is_none() {
        return Err(illegal("Null object is not allowed here."));
    }
    let backend = backend_of(parser);
    let script = argument(parser, "python-script");
    let model_prefix = require(parser, "model-prefix");
    let readable = |suffix: String| std::fs::File::open(format!("{model_prefix}{suffix}")).is_ok();
    let (snp_scorer, indel_scorer);
    match backend.as_str() {
        "JAVA_BGMM" => {
            if script.is_some() {
                return Err(illegal(
                    "Python script should not be provided when using JAVA_BGMM backend.",
                ));
            }
            // `deserializeScorerFromSerFiles`, SNP first; `deserialize` is the missing module.
            if readable(".snp.bgmmScorer.ser".to_string())
                || readable(".indel.bgmmScorer.ser".to_string())
            {
                return Err(Thrown::non_user(NOT_IMPLEMENTED, BGMM_NOT_IMPLEMENTED));
            }
            snp_scorer = false;
            indel_scorer = false;
        }
        "PYTHON_SCRIPT" => {
            can_read_file(script.as_deref())?;
            snp_scorer = readable(".snp.scorer.pkl".to_string());
            indel_scorer = readable(".indel.scorer.pkl".to_string());
        }
        _ => {
            if script.is_some() {
                return Err(illegal(
                    "Python script should not be provided when using PYTHON_IFOREST backend.",
                ));
            }
            return Err(import_stop("argparse"));
        }
    }
    if !snp_scorer && !indel_scorer {
        return Err(bad_input(format!(
            "At least one serialized scorer must be present in the model files with the prefix \
             {model_prefix}."
        )));
    }
    let wants = |kind: extract::VariantType| walker_arguments.modes.contains(&kind);
    if wants(extract::VariantType::Snp) && !snp_scorer {
        return Err(bad_input(
            "SNPs were indicated for extraction via the mode argument, but no serialized SNP \
             scorer was available in the model files with the prefix."
                .to_string(),
        ));
    }
    if wants(extract::VariantType::Indel) && !indel_scorer {
        return Err(bad_input(
            "INDELs were indicated for extraction via the mode argument, but no serialized INDEL \
             scorer was available in the model files with the prefix."
                .to_string(),
        ));
    }
    let snp_calibration = calibration_exists(&model_prefix, "snp")?;
    let indel_calibration = calibration_exists(&model_prefix, "indel")?;
    if !snp_calibration && scalar(parser, "snp-calibration-sensitivity-threshold").is_some() {
        return Err(bad_input(format!(
            "The snp-calibration-sensitivity-threshold argument was specified, but no SNP \
             calibration scores were provided in the model files with the prefix {model_prefix}."
        )));
    }
    if !indel_calibration && scalar(parser, "indel-calibration-sensitivity-threshold").is_some() {
        return Err(bad_input(format!(
            "The indel-calibration-sensitivity-threshold argument was specified, but no INDEL \
             calibration scores were provided in the model files with the prefix {model_prefix}."
        )));
    }
    validate_output_files(&[&format!("{prefix}.scores.hdf5")])?;

    // The first pass. A pass that extracts anything needs a script to score it, and the
    // annotations file is written before one is looked for.
    let traversed = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut extracted: Vec<Vec<extract::Row>> = Vec::new();
    for record in &traversed {
        let reduced_record = reduced(record);
        let found = metadata(&reduced_record, &resources, &walker_arguments, true)?;
        if !found.is_empty() {
            extracted.push(
                found
                    .iter()
                    .map(|e| extract::row(&reduced_record, e, &names, allele_specific))
                    .collect(),
            );
        }
    }
    if !extracted.is_empty() {
        write_file(
            &annotations_path,
            &matrix(
                &extracted,
                &names,
                &sorted_labels,
                allele_specific,
                flag(parser, "omit-alleles-in-hdf5"),
            ),
        )?;
        return Err(stop_at("python", "Running --python-script"));
    }

    // The second pass: nothing was extracted, so every record goes through as it was read.
    let written: Vec<htsjdk_vcf::variant::VariantContext> =
        traversed.into_iter().cloned().collect();
    let vcf = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    write_variant_output(parser, &vcf_path, &vcf)?;
    Ok(None)
}
