//! `CreateReadCountPanelOfNormals`, in a file of its own.
//!
//! The filters, the imputation, the truncation, the standardisation, the GC-bias correction, the
//! decomposition and the panel's HDF5 tree are [`gatk_tools::create_read_count_panel_of_normals`];
//! what is here is `runPipeline`: the argument checks, the read-counts files read into one matrix,
//! the annotated intervals validated against them, and the panel written.
use super::*;

/// One `SimpleCountCollection` TSV: its header block, its intervals and its counts.
struct CountsFile {
    header: String,
    intervals: Vec<gatk_tools::create_read_count_panel_of_normals::Interval>,
    counts: Vec<f64>,
}

/// `SimpleCountCollection.read`, for the TSV the collection writes. An HDF5 counts file is the
/// other format it reads, and this port has no HDF5 reader.
fn read_counts(path: &str) -> Result<CountsFile, Thrown> {
    use gatk_tools::create_read_count_panel_of_normals::Interval;
    can_read_file(path)?;
    if path.ends_with(".hdf5") || path.ends_with(".h5") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "an HDF5 read-counts file is a format this port does not read yet. This message is \
             the port's own and not GATK's.",
        ));
    }
    let text =
        std::fs::read_to_string(path).map_err(|error| Thrown::user(format!("{path}: {error}")))?;
    let header: String = text
        .lines()
        .take_while(|line| line.starts_with('@'))
        .map(|line| format!("{line}\n"))
        .collect();
    let mut intervals = Vec::new();
    let mut counts = Vec::new();
    for line in text
        .lines()
        .filter(|line| !line.starts_with('@') && !line.trim().is_empty())
        .skip(1)
    {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 4 {
            continue;
        }
        intervals.push(Interval {
            contig: columns[0].to_string(),
            start: columns[1].parse().unwrap_or_default(),
            end: columns[2].parse().unwrap_or_default(),
        });
        counts.push(columns[3].parse().unwrap_or(f64::NAN));
    }
    Ok(CountsFile {
        header,
        intervals,
        counts,
    })
}

/// `SAMTextHeaderCodec.encode(new SAMFileHeader(dictionary))`: a header of nothing but the
/// version and the dictionary. The `@SQ` lines are the counts file's own, which that codec wrote.
fn dictionary_text(header: &str) -> String {
    let mut text = String::from("@HD\tVN:1.6\n");
    for line in header.lines().filter(|line| line.starts_with("@SQ")) {
        text.push_str(line);
        text.push('\n');
    }
    text
}

/// `HDF5SVDReadCountPanelOfNormals.create`'s catch: whatever throws once the panel file is open is
/// logged with the output's path, the partial file deleted, and a `GATKException` thrown in its
/// place, so the run's own refusal is the log line rather than the wrapper's message. `exception`
/// is the caught one's `toString()`.
///
/// The log line carries the time of day as log4j prints it, in the container's zone; the coverage
/// harness reads the line with that stamp taken off, which is the only part no two runs share.
fn creation_failure(exception: &str, output: &str) -> Thrown {
    let since_midnight = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_millis() % 86_400_000)
        .unwrap_or_default();
    eprintln!(
        "{:02}:{:02}:{:02}.{:03} WARN  HDF5SVDReadCountPanelOfNormals - Exception encountered during \
         creation of panel of normals ({exception}).  Attempting to delete partial output in {}...",
        since_midnight / 3_600_000,
        since_midnight / 60_000 % 60,
        since_midnight / 1000 % 60,
        since_midnight % 1000,
        java_absolute_path(output)
    );
    let _ = std::fs::remove_file(output);
    Thrown::non_user(
        "org.broadinstitute.hellbender.exceptions.GATKException",
        "Could not create panel of normals.  It may be necessary to use stricter parameters for \
         filtering.  For example, use a larger value of minimum-interval-median-percentile.",
    )
}

/// `CreateReadCountPanelOfNormals.runPipeline`.
pub fn create_read_count_panel_of_normals(parser: &Parser) -> Outcome {
    use gatk_tools::create_read_count_panel_of_normals as pon;

    let inputs = arguments(parser, "input");
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let annotated = argument(parser, "annotated-intervals");
    let arguments_used = pon::Arguments {
        minimum_interval_median_percentile: double_or(
            parser,
            "minimum-interval-median-percentile",
            pon::DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE,
        ),
        maximum_zeros_in_sample_percentage: double_or(
            parser,
            "maximum-zeros-in-sample-percentage",
            pon::DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE,
        ),
        maximum_zeros_in_interval_percentage: double_or(
            parser,
            "maximum-zeros-in-interval-percentage",
            pon::DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE,
        ),
        extreme_sample_median_percentile: double_or(
            parser,
            "extreme-sample-median-percentile",
            pon::DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE,
        ),
        impute_zeros: scalar(parser, "do-impute-zeros").is_none_or(|value| value == "true"),
        extreme_outlier_truncation_percentile: double_or(
            parser,
            "extreme-outlier-truncation-percentile",
            pon::DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE,
        ),
        number_of_eigensamples: scalar(parser, "number-of-eigensamples")
            .and_then(|text| text.parse().ok())
            .unwrap_or(pon::DEFAULT_NUMBER_OF_EIGENSAMPLES),
    };
    let maximum_chunk_size: i64 = scalar(parser, "maximum-chunk-size")
        .and_then(|text| text.parse().ok())
        .unwrap_or(pon::DEFAULT_MAXIMUM_CHUNK_SIZE);

    // `validateArguments`: duplicates first, then each input that is a file must be readable (one
    // that does not exist passes here and is refused when it is read), then the output's folder.
    let mut seen = std::collections::HashSet::new();
    if !inputs.iter().all(|path| seen.insert(path.clone())) {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "List of input read-counts files cannot contain duplicates.",
        ));
    }
    for path in inputs.iter().chain(annotated.iter()) {
        if std::path::Path::new(path).is_file() {
            can_read_file(path)?;
        }
    }
    let output_path = std::path::Path::new(&output);
    let writable = if output_path.exists() {
        std::fs::metadata(output_path).is_ok_and(|m| !m.permissions().readonly())
    } else {
        let absolute = java_absolute_path(&output);
        std::path::Path::new(&absolute)
            .parent()
            .is_some_and(|parent| parent.is_dir())
    };
    if !writable {
        return Err(Thrown {
            failure: Failure::User,
            exception:
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile",
            message: Some(format!(
                "Couldn't write file {} because : The output file is not writeable.",
                java_absolute_path(&output)
            )),
        });
    }

    let sample_filenames: Vec<String> = inputs.iter().map(|p| java_absolute_path(p)).collect();
    let first = read_counts(&inputs[0])?;
    if first.intervals.len() as i64 > maximum_chunk_size {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            format!(
                "The number of intervals ({}) in each read-counts file cannot exceed the maximum \
                 chunk size ({maximum_chunk_size}).",
                first.intervals.len()
            ),
        ));
    }

    // `validateAnnotatedIntervals`: only the GC column is read, and the intervals must be the
    // counts' own, in their order.
    let gc_content: Option<Vec<f64>> = match &annotated {
        None => None,
        Some(path) => {
            can_read_file(path)?;
            let text = std::fs::read_to_string(path)
                .map_err(|error| Thrown::user(format!("{path}: {error}")))?;
            let mut rows = text
                .lines()
                .filter(|line| !line.starts_with('@') && !line.trim().is_empty());
            let columns: Vec<&str> = rows.next().unwrap_or_default().split('\t').collect();
            let gc_column = columns.iter().position(|c| *c == "GC_CONTENT");
            let mut intervals = Vec::new();
            let mut gc = Vec::new();
            for row in rows {
                let cells: Vec<&str> = row.split('\t').collect();
                intervals.push(pon::Interval {
                    contig: cells[0].to_string(),
                    start: cells[1].parse().unwrap_or_default(),
                    end: cells[2].parse().unwrap_or_default(),
                });
                gc.push(
                    gc_column
                        .and_then(|c| cells.get(c))
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(f64::NAN),
                );
            }
            if intervals != first.intervals {
                return Err(Thrown::non_user(
                    "java.lang.IllegalArgumentException",
                    "Annotated intervals do not match provided intervals.",
                ));
            }
            Some(gc)
        }
    };

    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(inputs.len());
    for (index, path) in inputs.iter().enumerate() {
        let counts = if index == 0 {
            CountsFile {
                header: String::new(),
                intervals: first.intervals.clone(),
                counts: first.counts.clone(),
            }
        } else {
            read_counts(path)?
        };
        if counts.intervals != first.intervals {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                pon::mismatched_intervals_message(path),
            ));
        }
        rows.push(counts.counts);
    }
    let original = pon::Matrix::new(&rows);

    let mut preprocessed =
        pon::preprocess_with_gc(&original, gc_content.as_deref(), &arguments_used).map_err(
            |message| {
                creation_failure(
            &format!("org.broadinstitute.hellbender.exceptions.UserException$BadInput: {message}"),
            &output,
        )
            },
        )?;
    let mut standardized = preprocessed.values.clone();
    pon::standardize(&mut standardized).map_err(|message| {
        creation_failure(
            &format!("java.lang.IllegalArgumentException: {message}"),
            &output,
        )
    })?;
    let panel_samples = standardized.samples;
    let k = pon::number_of_eigensamples(arguments_used.number_of_eigensamples, panel_samples);
    let decomposition = if panel_samples > 1 && k > 0 {
        let decomposition = pon::truncated_svd(&standardized, k);
        if !decomposition
            .singular_values
            .iter()
            .any(|s| *s > pon::SVD_EPSILON)
        {
            return Err(creation_failure(
                &format!(
                    "org.broadinstitute.hellbender.exceptions.UserException: {}",
                    pon::NO_NON_ZERO_SINGULAR_VALUES_MESSAGE
                ),
                &output,
            ));
        }
        Some(decomposition)
    } else {
        None
    };
    preprocessed.values = standardized;

    let command_line = crate::command_line::expanded("CreateReadCountPanelOfNormals", parser);
    let file = pon::write_panel(&pon::Panel {
        command_line: &command_line,
        sequence_dictionary: &dictionary_text(&first.header),
        original_counts: &original,
        sample_filenames: &sample_filenames,
        intervals: &first.intervals,
        gc_content: gc_content.as_deref(),
        preprocessed: &preprocessed,
        decomposition: decomposition.as_ref(),
        maximum_chunk_size,
    });
    write_file(&output, &file.to_bytes())?;
    Ok(None)
}
