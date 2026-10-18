//! The Spark tools (X.1), run as the `--spark-master local[1]` reference runs them.
//!
//! A `GATKSparkTool` is not a walker. `SparkCommandLineProgram.doWork` builds the Spark context
//! first, then `runPipeline` initialises the reference, the reads, the features and the intervals in
//! that order, validates the dictionaries, and hands the tool an RDD of reads. The port has no
//! context and no partitions: what it reproduces is what a single-partition run observes, which is
//! the claim X.1 makes first, and every refusal on the way in the order the reference reaches it.
//!
//! Four things separate this startup from [`read_walker_startup`]:
//!
//! * the Spark arguments are read before anything is opened: a `--conf` that is not
//!   `name=value` and a master URL Spark cannot parse are refused while the context is built;
//! * the reads header comes from Disq, which picks the format by EXTENSION and wraps whatever it
//!   throws in `Failed to read bam header from <path>`;
//! * the best available dictionary is the reference's, then the reads'; there is no master
//!   dictionary, because a Spark tool declares no `--sequence-dictionary`;
//! * one input only, for the tools that do not merge, and `--read-index` names nothing: Disq finds
//!   the index beside the BAM itself.
use super::*;

/// What `runPipeline` hands `runTool`: the reads, their header, the traversal intervals and the
/// resolved filter chain.
struct SparkStart {
    source: ReadsDataSource,
    header: SamHeader,
    intervals: Vec<gatk_engine::interval::SimpleInterval>,
    filters: Vec<gatk_tools::filter_resolution::ResolvedFilter>,
}

/// `SparkContext.parseMasterUrl`'s accepted shapes, for the masters that run in this process.
///
/// `local`, `local[N]`, `local[*]`, `local[N,F]` and `local[*,F]` all run the tool here, and the
/// port answers each the same way because a count does not depend on the thread count. A cluster
/// master would submit the job elsewhere, which the port does not do, and anything else is the
/// `SparkException` the context refuses with.
fn spark_master(parser: &Parser) -> Result<(), Thrown> {
    // `SparkContextFactory.DEFAULT_SPARK_MASTER`.
    let master = scalar(parser, "spark-master").unwrap_or_else(|| "local[*]".to_string());
    let threads = |text: &str| text == "*" || text.parse::<u32>().is_ok_and(|n| n > 0);
    let local = master == "local"
        || master
            .strip_prefix("local[")
            .and_then(|rest| rest.strip_suffix(']'))
            .is_some_and(|inside| match inside.split_once(',') {
                Some((n, failures)) => threads(n) && failures.parse::<u32>().is_ok(),
                None => threads(inside),
            });
    if local {
        return Ok(());
    }
    let cluster = ["spark://", "mesos://", "k8s://", "local-cluster["]
        .iter()
        .any(|prefix| master.starts_with(prefix))
        || master == "yarn";
    if cluster {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "--spark-master {master} submits the job to a cluster, which this port does not \
                 do: it runs a Spark tool as a local master does. This message is the port's own \
                 and not GATK's."
            ),
        ));
    }
    Err(Thrown::non_user(
        "org.apache.spark.SparkException",
        format!("Could not parse Master URL: '{master}'"),
    ))
}

/// `SparkCommandLineArgumentCollection.getSparkProperties`, for its refusal.
///
/// The split is Java's `String.split("=")`, which drops TRAILING empty strings: `a=` is one piece
/// and refused, `a=b=c` is three and refused, `=b` is two with an empty name and refused.
fn spark_properties(parser: &Parser) -> Result<(), Thrown> {
    for property in arguments(parser, "conf") {
        let mut pieces: Vec<&str> = property.split('=').collect();
        while pieces.len() > 1 && pieces.last() == Some(&"") {
            pieces.pop();
        }
        if pieces.len() != 2 || pieces[0].is_empty() || pieces[1].is_empty() {
            return Err(Thrown::command_line(format!(
                "Argument conf has a bad value: {property}. Expected a value of the form \
                 spark.property.name=value"
            )));
        }
    }
    Ok(())
}

/// `SparkContext.setLogLevel`, which upper-cases the level and requires it to be one of Log4j's.
///
/// The message lists `VALID_LOG_LEVELS` in the order the Scala `Set` iterates it, which is neither
/// alphabetical nor Log4j's; measured on the reference.
fn spark_verbosity(parser: &Parser) -> Result<(), Thrown> {
    const VALID: [&str; 8] = [
        "ALL", "DEBUG", "ERROR", "WARN", "INFO", "OFF", "FATAL", "TRACE",
    ];
    let Some(level) = scalar(parser, "spark-verbosity") else {
        return Ok(());
    };
    if VALID.contains(&level.to_uppercase().as_str()) {
        return Ok(());
    }
    Err(Thrown::non_user(
        "java.lang.IllegalArgumentException",
        format!(
            "requirement failed: Supplied level {level} did not match one of: {}",
            VALID.join(",")
        ),
    ))
}

/// `ReadsSparkSource.getHeader` for a local path: Disq's format by extension, then the file.
fn spark_reads_header(input: &str) -> Result<ReadsDataSource, Thrown> {
    let failed = |cause: String| {
        Thrown::user(format!(
            "Failed to read bam header from {input}\n Caused by:{cause}"
        ))
    };
    let lower = input.to_ascii_lowercase();
    if lower.ends_with(".cram") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "A CRAM input to a Spark tool is a GATK feature that this port does not carry yet. \
             This message is the port's own and not GATK's.",
        ));
    }
    if !(lower.ends_with(".bam") || lower.ends_with(".sam")) {
        return Err(failed(format!("Cannot find format extension for {input}")));
    }
    let path = std::path::Path::new(input);
    if !path.exists() {
        return Err(failed(format!("File {input} does not exist")));
    }
    ReadsDataSource::open_unindexed(path).map_err(|error| failed(format!("{error:?}")))
}

/// `SparkCommandLineProgram.doWork` up to `runTool`, for a tool that reads one input.
fn spark_startup(parser: &Parser, tool: &str) -> Result<SparkStart, Thrown> {
    // The plugin descriptor validates while the command line is parsed.
    let filters = resolve_read_filters(parser, tool)?;

    // `SparkContextFactory.getSparkContext(name, getSparkProperties(), master)`: the properties
    // are parsed as the call's argument, before the master is.
    spark_properties(parser)?;
    spark_master(parser)?;
    spark_verbosity(parser)?;

    // `initializeReference`, before the reads.
    let reference = reference_dictionary(parser)?;

    // `initializeReads`: `doNotMerge` is the default policy, and these tools keep it.
    let inputs = arguments(parser, "input");
    if inputs.len() > 1 {
        return Err(Thrown::user(
            "Sorry, we only support a single reads input for for this spark tool.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let source = spark_reads_header(&input)?;
    let header = source.header().clone();

    // `initializeIntervals`, against `getBestAvailableSequenceDictionary`: the reference's, then
    // the reads'.
    let best = reference.clone().unwrap_or_else(|| header.clone());
    let intervals = match interval_arguments(parser, &best)? {
        Some(parameters) => parameters.intervals,
        None => Vec::new(),
    };

    // `validateSequenceDictionaries`: the reference against the reads, and nothing else.
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(reference) = &reference {
            gatk_tools::sequence_dictionary::validate(
                "reference",
                &reference.sequences,
                "reads",
                &header.sequences,
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

    // `getUnfilteredReads` with intervals is a Disq traversal, which needs the index beside the
    // BAM. Without one the reference fails inside an executor task, and that failure is Spark's
    // log rather than a GATK refusal; the port says so in its own words.
    let intervals_specified = !arguments(parser, "intervals").is_empty()
        || !arguments(parser, "exclude-intervals").is_empty();
    let source = if intervals_specified {
        let path = std::path::Path::new(&input);
        match htsjdk_bam::sam_files::find_index(path) {
            Some(index) => ReadsDataSource::open(path, &index)
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
            None => {
                return Err(Thrown::non_user(
                    PORT_LIMITATION,
                    format!(
                        "Intervals set but no index file found for file:{input}: the reference \
                         fails this inside a Spark task, which this port does not reproduce. \
                         This message is the port's own and not GATK's."
                    ),
                ))
            }
        }
    } else {
        source
    };

    Ok(SparkStart {
        source,
        header,
        intervals,
        filters,
    })
}

/// `-O`, written with `print`: the text and no newline after it.
fn spark_output(parser: &Parser, text: &str) -> Result<(), Thrown> {
    if let Some(output) = argument(parser, "output") {
        std::fs::write(&output, text)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    }
    Ok(())
}

/// `CountReadsSpark.runTool`: `reads.count()`, PRINTED with `println` rather than returned.
///
/// `runTool` returns nothing and `doWork` returns null, so `handleResult` prints nothing: the
/// number on standard output is the tool's own line, with no `Tool returned:` above it.
pub fn count_reads_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "CountReadsSpark")?;
    let filter = read_filter(parser, &filters, &header)?;
    let count = gatk_tools::count_reads::count_reads(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    gatk_tools::tool_stdout::println(&count.to_string());
    spark_output(parser, &gatk_tools::count_reads::output(count))?;
    Ok(None)
}

/// `CountBasesSpark.runTool`: the read lengths folded into a sum, printed the same way.
pub fn count_bases_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "CountBasesSpark")?;
    let filter = read_filter(parser, &filters, &header)?;
    let count = gatk_tools::count_reads::count_bases(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    gatk_tools::tool_stdout::println(&count.to_string());
    spark_output(parser, &gatk_tools::count_reads::output(count))?;
    Ok(None)
}

/// `FlagStatSpark.runTool`: `FlagStatus` aggregated over the reads, and its `toString` printed.
///
/// One partition merges nothing, so the aggregate is the walker's sum record by record.
pub fn flag_stat_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "FlagStatSpark")?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let mut status = gatk_tools::counting_walkers::FlagStatus::default();
    for record in &records {
        status.add(
            record,
            contig_name(&header, record.reference_index),
            contig_name(&header, record.mate_reference_index),
        );
    }
    let text = status.to_text();
    gatk_tools::tool_stdout::println(&text);
    spark_output(parser, &text)?;
    Ok(None)
}
