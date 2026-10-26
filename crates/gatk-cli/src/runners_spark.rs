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

/// `SparkContextFactory.getSparkContext(name, getSparkProperties(), master)`: the properties are
/// parsed as the call's argument, before the master is, and the log level is set on the context.
fn spark_context(parser: &Parser) -> Result<(), Thrown> {
    spark_properties(parser)?;
    spark_master(parser)?;
    spark_verbosity(parser)
}

/// `-L` or `-XL` on the command line, which is what `intervalsSpecified` asks.
fn intervals_specified(parser: &Parser) -> bool {
    !arguments(parser, "intervals").is_empty() || !arguments(parser, "exclude-intervals").is_empty()
}

/// `validateSequenceDictionaries`: the reference against the reads, and nothing else.
fn validate_reference_against_reads(
    reference: &SamHeader,
    reads: &SamHeader,
) -> Result<(), Thrown> {
    gatk_tools::sequence_dictionary::validate(
        "reference",
        &reference.sequences,
        "reads",
        &reads.sequences,
        false,
        false,
    )
    .map_err(|refusal| Thrown {
        failure: Failure::User,
        exception: refusal.java_class(),
        message: Some(refusal.message()),
    })
}

/// `SparkCommandLineProgram.doWork` up to `runTool`, for a tool that reads one input.
fn spark_startup(parser: &Parser, tool: &str) -> Result<SparkStart, Thrown> {
    // The plugin descriptor validates while the command line is parsed.
    let filters = resolve_read_filters(parser, tool)?;
    spark_context(parser)?;

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

    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(reference) = &reference {
            validate_reference_against_reads(reference, &header)?;
        }
    }

    // `getUnfilteredReads` with intervals is a Disq traversal, which needs the index beside the
    // BAM. Without one the reference fails inside an executor task, and that failure is Spark's
    // log rather than a GATK refusal; the port says so in its own words.
    let source = if intervals_specified(parser) {
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
    // Disq reads from the record its guesser accepts, which is not always the first.
    let source = source
        .with_disq_split_start()
        .map_err(|error| Thrown::user(format!("{error:?}")))?;

    Ok(SparkStart {
        source,
        header,
        intervals,
        filters,
    })
}

/// `CountVariantsSpark.runTool`: Disq's variants over `-L`, counted and printed like the reads.
///
/// The tool replaces the plugin descriptors with none, so it has no read filter to resolve, and it
/// requires no reads: `-I` is opened only when given, for its dictionary. Disq queries a VCF over
/// intervals through its index and refuses one that has none, from inside the job, which surfaces
/// as the `IllegalArgumentException` itself rather than a GATK error.
pub fn count_variants_spark(parser: &Parser) -> Outcome {
    spark_context(parser)?;
    let reference = reference_dictionary(parser)?;
    let inputs = arguments(parser, "input");
    if inputs.len() > 1 {
        return Err(Thrown::user(
            "Sorry, we only support a single reads input for for this spark tool.",
        ));
    }
    let reads = match inputs.first() {
        Some(input) => Some(spark_reads_header(input)?.header().clone()),
        None => None,
    };

    // `getBestAvailableSequenceDictionary`: the reference's, then the reads', and with neither an
    // interval cannot be resolved at all.
    let intervals =
        match reference.as_ref().or(reads.as_ref()) {
            Some(best) => interval_arguments(parser, best)?.map(|parameters| parameters.intervals),
            None if intervals_specified(parser) => return Err(Thrown::user(
                "We require at least one input source that has a sequence dictionary (reference \
                 or reads) when intervals are specified",
            )),
            None => None,
        };
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let (Some(reference), Some(reads)) = (&reference, &reads) {
            validate_reference_against_reads(reference, reads)?;
        }
    }

    let variant = argument(parser, "variant").ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;
    if intervals.is_some() {
        let index = if variant.ends_with(".gz") || variant.ends_with(".bgz") {
            format!("{variant}.tbi")
        } else {
            format!("{variant}.idx")
        };
        if !std::path::Path::new(&index).exists() {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!("Intervals set but no index file found for {variant} at {index}"),
            ));
        }
    }
    let (codec, text) = open_feature_input(&variant)?;
    let features: Vec<Locus> = gatk_tools::feature_codec::features(&text, codec)
        .into_iter()
        .map(|(feature, _)| Locus {
            contig: feature.contig,
            start: feature.start,
            stop: feature.end,
        })
        .collect();
    let count = gatk_tools::count_variants::count(&features, intervals.as_deref(), true, &variant)
        .map_err(|error| Thrown {
            failure: Failure::User,
            exception: error.class(),
            message: Some(error.message()),
        })?;
    gatk_tools::tool_stdout::println(&count.to_string());
    spark_output(parser, &count.to_string())?;
    Ok(None)
}

/// `PileupSpark.processAlignments`: `Pileup`'s line per locus, without its deletion filter, saved
/// with `saveAsTextFile`.
///
/// `LocusWalkerSpark.getAlignments` cuts each interval (the whole dictionary when there is none)
/// into shards of `--read-shard-size` bases and runs a fresh locus iterator per shard, so a
/// downsampler starts over at every shard boundary. Without `--shuffle` the reads are joined to
/// the shards through `SparkSharder.joinOverlapping`, which fails inside a task when the filtered
/// reads are empty: every partition is, so the shards are given partition `-1`.
pub fn pileup_spark(parser: &Parser) -> Outcome {
    if !arguments(parser, "metadata").is_empty() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "PileupSpark's --metadata is a feature source this port does not open yet, and a run \
             that ignored it would answer without the column it asks for. This message is the \
             port's own and not GATK's.",
        ));
    }
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "PileupSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let mut reference = match argument(parser, "reference") {
        Some(path) => Some(
            gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
        ),
        None => None,
    };
    let dictionary = match reference.as_mut() {
        Some(source) => source.sequences().to_vec(),
        None => header
            .sequences
            .iter()
            .map(|sequence| (sequence.name.clone(), sequence.length as usize))
            .collect(),
    };

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    if records.is_empty() && !flag(parser, "shuffle") {
        eprintln!("ERROR Executor - Exception in task 0.0 in stage 1.0 (TID 1)");
        return Err(Thrown::non_user(
            "java.lang.ArrayIndexOutOfBoundsException",
            "Index -1 out of bounds for length 1",
        ));
    }

    let shard_size = number_or(parser, "read-shard-size", 10000).max(1);
    let whole: Vec<gatk_engine::interval::SimpleInterval> = if intervals.is_empty() {
        dictionary
            .iter()
            .map(|(name, length)| gatk_engine::interval::SimpleInterval {
                contig: name.clone(),
                start: 1,
                end: *length as i32,
            })
            .collect()
    } else {
        intervals.clone()
    };
    // `Shard.divideIntervalIntoShards` with no padding.
    let shards = whole.iter().flat_map(|interval| {
        (interval.start..=interval.end)
            .step_by(shard_size as usize)
            .map(move |start| gatk_engine::interval::SimpleInterval {
                contig: interval.contig.clone(),
                start,
                end: (start + shard_size - 1).min(interval.end),
            })
    });

    let mut bases: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    if let Some(source) = reference.as_mut() {
        for (name, length) in source.sequences().to_vec() {
            let contig = source
                .query(&name, 1, length as i32)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            bases.insert(name, contig);
        }
    }
    let output_insert_length = flag(parser, "output-insert-length");
    let show_verbose = flag(parser, "show-verbose");
    let mut text = String::new();
    for shard in shards {
        let applied = gatk_tools::locus_walker::traverse(
            &records,
            &header,
            None,
            Some(std::slice::from_ref(&shard)),
            gatk_tools::locus_walker::Options {
                max_depth_per_sample: number_or(parser, "max-depth-per-sample", 0),
                ..gatk_tools::locus_walker::Options::default()
            },
            &|_| true,
        )
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
        for one in &applied {
            let base = bases
                .get(&one.context.contig)
                .and_then(|contig| contig.get((one.context.position - 1) as usize))
                .map(|base| *base as char)
                .unwrap_or('N');
            // `saveAsTextFile` ends every record with a newline of its own, after the tool's.
            text.push_str(&gatk_tools::pileup::line_unfiltered(
                &one.context.pileup,
                base,
                &[],
                output_insert_length,
                show_verbose,
            ));
            text.push('\n');
        }
    }

    let directory = std::path::Path::new(&output);
    std::fs::create_dir_all(directory)
        .and_then(|()| std::fs::write(directory.join("part-00000"), &text))
        .and_then(|()| std::fs::write(directory.join("_SUCCESS"), b""))
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

/// `PrintReadsSpark.runTool`: the filtered reads, unsorted, through `ReadsSparkSink` and Disq.
///
/// The header reaches the writer through Kryo, which rebuilds `@HD` with a group order: `GO:none`
/// when the input named none, placed before `SO`. Disq then writes one of two layouts:
///
/// * `--sharded-output`: a directory holding `part-r-00000.bam`, header and records in one BGZF
///   stream and its terminator, beside an empty `_SUCCESS`. No index is written for a directory.
/// * otherwise one file assembled from parts: the header compressed ALONE, then the records'
///   own stream, then a single terminator; `.bai` and `.sbi` are named by APPENDING to the output.
///
/// Both deflate at GATKConfig's level with the deflater `--use-jdk-deflater` chooses.
pub fn print_reads_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "PrintReadsSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let mut header = header.clone();
    let sort_order = header.attributes.get("SO").map(str::to_string);
    let group_order = header.attributes.get("GO").unwrap_or("none").to_string();
    header.attributes.remove("SO");
    header.attributes.remove("GO");
    header.set_group_order(&group_order);
    if let Some(order) = sort_order {
        header.set_sort_order(&order);
    }

    let (level, deflater) = output_compression(parser);
    let failed = |what: &str, error: String| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {what}: {error}"))
    };
    if flag(parser, "sharded-output") {
        let mut writer =
            htsjdk_bam::writer::BamWriter::with_compression(Vec::new(), &header, level, deflater)
                .map_err(|error| failed(&output, error.to_string()))?;
        for record in &records {
            writer
                .write(record)
                .map_err(|error| failed(&output, error.to_string()))?;
        }
        let part = writer
            .finish()
            .map_err(|error| failed(&output, error.to_string()))?;
        let directory = std::path::Path::new(&output);
        std::fs::create_dir_all(directory)
            .and_then(|()| std::fs::write(directory.join("part-r-00000.bam"), &part))
            .and_then(|()| std::fs::write(directory.join("_SUCCESS"), b""))
            .map_err(|error| failed(&output, error.to_string()))?;
        return Ok(None);
    }

    // The header part: a header-only BAM without its terminator.
    let mut header_part =
        htsjdk_bam::writer::BamWriter::with_compression(Vec::new(), &header, level, deflater)
            .and_then(|writer| writer.finish())
            .map_err(|error| failed(&output, error.to_string()))?;
    header_part.truncate(header_part.len() - htsjdk_bgzf::EMPTY_GZIP_BLOCK.len());

    // The records' part, whose offsets the splitting index records relative to the whole file.
    let shift = (header_part.len() as u64) << 16;
    let mut stream = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
    let mut sbi = htsjdk_bam::sbi::SbiIndexWriter::new(
        number_or(parser, "splitting-index-granularity", 4096).max(1) as u64,
    );
    let mut encoded = Vec::new();
    for record in &records {
        encoded.clear();
        record
            .encode_into(&mut encoded)
            .map_err(|error| failed(&output, format!("{error:?}")))?;
        sbi.process_record(shift + stream.file_pointer())
            .map_err(|error| failed(&output, error))?;
        std::io::Write::write_all(&mut stream, &encoded)
            .map_err(|error| failed(&output, error.to_string()))?;
    }
    let records_part = stream
        .into_inner_without_terminator()
        .map_err(|error| failed(&output, error.to_string()))?;

    let mut bam = header_part;
    bam.extend_from_slice(&records_part);
    let terminator_offset = (bam.len() as u64) << 16;
    bam.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    std::fs::write(&output, &bam).map_err(|error| failed(&output, error.to_string()))?;

    if flag(parser, "create-output-bam-index") {
        let bai = htsjdk_bam::build_index::build_bam_index(&bam)
            .map_err(|error| failed(&format!("{output}.bai"), format!("{error:?}")))?;
        std::fs::write(format!("{output}.bai"), bai)
            .map_err(|error| failed(&format!("{output}.bai"), error.to_string()))?;
    }
    if flag(parser, "create-output-bam-splitting-index") {
        let index = sbi
            .finish(terminator_offset, bam.len() as u64, None, None)
            .map_err(|error| failed(&format!("{output}.sbi"), error))?;
        std::fs::write(format!("{output}.sbi"), index)
            .map_err(|error| failed(&format!("{output}.sbi"), error.to_string()))?;
    }
    Ok(None)
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
