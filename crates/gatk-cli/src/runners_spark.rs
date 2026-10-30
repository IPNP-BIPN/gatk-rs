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
    // `writeReads(ctx, output, reads, header, false)`: PrintReadsSpark asks for no sort.
    write_spark_reads(parser, &header, &records, &output, false)?;
    Ok(None)
}

/// `GATKSparkTool.writeReads` through `ReadsSparkSink` and Disq, as described above.
///
/// With `sort_to_header`, a coordinate-sorted header is sorted to first (`sortSamRecordsToMatchHeader`), over
/// `getRecommendedNumReducers()` partitions, and Disq writes ONE STREAM PER PARTITION: a part file
/// each for `--sharded-output`, and for a single file each partition's records compressed on
/// their own between the header and the terminator, the `.sbi` merged part by part.
fn write_spark_reads(
    parser: &Parser,
    header: &SamHeader,
    records: &[BamRecord],
    output: &str,
    sort_to_header: bool,
) -> Result<(), Thrown> {
    let mut header = header.clone();
    let sort_order = header.attributes.get("SO").map(str::to_string);
    let group_order = header.attributes.get("GO").unwrap_or("none").to_string();
    header.attributes.remove("SO");
    header.attributes.remove("GO");
    header.set_group_order(&group_order);
    if let Some(order) = &sort_order {
        header.set_sort_order(order);
    }
    let partitions = if sort_to_header {
        sorted_partitions(parser, sort_order.as_deref(), records)?
    } else {
        vec![records.to_vec()]
    };

    let (level, deflater) = output_compression(parser);
    let failed = |what: &str, error: String| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {what}: {error}"))
    };
    if flag(parser, "sharded-output") {
        let directory = std::path::Path::new(output);
        std::fs::create_dir_all(directory).map_err(|error| failed(output, error.to_string()))?;
        for (index, partition) in partitions.iter().enumerate() {
            let mut writer = htsjdk_bam::writer::BamWriter::with_compression(
                Vec::new(),
                &header,
                level,
                deflater,
            )
            .map_err(|error| failed(output, error.to_string()))?;
            for record in partition {
                writer
                    .write(record)
                    .map_err(|error| failed(output, error.to_string()))?;
            }
            let part = writer
                .finish()
                .map_err(|error| failed(output, error.to_string()))?;
            std::fs::write(directory.join(format!("part-r-{index:05}.bam")), &part)
                .map_err(|error| failed(output, error.to_string()))?;
        }
        std::fs::write(directory.join("_SUCCESS"), b"")
            .map_err(|error| failed(output, error.to_string()))?;
        return Ok(());
    }

    // The header part: a header-only BAM without its terminator.
    let mut header_part =
        htsjdk_bam::writer::BamWriter::with_compression(Vec::new(), &header, level, deflater)
            .and_then(|writer| writer.finish())
            .map_err(|error| failed(output, error.to_string()))?;
    header_part.truncate(header_part.len() - htsjdk_bgzf::EMPTY_GZIP_BLOCK.len());

    // Each partition's records, whose offsets the splitting index records relative to the whole
    // file. `SBIIndexMerger` keeps every part's own offsets, counted from that part's first
    // record, shifted to where the part starts.
    let granularity = number_or(parser, "splitting-index-granularity", 4096).max(1) as u64;
    let mut bam = header_part;
    let mut offsets: Vec<u64> = Vec::new();
    let mut encoded = Vec::new();
    for partition in &partitions {
        let shift = (bam.len() as u64) << 16;
        let mut stream = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
        for (index, record) in partition.iter().enumerate() {
            encoded.clear();
            record
                .encode_into(&mut encoded)
                .map_err(|error| failed(output, format!("{error:?}")))?;
            if (index as u64).is_multiple_of(granularity) {
                offsets.push(shift + stream.file_pointer());
            }
            std::io::Write::write_all(&mut stream, &encoded)
                .map_err(|error| failed(output, error.to_string()))?;
        }
        let part = stream
            .into_inner_without_terminator()
            .map_err(|error| failed(output, error.to_string()))?;
        bam.extend_from_slice(&part);
    }
    let terminator_offset = (bam.len() as u64) << 16;
    bam.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    std::fs::write(output, &bam).map_err(|error| failed(output, error.to_string()))?;

    if flag(parser, "create-output-bam-index") && sort_order.as_deref() == Some("coordinate") {
        // shortcut: the index is built from the whole file; Disq merges one per part with
        // `BAMIndexMerger`, which is the same index on every array row measured so far. Port the
        // merger if a multi-partition row ever disagrees.
        let bai = htsjdk_bam::build_index::build_bam_index(&bam)
            .map_err(|error| failed(&format!("{output}.bai"), format!("{error:?}")))?;
        std::fs::write(format!("{output}.bai"), bai)
            .map_err(|error| failed(&format!("{output}.bai"), error.to_string()))?;
    }
    if flag(parser, "create-output-bam-splitting-index") {
        let mut sbi = htsjdk_bam::sbi::SbiIndexWriter::new(1);
        for offset in &offsets {
            sbi.process_record(*offset)
                .map_err(|error| failed(&format!("{output}.sbi"), error))?;
        }
        let mut index = sbi
            .finish(terminator_offset, bam.len() as u64, None, None)
            .map_err(|error| failed(&format!("{output}.sbi"), error))?;
        // The merged header carries every record and the parts' granularity, where a writer fed
        // only the kept offsets counted those at a granularity of one.
        let records: usize = partitions.iter().map(Vec::len).sum();
        index[44..52].copy_from_slice(&(records as u64).to_le_bytes());
        index[52..60].copy_from_slice(&granularity.to_le_bytes());
        std::fs::write(format!("{output}.sbi"), index)
            .map_err(|error| failed(&format!("{output}.sbi"), error.to_string()))?;
    }
    Ok(())
}

/// `SparkUtils.sortUsingElementsAsKeys` for a coordinate-sorted header: `sortByKey` with
/// `HeaderlessSAMRecordCoordinateComparator` over a `RangePartitioner` of
/// `getRecommendedNumReducers()` partitions. Any other order is written as it was read.
///
/// The sort is stable, as Spark's `TimSort` is, so equal keys keep the traversal's order.
fn sorted_partitions(
    parser: &Parser,
    sort_order: Option<&str>,
    records: &[BamRecord],
) -> Result<Vec<Vec<BamRecord>>, Thrown> {
    if sort_order != Some("coordinate") {
        return Ok(vec![records.to_vec()]);
    }
    let mut sorted = records.to_vec();
    sorted.sort_by(htsjdk_bam::coordinate::compare);

    // `getRecommendedNumReducers`: the argument, else one per ten megabytes of input.
    let mut reducers = number_or(parser, "num-reducers", 0).max(0) as usize;
    if reducers == 0 {
        let size = argument(parser, "input")
            .and_then(|input| std::fs::metadata(input).ok())
            .map(|data| data.len())
            .unwrap_or(0);
        reducers = 1 + (size / (10 * 1024 * 1024)) as usize;
    }
    let bounds = range_bounds(&sorted, reducers)?;
    let mut partitions: Vec<Vec<BamRecord>> = vec![Vec::new(); bounds.len() + 1];
    for record in sorted {
        // `getPartition`: the first bound the key is not greater than.
        let partition = bounds
            .iter()
            .position(|bound| htsjdk_bam::coordinate::compare(&record, bound).is_le())
            .unwrap_or(bounds.len());
        partitions[partition].push(record);
    }
    Ok(partitions)
}

/// `RangePartitioner`'s bounds over one input partition, the sketch taking every key: weight one
/// each, then `determineBounds(candidates, min(partitions, candidates))`.
///
/// shortcut: one input partition holding at most `sampleSizePerPartition` keys, where the sketch is
/// the whole partition; more is a reservoir sample the port does not draw, and is refused by name.
fn range_bounds(sorted: &[BamRecord], partitions: usize) -> Result<Vec<BamRecord>, Thrown> {
    if partitions <= 1 || sorted.is_empty() {
        return Ok(Vec::new());
    }
    let sample_size = (20.0 * partitions as f64).min(1e6);
    let per_partition = (3.0 * sample_size).ceil() as usize;
    if sorted.len() > per_partition {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "Sorting {} reads into {partitions} partitions samples them, which this port does \
                 not reproduce yet. This message is the port's own and not GATK's.",
                sorted.len()
            ),
        ));
    }
    let partitions = partitions.min(sorted.len());
    let step = sorted.len() as f64 / partitions as f64;
    let mut cumulative = 0.0f64;
    let mut target = step;
    let mut bounds: Vec<BamRecord> = Vec::new();
    for key in sorted {
        if bounds.len() >= partitions - 1 {
            break;
        }
        cumulative += 1.0;
        if cumulative >= target
            && bounds
                .last()
                .is_none_or(|previous| htsjdk_bam::coordinate::compare(key, previous).is_gt())
        {
            bounds.push(key.clone());
            target += step;
        }
    }
    Ok(bounds)
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

/// `BaseRecalibratorSpark.runTool`: `BaseRecalibratorSparkFn` over the reads joined with the
/// known sites overlapping each one, written with `RecalUtils.outputRecalibrationReport`.
///
/// One partition is the walker's engine over the walker's reads, so the table is
/// `BaseRecalibrator`'s; what differs is the startup and the return. `runTool` returns nothing,
/// so there is no `Tool returned:` line where the walker prints `SUCCESS`.
pub fn base_recalibrator_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "BaseRecalibratorSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    spark_known_sites(parser)?;
    let table = super::recalibration_table(parser, &source, &header, &intervals, &filters)?;
    tables_combine(&header)?;
    std::fs::write(&output, table)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

/// `ApplyBQSRSpark.runTool`: the filtered reads through one `BQSRReadTransformer`, then
/// `writeReads`, which is PrintReadsSpark's writer.
pub fn apply_bqsr_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "ApplyBQSRSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    // The report is opened in `runTool`, after the reads are defined and before any is read.
    let recal = argument(parser, "bqsr-recal-file").unwrap_or_default();
    let recal_text = std::fs::read_to_string(&recal)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&recal, false)))?;
    let recalibrated = recalibrated(
        &header,
        &recal_text,
        &super::apply_bqsr_arguments(parser),
        &records,
    )?;
    write_spark_reads(parser, &header, &recalibrated, &output, true)?;
    Ok(None)
}

/// `BQSRPipelineSpark.runTool`: both passes in one job, the table never written.
///
/// The table is built from the initial reads filtered AGAIN by the six BQSR-specific filters,
/// which do not honour the command line; the table is applied to the initial reads, which carry
/// only the wellformed filter and whatever the command line added.
pub fn bqsr_pipeline_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "BQSRPipelineSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let mut for_table = filters.clone();
    for name in [
        "MappingQualityNotZeroReadFilter",
        "MappingQualityAvailableReadFilter",
        "MappedReadFilter",
        "NotSecondaryAlignmentReadFilter",
        "NotDuplicateReadFilter",
        "PassesVendorQualityCheckReadFilter",
    ] {
        for_table.push(gatk_tools::filter_resolution::ResolvedFilter {
            name: name.to_string(),
            negated: false,
        });
    }
    spark_known_sites(parser)?;
    let table = super::recalibration_table(parser, &source, &header, &intervals, &for_table)?;
    tables_combine(&header)?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    // `ApplyBQSRUniqueArgumentCollection.toApplyBQSRArgumentCollection(bqsrArgs)`: the
    // preserve threshold and the original-qualities switch are the table's own arguments.
    let recalibrated = recalibrated(
        &header,
        &table,
        &super::apply_bqsr_arguments(parser),
        &records,
    )?;
    write_spark_reads(parser, &header, &recalibrated, &output, true)?;
    Ok(None)
}

/// `addVCFsForSpark`: a known-sites file is shipped to the executors by its extension, so a BED,
/// which the walker reads, is refused before a read is.
fn spark_known_sites(parser: &Parser) -> Result<(), Thrown> {
    for path in arguments(parser, "known-sites") {
        if !(path.ends_with(".vcf") || path.ends_with(".vcf.gz")) {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                "Unrecognized known sites file extension. Must be .vcf or .vcf.gz",
            ));
        }
    }
    Ok(())
}

/// `treeAggregate`'s first combine: the empty table is sized by the read groups' DISTINCT
/// identifiers (the PU, else the ID), while `BaseRecalibrationEngine` sizes the partition's by
/// `getReadGroups().size()`. Two read groups sharing a PU make the two disagree, and
/// `RecalUtils.combineTables` refuses inside the executor whatever the reads were.
fn tables_combine(header: &SamHeader) -> Result<(), Thrown> {
    let mut identifiers: Vec<&str> = Vec::new();
    for group in &header.read_groups {
        let identifier = group.attributes.get("PU").unwrap_or(&group.id);
        if !identifiers.contains(&identifier) {
            identifiers.push(identifier);
        }
    }
    if identifiers.len() == header.read_groups.len() {
        return Ok(());
    }
    eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
    Err(Thrown::non_user(
        "java.lang.IllegalArgumentException",
        format!(
            "Table1 {},3 not equal to {},3",
            identifiers.len(),
            header.read_groups.len()
        ),
    ))
}

/// The transformer's refusals keep their class: a read group the table lacks is a
/// `GATKException`, exit 3.
fn recalibrated(
    header: &SamHeader,
    recal_text: &str,
    arguments: &gatk_engine::bqsr_transformer::ApplyBqsrArguments,
    records: &[BamRecord],
) -> Result<Vec<BamRecord>, Thrown> {
    gatk_tools::apply_bqsr::recalibrate(header, recal_text, arguments, records).map_err(|error| {
        Thrown {
            failure: if error.is_user() {
                Failure::User
            } else {
                Failure::Other
            },
            exception: error.java_class(),
            message: Some(error.message()),
        }
    })
}
