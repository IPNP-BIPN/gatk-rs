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

/// `LibraryIdGenerator.UNKNOWN_LIBRARY`.
const UNKNOWN_LIBRARY: &str = "Unknown Library";

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

    let source = spark_traversal_source(parser, &input, source)?;

    Ok(SparkStart {
        source,
        header,
        intervals,
        filters,
    })
}

/// What `getUnfilteredReads` traverses: the BAM with its index under intervals, from Disq's first
/// record.
fn spark_traversal_source(
    parser: &Parser,
    input: &str,
    source: ReadsDataSource,
) -> Result<ReadsDataSource, Thrown> {
    // `getUnfilteredReads` with intervals is a Disq traversal, which needs the index beside the
    // BAM. Without one the reference fails inside an executor task, and that failure is Spark's
    // log rather than a GATK refusal; the port says so in its own words.
    let source = if intervals_specified(parser) {
        let path = std::path::Path::new(input);
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
    source
        .with_disq_split_start()
        .map_err(|error| Thrown::user(format!("{error:?}")))
}

/// `CountVariantsSpark.runTool`: Disq's variants over `-L`, counted and printed like the reads.
///
/// The tool replaces the plugin descriptors with none, so it has no read filter to resolve, and it
/// requires no reads: `-I` is opened only when given, for its dictionary. Disq queries a VCF over
/// intervals through its index and refuses one that has none, from inside the job, which surfaces
/// as the `IllegalArgumentException` itself rather than a GATK error.
pub fn count_variants_spark(parser: &Parser) -> Outcome {
    let (variant, intervals) = variant_spark_startup(parser)?;
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

/// `GATKSparkTool`'s startup for a tool driven by `-V`, up to the variants: the context, the
/// reference, the optional reads for their dictionary, the intervals against the best dictionary,
/// and Disq's refusal of an interval query over a VCF with no index.
///
/// Returns the `-V` path and the traversal intervals.
fn variant_spark_startup(
    parser: &Parser,
) -> Result<(String, Option<Vec<gatk_engine::interval::SimpleInterval>>), Thrown> {
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
    Ok((variant, intervals))
}

/// `PrintVariantsSpark.processVariants`: the variants `-L` selects, unsorted, through
/// `VariantsSparkSink` and Disq.
///
/// The header is the input's as htsjdk reads it, with no command line added. Disq refuses BCF by
/// name before it writes anything, writes a plain `.vcf` as text with NO index whatever
/// `--create-output-variant-index` says, and writes a block-compressed name as the header
/// compressed alone, the records' own stream and one terminator, with a `.tbi` APPENDED when
/// indexing is on.
pub fn print_variants_spark(parser: &Parser) -> Outcome {
    let (variant, intervals) = variant_spark_startup(parser)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // Block compressed or not, as `CountVariantsSpark` opens it.
    let (_, text) = open_feature_input(&variant)?;
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let loci: Vec<Locus> = file
        .records
        .iter()
        .map(|record| Locus {
            contig: record.contig.clone(),
            start: record.start as i32,
            stop: record.stop as i32,
        })
        .collect();
    let kept: Vec<htsjdk_vcf::variant::VariantContext> =
        gatk_engine::variant_source::traverse(&loci, intervals.as_deref())
            .into_iter()
            .map(|locus| {
                // `traverse` hands back the loci themselves; their place is the record's.
                let index = (locus as *const Locus as usize - loci.as_ptr() as usize)
                    / std::mem::size_of::<Locus>();
                file.records[index].clone()
            })
            .collect();

    if output.ends_with(".bcf") || output.ends_with(".bcf.gz") {
        return Err(Thrown {
            failure: Failure::User,
            exception:
                "org.broadinstitute.hellbender.exceptions.UserException$UnimplementedFeature",
            message: Some(
                "It is currently not possible to write a BCF file on spark.  See \
                 https://github.com/broadinstitute/gatk/issues/4303 for more details ."
                    .to_string(),
            ),
        });
    }
    // A record with genotypes never reaches the writer: Kryo rebuilds its `LazyGenotypesContext`
    // without the parser, and the first `add` decodes through it. Without `--shuffle` that happens
    // where the driver collects the task's result; with it, inside the executor's join. Measured
    // over the whole covering array, the stage follows `--shuffle` alone.
    if !file.header.samples.is_empty() && !kept.is_empty() {
        const KRYO: &str = "com.esotericsoftware.kryo.KryoException: \
             java.lang.NullPointerException: Cannot invoke \
             \"htsjdk.variant.variantcontext.LazyGenotypesContext$LazyParser.parse(Object)\" \
             because \"this.parser\" is null\nSerialization trace:\n\
             genotypes (htsjdk.variant.variantcontext.VariantContext)\n";
        // shortcut: the executor's message names the lost task and its host, which no two runs
        // share; the stage's log line is what the covering array compares, and it is exact.
        let message = if flag(parser, "shuffle") {
            eprintln!("ERROR Executor - Exception in task 0.0 in stage 1.0 (TID 1)");
            format!("Job aborted due to stage failure: Task 0 in stage 1.0 failed 1 times: {KRYO}")
        } else {
            eprintln!("ERROR TaskResultGetter - Exception while getting task result");
            format!(
                "Job aborted due to stage failure: Exception while getting task result: {KRYO}\
                 interval (org.broadinstitute.hellbender.engine.spark.SparkSharder$PartitionLocatable)"
            )
        };
        return Err(Thrown::non_user("org.apache.spark.SparkException", message));
    }
    let failed = |error: String| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    };
    let header_text = file.header.write();
    let mut records_text = String::new();
    let encoder = htsjdk_vcf::encoder::VcfEncoder::new(&file.header);
    for record in &kept {
        encoder
            .encode_into(record, &mut records_text)
            .map_err(|error| failed(format!("{error:?}")))?;
        records_text.push('\n');
    }

    let block_compressed = output.ends_with(".gz") || output.ends_with(".bgz");
    if !block_compressed {
        std::fs::write(&output, format!("{header_text}{records_text}"))
            .map_err(|error| failed(error.to_string()))?;
        return Ok(None);
    }
    let (level, deflater) = output_compression(parser);
    let part = |text: &str| -> Result<Vec<u8>, Thrown> {
        let mut stream = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
        std::io::Write::write_all(&mut stream, text.as_bytes())
            .map_err(|error| failed(error.to_string()))?;
        stream
            .into_inner_without_terminator()
            .map_err(|error| failed(error.to_string()))
    };
    let mut bytes = part(&header_text)?;
    bytes.extend_from_slice(&part(&records_text)?);
    bytes.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    std::fs::write(&output, &bytes).map_err(|error| failed(error.to_string()))?;

    if flag(parser, "create-output-variant-index") {
        // Disq finalizes the index at the position the records end, BEFORE the terminator it
        // appends, where `IndexFeatureFile` reads the whole file and ends the last chunk after it.
        let source = index_feature_file::Source::new(&output);
        let records_end = bytes.len() - htsjdk_bgzf::EMPTY_GZIP_BLOCK.len();
        let index = index_feature_file::build_tabix(
            &bytes[..records_end],
            &source,
            &output,
            deflater,
            level,
        )
        .map_err(|refusal| Thrown::user(refusal.message()))?;
        let name = format!("{output}.tbi");
        std::fs::write(&name, index).map_err(|error| failed(error.to_string()))?;
    }
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

/// `SortSamSpark.runTool`: every read, the header's `SO` set to `--sort-order`, then `writeReads`
/// asked to sort to the header.
///
/// The tool sorts nothing itself. `ReadsSparkSink` sorts the `SAMRecord`s with the comparator the
/// header's order names, `HeaderlessSAMRecordCoordinateComparator` or htsjdk's
/// `SAMRecordQueryNameComparator`, and the reads keep `ALLOW_ALL_READS`, so nothing is dropped.
pub fn sort_sam_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "SortSamSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    // `SparkSortOrder`'s two constants, which Barclay has already checked; an enum, so read as a
    // scalar rather than as text.
    let order = scalar(parser, "sort-order").unwrap_or_else(|| "coordinate".to_string());
    let mut sorted = header.clone();
    sorted.set_sort_order(&order);
    write_spark_reads(parser, &sorted, &records, &output, true)?;
    Ok(None)
}

/// `getOutputMap`: where each read group goes, or the one file for all of them under `None`.
///
/// `customCommandLineValidation` runs `validateOutputParams` over the same arguments and never looks
/// at the list it returns, so a bad combination only fails when the tool reaches this.
fn revert_output_files(
    parser: &Parser,
    groups: &[String],
    extension: &str,
) -> Result<Vec<(Option<String>, String)>, Thrown> {
    let output = argument(parser, "output");
    if !flag(parser, "output-by-readgroup") {
        // `IOUtils.getPath(null)`, on a command line that named only `--output-map`.
        let output = output.ok_or_else(|| {
            Thrown::non_user(
                "java.lang.IllegalArgumentException",
                "Null object is not allowed here.",
            )
        })?;
        return Ok(vec![(None, output)]);
    }
    let Some(map) = argument(parser, "output-map") else {
        // Barclay's mutex guarantees one of the two, and `Paths.get(outputDir, fileName)` joins them.
        let directory = output.unwrap_or_default();
        return Ok(groups
            .iter()
            .map(|id| {
                let path = std::path::Path::new(&directory).join(format!("{id}{extension}"));
                (Some(id.clone()), path.display().to_string())
            })
            .collect());
    };

    // `createOutputMapFromFile` through `TableCodec(null)`: the first line names the columns, any
    // white space separates them, and `#` and `track` lines are no rows.
    let text = std::fs::read_to_string(&map).map_err(|_| {
        Thrown::non_user(
            "htsjdk.tribble.TribbleException$MalformedFeatureFile",
            format!(
                "Unable to create BasicFeatureReader using feature file , for input source: {map}"
            ),
        )
    })?;
    let mut lines = text.lines();
    let columns: Vec<&str> = lines.next().unwrap_or("").split_whitespace().collect();
    let column = |name: &str| {
        columns
            .iter()
            .position(|column| *column == name)
            .ok_or_else(|| {
                Thrown::non_user(
                    "java.lang.IllegalArgumentException",
                    format!("We don't have a column named {name}"),
                )
            })
    };
    let mut files: Vec<(Option<String>, String)> = Vec::new();
    for line in lines.filter(|line| !line.starts_with('#') && !line.starts_with("track")) {
        let (id_at, output_at) = (column("READ_GROUP_ID")?, column("OUTPUT")?);
        let fields: Vec<&str> = line.split_whitespace().collect();
        let (Some(id), Some(file)) = (fields.get(id_at), fields.get(output_at)) else {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "{map}: a row short of its columns fails the reference's TableFeature, which \
                     this port does not reproduce. This message is the port's own and not GATK's."
                ),
            ));
        };
        // A `HashMap`, so a later row of the same group wins.
        files.retain(|(group, _)| group.as_deref() != Some(*id));
        files.push((Some(id.to_string()), file.to_string()));
    }
    Ok(files)
}

/// `new SAMFileHeader()` with what `createOutHeader` carries over, and the read groups it is given.
fn revert_out_header(
    header: &SamHeader,
    sort_order: &str,
    groups: &[htsjdk_bam::header::ReadGroup],
    keep_alignment: bool,
) -> SamHeader {
    let mut out = SamHeader::default();
    out.set_sort_order(sort_order);
    if keep_alignment {
        out.sequences = header.sequences.clone();
        out.programs = header.programs.clone();
    }
    out.read_groups = groups.to_vec();
    out
}

fn read_group_of(record: &BamRecord) -> Option<&str> {
    match record.tags.get(htsjdk_bam::tag::Tag::new(b"RG")) {
        Some(htsjdk_bam::tag::TagValue::Str(id)) => Some(id),
        _ => None,
    }
}

/// `RevertSamSpark.revertReads`, one record.
///
/// Not Picard's `RevertSam`: `setMateIsUnplaced` marks every read paired on its way, and
/// supplementary reads are neither dropped nor unflagged.
fn revert_record(
    record: &mut BamRecord,
    restore_qualities: bool,
    keep_duplicates: bool,
    keep_alignment: bool,
    attributes: &[String],
) -> Result<(), Thrown> {
    use htsjdk_bam::tag::{Tag, TagValue};
    const PAIRED: u16 = 0x1;
    const PROPER_PAIR: u16 = 0x2;
    const UNMAPPED: u16 = 0x4;
    const MATE_UNMAPPED: u16 = 0x8;
    const REVERSE: u16 = 0x10;
    const MATE_REVERSE: u16 = 0x20;
    const SECONDARY: u16 = 0x100;
    const DUPLICATE: u16 = 0x400;

    if restore_qualities {
        // `getOriginalBaseQualities`: an empty `OQ` is none, and stays on the read.
        if let Some(TagValue::Str(original)) = record.tags.get(Tag::new(b"OQ")) {
            if !original.is_empty() {
                record.base_qualities = htsjdk_bam::fastq::fastq_to_phred(original);
                record.tags.remove(Tag::new(b"OQ"));
            }
        }
    }
    if !keep_duplicates {
        record.flags &= !DUPLICATE;
    }
    if keep_alignment {
        return Ok(());
    }
    if record.flags & REVERSE != 0 {
        // `SAMRecord.reverseComplement(true)`, whose default lists name these four tags.
        htsjdk_bam::sequence::reverse_complement(&mut record.read_bases);
        htsjdk_bam::sequence::reverse_qualities(&mut record.base_qualities);
        for name in [b"E2", b"SQ"] {
            if let Some(TagValue::Str(text)) = record.tags.get(Tag::new(name)) {
                let mut bytes = text.clone().into_bytes();
                htsjdk_bam::sequence::reverse_complement(&mut bytes);
                let text = String::from_utf8_lossy(&bytes).into_owned();
                record.tags.insert(Tag::new(name), TagValue::Str(text));
            }
        }
        for name in [b"OQ", b"U2"] {
            if let Some(TagValue::Str(text)) = record.tags.get(Tag::new(name)) {
                let text: String = text.chars().rev().collect();
                record.tags.insert(Tag::new(name), TagValue::Str(text));
            }
        }
        record.flags &= !REVERSE;
    }
    record.flags |= UNMAPPED | PAIRED | MATE_UNMAPPED;
    record.flags &= !(SECONDARY | PROPER_PAIR | MATE_REVERSE);
    record.reference_index = -1;
    record.alignment_start = 0;
    record.mapping_quality = 0;
    record.cigar = Default::default();
    record.inferred_insert_size = 0;
    record.mate_reference_index = -1;
    record.mate_alignment_start = 0;
    for name in attributes {
        // `ReadUtils.assertAttributeNameIsLegal`, inside a task.
        let bytes = name.as_bytes();
        let legal =
            bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric();
        if !legal {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "Read attribute {name} invalid: the reference fails this inside a Spark task, \
                     which this port does not reproduce. This message is the port's own and not \
                     GATK's."
                ),
            ));
        }
        record.tags.remove(Tag::new(&[bytes[0], bytes[1]]));
    }
    Ok(())
}

/// `QualityEncodingDetector.detect` over one read group's first 10000 reads, for a SAM file:
/// whether its qualities need the Illumina offset taken off.
///
/// `generateBestGuess` reduces to this because the schemes' ranges cannot leave `Standard` beside
/// another format: it needs an ASCII value below 59, which Solexa and Illumina both exclude. All
/// `I` qualities therefore read as Illumina, and that is what the reference does with them.
fn quality_offset_applies(records: &[&BamRecord]) -> Result<bool, Thrown> {
    if records.is_empty() {
        return Ok(false);
    }
    let (mut lowest, mut highest, mut standard) = (u8::MAX, 0u8, false);
    for record in records {
        // `getBaseQualityString`, where absent qualities are the one character `*`.
        let observed: Vec<u8> = if record.base_qualities.is_empty() {
            vec![b'*']
        } else {
            record.base_qualities.iter().map(|q| q + 33).collect()
        };
        for ascii in observed {
            lowest = lowest.min(ascii);
            highest = highest.max(ascii);
            standard |= (33..=58).contains(&ascii);
        }
    }
    let in_range = |low: u8| lowest >= low && highest <= 126;
    if standard && in_range(33) {
        return Ok(false);
    }
    if in_range(59) {
        return Ok(true);
    }
    Err(Thrown::non_user(
        "htsjdk.samtools.SAMException",
        "The quality format cannot be determined: all formats were excluded.",
    ))
}

/// `RevertSamSpark.sanitize`: the reads whose lengths agree, sorted by name unless the input's
/// header says they are, then each run of one name on a partition judged alone.
fn sanitize_reads(
    parser: &Parser,
    header: &SamHeader,
    records: Vec<BamRecord>,
    keep_first: bool,
) -> Result<Vec<Vec<BamRecord>>, Thrown> {
    const FIRST: u16 = 0x40;
    const SECOND: u16 = 0x80;
    const PAIRED: u16 = 0x1;

    // `createReadGroupFormatMap`, over the reads before any is dropped.
    let mut illumina = std::collections::HashSet::new();
    for group in &header.read_groups {
        let reads: Vec<&BamRecord> = records
            .iter()
            .filter(|record| read_group_of(record) == Some(group.id.as_str()))
            .take(10000)
            .collect();
        if quality_offset_applies(&reads)? {
            illumina.insert(group.id.clone());
        }
    }

    let kept: Vec<BamRecord> = records
        .into_iter()
        .filter(|record| record.read_bases.len() == record.base_qualities.len())
        .collect();
    let name_grouped = header.attributes.get("SO") == Some("queryname")
        || header.attributes.get("GO") == Some("query");
    let partitions = if name_grouped {
        vec![kept]
    } else {
        sorted_partitions(parser, Some("queryname"), &kept)?
    };

    let mut out = Vec::with_capacity(partitions.len());
    for partition in partitions {
        let mut admitted = Vec::new();
        for group in partition.chunk_by(|a, b| a.read_name == b.read_name) {
            let first = group
                .iter()
                .find(|r| r.flags & PAIRED != 0 && r.flags & FIRST != 0);
            let second = group
                .iter()
                .find(|r| r.flags & PAIRED != 0 && r.flags & SECOND != 0);
            let count = |flag: u16| {
                group
                    .iter()
                    .filter(|r| r.flags & PAIRED != 0 && r.flags & flag != 0)
                    .count()
            };
            let (firsts, seconds) = (count(FIRST), count(SECOND));
            let unpaired = group.iter().filter(|r| r.flags & PAIRED == 0).count();

            let chosen: Vec<&BamRecord> = if firsts > 0 || seconds > 0 {
                if firsts == 1 && seconds == 1 {
                    group.iter().collect()
                } else if keep_first && firsts >= 1 && seconds >= 1 {
                    vec![first.expect("counted"), second.expect("counted")]
                } else {
                    continue;
                }
            } else if unpaired > 1 {
                if !keep_first {
                    continue;
                }
                vec![group
                    .iter()
                    .find(|r| r.flags & PAIRED == 0)
                    .expect("counted")]
            } else {
                group.iter().collect()
            };

            for record in chosen {
                let mut record = record.clone();
                if read_group_of(&record).is_some_and(|id| illumina.contains(id)) {
                    // `ILLUMINA_TO_PHRED_SUBTRAHEND`. Solexa's lowest values go negative, which
                    // `setBaseQualities` refuses inside a task.
                    for quality in &mut record.base_qualities {
                        *quality = match quality.checked_sub(31) {
                            Some(lowered) => lowered,
                            None => {
                                let below = i32::from(*quality) - 31;
                                return Err(Thrown::non_user(
                                    "org.apache.spark.SparkException",
                                    format!("Base quality score {below} is invalid"),
                                ));
                            }
                        };
                    }
                }
                admitted.push(record);
            }
        }
        out.push(admitted);
    }
    Ok(out)
}

/// `RevertSamSpark.runTool`.
///
/// Its `--remove-duplicate-information` KEEPS the flags, as the field behind it is named, and
/// without `--remove-default-attributes-to-clear` no default tag is cleared. With it the list is
/// the unmodifiable defaults, and the `addAll` that follows always throws.
pub fn revert_sam_spark(parser: &Parser) -> Outcome {
    // `customCommandLineValidation`, while the command line is parsed.
    let sanitize = flag(parser, "sanitize");
    let keep_first = flag(parser, "keep-first-duplicate");
    if !sanitize && keep_first {
        return Err(Thrown::command_line(
            "Command Line Validation failed:'keepFirstDuplicate' cannot be used without 'sanitize'",
        ));
    }

    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "RevertSamSpark")?;
    let filter = read_filter(parser, &filters, &header)?;

    // `validateHeaderOverrides`, then the overrides on the driver's header.
    let sample = argument(parser, "sample-alias");
    let library = argument(parser, "library-name");
    if sample.is_some() || library.is_some() {
        let groups = &header.read_groups;
        let same = |value: fn(&htsjdk_bam::header::ReadGroup) -> Option<&str>| {
            groups
                .windows(2)
                .all(|pair| value(&pair[0]) == value(&pair[1]))
        };
        if sample.is_some() && !same(|group| group.attributes.get("SM")) {
            return Err(Thrown::user(
                "Read groups have multiple values for sample.  A value for sampleAlias cannot be \
                 supplied.",
            ));
        }
        if library.is_some() && !same(|group| group.attributes.get("LB")) {
            return Err(Thrown::user(
                "Read groups have multiple values for library name.  A value for library name \
                 cannot be supplied.",
            ));
        }
    }
    let mut driver = header.clone();
    for group in &mut driver.read_groups {
        if let Some(sample) = &sample {
            group.attributes.set("SM", sample);
        }
        if let Some(library) = &library {
            group.attributes.set("LB", library);
        }
    }

    // `getDefaultExtension`, by the first input's name.
    let extension = match scalar(parser, "output-by-readgroup-file-format").as_deref() {
        Some("sam") => ".sam",
        Some("bam") => ".bam",
        Some("cram") => ".cram",
        _ => match argument(parser, "input") {
            Some(input) if input.ends_with(".sam") => ".sam",
            Some(input) if input.ends_with(".cram") => ".cram",
            _ => ".bam",
        },
    };
    let ids: Vec<String> = driver.read_groups.iter().map(|g| g.id.clone()).collect();
    let files = revert_output_files(parser, &ids, extension)?;

    let sort_order = scalar(parser, "sort-order").unwrap_or_else(|| "queryname".to_string());
    let keep_alignment = flag(parser, "keep-alignment-information");
    if flag(parser, "output-by-readgroup") && driver.read_groups.is_empty() {
        return Err(Thrown::user("The header is missing its read group map"));
    }
    if flag(parser, "output-by-readgroup") {
        // `assertAllReadGroupsMapped`.
        for id in &ids {
            if !files.iter().any(|(group, _)| group.as_ref() == Some(id)) {
                let map: Vec<String> = files
                    .iter()
                    .map(|(group, file)| format!("{}={file}", group.as_deref().unwrap_or("null")))
                    .collect();
                return Err(Thrown::non_user(
                    "org.broadinstitute.hellbender.exceptions.GATKException",
                    format!(
                        "Read group id {id} not found in outputMap {{{}}}",
                        map.join(", ")
                    ),
                ));
            }
        }
    }
    if flag(parser, "remove-default-attributes-to-clear") {
        return Err(Thrown {
            failure: Failure::Other,
            exception: "java.lang.UnsupportedOperationException",
            message: None,
        });
    }

    let attributes = arguments(parser, "attributes-to-clear");
    let restore = !flag(parser, "dont-restore-original-qualities");
    let keep_duplicates = flag(parser, "remove-duplicate-information");
    let mut records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    for record in &mut records {
        revert_record(
            record,
            restore,
            keep_duplicates,
            keep_alignment,
            &attributes,
        )?;
    }
    // `r.getReadGroup().equals(key)` is the first thing to touch a read when a group is looked up:
    // `--attributes-to-clear RG`, or a read that never had one, fails the first task.
    let looks_up_group =
        (sanitize && !driver.read_groups.is_empty()) || flag(parser, "output-by-readgroup");
    if looks_up_group && records.iter().any(|record| read_group_of(record).is_none()) {
        eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
        return Err(Thrown::non_user(
            "org.apache.spark.SparkException",
            "Job aborted due to stage failure: Task 0 in stage 0.0 failed 1 times: \
             java.lang.NullPointerException: Cannot invoke \"String.equals(Object)\" because the \
             return value of \"org.broadinstitute.hellbender.utils.read.GATKRead.getReadGroup()\" \
             is null",
        ));
    }
    let partitions = if sanitize {
        sanitize_reads(parser, &driver, records, keep_first)?
    } else {
        vec![records]
    };

    for (group, file) in files {
        let Some(id) = group else {
            let out = revert_out_header(&driver, &sort_order, &driver.read_groups, keep_alignment);
            write_reverted(parser, &out, &sort_order, &partitions, &file)?;
            continue;
        };
        let groups: Vec<_> = driver
            .read_groups
            .iter()
            .filter(|group| group.id == id)
            .cloned()
            .collect();
        if groups.is_empty() {
            // `headerMap.get(id)` is null for a group only the map names.
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "The output map names read group {id}, which the header does not: the \
                     reference fails on a null header. This message is the port's own and not \
                     GATK's."
                ),
            ));
        }
        let out = revert_out_header(&driver, &sort_order, &groups, keep_alignment);
        let filtered: Vec<Vec<BamRecord>> = partitions
            .iter()
            .map(|partition| {
                partition
                    .iter()
                    .filter(|record| read_group_of(record) == Some(id.as_str()))
                    .cloned()
                    .collect()
            })
            .collect();
        write_reverted(parser, &out, &sort_order, &filtered, &file)?;
    }
    Ok(None)
}

/// `writeReads(ctx, path, reads, header, false)` for an output whose directory Hadoop would make.
fn write_reverted(
    parser: &Parser,
    header: &SamHeader,
    sort_order: &str,
    partitions: &[Vec<BamRecord>],
    output: &str,
) -> Result<(), Thrown> {
    if let Some(parent) = std::path::Path::new(output).parent() {
        std::fs::create_dir_all(parent).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
        })?;
    }
    // `--sharded-output` is BAM parts under whatever the directory is called.
    let sharded = flag(parser, "sharded-output");
    if output.ends_with(".sam") && !sharded {
        // Disq's SAM writer: the text, whatever the partitions, and no index of any kind.
        let records: Vec<BamRecord> = partitions.iter().flatten().cloned().collect();
        let text = htsjdk_bam::sam_file::write_sam(header, &records)
            .ok_or_else(|| Thrown::non_user(PORT_FAILURE, format!("could not write {output}")))?;
        return std::fs::write(output, text).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
        });
    }
    if !output.ends_with(".bam") && !sharded {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "{output}: Disq's CRAM writer is not ported yet. This message is the port's own \
                 and not GATK's."
            ),
        ));
    }
    write_spark_partitions(parser, header, Some(sort_order), partitions, output)
}

/// `ExtractOriginalAlignmentRecordsByNameSpark.runTool`: `getUnfilteredReads`, so `--read-filter`
/// is parsed and never applied, kept by name and written in the order read.
pub fn extract_original_alignment_records_by_name_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        ..
    } = spark_startup(parser, "ExtractOriginalAlignmentRecordsByNameSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let names_file = argument(parser, "read-name-file").ok_or_else(|| {
        Thrown::command_line(
            "Argument read-name-file was missing: Argument 'read-name-file' is required",
        )
    })?;
    // `parseReadNames`: `^@`, then `/1$`, then `/2$`, each off the line once.
    let text = std::fs::read_to_string(&names_file).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file {names_file}. Error was: {names_file} (No such file or directory)"
        ))
    })?;
    let names: std::collections::HashSet<String> = text
        .lines()
        .map(|line| {
            let line = line.strip_prefix('@').unwrap_or(line);
            let line = line.strip_suffix("/1").unwrap_or(line);
            line.strip_suffix("/2").unwrap_or(line).to_string()
        })
        .collect();
    let invert = flag(parser, "invert-match");
    let keep_all: Filter = Box::new(|_| true);
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &keep_all)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let kept: Vec<BamRecord> = records
        .into_iter()
        .filter(|record| names.contains(&record.read_name) != invert)
        .collect();
    write_spark_reads(parser, &header, &kept, &output, false)?;
    Ok(None)
}

/// `MatchType`, in declaration order: the order `--throw-on-diff` looks for them in.
#[derive(Clone, Copy, PartialEq)]
enum MatchType {
    Equal,
    SizeUnequal,
    ReadMismatch,
    DiffNumDupes,
    DifferentRepresentativeRead,
}

const MATCH_TYPES: [(MatchType, &str); 5] = [
    (MatchType::Equal, "EQUAL"),
    (MatchType::SizeUnequal, "SIZE_UNEQUAL"),
    (MatchType::ReadMismatch, "READ_MISMATCH"),
    (MatchType::DiffNumDupes, "DIFF_NUM_DUPES"),
    (
        MatchType::DifferentRepresentativeRead,
        "DIFFERENT_REPRESENTATIVE_READ",
    ),
];

/// `removeNonReadGroupAttributes`, then the read as the first BAM's header sees it: the tool
/// looks a read's contig up there by name, whichever BAM the read came from.
fn compared_read(mut record: BamRecord, own: &SamHeader, first: &SamHeader) -> BamRecord {
    let group = htsjdk_bam::tag::Tag::new(b"RG");
    let kept = record.tags.get(group).cloned();
    record.tags = htsjdk_bam::tag::Tags::new();
    if let Some(value) = kept {
        record.tags.insert(group, value);
    }
    let reindexed = |index: i32| {
        usize::try_from(index)
            .ok()
            .and_then(|index| own.sequences.get(index))
            .and_then(|contig| first.sequences.iter().position(|s| s.name == contig.name))
            .map_or(gatk_engine::read::NO_ALIGNMENT_REFERENCE_INDEX, |index| {
                index as i32
            })
    };
    record.reference_index = reindexed(record.reference_index);
    record.mate_reference_index = reindexed(record.mate_reference_index);
    record
}

/// `ReadUtils.getStrandedUnclippedStart`.
fn stranded_unclipped_start(record: &BamRecord) -> i32 {
    if gatk_engine::read::is_reverse_strand(record) {
        gatk_engine::read_utils::unclipped_end(record)
    } else {
        gatk_engine::read_utils::unclipped_start(record)
    }
}

/// `ReadsKey.longKeyForFragment`, sign extension included: its mask keeps an int negative when
/// the reference index is, so an unmapped read's key loses its start.
fn fragment_key(start: i32, reverse: bool, reference_index: i32, library: u8) -> i64 {
    (i64::from(start) << 32)
        | i64::from((reference_index << 16) & 0xFFFF_0000_u32 as i32)
        | i64::from((i32::from(library as i8) << 8) & 0xFF00)
        | i64::from(reverse)
}

/// The reads of one BAM by the MarkDuplicates key, then by start: `splitByStart` under `cogroup`.
///
/// A read whose group the first header cannot answer for fails `mapToPair`, inside `stage`.
fn duplicate_groups(
    records: &[BamRecord],
    header: &SamHeader,
    libraries: &[String],
    stage: u32,
) -> Result<std::collections::BTreeMap<(i64, i32), Vec<BamRecord>>, Thrown> {
    let mut groups: std::collections::BTreeMap<(i64, i32), Vec<BamRecord>> = Default::default();
    for record in records {
        let name = &record.read_name;
        let group = read_group_of(record);
        let known = group.and_then(|id| header.read_groups.iter().find(|rg| rg.id == id));
        let Some(known) = known else {
            eprintln!("ERROR Executor - Exception in task 0.0 in stage {stage}.0 (TID {stage})");
            return Err(Thrown::non_user(
                "org.apache.spark.SparkException",
                match group {
                    None => format!(
                        "SAM/BAM/CRAM file (unknown) is malformed: Read {name} is missing the read \
                         group (RG) tag, which is required by the GATK.  Please use \
                         http://gatkforums.broadinstitute.org/discussion/59/\
                         companion-utilities-replacereadgroups to fix this problem"
                    ),
                    Some(id) => format!(
                        "SAM/BAM/CRAM file (unknown) is malformed: Read {name} contains an (RG) \
                         tag with the group {id} which is not found in the file header."
                    ),
                },
            ));
        };
        let library = known.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY);
        let index = libraries.iter().position(|l| l == library).unwrap_or(0);
        let start = stranded_unclipped_start(record);
        let reference = if gatk_engine::read::is_unmapped(record) {
            gatk_engine::read::NO_ALIGNMENT_REFERENCE_INDEX
        } else {
            record.reference_index
        };
        let key = fragment_key(
            start,
            gatk_engine::read::is_reverse_strand(record),
            reference,
            index as u8,
        );
        groups.entry((key, start)).or_default().push(record.clone());
    }
    Ok(groups)
}

/// `CompareDuplicatesSpark.getDupes`.
fn match_type(first: &[BamRecord], second: &[BamRecord]) -> MatchType {
    if first.len() != second.len() {
        return MatchType::SizeUnequal;
    }
    let sorted = |reads: &[BamRecord]| {
        let mut reads = reads.to_vec();
        reads.sort_by(gatk_engine::read_utils::compare_read_coordinate);
        reads
    };
    let (first, second) = (sorted(first), sorted(second));
    if first
        .iter()
        .zip(&second)
        .any(|(one, other)| one.read_name != other.read_name)
    {
        return MatchType::ReadMismatch;
    }
    // Sets, so a read marked twice over counts once.
    let marked = |reads: &[BamRecord]| {
        let mut dupes: Vec<BamRecord> = Vec::new();
        for record in reads
            .iter()
            .filter(|record| gatk_engine::read::is_duplicate(record))
        {
            if !dupes.contains(record) {
                dupes.push(record.clone());
            }
        }
        dupes
    };
    let (first, second) = (marked(&first), marked(&second));
    if first.len() != second.len() {
        return MatchType::DiffNumDupes;
    }
    if first.iter().any(|record| !second.contains(record)) {
        return MatchType::DifferentRepresentativeRead;
    }
    MatchType::Equal
}

/// `CompareDuplicatesSpark.runTool`: the same reads under the same marks, or which way they differ.
///
/// The filter chain is `getReads()`'s, so it holds only for the first BAM; the second is read
/// through a `ReadsSparkSource` of its own, under the same intervals and no filter at all.
pub fn compare_duplicates_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "CompareDuplicatesSpark")?;
    let (output, output2) = (argument(parser, "output"), argument(parser, "output2"));
    if output.is_some() != output2.is_some() {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Arguments '--output' and '--output2' must both be specified together in order to \
             write mismatch bams or not at all",
        ));
    }
    let input2 = argument(parser, "input2").ok_or_else(|| {
        Thrown::command_line("Argument input2 was missing: Argument 'input2' is required")
    })?;
    let second_source = spark_reads_header(&input2)?;
    let second_header = second_source.header().clone();
    let second_source = spark_traversal_source(parser, &input2, second_source)?;

    let filter = read_filter(parser, &filters, &header)?;
    let keep_all: Filter = Box::new(|_| true);
    let traverse = |source, own: &SamHeader, filter| {
        gatk_tools::read_walker::traverse(source, &intervals, filter)
            .map(|records| {
                records
                    .into_iter()
                    .map(|record| compared_read(record, own, &header))
                    .collect::<Vec<_>>()
            })
            .map_err(|error| Thrown::user(format!("{error:?}")))
    };
    let first = traverse(&source, &header, &filter)?;
    let second = traverse(&second_source, &second_header, &keep_all)?;
    if first.len() != second.len() {
        return Err(Thrown::user(format!(
            "input bams have different numbers of mapped reads: {},{}",
            first.len(),
            second.len()
        )));
    }
    gatk_tools::tool_stdout::println(&format!(
        "processing bams with {} mapped reads",
        first.len()
    ));
    let duplicates = |reads: &[BamRecord]| {
        reads
            .iter()
            .filter(|r| gatk_engine::read::is_duplicate(r))
            .count()
    };
    let (first_dupes, second_dupes) = (duplicates(&first), duplicates(&second));
    if first_dupes != second_dupes {
        gatk_tools::tool_stdout::println(&format!(
            "BAMs have different number of total duplicates: {first_dupes},{second_dupes}"
        ));
    }
    gatk_tools::tool_stdout::println(&format!("first and second: {first_dupes},{second_dupes}"));

    // `constructLibraryIndex`.
    let mut libraries: Vec<String> = Vec::new();
    for group in &header.read_groups {
        let library = group.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY);
        if !libraries.iter().any(|known| known == library) {
            libraries.push(library.to_string());
        }
    }
    if libraries.len() > 255 {
        return Err(Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            format!(
                "Detected too many read libraries among read groups header, currently \
                 MarkDuplicatesSpark only supports up to 256 unique readgroup libraries but {} \
                 were found",
                libraries.len()
            ),
        ));
    }
    let first_groups = duplicate_groups(&first, &header, &libraries, 4)?;
    let second_groups = duplicate_groups(&second, &header, &libraries, 5)?;

    // A key the second BAM lacks hands `getDupes` a null, the first thing the cogroup's own stage
    // does with it. With no `--num-reducers` there is one partition, and otherwise `HashPartitioner`
    // says which task fails first: `Objects.hash` of the long key.
    let reducers = number_or(parser, "num-reducers", 0).max(1);
    let lost = first_groups
        .keys()
        .filter(|key| !second_groups.contains_key(*key))
        .map(|(key, _)| {
            let hash = 31i32.wrapping_add((*key ^ ((*key as u64) >> 32) as i64) as i32);
            i64::from(hash).rem_euclid(reducers as i64)
        })
        .min();
    if let Some(task) = lost {
        eprintln!(
            "ERROR Executor - Exception in task {task}.0 in stage 6.0 (TID {})",
            6 + task
        );
        return Err(Thrown {
            failure: Failure::Other,
            exception: "org.apache.spark.SparkException",
            message: None,
        });
    }
    let tagged: Vec<(MatchType, &[BamRecord], &[BamRecord])> = first_groups
        .iter()
        .map(|(key, one)| {
            let other = second_groups[key].as_slice();
            (match_type(one, other), one.as_slice(), other)
        })
        .collect();

    if let (Some(output), Some(output2)) = (&output, &output2) {
        let mut names = std::collections::HashSet::new();
        for (_, one, other) in tagged.iter().filter(|(kind, ..)| *kind != MatchType::Equal) {
            names.extend(one.iter().chain(*other).map(|r| r.read_name.as_str()));
        }
        // `setAttribute("in", ...)` lands on the `@HD` the reads were read under, behind `SO`.
        let (mut written, sort_order) = disq_header(&header);
        written.attributes.set("in", "original read file source");
        for (reads, path) in [(&first, output), (&second, output2)] {
            let kept: Vec<BamRecord> = reads
                .iter()
                .filter(|record| names.contains(record.read_name.as_str()))
                .cloned()
                .collect();
            let partitions = sorted_partitions(parser, sort_order.as_deref(), &kept)?;
            write_spark_partitions(parser, &written, sort_order.as_deref(), &partitions, path)?;
        }
    }

    let count = |kind: MatchType| tagged.iter().filter(|(found, ..)| *found == kind).count();
    if count(MatchType::SizeUnequal) != 0 {
        return Err(Thrown::user(
            "The number of reads by the MarkDuplicates key were unequal, indicating that the BAMs \
             are not the same",
        ));
    }
    if count(MatchType::ReadMismatch) != 0 {
        return Err(Thrown::user(
            "The reads grouped by the MarkDuplicates key were not the same, indicating that the \
             BAMs are not the same",
        ));
    }
    if flag(parser, "print-summary") {
        // The reference walks a `HashSet` of the constants, whose order is their identity hashes
        // and so another on every JVM. Theirs in declaration order is one of the orders it prints.
        gatk_tools::tool_stdout::println("##############################");
        for (kind, name) in MATCH_TYPES {
            gatk_tools::tool_stdout::println(&format!("{name}: {}", count(kind)));
        }
    }
    if flag(parser, "throw-on-diff") {
        for (kind, name) in MATCH_TYPES {
            if kind != MatchType::Equal && count(kind) != 0 {
                return Err(Thrown::user(format!(
                    "found difference between the two BAMs: {name} with count {}",
                    count(kind)
                )));
            }
        }
    }
    Ok(None)
}

/// The reads and headers a Spark metrics tool aggregates, and the `-O` it saves to.
fn spark_metric_reads(
    parser: &Parser,
    tool: &str,
) -> Result<(Vec<BamRecord>, Vec<String>), Thrown> {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    // `getMetricsFile`: the command line, then the time the run started.
    let headers = vec![
        crate::command_line::expanded(tool, parser),
        format!("Started on: {}", java_display_now()),
    ];
    Ok((records, headers))
}

/// `MetricsUtils.saveMetrics`. `-O` is optional and nothing checks it: without one,
/// `BucketUtils.createFile(null)` fails after the job has run. A `--chart` needs R, which the
/// pinned container does not have, and the port does not draw.
fn save_spark_metrics(parser: &Parser, text: &str) -> Result<(), Thrown> {
    let Some(output) = argument(parser, "output") else {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Null object is not allowed here.",
        ));
    };
    std::fs::write(&output, text).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })
}

fn refuse_spark_chart(parser: &Parser) -> Result<(), Thrown> {
    if argument(parser, "chart").is_some() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--chart runs R, which the port does not. This message is the port's own and not GATK's.",
        ));
    }
    Ok(())
}

/// `MeanQualityByCycleSpark.runTool`.
pub fn mean_quality_by_cycle_spark(parser: &Parser) -> Outcome {
    let (records, headers) = spark_metric_reads(parser, "MeanQualityByCycleSpark")?;
    refuse_spark_chart(parser)?;
    let (pf, aligned) = (
        flag(parser, "pfReadsOnly"),
        flag(parser, "alignedReadsOnly"),
    );
    let kept = records
        .iter()
        .filter(|record| gatk_tools::spark_metrics::metrics_read_filter(record, pf, aligned));
    let text = gatk_tools::spark_metrics::mean_quality_by_cycle(kept, &headers);
    save_spark_metrics(parser, &text)?;
    Ok(None)
}

/// `QualityScoreDistributionSpark.runTool`.
pub fn quality_score_distribution_spark(parser: &Parser) -> Outcome {
    let (records, headers) = spark_metric_reads(parser, "QualityScoreDistributionSpark")?;
    refuse_spark_chart(parser)?;
    let (pf, aligned) = (
        flag(parser, "pfReadsOnly"),
        flag(parser, "alignedReadsOnly"),
    );
    let kept = records
        .iter()
        .filter(|record| gatk_tools::spark_metrics::metrics_read_filter(record, pf, aligned));
    let text = gatk_tools::spark_metrics::quality_score_distribution(
        kept,
        flag(parser, "includeNoCalls"),
        &headers,
    );
    save_spark_metrics(parser, &text)?;
    Ok(None)
}

/// `CollectBaseDistributionByCycleSpark.runTool`.
pub fn collect_base_distribution_by_cycle_spark(parser: &Parser) -> Outcome {
    let (records, headers) = spark_metric_reads(parser, "CollectBaseDistributionByCycleSpark")?;
    refuse_spark_chart(parser)?;
    let (pf, aligned) = (
        flag(parser, "pf-reads-only"),
        flag(parser, "aligned-reads-only"),
    );
    let kept = records
        .iter()
        .filter(|record| gatk_tools::spark_metrics::metrics_read_filter(record, pf, aligned));
    let text = gatk_tools::spark_metrics::base_distribution_by_cycle(kept, &headers);
    save_spark_metrics(parser, &text)?;
    Ok(None)
}

/// `MetricsCollectorSparkTool.runTool` around `QualityYieldMetricsCollectorSpark`, which requires
/// a coordinate-sorted input because the default `getExpectedSortOrder` says so and nothing it
/// computes depends on it.
pub fn collect_quality_yield_metrics_spark(parser: &Parser) -> Outcome {
    let (records, headers) = spark_metric_reads(parser, "CollectQualityYieldMetricsSpark")?;
    let text = gatk_tools::spark_metrics::quality_yield(
        records.iter(),
        flag(parser, "use-original-qualities"),
        &headers,
    );
    save_spark_metrics(parser, &text)?;
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
    let (header, sort_order) = disq_header(header);
    let partitions = if sort_to_header {
        sorted_partitions(parser, sort_order.as_deref(), records)?
    } else {
        vec![records.to_vec()]
    };
    write_spark_partitions(parser, &header, sort_order.as_deref(), &partitions, output)
}

/// The header as Disq writes it, `GO` and then `SO` moved behind what else `@HD` holds, and its order.
fn disq_header(header: &SamHeader) -> (SamHeader, Option<String>) {
    let mut header = header.clone();
    let sort_order = header.attributes.get("SO").map(str::to_string);
    let group_order = header.attributes.get("GO").unwrap_or("none").to_string();
    header.attributes.remove("SO");
    header.attributes.remove("GO");
    header.set_group_order(&group_order);
    if let Some(order) = &sort_order {
        header.set_sort_order(order);
    }
    (header, sort_order)
}

/// [`write_spark_reads`] past the sort: the header as Disq writes it, and one stream per partition.
fn write_spark_partitions(
    parser: &Parser,
    header: &SamHeader,
    sort_order: Option<&str>,
    partitions: &[Vec<BamRecord>],
    output: &str,
) -> Result<(), Thrown> {
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
                header,
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
        htsjdk_bam::writer::BamWriter::with_compression(Vec::new(), header, level, deflater)
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
    for partition in partitions {
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

    if flag(parser, "create-output-bam-index") && sort_order == Some("coordinate") {
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

/// `SparkUtils.sortUsingElementsAsKeys` with `ReadsSparkSink.getSAMRecordComparator`: `sortByKey`
/// over a `RangePartitioner` of `getRecommendedNumReducers()` partitions, with
/// `HeaderlessSAMRecordCoordinateComparator` for a coordinate-sorted header and htsjdk's
/// `SAMRecordQueryNameComparator` for a queryname-sorted one. Any other order is written as it was
/// read.
///
/// The sort is stable, as Spark's `TimSort` is, so equal keys keep the traversal's order.
fn sorted_partitions(
    parser: &Parser,
    sort_order: Option<&str>,
    records: &[BamRecord],
) -> Result<Vec<Vec<BamRecord>>, Thrown> {
    let compare: fn(&BamRecord, &BamRecord) -> std::cmp::Ordering = match sort_order {
        Some("coordinate") => htsjdk_bam::coordinate::compare,
        Some("queryname") => htsjdk_bam::query_name::compare,
        _ => return Ok(vec![records.to_vec()]),
    };
    let mut sorted = records.to_vec();
    sorted.sort_by(compare);

    // `getRecommendedNumReducers`: the argument, else one per ten megabytes of input.
    let mut reducers = number_or(parser, "num-reducers", 0).max(0) as usize;
    if reducers == 0 {
        let size = argument(parser, "input")
            .and_then(|input| std::fs::metadata(input).ok())
            .map(|data| data.len())
            .unwrap_or(0);
        reducers = 1 + (size / (10 * 1024 * 1024)) as usize;
    }
    let bounds = range_bounds(&sorted, reducers, compare)?;
    let mut partitions: Vec<Vec<BamRecord>> = vec![Vec::new(); bounds.len() + 1];
    for record in sorted {
        // `getPartition`: the first bound the key is not greater than.
        let partition = bounds
            .iter()
            .position(|bound| compare(&record, bound).is_le())
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
fn range_bounds(
    sorted: &[BamRecord],
    partitions: usize,
    compare: fn(&BamRecord, &BamRecord) -> std::cmp::Ordering,
) -> Result<Vec<BamRecord>, Thrown> {
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
                .is_none_or(|previous| compare(key, previous).is_gt())
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
