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
pub(super) struct SparkStart {
    pub(super) source: ReadsDataSource,
    pub(super) header: SamHeader,
    pub(super) intervals: Vec<gatk_engine::interval::SimpleInterval>,
    pub(super) filters: Vec<gatk_tools::filter_resolution::ResolvedFilter>,
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
pub(super) fn spark_startup(parser: &Parser, tool: &str) -> Result<SparkStart, Thrown> {
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
    // Disq reads the header at `--read-validation-stringency`, and STRICT refuses this line.
    if let Some(group) = header.read_groups.iter().find(|group| {
        group.attributes.get("SM").is_none()
            && argument(parser, "read-validation-stringency").as_deref() == Some("STRICT")
    }) {
        return Err(Thrown::non_user(
            "htsjdk.samtools.SAMFormatException",
            format!(
                "Error parsing SAM header. @RG line missing SM tag. Line:\n@RG\tID:{}",
                group.id
            ),
        ));
    }

    // `initializeIntervals`, against `getBestAvailableSequenceDictionary`: the reference's, then
    // the reads'.
    let best = reference.clone().unwrap_or_else(|| header.clone());
    let intervals = match interval_arguments(parser, &best)? {
        Some(parameters) => parameters.intervals,
        None => Vec::new(),
    };

    // The aligners' `validateSequenceDictionaries`: unaligned reads need no dictionary.
    let aligning = matches!(tool, "BwaSpark" | "BwaAndMarkDuplicatesPipelineSpark")
        || (tool == "ReadsPipelineSpark" && flag(parser, "align"));
    if !aligning && !flag(parser, "disable-sequence-dictionary-validation") {
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
    // BAM. Without one the reference fails inside the first executor task.
    let source = if intervals_specified(parser) {
        let path = std::path::Path::new(input);
        match htsjdk_bam::sam_files::find_index(path) {
            Some(index) => ReadsDataSource::open(path, &index)
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
            None => {
                let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
                eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
                return Err(Thrown::non_user(
                    "org.apache.spark.SparkException",
                    format!(
                        "Intervals set but no index file found for file:{}",
                        absolute.display()
                    ),
                ));
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
    let (variant, intervals) = variant_spark_startup(parser, false)?;
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
///
/// A `VariantWalkerSpark` asks its driving variants first: `getBestAvailableSequenceDictionary` is
/// overridden there, so `PrintVariantsSpark`'s intervals resolve against the VCF's contigs whatever
/// `--reference` says, unless the header has none or only lengths an index synthesized.
fn variant_spark_startup(
    parser: &Parser,
    variant_walker: bool,
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
    let driving = match argument(parser, "variant") {
        Some(variant) if variant_walker && intervals_specified(parser) => {
            let (_, text) = open_feature_input(&variant)?;
            let dictionary = vcf_dictionary(&text);
            dictionary
                .sequences
                .iter()
                .any(|sequence| sequence.length > 0)
                .then_some(dictionary)
        }
        _ => None,
    };
    let intervals =
        match driving.as_ref().or(reference.as_ref()).or(reads.as_ref()) {
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
    let (variant, intervals) = variant_spark_startup(parser, true)?;
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

    // Nothing kept is not one answer but three, and which depends on the index. Disq keeps a split
    // only where a query touches one of the index's blocks: a tabix bin of 16384 bases, a Tribble
    // bin of 8000 (shortcut: the width `IndexFeatureFile` left this corpus's small VCF, which it
    // may optimise on a larger one). With no split at all the shard join without `--shuffle` asks
    // for no partitions, and a block-compressed output with an index has no part to merge.
    let touches_block = match &intervals {
        None => true,
        Some(intervals) => {
            let width = if variant.ends_with(".gz") || variant.ends_with(".bgz") {
                16384
            } else {
                8000
            };
            loci.iter().any(|locus| {
                intervals.iter().any(|interval| {
                    interval.contig == locus.contig
                        && (locus.start - 1) / width <= (interval.end - 1) / width
                        && (interval.start - 1) / width <= (locus.stop - 1) / width
                })
            })
        }
    };
    if kept.is_empty() && !touches_block && !flag(parser, "shuffle") {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "requirement failed: Number of partitions (0) must be positive.",
        ));
    }
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
    if kept.is_empty() && touches_block && !flag(parser, "shuffle") {
        eprintln!("ERROR Executor - Exception in task 0.0 in stage 1.0 (TID 1)");
        return Err(Thrown::non_user(
            "java.lang.ArrayIndexOutOfBoundsException",
            "Index -1 out of bounds for length 1",
        ));
    }
    if kept.is_empty()
        && !touches_block
        && (output.ends_with(".gz") || output.ends_with(".bgz"))
        && flag(parser, "create-output-variant-index")
    {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Cannot merge zero tabix files",
        ));
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

    write_spark_vcf(parser, &output, &header_text, &records_text)
}

/// `VariantsSparkSink.writeVariants`: the header and the records as two parts, block compressed by
/// name, with the index Disq finalizes.
pub(super) fn write_spark_vcf(
    parser: &Parser,
    output: &str,
    header_text: &str,
    records_text: &str,
) -> Outcome {
    let failed = |error: String| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    };
    let block_compressed = output.ends_with(".gz") || output.ends_with(".bgz");
    if !block_compressed {
        std::fs::write(output, format!("{header_text}{records_text}"))
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
    let mut bytes = part(header_text)?;
    bytes.extend_from_slice(&part(records_text)?);
    bytes.extend_from_slice(&htsjdk_bgzf::EMPTY_GZIP_BLOCK);
    std::fs::write(output, &bytes).map_err(|error| failed(error.to_string()))?;

    if flag(parser, "create-output-variant-index") {
        // Disq finalizes the index at the position the records end, BEFORE the terminator it
        // appends, where `IndexFeatureFile` reads the whole file and ends the last chunk after it.
        let source = index_feature_file::Source::new(output);
        let records_end = bytes.len() - htsjdk_bgzf::EMPTY_GZIP_BLOCK.len();
        let index = index_feature_file::build_tabix(
            &bytes[..records_end],
            &source,
            output,
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
    no_reads_refusal(parser, &records, &whole, shard_size)?;
    oversized_read_refusal(parser, &records, &header, shard_size)?;
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

/// `SparkSharder`'s overlaps partitioning, which `--shuffle` does not go through: no read leaves
/// `computePartitionReadExtents` nothing, so the first shard's partition index is -1 in whichever
/// task `parallelize` gave that shard. This happens while the RDD is built, before the tool runs.
pub(super) fn no_reads_refusal(
    parser: &Parser,
    records: &[BamRecord],
    intervals: &[gatk_engine::interval::SimpleInterval],
    shard_size: i32,
) -> Result<(), Thrown> {
    if records.is_empty() && !flag(parser, "shuffle") {
        let master = scalar(parser, "spark-master").unwrap_or_else(|| "local[*]".to_string());
        let slices = master
            .strip_prefix("local[")
            .and_then(|rest| rest.split([',', ']']).next())
            .map_or(1, |threads| {
                threads
                    .parse::<usize>()
                    .unwrap_or_else(|_| std::thread::available_parallelism().map_or(1, usize::from))
            });
        let shards: usize = intervals
            .iter()
            .map(|interval| ((interval.end - interval.start) / shard_size + 1) as usize)
            .sum();
        // `ParallelCollectionRDD.slice`: the slice `i` holds `i * n / p .. (i + 1) * n / p`.
        let task = (0..slices)
            .find(|i| i * shards / slices < (i + 1) * shards / slices)
            .unwrap_or(0);
        eprintln!(
            "ERROR Executor - Exception in task {task}.0 in stage 1.0 (TID {})",
            task + 1
        );
        return Err(Thrown::non_user(
            "java.lang.ArrayIndexOutOfBoundsException",
            "Index -1 out of bounds for length 1",
        ));
    }
    Ok(())
}

/// The same partitioning filling its shards, in coordinate order and inside a task, so after
/// anything the tool does before its first action: a read longer than a shard is refused.
pub(super) fn oversized_read_refusal(
    parser: &Parser,
    records: &[BamRecord],
    header: &SamHeader,
    shard_size: i32,
) -> Result<(), Thrown> {
    if flag(parser, "shuffle") {
        return Ok(());
    }
    for record in records.iter().filter(|record| !record.read_unmapped()) {
        let size = record.alignment_end() - record.alignment_start + 1;
        if size > shard_size {
            let contig = &header.sequences[record.reference_index as usize].name;
            return Err(Thrown::non_user(
                "org.apache.spark.SparkException",
                format!(
                    "Max size of locatable exceeded. Max size is {shard_size}, but locatable size \
                     is {size}. Try increasing shard size and/or padding. Locatable: {} \
                     {contig}:{}-{}",
                    record.read_name,
                    record.alignment_start,
                    record.alignment_end()
                ),
            ));
        }
    }
    Ok(())
}

/// `CalcMetadataSpark.runTool`: every read, with no read filter, into one partition's statistics.
pub fn calc_metadata_spark(parser: &Parser) -> Outcome {
    use gatk_tools::read_metadata::{MetadataFailure, ReadMetadata};
    let SparkStart {
        source,
        header,
        intervals,
        ..
    } = spark_startup(parser, "CalcMetadataSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // Disq cuts a BAM at the split size, where the port has one partition to report.
    let input = arguments(parser, "input").remove(0);
    let split = match number_or(parser, "bam-partition-size", 0) {
        0 => 32 * 1024 * 1024,
        size => i64::from(size),
    };
    if std::fs::metadata(&input).is_ok_and(|file| file.len() as i64 > split) {
        return Err(Thrown::user(format!(
            "{input} is larger than one partition, and this port reports a single partition's \
             statistics."
        )));
    }
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let metadata = ReadMetadata::new(&header, &records).map_err(|failure| match failure {
        MetadataFailure::Unpaired => {
            eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
            Thrown::non_user(
                "org.apache.spark.SparkException",
                "Job aborted due to stage failure: Task 0 in stage 0.0 failed 1 times: \
                 java.lang.IllegalStateException: Cannot get mate information for an unpaired read",
            )
        }
        MetadataFailure::NoReads => Thrown::non_user("java.lang.ArithmeticException", "/ by zero"),
        MetadataFailure::Order(error) => Thrown::user(format!("{error:?}")),
    })?;
    let bytes = if flag(parser, "write-as-binary") {
        metadata.flushed_binary()
    } else {
        metadata.text().into_bytes()
    };
    std::fs::write(&output, bytes).map_err(|error| {
        Thrown::user(format!(
            "Couldn't write file {output} because Could not create file at path:{output} due to \
             {output} ({})",
            java_io_reason(&error)
        ))
    })?;
    Ok(None)
}

/// `ExtractSVEvidenceSpark`: the metadata and the high-coverage regions first, whatever becomes of
/// the classification, then the evidence as a text file directory.
pub fn extract_sv_evidence_spark(parser: &Parser) -> Outcome {
    use gatk_tools::read_metadata::{MetadataArguments, MetadataFailure, ReadMetadata};
    use gatk_tools::sv_evidence;
    let SparkStart {
        source,
        header,
        intervals,
        ..
    } = spark_startup(parser, "ExtractSVEvidenceSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let input = arguments(parser, "input").remove(0);
    let split = match number_or(parser, "bam-partition-size", 0) {
        0 => 32 * 1024 * 1024,
        size => i64::from(size),
    };
    if std::fs::metadata(&input).is_ok_and(|file| file.len() as i64 > split) {
        return Err(Thrown::user(format!(
            "{input} is larger than one partition, and this port reports a single partition's \
             evidence."
        )));
    }
    let arguments_read = sv_evidence::EvidenceArguments {
        metadata: MetadataArguments {
            min_evidence_mapq: number_or(parser, "min-evidence-mapq", 20) as u8,
            min_evidence_match_length: i64::from(number_or(
                parser,
                "min-evidence-match-length",
                45,
            )),
            allowed_short_fragment_overhang: number_or(
                parser,
                "allowed-short-fragment-overhang",
                10,
            ),
            max_tracked_fragment_length: number_or(parser, "max-tracked-fragment-length", 2000)
                as usize,
        },
        high_depth_coverage_peak_factor: number_or(parser, "high-depth-coverage-peak-factor", 7),
        high_depth_coverage_factor: number_or(parser, "high-depth-coverage-factor", 3),
    };
    if header.attributes.get("SO") != Some("coordinate") {
        return Err(Thrown::non_user(
            "java.lang.IllegalStateException",
            "The reads must be coordinate sorted.",
        ));
    }
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let mut ignored_contigs = Vec::new();
    if let Some(path) = argument(parser, "cross-contigs-to-ignore") {
        let text = std::fs::read_to_string(&path).map_err(|error| {
            Thrown::user(format!(
                "Couldn't read file {path}. Error was: {path} ({})",
                java_io_reason(&error)
            ))
        })?;
        for line in text.lines() {
            match header
                .sequences
                .iter()
                .position(|sequence| sequence.name == line)
            {
                Some(contig) => ignored_contigs.push(contig as i32),
                None => {
                    return Err(Thrown::user(format!(
                        "crossContigToIgnoreFile contains an unrecognized contig name: {line}"
                    )))
                }
            }
        }
    }
    let metadata =
        ReadMetadata::with_arguments(&header, &records, &arguments_read.metadata, ignored_contigs)
            .map_err(|failure| {
                match failure {
        MetadataFailure::Unpaired => evidence_task_failure(
            "java.lang.IllegalStateException: Cannot get mate information for an unpaired read",
        ),
        MetadataFailure::NoReads => Thrown::non_user("java.lang.ArithmeticException", "/ by zero"),
        MetadataFailure::Order(error) => Thrown::user(format!("{error:?}")),
    }
            })?;
    let write_text = |path: &str, text: String| {
        std::fs::write(path, text).map_err(|error| {
            Thrown::user(format!(
                "Couldn't write file {path} because Could not create file at path:{path} due to \
                 {path} ({})",
                java_io_reason(&error)
            ))
        })
    };
    if let Some(path) = argument(parser, "read-metadata") {
        write_text(&path, metadata.text())?;
    }
    let lengths: Vec<i32> = header
        .sequences
        .iter()
        .map(|sequence| sequence.length)
        .collect();
    let ignored =
        sv_evidence::high_coverage_intervals(&records, &metadata, &lengths, &arguments_read);
    if let Some(path) = argument(parser, "high-coverage-intervals") {
        write_text(&path, sv_evidence::high_coverage_bed(&ignored, &metadata))?;
    }
    let lines = sv_evidence::classify(&records, &metadata, &ignored, &arguments_read.metadata)
        .map_err(|message| evidence_task_failure(&message))?;
    let mut text = String::new();
    for line in lines {
        text.push_str(&line);
        text.push('\n');
    }
    let directory = std::path::Path::new(&output);
    std::fs::create_dir_all(directory)
        .and_then(|()| std::fs::write(directory.join("part-00000"), &text))
        .and_then(|()| std::fs::write(directory.join("_SUCCESS"), b""))
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

fn evidence_task_failure(message: &str) -> Thrown {
    eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
    Thrown::non_user(
        "org.apache.spark.SparkException",
        format!("Job aborted due to stage failure: Task 0 in stage 0.0 failed 1 times: {message}"),
    )
}

/// A task's failure in an assembly tool, whose contig alignments are a cached RDD: the block
/// manager logs what stopped that one first, under the number the tool's RDDs have reached.
fn assembly_task_failure(cached_rdd: u32, message: &str) -> Thrown {
    if message.starts_with("no primary alignment for read") {
        eprintln!(
            "WARN  BlockManager - Putting block rdd_{cached_rdd}_0 failed due to exception org.broadinstitute.hellbender.exceptions.GATKException: {message}."
        );
    }
    evidence_task_failure(message)
}

/// `DiscoverVariantsFromContigAlignmentsSAMSpark.runTool`: the arguments, the external calls, the
/// output's path and the sample, then [`gatk_tools::sv_discovery`] and the VCF.
///
/// A BAM split over several partitions keeps every answer; what would change is only the order of
/// records that sort equal, which this port takes from one partition.
pub fn discover_variants_from_contig_alignments_sam_spark(parser: &Parser) -> Outcome {
    use gatk_tools::sv_discovery;
    let tool = "DiscoverVariantsFromContigAlignmentsSAMSpark";
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    // `getOutputPathInputs`' argument is required, and read after the arguments are validated.
    let prefix = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    validate_discover_stage(parser)?;

    let sample = sv_sample_id(&header)?;
    let cnv_calls = match argument(parser, "cnv-calls") {
        Some(path) => Some(load_cnv_calls(&path, &header, &sample)?),
        None => None,
    };
    // `getVcfOutputPath`.
    let output = match std::fs::metadata(&prefix) {
        Ok(found) if found.is_dir() => format!(
            "{prefix}{}{sample}_inv_del_ins.vcf",
            if prefix.ends_with('/') { "" } else { "/" }
        ),
        Ok(_) => {
            return Err(Thrown::user(format!(
                "Provided prefix for output is pointing to an existing file: {prefix}"
            )))
        }
        Err(_) => format!(
            "{prefix}{}{sample}_inv_del_ins.vcf",
            if prefix.ends_with('/') { "" } else { "_" }
        ),
    };

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let fasta = htsjdk_bam::fasta_index::IndexedFasta::open(std::path::Path::new(&reference_path))
        .map_err(|error| Thrown::user(error.message()))?;
    let mut reference = SparkReference { fasta };
    let adjacencies = sv_discovery::adjacencies(&records, &header)
        .map_err(|message| evidence_task_failure(&message))?;

    check_sv_truth(parser, &header, |padding| {
        adjacencies.check_padded_breakpoints(padding)
    })?;

    let variants = adjacencies
        .variants(&header, &mut reference, cnv_calls.as_deref())
        .map_err(|message| evidence_task_failure(&message))?;

    write_sv_vcf(parser, tool, &header, &variants, &output)?;
    Ok(None)
}

/// `SvDiscoveryUtils.evaluateIntervalsAndNarls`, which only logs, once it has read the truth and
/// made its intervals.
fn check_sv_truth(
    parser: &Parser,
    header: &SamHeader,
    check_breakpoints: impl FnOnce(i32) -> Result<(), String>,
) -> Result<(), Thrown> {
    let Some(truth) = argument(parser, "truth-vcf") else {
        return Ok(());
    };
    let padding = number_or(parser, "truth-interval-padding", 50);
    let refused = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);
    let file = sv_feature_file(&truth)?;
    for record in &file.records {
        let svtype = record
            .attributes
            .iter()
            .find(|(name, _)| name == "SVTYPE")
            .and_then(|(_, value)| value.format());
        // `getStructuralVariantType`: nothing without the attribute, and `valueOf` with it.
        let Some(svtype) = svtype else { continue };
        if !header
            .sequences
            .iter()
            .any(|sequence| sequence.name == record.contig)
        {
            return Err(Thrown::user(format!(
                "VCF contig {} does not appear in dictionary.",
                record.contig
            )));
        }
        if !["BND", "CNV", "DEL", "DUP", "INS", "INV"].contains(&svtype.as_str()) {
            return Err(refused(format!(
                "No enum constant htsjdk.variant.variantcontext.StructuralVariantType.{svtype}"
            )));
        }
        gatk_tools::sv_discovery::check_padded_start(record.start as i32 - padding)
            .map_err(refused)?;
    }
    check_breakpoints(padding).map_err(refused)
}

/// `SVVCFWriter.writeVariants`.
fn write_sv_vcf(
    parser: &Parser,
    tool: &str,
    header: &SamHeader,
    variants: &[htsjdk_vcf::variant::VariantContext],
    output: &str,
) -> Result<(), Thrown> {
    let mut vcf = htsjdk_vcf::header::VcfHeader::new();
    vcf.lines = gatk_tools::sv_discovery::header_lines();
    for (index, sequence) in header.sequences.iter().enumerate() {
        vcf.lines.push(htsjdk_vcf::header::HeaderLine::contig(
            &sequence.name,
            sequence.length as i64,
            index as i32,
        ));
    }
    vcf.lines
        .extend(default_tool_vcf_header_lines(parser, tool));
    let failed = |error: String| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    };
    let mut text = vcf.write();
    let encoder = htsjdk_vcf::encoder::VcfEncoder::new(&vcf);
    for variant in variants {
        encoder
            .encode_into(variant, &mut text)
            .map_err(|error| failed(format!("{error:?}")))?;
        text.push('\n');
    }
    std::fs::write(output, text).map_err(|error| {
        Thrown::user(format!(
            "Could not create file at path:{output} due to {output} ({})",
            java_io_reason(&error)
        ))
    })
}

/// `DiscoverVariantsFromContigAlignmentsArgumentCollection.validate`.
fn validate_discover_stage(parser: &Parser) -> Result<(), Thrown> {
    for (name, field, default) in [
        ("min-align-length", "minAlignLength", 50),
        (
            "assembly-imprecise-evidence-overlap-uncertainty",
            "assemblyImpreciseEvidenceOverlapUncertainty",
            100,
        ),
        (
            "imprecise-variant-evidence-threshold",
            "impreciseVariantEvidenceThreshold",
            7,
        ),
        ("truth-interval-padding", "truthIntervalPadding", 50),
        (
            "max-callable-imprecise-deletion-size",
            "maxCallableImpreciseVariantDeletionSize",
            15000,
        ),
    ] {
        let value = number_or(parser, name, default);
        if value <= 0 {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!("invalid value provided to {field}: {value}"),
            ));
        }
    }
    Ok(())
}

/// `SvDiscoverFromLocalAssemblyContigAlignmentsSpark.runTool`.
///
/// A VCF that was written stays written when a later step of the job fails, as in the tool.
pub fn sv_discover_from_local_assembly_contig_alignments_spark(parser: &Parser) -> Outcome {
    use gatk_tools::sv_discovery::assembly;
    let tool = "SvDiscoverFromLocalAssemblyContigAlignmentsSpark";
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    let prefix = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    validate_discover_stage(parser)?;
    // `broadcastCNVCalls`, whose sample comes from the same header.
    let cnv_calls = match argument(parser, "cnv-calls") {
        Some(path) => {
            let sample = sv_sample_id(&header)?;
            Some(load_cnv_calls(&path, &header, &sample)?)
        }
        None => None,
    };
    // `getOutputPrefix`.
    let prefix = match std::fs::metadata(&prefix) {
        Ok(found) if found.is_dir() => format!(
            "{prefix}{}{}_",
            if prefix.ends_with('/') { "" } else { "/" },
            sv_sample_id(&header)?
        ),
        Ok(_) => {
            return Err(Thrown::user(format!(
                "Provided prefix for output is pointing to an existing file: {prefix}"
            )))
        }
        Err(_) => format!(
            "{prefix}{}{}_",
            if prefix.ends_with('/') { "" } else { "_" },
            sv_sample_id(&header)?
        ),
    };
    let canonical = sv_canonical_chromosomes(parser, &header)?;
    sv_sample_id(&header)?;

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let fasta = htsjdk_bam::fasta_index::IndexedFasta::open(std::path::Path::new(&reference_path))
        .map_err(|error| Thrown::user(error.message()))?;
    let mut reference = SparkReference { fasta };
    let in_task = |message: String| assembly_task_failure(23, &message);

    let classified = assembly::preprocess(&records, &header, &canonical).map_err(in_task)?;

    // `extractSimpleVariants`.
    let simple_calls =
        assembly::simple_calls(&classified, &header, &mut reference).map_err(in_task)?;
    check_sv_truth(parser, &header, |padding| {
        simple_calls.check_padded_breakpoints(padding)
    })?;
    let simple = simple_calls
        .variants(&header, cnv_calls.as_deref())
        .map_err(in_task)?;
    let write = |name: &str, variants: &[htsjdk_vcf::variant::VariantContext]| {
        write_sv_vcf(
            parser,
            tool,
            &header,
            &assembly::sorted(variants, &header),
            &format!("{prefix}{name}"),
        )
    };
    write("NonComplex.vcf", &simple)?;

    // `extractCpxVariants`.
    let complex = assembly::complex_calls(&classified, &header, &mut reference).map_err(in_task)?;
    write("Complex.vcf", &assembly::complex_variants(&complex))?;
    let one_segment =
        assembly::zero_or_one_segment_variants(&complex, &mut reference).map_err(in_task)?;
    let from_segments =
        assembly::multi_segment_variants(&complex, &mut reference).map_err(in_task)?;
    let pairs = assembly::pair_iteration_adjacencies(&complex, &records, &header, &canonical)
        .map_err(in_task)?;
    check_sv_truth(parser, &header, |padding| {
        pairs.check_padded_breakpoints(padding)
    })?;
    let from_pairs = assembly::pair_iteration_variants(
        &pairs,
        &complex,
        &header,
        &mut reference,
        cnv_calls.as_deref(),
    )
    .map_err(in_task)?;
    let multi_segment = assembly::remove_duplicates(from_segments, from_pairs).map_err(in_task)?;
    write("cpx_reinterpreted_simple_1_seg.vcf", &one_segment)?;
    write("cpx_reinterpreted_simple_multi_seg.vcf", &multi_segment)?;

    // `writeSAMfilesForUnknown`.
    if flag(parser, "write-sam") {
        let sets = classified.unknown_reads(&records).map_err(in_task)?;
        // What reaches the writer has no other `@HD` tag, whatever the input's line carried.
        let mut sorted_header = header.clone();
        sorted_header.attributes = SamHeader::default().attributes;
        sorted_header.set_group_order("none");
        sorted_header.set_sort_order("queryname");
        for (name, set) in ["AMBIGUOUS", "INCOMPLETE", "UNINFORMATIVE"]
            .iter()
            .zip(sets)
        {
            let path = format!("{prefix}{name}.bam");
            let failed = |error: String| {
                Thrown::non_user(PORT_FAILURE, format!("could not write {path}: {error}"))
            };
            let (level, deflater) = output_compression(parser);
            let mut writer = htsjdk_bam::writer::BamWriter::with_compression(
                Vec::new(),
                &sorted_header,
                level,
                deflater,
            )
            .map_err(|error| failed(error.to_string()))?;
            for record in set {
                writer
                    .write(record)
                    .map_err(|error| failed(error.to_string()))?;
            }
            let bytes = writer.finish().map_err(|error| failed(error.to_string()))?;
            std::fs::write(&path, bytes).map_err(|error| failed(error.to_string()))?;
        }
    }

    let mut merged = simple;
    merged.extend(one_segment);
    merged.extend(multi_segment);
    let filtered = assembly::filter_merged(
        merged,
        number_or(parser, "min-mq", 30),
        number_or(parser, "min-align-length", 50),
    );
    write("merged_simple.vcf", &filtered)?;
    Ok(None)
}

/// `SVUtils.getCanonicalChromosomes`.
fn sv_canonical_chromosomes(parser: &Parser, header: &SamHeader) -> Result<Vec<String>, Thrown> {
    let mut canonical: Vec<String> = header
        .sequences
        .iter()
        .map(|sequence| sequence.name.clone())
        .collect();
    if let Some(path) = argument(parser, "non-canonical-contig-names-file") {
        let text = std::fs::read_to_string(&path).map_err(|_| {
            Thrown::user(format!(
                "Can't read nonCanonicalContigNamesFile file {path}"
            ))
        })?;
        for line in text.lines() {
            canonical.retain(|name| name != line);
        }
    }
    Ok(canonical)
}

/// `CpxVariantReInterpreterSpark.runTool`: the complex calls of a VCF, the reads of the contigs that
/// made them, and the two files of simple variants the assembly step derives from its own.
///
/// The header's dictionary is the VCF's, not the reads', and nothing here looks at a sample.
pub fn cpx_variant_re_interpreter_spark(parser: &Parser) -> Outcome {
    use gatk_tools::sv_discovery::assembly;
    let tool = "CpxVariantReInterpreterSpark";
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    let canonical = sv_canonical_chromosomes(parser, &header)?;
    let complex_vcf = argument(parser, "cpx-vcf").ok_or_else(|| {
        Thrown::command_line("Argument cpx-vcf was missing: Argument 'cpx-vcf' is required")
    })?;
    let prefix = argument(parser, "prefix-out-vcf").ok_or_else(|| {
        Thrown::command_line(
            "Argument prefix-out-vcf was missing: Argument 'prefix-out-vcf' is required",
        )
    })?;
    // Hadoop-BAM's input format, which reads the path before anything else does.
    if !std::path::Path::new(&complex_vcf).exists() {
        return Err(Thrown::non_user(
            "java.lang.RuntimeException",
            format!("File {complex_vcf} does not exist"),
        ));
    }
    let (_, text) = open_feature_input(&complex_vcf)?;
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let variants: Vec<_> = if intervals.is_empty() {
        file.records
    } else {
        let index = if complex_vcf.ends_with(".gz") || complex_vcf.ends_with(".bgz") {
            format!("{complex_vcf}.tbi")
        } else {
            format!("{complex_vcf}.idx")
        };
        if !std::path::Path::new(&index).exists() {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!("Intervals set but no index file found for {complex_vcf} at {index}"),
            ));
        }
        file.records
            .into_iter()
            .filter(|record| {
                intervals.iter().any(|interval| {
                    interval.contig == record.contig
                        && i64::from(interval.start) <= record.stop
                        && record.start <= i64::from(interval.end)
                })
            })
            .collect()
    };

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let fasta = htsjdk_bam::fasta_index::IndexedFasta::open(std::path::Path::new(&reference_path))
        .map_err(|error| Thrown::user(error.message()))?;
    let mut reference = SparkReference { fasta };
    let in_task = |message: String| assembly_task_failure(32, &message);

    let complex = assembly::complex_calls_from_variants(variants).map_err(in_task)?;
    let one_segment =
        assembly::zero_or_one_segment_variants(&complex, &mut reference).map_err(in_task)?;
    let from_segments =
        assembly::multi_segment_variants(&complex, &mut reference).map_err(in_task)?;
    let pairs = assembly::pair_iteration_adjacencies(&complex, &records, &header, &canonical)
        .map_err(in_task)?;
    // The pairs go through the same interpreter as the discovery's, truth and all.
    check_sv_truth(parser, &header, |padding| {
        pairs.check_padded_breakpoints(padding)
    })?;
    let from_pairs =
        assembly::pair_iteration_variants(&pairs, &complex, &header, &mut reference, None)
            .map_err(in_task)?;
    let multi_segment = assembly::remove_duplicates(from_segments, from_pairs).map_err(in_task)?;

    let dictionary = vcf_dictionary(&text);
    // `getSequenceDictionary` is null with no `##contig`, and the writer refuses it however few
    // variants there are to sort.
    if dictionary.sequences.is_empty() {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Null object is not allowed here.",
        ));
    }
    for (suffix, variants) in [
        ("_1_seg.vcf", &one_segment),
        ("_multi_seg.vcf", &multi_segment),
    ] {
        write_sv_vcf(
            parser,
            tool,
            &dictionary,
            &assembly::sorted(variants, &dictionary),
            &format!("{prefix}{suffix}"),
        )?;
    }
    Ok(None)
}

/// `new FeatureDataSource<>(path, ...)` on a VCF, whose refusal of a missing file names it as a URI.
fn sv_feature_file(path: &str) -> Result<htsjdk_vcf::reader::VcfFile, Thrown> {
    if !std::path::Path::new(path).exists() {
        let absolute = std::path::absolute(path).unwrap_or_else(|_| std::path::PathBuf::from(path));
        return Err(Thrown::user(format!(
            "Couldn't read file file://{}. Error was: It doesn't exist.",
            absolute.display()
        )));
    }
    let (_, text) = open_feature_input(path)?;
    htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })
}

/// `SVUtils.getSampleId`.
fn sv_sample_id(header: &SamHeader) -> Result<String, Thrown> {
    let mut samples: Vec<String> = Vec::new();
    for group in &header.read_groups {
        // A read group with no sample is `null` in the set, and in the path it builds.
        let sample = group
            .attributes
            .get("SM")
            .map_or_else(|| "null".to_string(), str::to_string);
        if !samples.contains(&sample) {
            samples.push(sample);
        }
    }
    if samples.len() == 1 {
        return Ok(samples.remove(0));
    }
    let shown = gatk_engine::java_hash::hash_set_order(&samples).unwrap_or(samples);
    Err(Thrown::non_user(
        "java.lang.IllegalStateException",
        format!(
            "Read groups must contain reads from one and only one sample, but we are finding the \
             following ones in the given header: \t[{}]",
            shown.join(", ")
        ),
    ))
}

/// `CNVInputReader.loadCNVCalls`, the tree's order being its intervals'.
fn load_cnv_calls(
    path: &str,
    header: &SamHeader,
    sample: &str,
) -> Result<Vec<gatk_tools::sv_discovery::CnvCall>, Thrown> {
    let state = |message: String| Thrown::non_user("java.lang.IllegalStateException", message);
    let file = sv_feature_file(path)?;
    if file.header.samples.len() != 1 {
        return Err(state("CNV call VCF should be single sample".to_string()));
    }
    if !file.header.samples.iter().any(|name| name == sample) {
        return Err(state(format!(
            "CNV call VCF does not contain calls for sample {sample}"
        )));
    }
    let contigs: Vec<(String, i64)> = file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Contig { fields, .. } => {
                let field = |key: &str| {
                    fields
                        .iter()
                        .find(|(name, _)| name == key)
                        .map(|(_, value)| value.clone())
                };
                Some((
                    field("ID")?,
                    field("length")
                        .and_then(|text| text.parse().ok())
                        .unwrap_or(0),
                ))
            }
            _ => None,
        })
        .collect();
    if contigs.is_empty() {
        return Err(state(
            "CNV calls file does not have a valid sequence dictionary".to_string(),
        ));
    }
    let same = contigs.len() == header.sequences.len()
        && contigs
            .iter()
            .zip(&header.sequences)
            .all(|(one, two)| one.0 == two.name && one.1 == two.length as i64);
    if !same {
        return Err(state(
            "CNV calls file does not have the same sequence dictionary as the read evidence"
                .to_string(),
        ));
    }
    let mut calls: Vec<gatk_tools::sv_discovery::CnvCall> = Vec::new();
    for record in &file.records {
        let genotype = record
            .genotypes
            .iter()
            .find(|genotype| genotype.sample_name == sample);
        let field = |key: &str| {
            genotype
                .and_then(|genotype| genotype.extended.iter().find(|(name, _)| name == key))
                .and_then(|(_, value)| value.format())
                .unwrap_or_else(|| "null".to_string())
        };
        let call = gatk_tools::sv_discovery::CnvCall {
            contig: header
                .sequences
                .iter()
                .position(|sequence| sequence.name == record.contig)
                .map_or(-1, |index| index as i32),
            start: record.start as i32,
            end: record.stop as i32,
            annotation: format!("{}:{}:{}", record.id, field("CN"), field("CNQ")),
        };
        // `SVIntervalTree.put` replaces the value at an interval it already holds.
        match calls
            .iter_mut()
            .find(|held| (held.contig, held.start, held.end) == (call.contig, call.start, call.end))
        {
            Some(held) => *held = call,
            None => calls.push(call),
        }
    }
    calls.sort_by_key(|call| (call.contig, call.start, call.end));
    Ok(calls)
}

/// `ReferenceMultiSparkSource.getBases`.
struct SparkReference {
    fasta: htsjdk_bam::fasta_index::IndexedFasta<htsjdk_bam::fasta_index::PositionalFile>,
}

impl gatk_tools::sv_discovery::Reference for SparkReference {
    fn bases(&mut self, contig: &str, start: i32, end: i32) -> Result<Vec<u8>, String> {
        self.fasta
            .query(contig, i64::from(start), i64::from(end))
            .map_err(|error| error.message())
    }
}

/// `PathSeqScoreSpark`, whose reads come from its own two arguments.
pub fn pathseq_score_spark(parser: &Parser) -> Outcome {
    use gatk_tools::pathseq_score as score;
    let tool = "PathSeqScoreSpark";
    resolve_read_filters(parser, tool)?;
    spark_context(parser)?;
    if !arguments(parser, "input").is_empty() {
        return Err(bad_input(
            "Please use --paired-input or --unpaired-input instead of --input".to_string(),
        ));
    }
    // `readInputWithHeader`, which skips a path that is not there.
    let read_input = |name: &str| -> Result<Option<(Vec<BamRecord>, SamHeader)>, Thrown> {
        let Some(path) = argument(parser, name) else {
            return Ok(None);
        };
        if !std::path::Path::new(&path).exists() {
            return Ok(None);
        }
        let source = spark_reads_header(&path)?;
        let header = source.header().clone();
        let records = gatk_tools::read_walker::traverse(&source, &[], &|_| true)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        Ok(Some((score::primary_reads(records), header)))
    };
    let paired = read_input("paired-input")?;
    let unpaired = read_input("unpaired-input")?;
    if let Some((_, header)) = &paired {
        if header.attributes.get("SO") != Some("queryname") {
            return Err(bad_input(
                "Paired input BAM must be sorted by queryname".to_string(),
            ));
        }
    }
    // `joinBamHeaders`: what the unpaired input adds, unless the two are equal.
    let header = match (&paired, &unpaired) {
        (Some((_, first)), Some((_, second))) => {
            let mut joined = first.clone();
            for sequence in &second.sequences {
                if !joined
                    .sequences
                    .iter()
                    .any(|known| known.name == sequence.name)
                {
                    joined.sequences.push(sequence.clone());
                }
            }
            for group in &second.read_groups {
                if !joined.read_groups.iter().any(|known| known.id == group.id) {
                    joined.read_groups.push(group.clone());
                }
            }
            joined
        }
        (Some((_, only)), None) | (None, Some((_, only))) => only.clone(),
        (None, None) => return Err(bad_input("No headers were loaded".to_string())),
    };
    let contigs: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let database = pathseq_taxonomy(parser)?;
    let write = |path: &str, text: &str| {
        std::fs::write(path, text).map_err(|error| Thrown::user(format!("{path}: {error}")))
    };
    if let Some(path) = argument(parser, "score-warnings") {
        write(&path, &score::missing_accessions(&contigs, &database))?;
    }
    // A SAMRecord finds its contig by name, so the second input's indices move to the join's.
    let unpaired = unpaired.map(|(mut records, own)| {
        let moved = |index: i32| {
            usize::try_from(index)
                .ok()
                .and_then(|index| own.sequences.get(index))
                .and_then(|sequence| contigs.iter().position(|name| *name == sequence.name))
                .map_or(-1, |at| at as i32)
        };
        for record in &mut records {
            record.reference_index = moved(record.reference_index);
            record.mate_reference_index = moved(record.mate_reference_index);
        }
        records
    });
    let scored = score::score_reads(
        paired.map(|(records, _)| records),
        unpaired,
        &contigs,
        &database,
        double_or(parser, "min-score-identity", 0.90),
        double_or(parser, "identity-margin", 0.02),
        flag(parser, "divide-by-genome-length"),
        flag(parser, "not-normalized-by-kingdom"),
    )
    .map_err(|error| match error {
        score::ScoreError::BadInput(message) => {
            pathseq_task_failure(format!("Bad input: {message}"))
        }
        score::ScoreError::Other { message, .. } => pathseq_task_failure(message),
    })?;
    write(
        &argument(parser, "scores-output").expect("a required argument"),
        &scored.scores,
    )?;
    if let Some(path) = argument(parser, "score-metrics") {
        let headers = [
            crate::command_line::expanded(tool, parser),
            format!("Started on: {}", java_display_now()),
        ];
        write(&path, &score::metrics_file(&scored.metrics, &headers))?;
    }
    // Not `writeReads`: the reducers are each input's own added up, this overload of the sink writes
    // both indices whatever the arguments say, and with no header merger there is no `GO` to add.
    if let Some(output) = argument(parser, "output") {
        let mut reducers = 0;
        for input in ["paired-input", "unpaired-input"] {
            if let Some(path) = argument(parser, input) {
                if std::path::Path::new(&path).exists() {
                    reducers += pathseq_reducers(parser, &path);
                }
            }
        }
        let unpaired = scored.reads.len() - scored.paired;
        write_pathseq_bam(
            parser,
            &header,
            &scored.reads,
            &output,
            reducers,
            &[scored.paired, unpaired],
        )?;
    }
    Ok(None)
}

/// `PSUtils.pathseqGetRecommendedNumReducers`.
fn pathseq_reducers(parser: &Parser, path: &str) -> usize {
    let asked = number_or(parser, "num-reducers", 0).max(0) as usize;
    if asked != 0 {
        return asked;
    }
    let size = std::fs::metadata(path).map_or(0, |data| data.len());
    1 + (size / (10 * 1024 * 1024)) as usize
}

/// `ReadsSparkSink.writeReads`, sorting to the header's order and writing both indices.
fn write_pathseq_bam(
    parser: &Parser,
    header: &SamHeader,
    reads: &[BamRecord],
    output: &str,
    reducers: usize,
    input_partitions: &[usize],
) -> Result<(), Thrown> {
    let sort_order = header.attributes.get("SO").map(str::to_string);
    let compare: Option<fn(&BamRecord, &BamRecord) -> std::cmp::Ordering> =
        match sort_order.as_deref() {
            Some("coordinate") => Some(htsjdk_bam::coordinate::compare),
            Some("queryname") => Some(htsjdk_bam::query_name::compare),
            _ => None,
        };
    let partitions = match compare {
        Some(compare) => range_partitions(reads, reducers, input_partitions, compare)?,
        None => vec![reads.to_vec()],
    };
    write_spark_indexed(
        parser,
        header,
        sort_order.as_deref(),
        &partitions,
        output,
        (true, true),
    )
}

/// `new PSFilter`'s arguments, over a dictionary `checkAndClearHeaderSequences` emptied.
fn pathseq_filter_arguments(
    parser: &Parser,
) -> Result<gatk_tools::pathseq_filter::FilterArguments, Thrown> {
    let contigs_to_ignore = arguments(parser, "ignore-alignment-contigs");
    let aligned_input = flag(parser, "is-host-aligned");
    if aligned_input {
        if let Some(contig) = contigs_to_ignore.first() {
            return Err(bad_input(format!(
                "Ignored sequence {contig} not found in input header."
            )));
        }
    }
    Ok(gatk_tools::pathseq_filter::FilterArguments {
        aligned_input,
        contigs_to_ignore,
        min_identity: number_or(parser, "host-min-identity", 30),
        skip_filters: flag(parser, "skip-quality-filters"),
        max_adapter_mismatches: number_or(parser, "max-adapter-mismatches", 1),
        min_adapter_length: number_or(parser, "min-adapter-length", 12),
        dust_mask: number_or(parser, "dust-mask-quality", 2),
        dust_window: number_or(parser, "dust-window", 64),
        dust_t: double_or(parser, "dust-t", 20.0),
        read_trim_threshold: number_or(parser, "quality-threshold", 15),
        min_read_length: number_or(parser, "min-clipped-read-length", 31),
        quality_threshold: number_or(parser, "min-base-quality", 15),
        max_ambiguous_bases: number_or(parser, "max-masked-bases", 2),
        host_kmer_threshold: number_or(parser, "host-kmer-thresh", 1),
        filter_duplicates: flag(parser, "filter-duplicates"),
        reads_per_partition: number_or(parser, "filter-reads-per-partition", 200000),
    })
}

/// `PSFilter.doFilter` and `PSFilterFileLogger.close`, which the try-with-resources runs before
/// anything is written.
fn pathseq_run_filter(
    parser: &Parser,
    tool: &str,
    records: Vec<BamRecord>,
    contigs: &[String],
    arguments: &gatk_tools::pathseq_filter::FilterArguments,
) -> Result<gatk_tools::pathseq_filter::Filtered, Thrown> {
    use gatk_tools::pathseq_filter as filter;
    let kmer_file = argument(parser, "kmer-file");
    let load_kmers = kmer_file.as_ref().map(|path| {
        let path = path.clone();
        move || -> Result<filter::KmerLibrary, String> {
            if !path.ends_with(".hss") && !path.ends_with(".bfi") {
                return Err(format!(
                    "Bad input: Unknown kmer set extension in file name {path}"
                ));
            }
            let bytes =
                std::fs::read(&path).map_err(|_| format!("{path} (No such file or directory)"))?;
            filter::read_kmer_library(&path, &bytes)
        }
    });
    let load_kmers = load_kmers
        .as_ref()
        .map(|load| load as &dyn Fn() -> Result<filter::KmerLibrary, String>);
    let host_image = argument(parser, "filter-bwa-image");
    let min_identity = arguments.min_identity;
    let seed_length = number_or(parser, "filter-bwa-seed-length", 19);
    let host_filter = host_image.as_ref().map(|image| {
        move |reads: Vec<BamRecord>| -> Result<Vec<BamRecord>, String> {
            // `new PSBwaFilter`, which every partition makes whatever it holds.
            let bytes = std::fs::read(image)
                .ok()
                .filter(|bytes| !bytes.is_empty())
                .ok_or_else(|| {
                    format!("cannot read the image file '{image}': is empty or is not readable")
                })?;
            if reads.is_empty() {
                return Ok(reads);
            }
            let threads = std::thread::available_parallelism().map_or(1, |count| count.get());
            // The lines name the image's contigs, which the header no longer holds.
            let contigs = gatk_tools::bwa_aligner::image_contigs(&bytes)?;
            let names = SamHeader {
                sequences: contigs
                    .iter()
                    .map(|name| htsjdk_bam::header::SequenceRecord::new(name, 1))
                    .collect(),
                ..SamHeader::default()
            };
            gatk_tools::bwa_aligner::host_filter(
                reads,
                &bytes,
                &names,
                seed_length,
                |line| filter::host_alignment_keeps(line, min_identity),
                threads,
            )
        }
    });
    let host_filter = host_filter
        .as_ref()
        .map(|run| run as &filter::HostFilter<'_>);
    let filtered = filter::do_filter(records, contigs, arguments, load_kmers, host_filter)
        .map_err(pathseq_task_failure)?;
    if let Some(path) = argument(parser, "filter-metrics") {
        let headers = [
            crate::command_line::expanded(tool, parser),
            format!("Started on: {}", java_display_now()),
        ];
        let text = filter::metrics_file(&filtered.metrics, &headers);
        std::fs::write(&path, text).map_err(|error| Thrown::user(format!("{path}: {error}")))?;
    }
    Ok(filtered)
}

/// `PathSeqPipelineSpark.runTool`: the filter, the microbe alignment and the scores in one run.
pub fn pathseq_pipeline_spark(parser: &Parser) -> Outcome {
    use gatk_tools::pathseq_filter::spark_repartition;
    use gatk_tools::pathseq_score as score;
    let tool = "PathSeqPipelineSpark";
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    if number_or(parser, "num-reducers", 0) > 0 {
        return Err(bad_input(
            "Use --readsPerPartitionOutput instead of --num-reducers.".to_string(),
        ));
    }
    let filter_arguments = pathseq_filter_arguments(parser)?;
    let filter_reads = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter_reads)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let contigs: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let filtered = pathseq_run_filter(parser, tool, records, &contigs, &filter_arguments)?;
    let total = filtered.paired.len() + filtered.unpaired.len();

    // `repartition`, which moves pairs as a whole.
    let per_partition = number_or(parser, "pipeline-reads-per-partition", 5000) as usize;
    let paired_partitions = 1 + filtered.paired.len() / per_partition;
    let unpaired_partitions = 1 + filtered.unpaired.len() / per_partition;
    let pairs: Vec<Vec<BamRecord>> = filtered.paired.chunks(2).map(<[_]>::to_vec).collect();
    let paired_batches: Vec<Vec<BamRecord>> = spark_repartition(pairs, paired_partitions)
        .into_iter()
        .map(|pairs| pairs.into_iter().flatten().collect())
        .collect();
    let unpaired_batches = spark_repartition(filtered.unpaired, unpaired_partitions);

    let mut microbe = header.clone();
    microbe.sequences = pathseq_microbe_dictionary(parser)?;
    let names: Vec<String> = microbe.sequences.iter().map(|s| s.name.clone()).collect();

    let database = pathseq_taxonomy(parser)?;
    let write = |path: &str, text: &str| {
        std::fs::write(path, text).map_err(|error| Thrown::user(format!("{path}: {error}")))
    };
    if let Some(path) = argument(parser, "score-warnings") {
        write(&path, &score::missing_accessions(&names, &database))?;
    }

    let image = argument(parser, "microbe-bwa-image").expect("a required argument");
    let bytes = bwa_image(&image)?;
    let options = pathseq_bwa_options(parser);
    let threads = std::thread::available_parallelism().map_or(1, |count| count.get());
    let align = |batches: Vec<Vec<BamRecord>>, paired: bool| -> Result<Vec<BamRecord>, Thrown> {
        let mut aligned = Vec::new();
        for batch in &batches {
            aligned.extend(
                gatk_tools::bwa_aligner::align_for_pathseq(
                    batch, paired, &bytes, &microbe, &options, threads,
                )
                .map_err(pathseq_task_failure)?,
            );
        }
        Ok(aligned)
    };
    let aligned_paired = align(paired_batches, true)?;
    let aligned_unpaired = align(unpaired_batches, false)?;

    let scored = score::score_reads(
        Some(aligned_paired),
        Some(aligned_unpaired),
        &names,
        &database,
        double_or(parser, "min-score-identity", 0.90),
        double_or(parser, "identity-margin", 0.02),
        flag(parser, "divide-by-genome-length"),
        flag(parser, "not-normalized-by-kingdom"),
    )
    .map_err(|error| match error {
        score::ScoreError::BadInput(message) => {
            pathseq_task_failure(format!("Bad input: {message}"))
        }
        score::ScoreError::Other { message, .. } => pathseq_task_failure(message),
    })?;
    write(
        &argument(parser, "scores-output").expect("a required argument"),
        &scored.scores,
    )?;
    let mut reads = scored.reads;
    pathseq_used_sequences(
        &mut microbe,
        &mut reads,
        paired_partitions + unpaired_partitions,
    );
    if let Some(path) = argument(parser, "score-metrics") {
        let headers = [
            crate::command_line::expanded(tool, parser),
            format!("Started on: {}", java_display_now()),
        ];
        write(&path, &score::metrics_file(&scored.metrics, &headers))?;
    }
    if let Some(output) = argument(parser, "output") {
        let per_output = number_or(parser, "readsPerPartitionOutput", 1000000) as usize;
        let reducers = (total / per_output).max(1);
        // The header went through a broadcast, and Kryo gives it a group order.
        let (microbe, _) = disq_header(&microbe);
        write_pathseq_bam(parser, &microbe, &reads, &output, reducers, &[reads.len()])?;
    }
    Ok(None)
}

/// `PSBwaUtils.addReferenceSequencesToHeader`'s dictionary.
fn pathseq_microbe_dictionary(
    parser: &Parser,
) -> Result<Vec<htsjdk_bam::header::SequenceRecord>, Thrown> {
    let dictionary = argument(parser, "microbe-dict").expect("a required argument");
    // `ReferenceUtils.loadFastaDictionary` over a `GATKPath`, whose stream names the absolute path.
    let text = std::fs::read_to_string(&dictionary).map_err(|_| {
        let absolute = std::path::absolute(&dictionary)
            .map_or(dictionary.clone(), |path| path.display().to_string());
        Thrown::non_user("java.lang.RuntimeException", absolute)
    })?;
    Ok(htsjdk_bam::reader::parse_header_text(&text).sequences)
}

/// `PSScorer.readTaxonomyDatabase`.
fn pathseq_taxonomy(
    parser: &Parser,
) -> Result<gatk_tools::pathseq_score::TaxonomyDatabase, Thrown> {
    let database = argument(parser, "taxonomy-file").expect("a required argument");
    let bytes = std::fs::read(&database).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file. Error was: {database} (No such file or directory)"
        ))
    })?;
    gatk_tools::pathseq_score::read_taxonomy_database(&bytes)
        .map_err(|message| Thrown::non_user("com.esotericsoftware.kryo.KryoException", message))
}

/// `PathSeqFilterSpark.runTool`.
pub fn pathseq_filter_spark(parser: &Parser) -> Outcome {
    let tool = "PathSeqFilterSpark";
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    let arguments = pathseq_filter_arguments(parser)?;
    let filter_reads = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter_reads)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let contigs: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let filtered = pathseq_run_filter(parser, tool, records, &contigs, &arguments)?;

    let mut cleared = header.clone();
    cleared.sequences.clear();
    if !filtered.paired.is_empty() {
        cleared.set_sort_order("queryname");
        let output = argument(parser, "paired-output").ok_or_else(null_pathseq_output)?;
        write_spark_reads(parser, &cleared, &filtered.paired, &output, true)?;
    }
    if !filtered.unpaired.is_empty() {
        cleared.set_sort_order("unsorted");
        let output = argument(parser, "unpaired-output").ok_or_else(null_pathseq_output)?;
        write_spark_reads(parser, &cleared, &filtered.unpaired, &output, true)?;
    }
    Ok(None)
}

/// A refusal inside a task, which Spark's exception carries by message.
fn pathseq_task_failure(message: String) -> Thrown {
    eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
    Thrown::non_user("org.apache.spark.SparkException", message)
}

/// `writeReads` handed the output nobody named.
fn null_pathseq_output() -> Thrown {
    Thrown::non_user(
        "java.lang.IllegalArgumentException",
        "Null object is not allowed here.".to_string(),
    )
}

/// `FindBadGenomicKmersSpark.runTool`: no reads are needed, and the reference is read whole.
pub fn find_bad_genomic_kmers_spark(parser: &Parser) -> Outcome {
    use gatk_tools::find_bad_genomic_kmers as kmers;
    let tool = "FindBadGenomicKmersSpark";
    resolve_read_filters(parser, tool)?;
    spark_context(parser)?;
    let dictionary = reference_dictionary(parser)?.ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
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
    interval_arguments(parser, &dictionary)?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(reads) = &reads {
            validate_reference_against_reads(&dictionary, reads)?;
        }
    }
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let k = number_or(parser, "k-size", 51);
    let max_dust = number_or(parser, "kmer-max-dust-score", 49);

    let illegal =
        |message: &str| Thrown::non_user("java.lang.IllegalArgumentException", message.to_string());
    if k == 0 {
        return Err(illegal("provided kmer size is zero"));
    }
    if i64::from(k) > kmers::REF_RECORD_LEN as i64 {
        // `start += effectiveRecLen` does not advance: the reference fills the heap.
        return Err(Thrown::non_user(
            "java.lang.OutOfMemoryError",
            "Java heap space",
        ));
    }
    let reference_path = argument(parser, "reference").expect("the dictionary came from it");
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let mut chunks = Vec::new();
    for sequence in &dictionary.sequences {
        let bases = reference
            .query(&sequence.name, 1, sequence.length)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        kmers::chunk_contig(&bases, k, &mut chunks);
    }
    if k < 0 {
        return Err(illegal("provided kmer size is non positive"));
    }
    if max_dust <= 0 {
        return Err(illegal("provided DUST filter score is non positive"));
    }
    let task_failure = |message: String| {
        eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
        Thrown::non_user(
            "org.apache.spark.SparkException",
            format!(
                "Job aborted due to stage failure: Task 0 in stage 0.0 failed 1 times: {message}"
            ),
        )
    };
    if k < 4 {
        return Err(task_failure(
            "org.broadinstitute.hellbender.exceptions.GATKException: kmer size must be at least 4 \
             for this filter to work properly."
                .to_string(),
        ));
    }
    const EVEN: &str = "Kmer length must be odd to canonicalize.";
    let mut kill_list = kmers::collect_ubiquitous_kmers(k, max_dust, kmers::MAX_KMER_FREQ, &chunks)
        .map_err(|failure| match failure {
            kmers::KmerFailure::EvenSize => {
                task_failure(format!("java.lang.IllegalArgumentException: {EVEN}"))
            }
            other => task_failure(format!("java.lang.IllegalStateException: {other:?}")),
        })?;
    if let Some(fasta) = argument(parser, "high-copy-fasta") {
        let lines = high_copy_lines(&fasta)?;
        let extra = kmers::process_fasta(&lines, k, max_dust).map_err(|failure| match failure {
            kmers::KmerFailure::EvenSize => illegal(EVEN),
            kmers::KmerFailure::EmptyLine => Thrown::non_user(
                "java.lang.StringIndexOutOfBoundsException",
                "String index out of range: 0",
            ),
            kmers::KmerFailure::Table(failed) => {
                Thrown::non_user("java.lang.IllegalStateException", failed.0)
            }
        })?;
        kill_list = kmers::uniquify(&kill_list, &extra)
            .map_err(|failed| Thrown::non_user("java.lang.IllegalStateException", failed.0))?;
    }
    std::fs::write(&output, kmers::kmers_file(&kill_list, k)).map_err(|error| {
        Thrown::user(format!(
            "Couldn't write file {output} because Could not create file at path:{output} due to \
             {output} ({})",
            java_io_reason(&error)
        ))
    })?;
    Ok(None)
}

/// `BucketUtils.openFile` then `readLine`: a block compressed name is inflated, and each UTF-16 unit
/// outside ASCII is one character that is not a base.
fn high_copy_lines(path: &str) -> Result<Vec<Vec<u8>>, Thrown> {
    let bytes = std::fs::read(path).map_err(|error| {
        Thrown::user(format!(
            "Couldn't read file {path}. Error was: {path} ({})",
            java_io_reason(&error)
        ))
    })?;
    let bytes = if path.ends_with(".gz") || path.ends_with(".bgz") || path.ends_with(".bgzf") {
        let mut inflated = Vec::new();
        std::io::Read::read_to_end(
            &mut flate2::read::MultiGzDecoder::new(bytes.as_slice()),
            &mut inflated,
        )
        .map_err(|error| Thrown::user(format!("Couldn't read file {path}. Error was: {error}")))?;
        inflated
    } else {
        bytes
    };
    let text = String::from_utf8_lossy(&bytes);
    let mut lines = Vec::new();
    let mut line = Vec::new();
    let mut chars = text.encode_utf16().peekable();
    let mut pending = false;
    while let Some(unit) = chars.next() {
        pending = true;
        match unit {
            0x0a => {
                lines.push(std::mem::take(&mut line));
                pending = false;
            }
            0x0d => {
                if chars.peek() == Some(&0x0a) {
                    chars.next();
                }
                lines.push(std::mem::take(&mut line));
                pending = false;
            }
            unit if unit < 0x80 => line.push(unit as u8),
            _ => line.push(b'?'),
        }
    }
    if pending {
        lines.push(line);
    }
    Ok(lines)
}

/// What `FileInputStream` puts in the parentheses of its `FileNotFoundException`.
fn java_io_reason(error: &std::io::Error) -> &'static str {
    match error.kind() {
        std::io::ErrorKind::NotFound => "No such file or directory",
        std::io::ErrorKind::PermissionDenied => "Permission denied",
        _ => "Is a directory",
    }
}

/// `CollectAllelicCountsSpark.processAlignments`: `CollectAllelicCounts`' table, every locus of
/// every interval, over `LocusWalkerSpark`'s shards as in [`pileup_spark`].
pub fn collect_allelic_counts_spark(parser: &Parser) -> Outcome {
    let tool = "CollectAllelicCountsSpark";
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, tool)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let shard_size = number_or(parser, "read-shard-size", 10000).max(1);
    no_reads_refusal(parser, &records, &intervals, shard_size)?;
    let sample = sample_name(&header)?;
    oversized_read_refusal(parser, &records, &header, shard_size)?;
    let mut applied = Vec::new();
    for interval in &intervals {
        for start in (interval.start..=interval.end).step_by(shard_size as usize) {
            let shard = gatk_engine::interval::SimpleInterval {
                contig: interval.contig.clone(),
                start,
                end: (start + shard_size - 1).min(interval.end),
            };
            // The shuffle groups reads by shard, so a shard no read reaches never exists.
            if flag(parser, "shuffle")
                && gatk_tools::read_walker::traverse(&source, std::slice::from_ref(&shard), &filter)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?
                    .is_empty()
            {
                continue;
            }
            applied.extend(
                gatk_tools::locus_walker::traverse(
                    &records,
                    &header,
                    None,
                    Some(std::slice::from_ref(&shard)),
                    gatk_tools::locus_walker::Options {
                        max_depth_per_sample: number_or(parser, "max-depth-per-sample", 0),
                        emit_empty_loci: true,
                        ..gatk_tools::locus_walker::Options::default()
                    },
                    &|_| true,
                )
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
            );
        }
    }
    allelic_counts_table(parser, &header, &mut reference, &applied, &sample, &output)
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

/// `BwaSparkEngine`: its header, and BWA over one partition of reads.
///
/// `header` is the reads' when it has a dictionary, otherwise the reference's, and the image is
/// the one named, otherwise the reference's own. A paired run takes every read and refuses any
/// order but `queryname`.
fn bwa_engine(parser: &Parser, header: &SamHeader) -> Result<(SamHeader, String), Thrown> {
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut aligned_header = header.clone();
    if aligned_header.sequences.is_empty() {
        if let Some(dictionary) = reference_dictionary(parser)? {
            aligned_header.sequences = dictionary.sequences;
        }
    }
    let image = match argument(parser, "bwa-mem-index-image") {
        Some(image) => image,
        None => {
            // `ctx.addFile`, which looks before anything runs.
            let image = format!("{reference}.img");
            if !std::path::Path::new(&image).exists() {
                let absolute = std::path::absolute(&image)
                    .map_or(image.clone(), |path| path.display().to_string());
                return Err(Thrown::non_user(
                    "java.io.FileNotFoundException",
                    format!("File file:{absolute} does not exist"),
                ));
            }
            image
        }
    };
    Ok((aligned_header, image))
}

/// `ReadAligner`, past its construction.
fn bwa_align(
    reads: &[BamRecord],
    paired: bool,
    image: &str,
    header: &SamHeader,
) -> Result<Vec<BamRecord>, Thrown> {
    // `BwaMemIndexCache.getInstance` comes first in the constructor, in the executor.
    let bytes = std::fs::read(image)
        .ok()
        .filter(|bytes| !bytes.is_empty())
        .ok_or_else(|| {
            evidence_task_failure(&format!(
                "cannot read the image file '{image}': is empty or is not readable"
            ))
        })?;
    if paired && header.attributes.get("SO") != Some("queryname") {
        return Err(evidence_task_failure(
            "Input must be queryname sorted unless you use single-ended alignment mode.",
        ));
    }
    let threads = std::thread::available_parallelism().map_or(1, |count| count.get());
    gatk_tools::bwa_aligner::align(reads, paired, &bytes, header, threads)
        .map_err(|message| evidence_task_failure(&message))
}

/// `BwaSpark.runTool`: BWA over the reads, in one partition, sorted back to the header.
///
/// Single-end alignment filters first. Paired alignment filters after, since a filter does not
/// know about pairs, which is how a supplementary line of a pair never reaches the output.
pub fn bwa_spark(parser: &Parser) -> Outcome {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "BwaSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let (aligned_header, image) = bwa_engine(parser, &header)?;
    let filter = read_filter(parser, &filters, &aligned_header)?;
    let paired = !flag(parser, "single-end-alignment");
    let reads = if paired {
        gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
    } else {
        gatk_tools::read_walker::traverse(&source, &intervals, &filter)
    }
    .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let aligned = bwa_align(&reads, paired, &image, &aligned_header)?;
    let kept: Vec<_> = if paired {
        aligned.into_iter().filter(|read| filter(read)).collect()
    } else {
        aligned
    };
    // The sink's argument after the parts directory is `shardedOutput`: a single file is not sorted,
    // which leaves each read's primary line ahead of the supplementary one the comparator would move.
    write_spark_reads(
        parser,
        &aligned_header,
        &kept,
        &output,
        flag(parser, "sharded-output"),
    )?;
    Ok(None)
}

/// `BwaAndMarkDuplicatesPipelineSpark.runTool`: paired alignment, the tool's filter, and
/// `MarkDuplicatesSpark.mark` with a default `OpticalDuplicateFinder`, written in the header's own
/// order: this pipeline never sets it to coordinate, and removes nothing whatever it is asked.
pub fn bwa_and_mark_duplicates_pipeline_spark(parser: &Parser) -> Outcome {
    use gatk_tools::mark_duplicates_spark as md;
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = spark_startup(parser, "BwaAndMarkDuplicatesPipelineSpark")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let (aligned_header, image) = bwa_engine(parser, &header)?;
    let filter = read_filter(parser, &filters, &aligned_header)?;
    let reads = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let aligned: Vec<_> = bwa_align(&reads, true, &image, &aligned_header)?
        .into_iter()
        .filter(|read| filter(read))
        .collect();
    drop(filter);
    let settings = mark_arguments(parser, parser, false, None, 100)?;
    if aligned_header.read_groups.is_empty() {
        return Err(Thrown::user(
            "Bad input: Sam file header missing Read Group fields. MarkDuplicatesSpark currently \
             requires reads to be labeled with read group tags, please add read groups tags to \
             your reads",
        ));
    }
    let marked = md::mark(vec![aligned], &aligned_header, &settings)
        .map_err(|failure| evidence_task_failure(&failure.message))?;
    write_spark_reads(
        parser,
        &aligned_header,
        &marked.reads,
        &output,
        !flag(parser, "sharded-output"),
    )?;
    Ok(None)
}

/// `PSBwaAligner`'s setters, as `bwa mem` spells them.
fn pathseq_bwa_options(parser: &Parser) -> Vec<String> {
    let hits = number_or(parser, "max-alternate-hits", 5000);
    vec![
        "-h".to_string(),
        format!("{hits},{hits}"),
        "-k".to_string(),
        number_or(parser, "microbe-min-seed-length", 19).to_string(),
        "-T".to_string(),
        number_or(parser, "bwa-score-threshold", 30).to_string(),
    ]
}

/// `BwaMemIndexCache.getInstance`, in the executor.
fn bwa_image(image: &str) -> Result<Vec<u8>, Thrown> {
    std::fs::read(image)
        .ok()
        .filter(|bytes| !bytes.is_empty())
        .ok_or_else(|| {
            pathseq_task_failure(format!(
                "cannot read the image file '{image}': is empty or is not readable"
            ))
        })
}

/// `PSBwaUtils.removeUnmappedHeaderSequences`: the contigs some read landed on, in the order
/// `distinct` over `partitions` partitions collects them, and the reads' indices moved to the
/// shorter dictionary.
fn pathseq_used_sequences(header: &mut SamHeader, reads: &mut [BamRecord], partitions: usize) {
    let hash = |name: &String| gatk_engine::java_hash::string_hash_code(name);
    let mut buckets: Vec<Vec<(String, ())>> = vec![Vec::new(); partitions];
    for read in reads.iter() {
        if gatk_engine::read::is_unmapped(read) {
            continue;
        }
        let Some(sequence) = usize::try_from(read.reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
        else {
            continue;
        };
        let bucket = hash(&sequence.name).rem_euclid(partitions as i32) as usize;
        buckets[bucket].push((sequence.name.clone(), ()));
    }
    let reduced: Vec<String> = buckets
        .into_iter()
        .flat_map(|names| {
            gatk_tools::sv_discovery::append_only_group(names, hash)
                .into_iter()
                .map(|(name, _)| name)
        })
        .collect();
    let used: Vec<_> = reduced
        .into_iter()
        .filter_map(|name| {
            header
                .sequences
                .iter()
                .find(|sequence| sequence.name == name)
                .cloned()
        })
        .collect();
    let moved = |index: i32, from: &[htsjdk_bam::header::SequenceRecord]| {
        usize::try_from(index)
            .ok()
            .and_then(|index| from.get(index))
            .and_then(|sequence| used.iter().position(|kept| kept.name == sequence.name))
            .map_or(-1, |at| at as i32)
    };
    for read in reads.iter_mut() {
        read.reference_index = moved(read.reference_index, &header.sequences);
        read.mate_reference_index = moved(read.mate_reference_index, &header.sequences);
    }
    header.sequences = used;
}

/// `ParallelCopyGCSDirectoryIntoHDFSSpark.runTool`: what a machine with neither service can see.
pub fn parallel_copy_gcs_directory_into_hdfs_spark(parser: &Parser) -> Outcome {
    resolve_read_filters(parser, "ParallelCopyGCSDirectoryIntoHDFSSpark")?;
    spark_context(parser)?;
    let input = argument(parser, "input-gcs-path").expect("a required argument");
    if !input.starts_with("gs://") {
        return Err(Thrown::user(format!("Input path {input} is not a GCS URI")));
    }
    let output = argument(parser, "output-hdfs-directory").expect("a required argument");
    if !output.starts_with("hdfs://") {
        return Err(Thrown::user(format!(
            "Output directory {output} is not an HDFS URI"
        )));
    }
    Err(Thrown::non_user(
        PORT_LIMITATION,
        "ParallelCopyGCSDirectoryIntoHDFSSpark needs a Google Cloud Storage and an HDFS, which \
         this port does not reach. This message is the port's own and not GATK's.",
    ))
}

/// `PathSeqBwaSpark.runTool`: each input aligned on its own and written with only the contigs its
/// reads landed on.
pub fn pathseq_bwa_spark(parser: &Parser) -> Outcome {
    let tool = "PathSeqBwaSpark";
    resolve_read_filters(parser, tool)?;
    spark_context(parser)?;
    if !arguments(parser, "input").is_empty() {
        return Err(bad_input(
            "Please use --paired-input or --unpaired-input instead of --input".to_string(),
        ));
    }
    for output in ["paired-output", "unpaired-output"] {
        if argument(parser, output).is_some_and(|path| !path.ends_with(".bam")) {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                "Only BAM output is supported",
            ));
        }
    }
    let image = argument(parser, "microbe-bwa-image").expect("a required argument");
    let options = pathseq_bwa_options(parser);
    let threads = std::thread::available_parallelism().map_or(1, |count| count.get());
    let mut loaded = false;
    for (input, output, paired) in [
        ("paired-input", "paired-output", true),
        ("unpaired-input", "unpaired-output", false),
    ] {
        let Some(path) = argument(parser, input) else {
            continue;
        };
        if !std::path::Path::new(&path).exists() {
            continue;
        }
        loaded = true;
        let source = spark_reads_header(&path)?;
        let mut header = source.header().clone();
        if !header.sequences.is_empty() {
            return Err(bad_input(
                "Input BAM should be unaligned, but found one or more sequences in the header."
                    .to_string(),
            ));
        }
        header.sequences = pathseq_microbe_dictionary(parser)?;
        if paired && header.attributes.get("SO") != Some("queryname") {
            return Err(bad_input(
                "Paired input BAM must be sorted by queryname".to_string(),
            ));
        }
        let reads = gatk_tools::read_walker::traverse(&source, &[], &|_| true)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let bytes = bwa_image(&image)?;
        if paired && reads.len() % 2 != 0 {
            return Err(pathseq_task_failure(
                "Bad input: Expected paired reads but there are an odd number".to_string(),
            ));
        }
        let mut aligned = gatk_tools::bwa_aligner::align_for_pathseq(
            &reads, paired, &bytes, &header, &options, threads,
        )
        .map_err(pathseq_task_failure)?;
        pathseq_used_sequences(&mut header, &mut aligned, 1);
        let output = argument(parser, output).ok_or_else(null_pathseq_output)?;
        let reducers = pathseq_reducers(parser, &path);
        write_pathseq_bam(
            parser,
            &header,
            &aligned,
            &output,
            reducers,
            &[aligned.len()],
        )?;
    }
    if !loaded {
        return Err(bad_input(
            "No reads were loaded. Ensure --paired-input and/or --unpaired-input are set and valid."
                .to_string(),
        ));
    }
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
    if !output.ends_with(".bam") && !output.ends_with(".sam") && !sharded {
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
    let start = spark_startup(parser, tool)?;
    spark_metric_reads_of(parser, tool, &start)
}

/// [`spark_metric_reads`] past the startup, for a tool with a check of its own between the two.
fn spark_metric_reads_of(
    parser: &Parser,
    tool: &str,
    start: &SparkStart,
) -> Result<(Vec<BamRecord>, Vec<String>), Thrown> {
    let SparkStart {
        source,
        header,
        intervals,
        filters,
    } = start;
    let filter = read_filter(parser, filters, header)?;
    let records = gatk_tools::read_walker::traverse(source, intervals, &filter)
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

/// The arguments of `InsertSizeMetricsArgumentCollection`, and `ALL_READS` without a level.
fn insert_size_settings(parser: &Parser) -> gatk_tools::spark_metrics::InsertSizeArguments {
    let mut levels = arguments(parser, "metric-accumulation-level");
    if levels.is_empty() {
        levels.push("ALL_READS".to_string());
    }
    gatk_tools::spark_metrics::InsertSizeArguments {
        minimum_pct: scalar(parser, "min-category-reads-percentage")
            .and_then(|text| text.parse().ok())
            .unwrap_or(0.05),
        deviations: double_or(parser, "histogram-plot-deviations-tolerance", 10.0),
        width: scalar(parser, "width").and_then(|text| text.parse().ok()),
        levels,
    }
}

/// `CollectMultipleMetricsSpark.runTool`: the collectors in the order asked for, each over the
/// reads through its own default filters. The command line's filters are parsed and validated and
/// then bypassed, and each collector's arguments are its defaults but for the accumulation levels,
/// which only the insert sizes read.
pub fn collect_multiple_metrics_spark(parser: &Parser) -> Outcome {
    let tool = "CollectMultipleMetricsSpark";
    let start = spark_startup(parser, tool)?;
    let base = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let mut collectors = arguments(parser, "collectors");
    if collectors.is_empty() {
        collectors = vec![
            "CollectInsertSizeMetrics".to_string(),
            "CollectQualityYieldMetrics".to_string(),
        ];
    }

    let allow_all = |_: &BamRecord| true;
    let mut reads = None;
    // The stage and task the next job takes: an aggregate is one, and `isEmpty` one more.
    let mut stage = 0;
    for collector in &collectors {
        let insert_size = collector == "CollectInsertSizeMetrics";
        // `validateCollector`, naming the enum constant's own class, and with the output's name
        // chosen first.
        let order = start.header.attributes.get("SO").unwrap_or("unsorted");
        if order != "coordinate" {
            let constant = if insert_size { 1 } else { 2 };
            return Err(Thrown::user(format!(
                "Input \"org.broadinstitute.hellbender.tools.spark.pipelines.metrics.\
                 CollectMultipleMetricsSpark$SparkCollectors${constant}\" has sort order \
                 \"{order}\" but \"coordinate\" is required.If you believe the file to be \
                 sorted correctly, use assume-sorted=true"
            )));
        }
        if reads.is_none() {
            let records =
                gatk_tools::read_walker::traverse(&start.source, &start.intervals, &allow_all)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?;
            let headers = vec![
                crate::command_line::expanded(tool, parser),
                format!("Started on: {}", java_display_now()),
            ];
            reads = Some((records, headers));
        }
        let Some((records, headers)) = &reads else {
            unreachable!("the reads were just traversed");
        };
        if insert_size {
            // `InsertSizeMetricsCollector.getDefaultReadFilters`; its `MappingQualityReadFilter`
            // is parameterised at nought and keeps every read.
            let chain =
                gatk_tools::plugin_ownership::default_filters("CollectInsertSizeMetricsSpark")
                    .unwrap_or_default()
                    .iter()
                    .filter(|name| **name != "MappingQualityReadFilter")
                    .map(|name| gatk_tools::filter_resolution::ResolvedFilter {
                        name: name.to_string(),
                        negated: false,
                    })
                    .collect::<Vec<_>>();
            let filter = read_filter(parser, &chain, &start.header)?;
            let kept = records
                .iter()
                .filter(|read| filter(read))
                .cloned()
                .collect::<Vec<_>>();
            if kept.is_empty() {
                return Err(Thrown::non_user(
                    "org.broadinstitute.hellbender.exceptions.GATKException",
                    "No valid reads found in input file.",
                ));
            }
            let text = gatk_tools::spark_metrics::insert_size(
                &kept,
                &start.header,
                &insert_size_settings(parser),
                headers,
            )
            .map_err(|message| {
                let task = stage + 1;
                eprintln!("ERROR Executor - Exception in task 0.0 in stage {task}.0 (TID {task})");
                Thrown::non_user("org.apache.spark.SparkException", message)
            })?;
            if let Some(text) = text {
                write_multiple_metrics(&format!("{base}.insertSizeMetrics.txt"), &text)?;
            }
            stage += 2;
        } else {
            let text = gatk_tools::spark_metrics::quality_yield(records.iter(), false, headers);
            write_multiple_metrics(&format!("{base}.qualityYieldMetrics.txt"), &text)?;
            stage += 1;
        }
    }
    Ok(None)
}

fn write_multiple_metrics(output: &str, text: &str) -> Result<(), Thrown> {
    std::fs::write(output, text).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })
}

/// `CollectInsertSizeMetricsSpark.runTool`: `MetricsCollectorSparkTool`'s lifecycle over
/// `InsertSizeMetricsCollector`, which is GATK's own and not Picard's.
pub fn collect_insert_size_metrics_spark(parser: &Parser) -> Outcome {
    let tool = "CollectInsertSizeMetricsSpark";
    let start = spark_startup(parser, tool)?;
    // `validateExpectedSortOrder`, with `assumeSorted` false.
    let order = start.header.attributes.get("SO").unwrap_or("unsorted");
    if order != "coordinate" {
        let input = arguments(parser, "input").join(",");
        return Err(Thrown::user(format!(
            "Input \"{input}\" has sort order \"{order}\" but \"coordinate\" is required.If you \
             believe the file to be sorted correctly, use assume-sorted=true"
        )));
    }
    let (records, headers) = spark_metric_reads_of(parser, tool, &start)?;
    if records.is_empty() {
        return Err(Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            "No valid reads found in input file.",
        ));
    }
    let settings = insert_size_settings(parser);
    let text = gatk_tools::spark_metrics::insert_size(&records, &start.header, &settings, &headers)
        .map_err(|message| {
            // `isEmpty` is stage 0, and `mapPartitions` the stage after it.
            eprintln!("ERROR Executor - Exception in task 0.0 in stage 1.0 (TID 1)");
            Thrown::non_user("org.apache.spark.SparkException", message)
        })?;
    if let Some(text) = text {
        save_spark_metrics(parser, &text)?;
        if flag(parser, "produce-plot") {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                "--produce-plot runs R, which the port does not. This message is the port's own \
                 and not GATK's.",
            ));
        }
    }
    Ok(None)
}

/// `MarkDuplicatesSpark.runTool`.
///
/// One input, a queryname sort when the header does not promise the names are grouped, the marking,
/// the metrics, and the file coordinate-sorted whatever came in. Reads are one partition here: the
/// reducers only move the partition boundaries, which a name never spans.
pub fn mark_duplicates_spark(parser: &Parser, mode: &Parser) -> Outcome {
    use gatk_tools::mark_duplicates_spark as md;
    let tool = "MarkDuplicatesSpark";
    // `initializeReads`: `concatMerge`. The inputs are a map's keys, so a path given twice is read
    // once, and the checks below still see two.
    let filters = resolve_read_filters(parser, tool)?;
    spark_context(parser)?;
    let reference = reference_dictionary(parser)?;
    let given = arguments(parser, "input");
    let mut inputs: Vec<String> = Vec::new();
    for input in &given {
        if !inputs.contains(input) {
            inputs.push(input.clone());
        }
    }
    if inputs.is_empty() {
        return Err(Thrown::command_line(
            "Argument input was missing: Argument 'input' is required",
        ));
    }
    let sources = inputs
        .iter()
        .map(|input| spark_reads_header(input))
        .collect::<Result<Vec<ReadsDataSource>, Thrown>>()?;
    let own_headers: Vec<SamHeader> = sources.iter().map(|s| s.header().clone()).collect();
    let merged = gatk_tools::sam_file_header_merger::merge(&own_headers)
        .map_err(|message| Thrown::non_user("htsjdk.samtools.SAMException", message))?;
    let mut header = merged.header;
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
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let treat_unsorted = flag(parser, "treat-unsorted-as-querygroup-ordered");
    let unsorted = |input: &SamHeader| {
        matches!(
            input.attributes.get("SO").unwrap_or("unsorted"),
            "unknown" | "unsorted"
        )
    };
    if given.len() > 1 {
        let allow = flag(parser, "allow-multiple-sort-orders-in-input");
        // `treatAsReadGroupOrdered`.
        let badly_sorted = given.iter().find(|input| {
            let at = inputs.iter().position(|known| known == *input).unwrap_or(0);
            let own = &own_headers[at];
            let grouped = own.attributes.get("SO") == Some("queryname")
                || own.attributes.get("GO") == Some("query");
            !(grouped || (treat_unsorted && unsorted(own)))
        });
        match badly_sorted {
            Some(input) if !allow => {
                let at = inputs.iter().position(|known| known == input).unwrap_or(0);
                return Err(Thrown::user(format!(
                    "Multiple inputs to MarkDuplicatesSpark detected. MarkDuplicatesSpark requires \
                     all inputs to be queryname sorted or querygroup-sorted for multi-input \
                     processing but input {input} was sorted in {} order",
                    gatk_tools::sam_file_header_merger::header_to_string(&own_headers[at])
                )));
            }
            None if !allow => header.set_group_order("query"),
            _ => {}
        }
    } else if treat_unsorted && unsorted(&header) {
        header.set_group_order("query");
    }

    let sources = inputs
        .iter()
        .zip(sources)
        .map(|(input, source)| spark_traversal_source(parser, input, source))
        .collect::<Result<Vec<ReadsDataSource>, Thrown>>()?;
    let filter = read_filter(parser, &filters, &header)?;
    // `getReads`: a union, one partition an input, of the reads under the merged header. A read
    // holds its contig by name, which is what makes its index the merged dictionary's.
    let contigs: std::collections::HashMap<&str, i32> = header
        .sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.name.as_str(), index as i32))
        .collect();
    let mut partitions: Vec<Vec<BamRecord>> = Vec::new();
    for (source, own) in sources.iter().zip(&own_headers) {
        let translate = |record: &mut BamRecord| {
            let index = |old: i32| match usize::try_from(old) {
                Ok(at) => own
                    .sequences
                    .get(at)
                    .and_then(|sequence| contigs.get(sequence.name.as_str()).copied())
                    .unwrap_or(old),
                Err(_) => old,
            };
            record.reference_index = index(record.reference_index);
            record.mate_reference_index = index(record.mate_reference_index);
        };
        let mut records = if merged.merged_dictionaries {
            gatk_tools::read_walker::traverse(source, &intervals, &|record: &BamRecord| {
                let mut record = record.clone();
                translate(&mut record);
                filter(&record)
            })
        } else {
            gatk_tools::read_walker::traverse(source, &intervals, &filter)
        }
        .map_err(|error| match error {
            // Disq's query, in the executor, over a file that lacks a contig the merge has.
            gatk_engine::reads::ReadsError::ContigNotInDictionary(contig) => {
                eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
                Thrown::non_user(
                    "org.apache.spark.SparkException",
                    format!("Contig {contig} not present in reads sequence dictionary"),
                )
            }
            error => Thrown::user(format!("{error:?}")),
        })?;
        if merged.merged_dictionaries {
            records.iter_mut().for_each(translate);
        }
        partitions.push(records);
    }
    drop(filter);

    let (remove_all, remove_sequencing) = (
        flag(parser, "remove-all-duplicates"),
        flag(parser, "remove-sequencing-duplicates"),
    );
    let settings = mark_arguments(
        parser,
        mode,
        remove_sequencing,
        argument(parser, "read-name-regex").as_deref(),
        number_or(parser, "optical-duplicate-pixel-distance", 100),
    )?;

    // `getHeaderReadGroupIndexMap`, on the driver.
    if header.read_groups.is_empty() {
        return Err(Thrown::user(
            "Bad input: Sam file header missing Read Group fields. MarkDuplicatesSpark currently \
             requires reads to be labeled with read group tags, please add read groups tags to \
             your reads",
        ));
    }
    // `querynameSortReadsIfNecessary`.
    let grouped = header.attributes.get("SO") == Some("queryname")
        || header.attributes.get("GO") == Some("query");
    if !grouped {
        let mut records: Vec<BamRecord> = partitions.into_iter().flatten().collect();
        records.sort_by(htsjdk_bam::query_name::compare);
        partitions = vec![records];
    }
    let marked = md::mark(partitions, &header, &settings).map_err(|failure| {
        // Measured on one partition: a queryname sort of its own takes two stages and a task,
        // and the zip comes two stages and a task after the grouping.
        let (stage, task) = match (grouped, failure.late) {
            (true, false) => (1, 1),
            (true, true) => (3, 3),
            (false, false) => (3, 2),
            (false, true) => (5, 4),
        };
        eprintln!("ERROR Executor - Exception in task 0.0 in stage {stage}.0 (TID {task})");
        Thrown::non_user("org.apache.spark.SparkException", failure.message)
    })?;
    if let Some(metrics) = argument(parser, "metrics-file") {
        let headers = vec![
            crate::command_line::expanded(tool, parser),
            format!("Started on: {}", java_display_now()),
        ];
        std::fs::write(&metrics, md::metrics_file(&marked, &header, &headers)).map_err(
            |error| Thrown::non_user(PORT_FAILURE, format!("could not write {metrics}: {error}")),
        )?;
    }
    let kept: Vec<BamRecord> = marked
        .reads
        .into_iter()
        .filter(|record| {
            if remove_all {
                !gatk_engine::read::is_duplicate(record)
            } else if remove_sequencing {
                !matches!(
                    record.tags.get(htsjdk_bam::tag::Tag::new(b"DT")),
                    Some(htsjdk_bam::tag::TagValue::Str(kind)) if kind == "SQ"
                )
            } else {
                true
            }
        })
        .collect();
    header.set_sort_order("coordinate");
    write_spark_reads(parser, &header, &kept, &output, true)?;
    Ok(None)
}

/// `MarkDuplicatesSparkArgumentCollection`, with the optical finder's two settings as the caller's
/// tool has them.
fn mark_arguments(
    parser: &Parser,
    mode: &Parser,
    remove_sequencing: bool,
    read_name_regex: Option<&str>,
    pixel_distance: i32,
) -> Result<gatk_tools::mark_duplicates_spark::MarkArguments, Thrown> {
    use gatk_tools::mark_duplicates_spark as md;
    let mut tagging = match scalar(parser, "duplicate-tagging-policy").as_deref() {
        Some("All") => md::Tagging::All,
        Some("OpticalOnly") => md::Tagging::OpticalOnly,
        _ => md::Tagging::DontTag,
    };
    if remove_sequencing && tagging == md::Tagging::DontTag {
        tagging = md::Tagging::OpticalOnly;
    }
    let regex = md::ReadNameRegex::from_arg(read_name_regex)
        .map_err(|error| Thrown::non_user("java.util.regex.PatternSyntaxException", error))?;
    Ok(md::MarkArguments {
        scoring: match scalar(parser, "duplicate-scoring-strategy").as_deref() {
            Some("TOTAL_MAPPED_REFERENCE_LENGTH") => md::Scoring::TotalMappedReferenceLength,
            Some("FLOW_SUM_OF_BASE_QUALITIES") => md::Scoring::FlowSumOfBaseQualities,
            _ => md::Scoring::SumOfBaseQualities,
        },
        dont_mark_unmapped_mates: flag(parser, "do-not-mark-unmapped-mates"),
        tagging,
        read_name_regex: regex,
        pixel_distance,
        flow: md::FlowArguments {
            quality_sum_strategy: flag(mode, "flow-quality-sum-strategy"),
            end_location_significant: flag(mode, "single-end-reads-end-position-significant"),
            end_uncertainty: number_or(mode, "flow-end-pos-uncertainty", 0),
            use_clipped_locations: flag(mode, "single-end-reads-clipping-is-end"),
            skip_start_homopolymers: number_or(mode, "flow-skip-start-homopolymers", 0),
            q_is_known_end: flag(mode, "flow-q-is-known-end"),
        },
    })
}

/// `ReadsPipelineSpark.runTool` up to the caller: the reads marked, sorted, recalibrated from a
/// table of their own, written if asked, and put through the caller's filters at their defaults.
/// The duplicates stay: only `MarkDuplicatesSpark.runTool` removes any, and `new
/// OpticalDuplicateFinder()` is the default regex and distance whatever the tool was given.
pub(super) fn reads_pipeline_reads(
    parser: &Parser,
    source: &ReadsDataSource,
    header: &SamHeader,
    intervals: &[gatk_engine::interval::SimpleInterval],
    filters: &[gatk_tools::filter_resolution::ResolvedFilter],
) -> Result<(Vec<BamRecord>, SamHeader), Thrown> {
    use gatk_tools::mark_duplicates_spark as md;
    spark_known_sites(parser)?;
    // `BwaSparkEngine`'s header is what follows the alignment.
    let engine = if flag(parser, "align") {
        Some(bwa_engine(parser, header)?)
    } else {
        None
    };
    let header = engine.as_ref().map_or(header, |(aligned, _)| aligned);
    let filter = read_filter(parser, filters, header)?;
    let mut records = match &engine {
        Some((aligned, image)) if !flag(parser, "single-end-alignment") => {
            // After alignment, since filtering does not know about pairs.
            let reads = gatk_tools::read_walker::traverse(source, intervals, &|_| true)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            bwa_align(&reads, true, image, aligned)?
                .into_iter()
                .filter(|read| filter(read))
                .collect()
        }
        Some((aligned, image)) => {
            let reads = gatk_tools::read_walker::traverse(source, intervals, &filter)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            bwa_align(&reads, false, image, aligned)?
        }
        None => gatk_tools::read_walker::traverse(source, intervals, &filter)
            .map_err(|error| Thrown::user(format!("{error:?}")))?,
    };
    drop(filter);
    let settings = mark_arguments(
        parser,
        parser,
        false,
        Some("<optimized capture of last three ':' separated fields as numeric values>"),
        100,
    )?;
    if header.read_groups.is_empty() {
        return Err(Thrown::user(
            "Bad input: Sam file header missing Read Group fields. MarkDuplicatesSpark currently \
             requires reads to be labeled with read group tags, please add read groups tags to \
             your reads",
        ));
    }
    let grouped = header.attributes.get("SO") == Some("queryname")
        || header.attributes.get("GO") == Some("query");
    if !grouped {
        records.sort_by(htsjdk_bam::query_name::compare);
    }
    let marked = md::mark(vec![records], header, &settings).map_err(|failure| {
        eprintln!("ERROR Executor - Exception in task 0.0 in stage 1.0 (TID 1)");
        Thrown::non_user("org.apache.spark.SparkException", failure.message)
    })?;
    let mut sorted = marked.reads;
    sorted.sort_by(htsjdk_bam::coordinate::compare);

    let mut for_table = Vec::new();
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
    let table_filter = read_filter(parser, &for_table, header)?;
    let table_reads: Vec<BamRecord> = sorted
        .iter()
        .filter(|record| table_filter(record))
        .cloned()
        .collect();
    drop(table_filter);
    let table = super::recalibration_table_over(
        parser,
        header,
        super::RecalibrationReads::Chosen(&table_reads),
    )
    .map_err(table_task_failure)?;
    tables_combine(header)?;
    let finished = recalibrated(
        header,
        &table,
        &super::apply_bqsr_arguments(parser),
        &sorted,
    )?;
    if let Some(output_bam) = argument(parser, "output-bam") {
        write_spark_reads(parser, header, &finished, &output_bam, true)?;
    }
    let caller_filter = super::standard_hc_read_filter(header)?;
    let for_caller: Vec<BamRecord> = finished
        .into_iter()
        .filter(|record| caller_filter(record))
        .collect();
    let mut reads_header = header.clone();
    reads_header.set_sort_order("coordinate");
    Ok((for_caller, reads_header))
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
    let indices = (
        flag(parser, "create-output-bam-index"),
        flag(parser, "create-output-bam-splitting-index"),
    );
    write_spark_indexed(parser, header, sort_order, partitions, output, indices)
}

/// [`write_spark_partitions`], told which of the `.bai` and the `.sbi` to write.
fn write_spark_indexed(
    parser: &Parser,
    header: &SamHeader,
    sort_order: Option<&str>,
    partitions: &[Vec<BamRecord>],
    output: &str,
    (write_bai, write_sbi): (bool, bool),
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

    if output.ends_with(".sam") {
        // Disq's SAM writer: the text, whatever the partitions, and no index of any kind.
        let records: Vec<BamRecord> = partitions.iter().flatten().cloned().collect();
        let text = htsjdk_bam::sam_file::write_sam(header, &records)
            .ok_or_else(|| failed(output, "the SAM text could not be built".to_string()))?;
        return std::fs::write(output, text).map_err(|error| failed(output, error.to_string()));
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
                let pointer = stream
                    .file_pointer()
                    .map_err(|error| failed(output, error.to_string()))?;
                offsets.push(shift + pointer);
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

    if write_bai && sort_order == Some("coordinate") {
        // shortcut: the index is built from the whole file; Disq merges one per part with
        // `BAMIndexMerger`, which is the same index on every array row measured so far. Port the
        // merger if a multi-partition row ever disagrees.
        let bai = htsjdk_bam::build_index::build_bam_index(&bam)
            .map_err(|error| failed(&format!("{output}.bai"), format!("{error:?}")))?;
        std::fs::write(format!("{output}.bai"), bai)
            .map_err(|error| failed(&format!("{output}.bai"), error.to_string()))?;
    }
    if write_sbi {
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
    // `getRecommendedNumReducers`: the argument, else one per ten megabytes of input.
    let mut reducers = number_or(parser, "num-reducers", 0).max(0) as usize;
    if reducers == 0 {
        let mut inputs = arguments(parser, "input");
        inputs.sort();
        inputs.dedup();
        let size: u64 = inputs
            .iter()
            .filter_map(|input| std::fs::metadata(input).ok())
            .map(|data| data.len())
            .sum();
        reducers = 1 + (size / (10 * 1024 * 1024)) as usize;
    }
    range_partitions(records, reducers, &[records.len()], compare)
}

/// The sort itself: `reducers` partitions over input partitions of the sizes given.
fn range_partitions(
    records: &[BamRecord],
    reducers: usize,
    input_partitions: &[usize],
    compare: fn(&BamRecord, &BamRecord) -> std::cmp::Ordering,
) -> Result<Vec<Vec<BamRecord>>, Thrown> {
    let mut sorted = records.to_vec();
    sorted.sort_by(compare);
    let bounds = range_bounds(&sorted, reducers, input_partitions, compare)?;
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

/// `RangePartitioner`'s bounds, the sketch taking every key: weight one each, then
/// `determineBounds(candidates, min(partitions, candidates))`.
///
/// shortcut: input partitions each holding at most `sampleSizePerPartition` keys, where the sketch
/// is the whole partition; more is a reservoir sample the port does not draw, and is refused by name.
fn range_bounds(
    sorted: &[BamRecord],
    partitions: usize,
    input_partitions: &[usize],
    compare: fn(&BamRecord, &BamRecord) -> std::cmp::Ordering,
) -> Result<Vec<BamRecord>, Thrown> {
    if partitions <= 1 || sorted.is_empty() {
        return Ok(Vec::new());
    }
    let sample_size = (20.0 * partitions as f64).min(1e6);
    let per_partition = (3.0 * sample_size / input_partitions.len() as f64).ceil() as usize;
    if input_partitions.iter().any(|size| *size > per_partition) {
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
    let table = super::recalibration_table(parser, &source, &header, &intervals, &filters)
        .map_err(table_task_failure)?;
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
    let table = super::recalibration_table(parser, &source, &header, &intervals, &for_table)
        .map_err(table_task_failure)?;
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

/// `updatePlatformForRead` throws inside the executor, so its `MalformedRead` reaches the driver as
/// the job's failure.
fn table_task_failure(thrown: Thrown) -> Thrown {
    match &thrown.message {
        Some(message) if message.contains("reads with no platform information") => {
            eprintln!("ERROR Executor - Exception in task 0.0 in stage 0.0 (TID 0)");
            Thrown::non_user("org.apache.spark.SparkException", message.clone())
        }
        _ => thrown,
    }
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
