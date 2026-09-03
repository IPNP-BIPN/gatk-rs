//! The tools this port can run, and the file plumbing under them.
//!
//! A port's own API takes and returns bytes: `index_feature_file::build` is handed a whole file as
//! a `&str` and answers with the index as a `Vec<u8>`, which is what a conformance suite comparing
//! whole outputs against a golden wants. A command line wants neither: it names paths, and
//! something has to read one and write the other.
//!
//! That something is here rather than in the port, deliberately. The suite's claim is about the
//! bytes, and a `std::fs::read` in the middle of the ported function would put the filesystem
//! inside the thing being compared.
//!
//! Ported from `org.broadinstitute.hellbender.tools.IndexFeatureFile`,
//! `org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation` and
//! `org.broadinstitute.hellbender.tools.CountReads`.

use gatk_barclay::{Parser, Value};
use gatk_engine::interval::MergingRule;
use gatk_engine::interval_arguments::{SetRule, TraversalParameters};
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::index_feature_file::{self, Refusal, Source};
use gatk_tools::main_entry::{Failure, Thrown, PORT_FAILURE, PORT_LIMITATION};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// What a runner answers: what the tool returned, or what it threw.
///
/// A [`Thrown`] rather than a message, because the two handlers `mainEntry` calls write different
/// things: `handleUserException` decorates the message with a banner, and `handleNonUserException`
/// prints the exception's own CLASS in front of it. A runner that only reported a message left the
/// dispatcher with one banner to print for both (`main-non-user`).
pub type Outcome = Result<Option<String>, Thrown>;

/// The value of one named argument, as the parser left it.
///
/// A path argument holds a `Tagged` value, whose tag is `None` when nobody wrote one; a plain
/// string argument holds a `Str`. Both are read the same way here because the tool wants the path
/// and not the tag.
pub fn argument(parser: &Parser, long_name: &str) -> Option<String> {
    parser
        .definitions()
        .iter()
        .find(|definition| definition.long_name() == long_name)
        .and_then(|definition| match &definition.value {
            Value::Tagged { value, .. } => Some(value.clone()),
            Value::Str(text) => Some(text.clone()),
            _ => None,
        })
}

/// `IndexFeatureFile.doWork`, with the two paths read and written.
///
/// The reference returns the index's path, which `handleResult` then prints, so this returns it
/// too. Every refusal is the port's own [`Refusal`], and each of them is a `UserException`, which
/// is status two.
pub fn index_feature_file(parser: &Parser) -> Outcome {
    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let output = argument(parser, "output");
    // Almost every refusal here is a `UserException`, which is status two; the one that is not
    // names its own class and takes the other handler.
    let refused = |refusal: Refusal| {
        if refusal.is_user() {
            Thrown::user(refusal.message())
        } else {
            Thrown::non_user(refusal.java_class(), refusal.message())
        }
    };
    // The reference reads the file to find a codec for it, so a file that is not there is refused
    // before anything else is asked of it.
    let bytes = std::fs::read(&input).map_err(|_| {
        refused(Refusal::CouldNotReadInputFile {
            path: input.clone(),
        })
    })?;
    let text = decode(&bytes, &input).map_err(refused)?;
    let name = output
        .clone()
        .unwrap_or_else(|| index_feature_file::default_output(&input));
    index_feature_file::check_output(&input, &name, &input).map_err(refused)?;
    // The header records the file's own identity, and its timestamp with it: a caller that wants
    // the reference's bytes has to supply the real one, which is what this does.
    let mut source = Source::new(&input);
    source.timestamp = modified_millis(&input);
    // A block compressed input gets a tabix index, whose positions are the pointers a BGZF reader
    // reports rather than offsets into the text, so it is handed the FILE and not what is in it.
    let index = match index_feature_file::index_kind(&input) {
        index_feature_file::IndexKind::Tabix => {
            // A `.tbi` is a BGZF file, and GATK replaces htsjdk's static deflater factory: the
            // reference's bytes are Intel's GKL unless `--use-jdk-deflater` says otherwise, which
            // is the argument the tool declares for exactly this.
            let deflater = if flag(parser, "use-jdk-deflater") {
                htsjdk_bgzf::Deflater::Jdk
            } else {
                htsjdk_bgzf::Deflater::Gkl
            };
            // And the LEVEL is GATKConfig's, which is two rather than htsjdk's five: `Main`
            // installs it as a system property before any tool runs, so every block-compressed
            // file a real invocation writes uses it (#1032).
            let level = gatk_tools::gatk_config::compression_level(
                std::env::var(gatk_tools::gatk_config::COMPRESSION_LEVEL)
                    .ok()
                    .as_deref(),
            );
            index_feature_file::build_tabix(&bytes, &source, &input, deflater, level)
                .map_err(refused)?
        }
        _ => index_feature_file::build(&text, &source, &input).map_err(refused)?,
    };
    std::fs::write(&name, index).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {name}: {error}"))
    })?;
    Ok(Some(name))
}

/// `PrintBGZFBlockInformation.doWork`, with the file read and the report written.
///
/// The tool prints to standard output when it is given no `--output`, which is the one place a
/// runner's answer is not a file: what it returns is the report itself, and `handleResult` prints
/// what a tool returns.
pub fn print_bgzf_block_information(parser: &Parser) -> Outcome {
    let input = argument(parser, "bgzf-file").ok_or_else(|| {
        Thrown::command_line("Argument bgzf-file was missing: Argument 'bgzf-file' is required")
    })?;
    let bytes = std::fs::read(&input).map_err(|_| {
        Thrown::user(
            gatk_tools::print_bgzf_block_information::Refusal::DoesNotExist {
                path: input.clone(),
            }
            .message(),
        )
    })?;
    let name = std::path::Path::new(&input)
        .file_name()
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_else(|| input.clone());
    let (report, refusal) = gatk_tools::print_bgzf_block_information::report(&bytes, &name, &input);
    match argument(parser, "output") {
        Some(path) => std::fs::write(&path, report).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {path}: {error}"))
        })?,
        // With no output the report goes to standard output, which the dispatcher prints as the
        // tool's own return value.
        None => {
            if let Some(refusal) = refusal {
                return Err(Thrown::user(refusal.message()));
            }
            return Ok(Some(report));
        }
    }
    match refusal {
        Some(refusal) => Err(Thrown::user(refusal.message())),
        None => Ok(None),
    }
}

/// The file's text, decompressed when the name says it is block compressed.
fn decode(bytes: &[u8], path: &str) -> Result<String, Refusal> {
    let raw = if path.ends_with(".gz") {
        htsjdk_bgzf::read::decompress_all(bytes).map_err(|_| Refusal::NoSuitableCodecs {
            path: path.to_string(),
        })?
    } else {
        bytes.to_vec()
    };
    String::from_utf8(raw).map_err(|_| Refusal::NoSuitableCodecs {
        path: path.to_string(),
    })
}

/// `File.lastModified()` in milliseconds, which the index header carries.
fn modified_millis(path: &str) -> i64 {
    std::fs::metadata(path)
        .and_then(|data| data.modified())
        .ok()
        .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|since| since.as_millis() as i64)
        .unwrap_or(0)
}

/// The values of a collection argument, as the parser left them.
pub fn arguments(parser: &Parser, long_name: &str) -> Vec<String> {
    parser
        .definitions()
        .iter()
        .find(|definition| definition.long_name() == long_name)
        .map(|definition| match &definition.value {
            Value::List(values) => values
                .iter()
                .map(|value| match value {
                    Value::Tagged { value, .. } => value.clone(),
                    other => other.to_java_string(),
                })
                .collect(),
            _ => Vec::new(),
        })
        .unwrap_or_default()
}

/// A scalar argument's value as text, whatever class it holds.
///
/// [`argument`] answers for the two classes a path arrives in; a number arrives as an `Int`, and
/// reading it through that function would silently answer `None` and leave a default in place. The
/// `count-reads-plumbing` golden is what caught it: `--minimum-mapping-quality 70` counted eight
/// reads instead of none.
pub fn scalar(parser: &Parser, long_name: &str) -> Option<String> {
    parser
        .definitions()
        .iter()
        .find(|definition| definition.long_name() == long_name)
        .and_then(|definition| match &definition.value {
            Value::Null => None,
            Value::Tagged { value, .. } => Some(value.clone()),
            other => Some(other.to_java_string()),
        })
}

/// Whether a flag argument was set, which is a boolean whose default is false.
pub fn flag(parser: &Parser, long_name: &str) -> bool {
    matches!(
        parser
            .definitions()
            .iter()
            .find(|definition| definition.long_name() == long_name)
            .map(|definition| &definition.value),
        Some(Value::Bool(true))
    )
}

/// The conjunction a command line's read filters make, over records the walker hands it.
pub type Filter<'a> = Box<dyn Fn(&BamRecord) -> bool + 'a>;

/// The read filter a command line asks for, which is a conjunction and not a choice.
///
/// `--read-filter` ADDS to the tool's defaults; `--disable-tool-default-read-filters` is what
/// replaces them. Both are in the `count-reads-plumbing` golden, one case each, and the filter
/// order follows the reference's: the defaults first, then the named ones in the order they were
/// named.
///
/// A filter this port does not carry is refused rather than ignored, because ignoring one would
/// count reads the reference filtered out and answer with a number that looks right.
/// `validateAndResolvePlugins`, which the reference runs while it PARSES the command line.
///
/// That is the whole reason this is a function of its own: the refusals happen before a file is
/// opened, an interval is resolved or a dictionary is compared, so a port that resolved its
/// filters at the point of use answered a later question first. A covering-array row naming a
/// filter both enabled and disabled read `Dictionary cannot have size zero` in the port and
/// `are both enabled and disabled` in the reference (#69).
///
/// A tool that applies no read filter still validates them: the descriptor is the command line's,
/// not the traversal's.
fn resolve_read_filters(
    parser: &Parser,
    tool: &str,
) -> Result<Vec<gatk_tools::filter_resolution::ResolvedFilter>, Thrown> {
    // The descriptor owns FOUR arguments, and reading two of them was a port that ignored
    // `--disable-read-filter` and `--inverted-read-filter` entirely. `filter-resolution` measured
    // what all four decide, including the order and the six refusals.
    gatk_tools::filter_resolution::resolve(
        gatk_tools::plugin_ownership::default_filters(tool).unwrap_or(&[]),
        &gatk_tools::plugin_ownership::CATALOGUE,
        &arguments(parser, "read-filter"),
        &arguments(parser, "disable-read-filter"),
        &arguments(parser, "inverted-read-filter"),
        flag(parser, "disable-tool-default-read-filters"),
    )
    .map_err(|error| Thrown {
        // Every one of them is a `CommandLineException`, which is status ONE.
        failure: Failure::CommandLine,
        exception: error.java_class(),
        message: Some(error.message()),
    })
}

/// The predicate the resolved list becomes, once a header exists to read a filter against.
fn read_filter<'a>(
    parser: &'a Parser,
    resolved: &[gatk_tools::filter_resolution::ResolvedFilter],
    header: &'a SamHeader,
) -> Result<Filter<'a>, Thrown> {
    let mut plain: Vec<(gatk_readfilter::ReadFilter, bool)> = Vec::new();
    let mut wellformed: Option<bool> = None;
    let mut parameterized: Vec<(gatk_readfilter::Parameterized, bool)> = Vec::new();
    for filter in resolved {
        let name = filter.name.as_str();
        if name == "WellformedReadFilter" {
            wellformed = Some(filter.negated);
        } else if let Some(plain_filter) = gatk_readfilter::by_name(name) {
            plain.push((plain_filter, filter.negated));
        } else if name == "MappingQualityReadFilter" {
            let minimum = scalar(parser, "minimum-mapping-quality")
                .and_then(|text| text.parse::<i32>().ok())
                .unwrap_or(10);
            let maximum =
                scalar(parser, "maximum-mapping-quality").and_then(|text| text.parse::<i32>().ok());
            parameterized.push((
                gatk_readfilter::Parameterized::MappingQuality {
                    min: minimum,
                    max: maximum,
                },
                filter.negated,
            ));
        } else {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "{name} is a GATK read filter that this port does not carry yet. This message is the port's own and not GATK's."
                ),
            ));
        }
    }
    Ok(Box::new(move |read: &BamRecord| {
        // `ReadFilterNegate` wraps the filter rather than replacing it, so a negated one answers
        // the opposite of what the filter itself answers, on the same read.
        if let Some(negated) = wellformed {
            if gatk_readfilter::with_header::wellformed(read, header) == negated {
                return false;
            }
        }
        if !plain
            .iter()
            .all(|(filter, negated)| filter(read) != *negated)
        {
            return false;
        }
        parameterized
            .iter()
            .all(|(filter, negated)| filter.decide(read).unwrap_or(false) != *negated)
    }))
}

/// `--read-index`: the index to use for each reads input, in place of the one the name implies.
///
/// Two rules, and both fire while the reads are OPENED, before the intervals are resolved or any
/// dictionary is compared.
///
/// The count has to match: `ReadsPathDataSource` refuses a command line with a different number of
/// indices and inputs, and the message counts both.
///
/// And an index is refused only by the LAST branch of `SamReaderFactory`, which is plain text. The
/// order there is BAM, then block compressed, then gzip, then CRAM, then SRA, then text -- and the
/// two compressed branches build a `SAMTextReader` over the decompressed stream without looking at
/// the index at all. So a `.vcf.gz` handed to `--input` with an index beside it is accepted and the
/// index ignored, where the same file uncompressed is a `RuntimeException` and exits three.
///
/// `Ok(None)` where the argument was not given, which is what tells the caller to look for the
/// index the file's own name implies (`SamFiles.findIndex`).
fn read_index(
    parser: &Parser,
    inputs: usize,
    is_binary: bool,
    is_compressed: bool,
) -> Result<Option<std::path::PathBuf>, Thrown> {
    let indices = arguments(parser, "read-index");
    if indices.is_empty() {
        return Ok(None);
    }
    if indices.len() != inputs {
        return Err(Thrown::user(format!(
            "Must have the same number of BAM/CRAM/SAM paths and indices. Saw {inputs} \
             BAM/CRAM/SAMs but {} indices",
            indices.len()
        )));
    }
    if !is_binary && !is_compressed {
        return Err(Thrown::non_user(
            "java.lang.RuntimeException",
            "Cannot use index file with textual SAM file",
        ));
    }
    // A compressed text stream takes a reader that never asks for one, so the index is accepted
    // and dropped rather than used.
    Ok(is_binary.then(|| std::path::PathBuf::from(&indices[0])))
}

/// `--reference`: the reference dictionary, read off the `.dict` beside the FASTA.
///
/// GATK requires a reference to carry both an index and a dictionary, and the dictionary's name
/// REPLACES the FASTA's extension rather than appending to it, which is the opposite of a feature
/// file's index and the same rule `SamFiles.findIndex` follows for a BAM.
///
/// What the port needs from a reference is that dictionary: it takes part in the validation, and
/// it is the best available one when no `--sequence-dictionary` was given.
fn reference_dictionary(parser: &Parser) -> Result<Option<SamHeader>, Thrown> {
    let Some(path) = argument(parser, "reference") else {
        return Ok(None);
    };
    let dictionary = std::path::Path::new(&path).with_extension("dict");
    let text = std::fs::read_to_string(&dictionary)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&path, false)))?;
    Ok(Some(htsjdk_bam::reader::parse_header_text(&text)))
}

/// `--sequence-dictionary`: the MASTER dictionary, read off a `.dict` file.
///
/// A `.dict` is a SAM header with `@SQ` lines and no records, so it parses as one. `None` where
/// the argument was not given, which is what `masterSequenceDictionary == null` means.
fn master_dictionary(parser: &Parser) -> Result<Option<SamHeader>, Thrown> {
    let Some(path) = argument(parser, "sequence-dictionary") else {
        return Ok(None);
    };
    let text = std::fs::read_to_string(&path)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&path, false)))?;
    Ok(Some(htsjdk_bam::reader::parse_header_text(&text)))
}

/// `validateDictionaries` against the master, which runs before the pairs that do not involve it.
///
/// `requireSuperset` is `hasCramInput()`, which is false for every input this port opens, and the
/// contig ordering is not checked. The name in the message is the reference's own.
fn validate_against_master(
    master: &SamHeader,
    other_name: &str,
    other: &[htsjdk_bam::header::SequenceRecord],
) -> Result<(), Thrown> {
    gatk_tools::sequence_dictionary::validate(
        "master sequence dictionary",
        &master.sequences,
        other_name,
        other,
        false,
        false,
    )
    .map_err(|refusal| Thrown {
        failure: Failure::User,
        exception: refusal.java_class(),
        message: Some(refusal.message()),
    })
}

/// The five arguments of `IntervalArgumentCollection`, resolved against a dictionary.
///
/// `--intervals` was the only one the runners read: `--exclude-intervals`, `--interval-set-rule`,
/// `--interval-merging-rule`, `--interval-padding` and `--interval-exclusion-padding` changed
/// nothing here and change the answer in the reference. They are measured in `interval-arguments`
/// and ported in [`gatk_engine::interval_arguments`]; this is the layer that reads them off a
/// command line.
///
/// `None` where the collection was not specified at all, which is what a walker traverses
/// everything for.
fn interval_arguments(
    parser: &Parser,
    header: &SamHeader,
) -> Result<Option<TraversalParameters>, Thrown> {
    let include = arguments(parser, "intervals");
    let exclude = arguments(parser, "exclude-intervals");
    if include.is_empty() && exclude.is_empty() {
        return Ok(None);
    }

    // The two enums arrive as their constant names, which is what the declaration's domain holds.
    let set_rule = match scalar(parser, "interval-set-rule").as_deref() {
        Some("INTERSECTION") => SetRule::Intersection,
        _ => SetRule::Union,
    };
    let merging_rule = match scalar(parser, "interval-merging-rule").as_deref() {
        Some("OVERLAPPING_ONLY") => MergingRule::OverlappingOnly,
        _ => MergingRule::All,
    };
    let padding = number(parser, "interval-padding");
    let exclusion_padding = number(parser, "interval-exclusion-padding");

    let parameters = gatk_engine::interval_arguments::traversal_parameters(
        &include,
        &exclude,
        header,
        set_rule,
        merging_rule,
        padding,
        exclusion_padding,
    )
    .map_err(|error| Thrown {
        failure: match error {
            // A bad argument value is a `CommandLineException`, which is status ONE, and the two
            // interval refusals are exactly that where the parse failures are status two.
            gatk_engine::interval_arguments::IntervalArgumentError::EmptyIntersection {
                ..
            }
            | gatk_engine::interval_arguments::IntervalArgumentError::ExcludedEverything {
                ..
            } => Failure::CommandLine,
            _ => Failure::User,
        },
        exception: error.java_class(),
        message: Some(error.message()),
    })?;

    if parameters.traverse_unmapped {
        // `-L unmapped` asks the traversal for the records with no position, which neither of
        // these tools' ported traversals can produce. Refusing is the port's own answer and says
        // so; counting the mapped ones and calling it the total would not.
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "-L unmapped asks for a traversal of unmapped records that this port does not carry \
             yet. This message is the port's own and not GATK's.",
        ));
    }
    Ok(Some(parameters))
}

/// An integer argument, or zero where it was not given.
pub fn number(parser: &Parser, long_name: &str) -> i32 {
    scalar(parser, long_name)
        .and_then(|text| text.parse().ok())
        .unwrap_or(0)
}

/// `CountReads.doWork`, with the input read and the output written.
///
/// Three things the `count-reads-plumbing` golden pins and this reproduces: the tool RETURNS the
/// count, so `handleResult` prints a number; `-O` receives that number and nothing else, with no
/// trailing newline, because the reference writes it with `print`; and `-O` does not suppress the
/// return, so the file is written AND the value comes back.
pub fn count_reads(parser: &Parser) -> Outcome {
    // The plugin descriptor is validated while the command line is PARSED, so its refusals come
    // before the input is even opened.
    let resolved_filters = resolve_read_filters(parser, "CountReads")?;
    // `--input` is a COLLECTION on a read walker, not a scalar: the reference takes more than one
    // BAM and merges their headers. This port reads one, which is what every case of the golden
    // hands it, and refuses the rest rather than silently counting the first.
    let inputs = arguments(parser, "input");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --input is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let path = std::path::Path::new(&input);
    let intervals_given = !arguments(parser, "intervals").is_empty();

    // `GATKTool.onStartup` fixes the order, and the order is most of what a covering-array row
    // over this tool measures:
    //
    //   loadMasterSequenceDictionary, initializeReads, initializeIntervals,
    //   validateSequenceDictionaries, then the traversal -- which sets its bounds before it reads
    //   a record.
    //
    // So a file that is not a BAM is refused LAST, by the record parse, and everything the
    // arguments decide is refused before it. A port that asked the reader first answered a later
    // question first (#69).
    let master = master_dictionary(parser)?;
    // `initializeReference` comes before `initializeReads`, and its dictionary outranks the reads'
    // in `getBestAvailableSequenceDictionary`.
    let reference = reference_dictionary(parser)?;

    // `initializeReads`: the file itself, which is a refusal only when it cannot be read at all.
    let bytes = std::fs::read(path).ok();
    let compressed = bytes
        .as_deref()
        .map(gatk_tools::read_walker_refusal::is_block_compressed)
        .unwrap_or(false);
    let decompressed = match (&bytes, compressed) {
        (None, _) => None,
        (Some(bytes), false) => Some(bytes.clone()),
        (Some(bytes), true) => htsjdk_bgzf::read::decompress_all(bytes).ok(),
    };
    let reader_refusal = gatk_tools::read_walker_refusal::refusal(
        &input,
        path.exists(),
        path.is_dir(),
        decompressed.as_deref(),
        compressed,
        // The empty-dictionary refusal belongs to the INTERVAL step, and the dictionary it asks
        // for is the best available one, so a `--sequence-dictionary` supplies it.
        intervals_given && master.is_none(),
    );
    if let Some(refusal) = &reader_refusal {
        if !matches!(
            refusal,
            gatk_tools::read_walker_refusal::Refusal::NotSamText { .. }
        ) {
            // The refusal names the exception the reference throws, and the non-user handler
            // PRINTS that class: a walker refusing an interval over a stream with no dictionary
            // answers `java.lang.IllegalArgumentException: ...` and not a banner (#1020).
            return Err(if refusal.is_user() {
                Thrown::user(refusal.message())
            } else {
                Thrown::non_user(refusal.exception(), refusal.message())
            });
        }
    }
    // A stream that is not SAM text has a header of no sequences, which is what the checks below
    // see; the refusal itself waits for the traversal.
    let deferred_parse_refusal = reader_refusal;

    // The index is the one htsjdk's own search finds, not the one a single `with_extension` call
    // guesses. `reads.bam.bai` was the only name asked for here, and htsjdk writes `reads.bai` at
    // least as often and looks for it FIRST, so an interval query over a file indexed the other
    // way found no index and answered zero rather than refusing (#1020).
    // `--read-index` names the index outright, and its two refusals fire while the reads are
    // opened. Without it the index is the one htsjdk's own search finds.
    let is_binary = decompressed
        .as_deref()
        .is_some_and(|bytes| bytes.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC));
    let named_index = read_index(parser, 1, is_binary, compressed)?;
    // Only a BAM has an index at all: the compressed and plain text branches build a reader that
    // has none, so `indicesAvailable` is false there whatever the command line named.
    let index = if is_binary {
        named_index.or_else(|| htsjdk_bam::sam_files::find_index(path))
    } else {
        None
    };
    let source = if deferred_parse_refusal.is_some() {
        None
    } else {
        Some(match &index {
            // An index that does not parse is not refused here: htsjdk opens it lazily, so the
            // checks after this one run first and one of them is usually what refuses.
            Some(index) => match ReadsDataSource::open(path, index) {
                Ok(source) => source,
                Err(_) => ReadsDataSource::open_unindexed(path)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?,
            },
            None => ReadsDataSource::open_unindexed(path)
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
        })
    };
    let header = source
        .as_ref()
        .map(|source| source.header().clone())
        .unwrap_or_default();

    // `initializeIntervals`, against the best available dictionary: master, then reference, then
    // reads, which is `getBestAvailableSequenceDictionary`'s own order.
    let best = master
        .clone()
        .or_else(|| reference.clone())
        .unwrap_or_else(|| header.clone());
    let intervals = interval_arguments(parser, &best)?
        .map(|parameters| parameters.intervals)
        .unwrap_or_default();

    // `validateSequenceDictionaries`, which the argument turns off wholesale. The master block
    // runs first and checks the reference before the reads; then the reference is checked against
    // the reads on its own.
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            validate_against_master(master, "reads", &header.sequences)?;
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
        }
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

    // `setTraversalBounds`, which the traversal calls before it reads anything.
    if intervals_given && index.is_none() {
        return Err(Thrown::user(
            "Traversal by intervals was requested but some input files are not indexed.",
        ));
    }

    // And only now the record parse.
    if let Some(refusal) = deferred_parse_refusal {
        return Err(Thrown::non_user(refusal.exception(), refusal.message()));
    }
    let source = source.expect("a source, since the parse refusal was not taken");

    let filter = read_filter(parser, &resolved_filters, &header)?;
    let count = gatk_tools::count_reads::count_reads(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;

    if let Some(output) = argument(parser, "output") {
        // `print`, not `println`: the file is the number's digits and nothing else.
        std::fs::write(&output, gatk_tools::count_reads::output(count))
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    }
    // The tool returns the count itself, which is what `handleResult` prints.
    Ok(Some(count.to_string()))
}

/// The sequence dictionary a VCF's own header declares, which is what `-L` resolves against.
///
/// `##contig=<ID=chr1,length=100000>`, in the order the header writes them. Anything else on the
/// line is ignored: `assembly` and `URL` are carried by `SAMSequenceRecord` and no interval query
/// reads them.
fn vcf_dictionary(text: &str) -> SamHeader {
    let mut header = SamHeader::default();
    for line in text.lines() {
        let Some(body) = line.strip_prefix("##contig=<") else {
            if line.starts_with("#CHROM") {
                // The header ends here; a `##contig` after it is not a header line.
                break;
            }
            continue;
        };
        let body = body.trim_end_matches('>');
        let mut name = None;
        let mut length = None;
        for field in body.split(',') {
            match field.split_once('=') {
                Some(("ID", value)) => name = Some(value.to_string()),
                Some(("length", value)) => length = value.parse::<i32>().ok(),
                _ => {}
            }
        }
        if let (Some(name), Some(length)) = (name, length) {
            header
                .sequences
                .push(htsjdk_bam::header::SequenceRecord::new(&name, length));
        }
    }
    header
}

/// Whether the input supports random access, which is what `-L` needs before any record is read.
///
/// `FeatureDataSource` asks the codec for an index beside the file: `.idx` for a plain feature
/// file and `.tbi` for a block compressed one, both APPENDED to the whole name rather than
/// replacing anything. That is `Tribble.indexPath` and `Tribble.tabixIndexPath`, and it is not
/// `SamFiles.findIndex`'s rule: a feature file's index is never named by replacing its extension.
fn has_feature_index(path: &str) -> bool {
    std::path::Path::new(&index_feature_file::default_output(path)).is_file()
}

/// Everything `GATKTool.onStartup` does for a VARIANT walker, up to and including the intervals.
///
/// It is shared rather than copied because the ORDER is the part that is measured: a port that
/// asks the reader before it validates a dictionary answers a later question first, and a
/// covering-array row shows that as a wrong message rather than a wrong number. The sequence is
/// the reference's own -- the read filters while the command line is parsed, then the driving
/// variants, then the reads, then the master dictionary, then the reference, then the pairwise
/// validation, and only then `-L`.
struct VariantStart {
    /// The `--variant` path as it was given, which several refusals quote.
    input: String,
    /// The driving variants as text, decompressed if the file was block compressed.
    text: String,
    codec: gatk_tools::feature_codec::Codec,
    intervals: Option<Vec<gatk_engine::interval::SimpleInterval>>,
}

fn variant_walker_startup(parser: &Parser, tool: &str) -> Result<VariantStart, Thrown> {
    let _ = resolve_read_filters(parser, tool)?;
    // `--variant` is a SCALAR on this tool, where a read walker's `--input` is a collection: the
    // declaration says `collection: false`, and reading it as a list finds nothing at all.
    let input = argument(parser, "variant").ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;

    let codec = gatk_tools::feature_codec::codec_for(&input).ok_or_else(|| {
        Thrown::user(
            index_feature_file::Refusal::NoSuitableCodecs {
                path: input.clone(),
            }
            .message(),
        )
    })?;
    let bytes = std::fs::read(&input).map_err(|_| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: input.clone(),
            }
            .message(),
        )
    })?;
    // A block compressed feature file is read through its own decompression, and its index is a
    // `.tbi` rather than a `.idx`; neither changes what the traversal counts.
    let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .ok()
            .and_then(|bytes| String::from_utf8(bytes).ok())
            .ok_or_else(|| {
                Thrown::non_user(
                    gatk_tools::read_walker_refusal::SAM_FORMAT,
                    format!("{input} is not a block compressed file"),
                )
            })?
    } else {
        String::from_utf8_lossy(&bytes).into_owned()
    };

    let header = vcf_dictionary(&text);

    // `--read-index` is counted against the READS inputs, and both its refusals fire while they
    // are opened -- which a variant walker does whenever a command line names any, and before
    // anything decides whether the dictionaries are compared at all.
    let reads_inputs = arguments(parser, "input").len();
    let reads_dictionaries: Vec<Vec<htsjdk_bam::header::SequenceRecord>> =
        arguments(parser, "input")
            .iter()
            .map(|reads| reads_dictionary(parser, reads, reads_inputs))
            .collect::<Result<_, _>>()?;

    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    // `validateSequenceDictionaries` is ONE method and the argument turns all of it off, the
    // master block included: a guard around part of it refuses command lines the reference runs.
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            // The master block runs before the reference/reads/features loop, and inside it the
            // READS come before the reference and the reference before the features.
            for reads in &reads_dictionaries {
                validate_against_master(master, "reads", reads)?;
            }
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
            validate_against_master(master, "features", &header.sequences)?;
        }
        if let Some(reference) = &reference {
            // The reference against the reads, then the reference against the features.
            for reads in &reads_dictionaries {
                gatk_tools::sequence_dictionary::validate(
                    "reference",
                    &reference.sequences,
                    "reads",
                    reads,
                    false,
                    false,
                )
                .map_err(|refusal| Thrown {
                    failure: Failure::User,
                    exception: refusal.java_class(),
                    message: Some(refusal.message()),
                })?;
            }
            gatk_tools::sequence_dictionary::validate(
                "reference",
                &reference.sequences,
                "features",
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
    // `GATKTool.onStartup` validates the dictionaries against each other BEFORE the traversal, and
    // a variant walker still opens the reads when a command line names any: the pair goes through
    // `validateDictionaries("reads", readDict, "features", featureDict)`, whose four-argument
    // overload requires no superset and does not check the contig ordering. A corpus whose BAM and
    // VCF share no contig is therefore refused whatever the intervals say (#1038).
    if !flag(parser, "disable-sequence-dictionary-validation") {
        for reads in &reads_dictionaries {
            gatk_tools::sequence_dictionary::validate(
                "reads",
                reads,
                "features",
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

    // The master dictionary is validated against the features, and then again as `best available`.
    // What `-L` resolves against is NOT the master here: a variant walker prefers the DRIVING
    // VARIANTS' dictionary unless that one was synthesized from an index, and a VCF carrying
    // `##contig` lines gives a real one.
    let best = if header.sequences.is_empty() {
        master.clone().unwrap_or_else(|| header.clone())
    } else {
        header.clone()
    };
    let intervals = interval_arguments(parser, &best)?.map(|parameters| parameters.intervals);

    Ok(VariantStart {
        input,
        text,
        codec,
        intervals,
    })
}
/// `CountVariants.doWork`, with the input read, the intervals resolved and the count written.
///
/// Four things the `count-variants` golden pins and this reproduces: the count reaches no stream
/// without `-O`, whatever the class documentation says; a record is selected by its whole SPAN,
/// `END` or the length of `REF`, so an interval reaches a record whose position it does not hold;
/// `-L` against an input with no index is refused BEFORE any record is read; and the refusal for
/// an unwritable `-O` carries the path and nothing else.
pub fn count_variants(parser: &Parser) -> Outcome {
    // A variant walker applies no read filter and still VALIDATES the ones a command line names:
    // the descriptor belongs to the command line rather than to the traversal, so `--read-filter`
    // and its three companions are refused here exactly as they are on a read walker.
    let VariantStart {
        input,
        text,
        codec,
        intervals,
    } = variant_walker_startup(parser, "CountVariants")?;

    let features: Vec<Locus> = gatk_tools::feature_codec::features(&text, codec)
        .into_iter()
        .map(|(feature, _)| Locus {
            contig: feature.contig,
            start: feature.start,
            stop: feature.end,
        })
        .collect();

    let count = gatk_tools::count_variants::count(
        &features,
        intervals.as_deref(),
        has_feature_index(&input),
        &input,
    )
    .map_err(|error| Thrown {
        failure: Failure::User,
        exception: error.class(),
        message: Some(error.message()),
    })?;

    let output = argument(parser, "output");
    gatk_tools::count_variants::write_output(output.as_deref().map(std::path::Path::new), count)
        .map_err(|error| Thrown {
            failure: Failure::User,
            exception: error.class(),
            message: Some(error.message()),
        })?;
    // The tool returns the count, which `handleResult` prints.
    Ok(Some(count.to_string()))
}

/// The sequence dictionary a reads input carries, which for a file that is not a BAM is EMPTY.
///
/// `ReadsPathDataSource` opens whatever it is given and asks for the header before it reads a
/// record, so a VCF or a BED handed to `--input` is a stream with no `@SQ` line rather than a
/// refusal: the dictionary comes back empty and the comparison against the features' dictionary
/// finds no common contigs. The refusal a read WALKER makes for the same file is a later one, and
/// it is the record parse rather than the header that makes it (`read-walker-refusals`).
fn reads_dictionary(
    parser: &Parser,
    path: &str,
    inputs: usize,
) -> Result<Vec<htsjdk_bam::header::SequenceRecord>, Thrown> {
    let bytes = std::fs::read(path)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(path, false)))?;
    let compressed = gatk_tools::read_walker_refusal::is_block_compressed(&bytes);
    let decompressed = if compressed {
        htsjdk_bgzf::read::decompress_all(&bytes).unwrap_or_default()
    } else {
        bytes
    };
    let is_binary = decompressed.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC);
    // `--read-index` is refused while the reads are OPENED, so it refuses here too: a walker that
    // only wants the dictionary still opens them.
    let _ = read_index(parser, inputs, is_binary, compressed)?;
    if !is_binary {
        return Ok(Vec::new());
    }
    let source = ReadsDataSource::open_unindexed(std::path::Path::new(path))
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    Ok(source.header().sequences.clone())
}

/// One decoded locus, which is all the traversal looks at.
struct Locus {
    contig: String,
    start: i32,
    stop: i32,
}

impl gatk_engine::variant_source::Located for Locus {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn stop(&self) -> i32 {
        self.stop
    }
}

/// `UserException$BadInput`, whose constructor puts `Bad input: ` in front of the message.
fn bad_input(message: String) -> Thrown {
    Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(format!("Bad input: {message}")),
    }
}

/// Whether a BAM's header says `SO:coordinate`, which is what a `.bai` needs.
fn is_coordinate_sorted(bam: &[u8]) -> bool {
    htsjdk_bgzf::read::decompress_all(bam)
        .ok()
        .map(|bytes| {
            // The header text starts at byte 8, after the magic and its length.
            let length = bytes
                .get(4..8)
                .map(|four| i32::from_le_bytes(four.try_into().unwrap_or_default()) as usize)
                .unwrap_or(0);
            String::from_utf8_lossy(bytes.get(8..8 + length).unwrap_or_default()).into_owned()
        })
        .is_some_and(|text| {
            text.lines()
                .find(|line| line.starts_with("@HD"))
                .is_some_and(|line| line.contains("SO:coordinate"))
        })
}

/// Every record's virtual offset in a BAM, plus the offset after the last one and the file's size.
///
/// `SBIIndexWriter` is fed one pointer per record and closed with the position the next record
/// would have gone to, so a walker is what the tool needs from the file and not its records: the
/// bytes are skipped by their own length prefix and never decoded.
fn record_offsets(bam: &[u8]) -> Result<(Vec<u64>, u64), Thrown> {
    use std::io::Read;

    let truncated = || Thrown::user("The BAM ends inside a record".to_string());
    let mut reader = htsjdk_bgzf::BgzfReader::new(bam);

    // The header, which is skipped whole: its text and its reference names decide nothing here.
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic).map_err(|_| truncated())?;
    if magic != htsjdk_bam::writer::BAM_MAGIC {
        return Err(Thrown::user("The file is not a BAM".to_string()));
    }
    let read_i32 = |reader: &mut htsjdk_bgzf::BgzfReader<&[u8]>| -> Result<i32, Thrown> {
        let mut bytes = [0u8; 4];
        reader.read_exact(&mut bytes).map_err(|_| truncated())?;
        Ok(i32::from_le_bytes(bytes))
    };
    let text_length = read_i32(&mut reader)?;
    std::io::copy(
        &mut std::io::Read::by_ref(&mut reader).take(text_length as u64),
        &mut std::io::sink(),
    )
    .map_err(|_| truncated())?;
    let references = read_i32(&mut reader)?;
    for _ in 0..references {
        let name_length = read_i32(&mut reader)?;
        std::io::copy(
            &mut std::io::Read::by_ref(&mut reader).take(name_length as u64),
            &mut std::io::sink(),
        )
        .map_err(|_| truncated())?;
        let _length = read_i32(&mut reader)?;
    }

    let mut offsets = Vec::new();
    let mut next_start = reader.virtual_pos();
    loop {
        let start = reader.virtual_pos();
        let mut size_bytes = [0u8; 4];
        match reader.read_exact(&mut size_bytes) {
            Ok(()) => {}
            Err(_) => break,
        }
        let block_size = i32::from_le_bytes(size_bytes);
        std::io::copy(
            &mut std::io::Read::by_ref(&mut reader).take(block_size as u64),
            &mut std::io::sink(),
        )
        .map_err(|_| truncated())?;
        offsets.push(start);
        next_start = reader.virtual_pos();
    }
    Ok((offsets, next_start))
}

/// `CreateHadoopBamSplittingIndex.doWork`, with the BAM read and the index written.
///
/// Four things the `splitting-index` golden pins and this reproduces: the granularity is refused
/// BEFORE anything is opened; the input's extension is refused next, by its extension and not by
/// its contents; the default output APPENDS `.sbi` where the `.bai` companion REPLACES an
/// extension; and the last entry is where the next record would have gone, which for an empty BAM
/// is the file's own length rather than anything inside it.
pub fn create_hadoop_bam_splitting_index(parser: &Parser) -> Outcome {
    use gatk_tools::create_hadoop_bam_splitting_index as sbi;

    let granularity = scalar(parser, "splitting-index-granularity")
        .and_then(|text| text.parse::<i64>().ok())
        .unwrap_or(sbi::DEFAULT_GRANULARITY as i64);
    // `doWork`'s first line: the argument is judged before the input is looked at, and it is the
    // PARSER that refuses it, so the message names the argument and the status is one.
    sbi::assert_granularity(granularity).map_err(|message| {
        Thrown::command_line(format!(
            "Argument splitting-index-granularity has a bad value: {granularity}. {message}"
        ))
    })?;

    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    // `UserException$BadInput`, whose own constructor puts `Bad input: ` in front of whatever it
    // is handed. The golden carries the prefix and the class both.
    sbi::assert_is_bam(&input).map_err(bad_input)?;

    let bytes = std::fs::read(&input)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&input, false)))?;
    let (offsets, next_start) = record_offsets(&bytes)?;
    let entries = sbi::offsets(&offsets, granularity as u64, next_start);
    let index = sbi::write(
        bytes.len() as u64,
        offsets.len() as u64,
        granularity as u64,
        &entries,
    );

    let output = argument(parser, "output").unwrap_or_else(|| sbi::default_output(&input));
    std::fs::write(&output, index).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;

    if flag(parser, "create-bai") {
        // Only the `.bai` path reads the records, so only it cares how they are sorted.
        if !is_coordinate_sorted(&bytes) {
            return Err(bad_input(
                "Cannot create a .bai index for a file that isn't coordinate sorted.".to_string(),
            ));
        }
        // The companion's name REPLACES the index's extension, so `reads.bam.sbi` becomes
        // `reads.bam.bai` and an output named `elsewhere.idx` becomes `elsewhere.bai`.
        let companion = sbi::bai_companion(&output);
        let bai = htsjdk_bam::build_index::build_bam_index(&bytes)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        std::fs::write(&companion, bai).map_err(|error| {
            Thrown::non_user(
                PORT_FAILURE,
                format!("could not write {companion}: {error}"),
            )
        })?;
    }
    // The tool returns nothing, so `handleResult` prints nothing.
    Ok(None)
}

/// The BGZF compression a tool run writes at: GATK's deflater and GATK's level.
///
/// `--use-jdk-deflater` chooses the first and `GATKConfig`'s `samjdk.compression_level` the
/// second, and neither is htsjdk's own default. Every file a tool writes block compressed depends
/// on both (#1032).
fn output_compression(parser: &Parser) -> (u32, htsjdk_bgzf::Deflater) {
    let level = gatk_tools::gatk_config::compression_level(
        std::env::var(gatk_tools::gatk_config::COMPRESSION_LEVEL)
            .ok()
            .as_deref(),
    );
    let deflater = if flag(parser, "use-jdk-deflater") {
        htsjdk_bgzf::Deflater::Jdk
    } else {
        htsjdk_bgzf::Deflater::Gkl
    };
    (level, deflater)
}

/// `PrintReads.doWork`: the reads that survive the traversal, written back out.
///
/// It is `CountReads` with a writer at the end, and the writer is where the arguments this tool
/// has and that one does not finally reach something: `--create-output-bam-index` decides whether
/// a `.bai` is written beside the BAM, and `--add-output-sam-program-record` whether an `@PG` line
/// is added at all. The `CL` that line carries is the expanded command line, which is why this
/// tool needed [`crate::command_line::expanded`] before it could have a runner.
pub fn print_reads(parser: &Parser) -> Outcome {
    let resolved_filters = resolve_read_filters(parser, "PrintReads")?;

    let inputs = arguments(parser, "input");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --input is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let path = std::path::Path::new(&input);

    let bytes = std::fs::read(path)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&input, false)))?;
    let is_binary = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .map(|inflated| inflated.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC))
            .unwrap_or(false)
    } else {
        bytes.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC)
    };
    let named_index = read_index(
        parser,
        1,
        is_binary,
        gatk_tools::read_walker_refusal::is_block_compressed(&bytes),
    )?;
    let index = if is_binary {
        named_index.or_else(|| htsjdk_bam::sam_files::find_index(path))
    } else {
        None
    };
    // An index that is not one is refused by its MAGIC, and the refusal WAITS: `initializeReads`
    // opens the file, `initializeIntervals` runs next, and a query that does not resolve is
    // refused before anything asks the index what it is.
    let (source, bad_index) = match &index {
        Some(index) => match ReadsDataSource::open(path, index) {
            Ok(source) => (source, None),
            Err(_) => (
                ReadsDataSource::open_unindexed(path)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?,
                Some(
                    index
                        .file_name()
                        .map(|name| name.to_string_lossy().into_owned())
                        .unwrap_or_default(),
                ),
            ),
        },
        None => (
            ReadsDataSource::open_unindexed(path)
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
            None,
        ),
    };

    let header = source.header().clone();
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            validate_against_master(master, "reads", &header.sequences)?;
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
        }
    }
    let best = master
        .clone()
        .or_else(|| reference.clone())
        .unwrap_or_else(|| header.clone());
    let intervals = interval_arguments(parser, &best)?
        .map(|parameters| parameters.intervals)
        .unwrap_or_default();
    if !intervals.is_empty() && index.is_none() {
        return Err(Thrown::user(
            "Traversal by intervals was requested but some input files are not indexed.",
        ));
    }
    if let Some(name) = bad_index {
        return Err(Thrown::non_user(
            gatk_tools::read_walker_refusal::SAM_FORMAT,
            format!("Unknown BAM index file type: {name}"),
        ));
    }
    let filter = read_filter(parser, &resolved_filters, &header)?;

    let command_line = crate::command_line::expanded("PrintReads", parser);
    let options = gatk_tools::print_reads::Options {
        intervals,
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let (bam, bai) =
        gatk_tools::print_reads::print_reads_with(&source, &options, &filter, level, deflater)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let written = bam;
    std::fs::write(&output, &written).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    if let Some(bai) = bai {
        // The index REPLACES the output's extension: `out.bam` is indexed by `out.bai`.
        let companion = std::path::Path::new(&output).with_extension("bai");
        std::fs::write(&companion, bai).map_err(|error| {
            Thrown::non_user(
                PORT_FAILURE,
                format!("could not write {}: {error}", companion.display()),
            )
        })?;
    }
    if flag(parser, "create-output-bam-md5") {
        // The digest APPENDS where the index replaces: `out.bam` is checksummed by `out.bam.md5`,
        // and the file is the thirty-two hex characters and nothing else.
        let digest = format!("{output}.md5");
        std::fs::write(&digest, gatk_tools::gather_bam_files::md5_file(&written)).map_err(
            |error| Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}")),
        )?;
    }
    // The tool returns nothing, so `handleResult` prints nothing.
    Ok(None)
}

/// The `.idx` a VCF WRITER builds, which is not the one `IndexFeatureFile` builds from the file.
///
/// `IndexingVariantContextWriter` hands the creator the position BEFORE each record, absolute in
/// the output stream so the header is counted, and closes it with the whole file's length. Then
/// `setIndexSequenceDictionary` puts the dictionary in the creator's own property map, before
/// `finalizeIndex` appends its statistics -- so the dictionary is `DICT:` properties and the flag
/// that used to carry it is zero. And nothing ever stats the file being written, so its size, its
/// timestamp and its md5 are left at zero where `IndexFeatureFile` fills them in.
fn on_the_fly_index(
    text: &str,
    dictionary: &[(String, i32)],
    path: &str,
    size: i64,
    timestamp: i64,
) -> Vec<u8> {
    use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR, VERSION};
    use htsjdk_tribble::index_write::{BalanceApproach, BuiltIndex, DynamicIndexCreator, Feature};

    let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
    let mut at: i64 = 0;
    for line in text.split_inclusive('\n') {
        let body = line.trim_end_matches('\n');
        if !body.starts_with('#') {
            let columns: Vec<&str> = body.split('\t').collect();
            if columns.len() >= 8 {
                if let Ok(start) = columns[1].parse::<i32>() {
                    let end = columns[7]
                        .split(';')
                        .filter_map(|field| field.split_once('='))
                        .find(|(key, _)| *key == "END")
                        .and_then(|(_, value)| value.parse::<i32>().ok())
                        .unwrap_or(start + columns[3].len() as i32 - 1);
                    creator.add_feature(
                        &Feature {
                            contig: columns[0].to_string(),
                            start,
                            end,
                        },
                        at,
                    );
                }
            }
        }
        at += line.len() as i64;
    }

    let mut properties: Vec<(String, String)> = dictionary
        .iter()
        .map(|(name, length)| (format!("DICT:{name}"), length.to_string()))
        .collect();
    properties.extend(creator.properties());

    let (index_type, contigs, interval_contigs) = match creator.finalize(text.len() as i64) {
        Ok(BuiltIndex::Linear(contigs)) => (LINEAR, contigs, Vec::new()),
        Ok(BuiltIndex::IntervalTree(intervals)) => (INTERVAL_TREE, Vec::new(), intervals),
        Err(_) => (LINEAR, Vec::new(), Vec::new()),
    };

    TribbleIndex {
        index_type,
        version: VERSION,
        indexed_path: format!("file://{path}"),
        // `close()` writes the index with `writeBasedOnFeaturePath`, which STATS the file it has
        // just finished: the size and the timestamp are the written file's. The md5 is not
        // computed and stays empty.
        indexed_file_size: size,
        indexed_file_timestamp: timestamp,
        indexed_file_md5: String::new(),
        // Zero, not the dictionary flag: from version 3 the dictionary is properties.
        flags: 0,
        properties,
        contigs,
        interval_contigs,
    }
    .write()
    .unwrap_or_default()
}

/// One input VCF, read as far as the gather looks at it.
///
/// The gather compares dictionaries, sample lists and record positions, and copies lines: nothing
/// it decides needs an allele parsed, so the file is read as text and the header kept whole.
fn gather_shard(name: &str) -> Result<(gatk_tools::gather_vcfs::Shard, Vec<String>, bool), Thrown> {
    let bytes = std::fs::read(name)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(name, false)))?;
    let compressed = gatk_tools::read_walker_refusal::is_block_compressed(&bytes);
    let text = if compressed {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .ok()
            .and_then(|bytes| String::from_utf8(bytes).ok())
            .ok_or_else(|| {
                Thrown::non_user(
                    gatk_tools::read_walker_refusal::SAM_FORMAT,
                    format!("{name} is not a block compressed file"),
                )
            })?
    } else {
        String::from_utf8_lossy(&bytes).into_owned()
    };

    let mut header = Vec::new();
    let mut dictionary = Vec::new();
    let mut samples = Vec::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if line.starts_with("##") {
            header.push(line.to_string());
            if let Some(body) = line.strip_prefix("##contig=<") {
                if let Some(id) = body
                    .trim_end_matches('>')
                    .split(',')
                    .find_map(|field| field.strip_prefix("ID="))
                {
                    dictionary.push(id.to_string());
                }
            }
        } else if line.starts_with("#CHROM") {
            header.push(line.to_string());
            // The samples are every column past FORMAT, which is the ninth.
            samples = line.split('\t').skip(9).map(str::to_string).collect();
        } else if !line.is_empty() {
            let mut fields = line.split('\t');
            let contig = fields.next().unwrap_or_default().to_string();
            let position = fields
                .next()
                .and_then(|text| text.parse::<i32>().ok())
                .unwrap_or(0);
            records.push((contig, position));
        }
    }
    Ok((
        gatk_tools::gather_vcfs::Shard {
            name: name.to_string(),
            dictionary,
            samples,
            records,
        },
        header,
        compressed,
    ))
}

/// `GatherVcfsCloud.doWork`: the shards' records, in order, under the first shard's header.
///
/// The tool has two paths and this carries one. CONVENTIONAL re-reads and re-writes the records;
/// BLOCK copies the compressed BLOCKS of each input, which is a different set of bytes for the
/// same records and is not something a text writer can produce. The port refuses that path in its
/// own words rather than writing the conventional bytes under its name.
pub fn gather_vcfs_cloud(parser: &Parser) -> Outcome {
    use gatk_tools::gather_vcfs::{Arguments, GatherType};

    let inputs = arguments(parser, "input");
    if inputs.is_empty() {
        return Err(Thrown::command_line(
            "Argument input was missing: Argument 'input' is required",
        ));
    }
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let mut shards = Vec::new();
    let mut lines = Vec::new();
    let mut all_compressed = true;
    for input in &inputs {
        let (shard, header, compressed) = gather_shard(input)?;
        all_compressed &= compressed;
        lines.push((header, input.clone()));
        shards.push(shard);
    }

    let output_is_block_compressed = output.ends_with(".gz") || output.ends_with(".bgz");
    let gather_type = match scalar(parser, "gather-type").as_deref() {
        Some("BLOCK") => GatherType::Block,
        Some("CONVENTIONAL") => GatherType::Conventional,
        _ => GatherType::Automatic,
    };
    let arguments_for_gather = Arguments {
        gather_type,
        ignore_safety_checks: flag(parser, "ignore-safety-checks"),
        disable_contig_ordering_check: flag(parser, "disable-contig-ordering-check"),
        output_is_block_compressed,
        inputs_are_block_compressed: all_compressed,
    };

    let written =
        gatk_tools::gather_vcfs::gather(&shards, &arguments_for_gather).map_err(|error| {
            let class = error.java_class();
            Thrown {
                failure: Failure::User,
                exception: class,
                // `UserException$BadInput`'s constructor puts `Bad input: ` in front of whatever
                // it is handed, and the port's message is what it was handed.
                message: Some(if class.ends_with("$BadInput") {
                    format!("Bad input: {}", error.message())
                } else {
                    error.message()
                }),
            }
        })?;

    // `AUTOMATIC` resolves to BLOCK when everything in sight is block compressed, and BLOCK copies
    // bytes rather than records.
    let effective_block = gather_type == GatherType::Block
        || (gather_type == GatherType::Automatic && all_compressed && output_is_block_compressed);
    if effective_block {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "Block gathering copies the inputs' compressed blocks, which this port does not carry \
             yet. This message is the port's own and not GATK's.",
        ));
    }

    // The dictionary with its LENGTHS, which the writer's index carries as properties where the
    // gather only needed the names.
    let dictionary_lengths: Vec<(String, i32)> = lines[0]
        .0
        .iter()
        .filter_map(|line| line.strip_prefix("##contig=<"))
        .filter_map(|body| {
            let body = body.trim_end_matches('>');
            let name = body.split(',').find_map(|f| f.strip_prefix("ID="))?;
            let length = body
                .split(',')
                .find_map(|f| f.strip_prefix("length="))
                .and_then(|value| value.parse().ok())
                .unwrap_or(0);
            Some((name.to_string(), length))
        })
        .collect();

    // The header is the FIRST shard's, whole, and then every record the gather selected.
    let mut text = String::new();
    for line in &lines[0].0 {
        text.push_str(line);
        text.push('\n');
    }
    let bodies: Vec<Vec<String>> = inputs
        .iter()
        .map(|input| {
            std::fs::read(input)
                .ok()
                .map(|bytes| {
                    let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
                        htsjdk_bgzf::read::decompress_all(&bytes).unwrap_or_default()
                    } else {
                        bytes
                    };
                    String::from_utf8_lossy(&text)
                        .lines()
                        .filter(|line| !line.starts_with('#') && !line.is_empty())
                        .map(str::to_string)
                        .collect()
                })
                .unwrap_or_default()
        })
        .collect();
    for (shard, record) in &written {
        text.push_str(&bodies[*shard][*record]);
        text.push('\n');
    }

    // The bytes on disk, which for a plain output are the text itself and for a `.gz` are that
    // text block compressed.
    let bytes = if output_is_block_compressed {
        let (level, deflater) = output_compression(parser);
        let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
        std::io::Write::write_all(&mut writer, text.as_bytes())
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
        writer
            .into_inner()
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?
    } else {
        text.clone().into_bytes()
    };
    std::fs::write(&output, &bytes).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;

    if flag(parser, "create-output-variant-index") {
        // The index a feature file's name implies: a `.tbi` for a block compressed output and a
        // Tribble `.idx` for a plain one, both APPENDED to the whole name.
        //
        // The Tribble header records the file's URI, its SIZE and its lastModified, so an index
        // built with the zero `Source` defaults to differs from the reference's in bytes it never
        // looks at again. The file has just been written, so its mtime is there to be read.
        let mut source = index_feature_file::Source::new(&output);
        source.timestamp = modified_millis(&output);
        let index = match index_feature_file::index_kind(&output) {
            index_feature_file::IndexKind::Tabix => {
                let (level, deflater) = output_compression(parser);
                index_feature_file::build_tabix(&bytes, &source, &output, deflater, level)
                    .map_err(|refusal| Thrown::user(refusal.message()))?
            }
            // A writer's `.idx` is not the one `IndexFeatureFile` builds from the same file.
            _ => on_the_fly_index(
                &text,
                &dictionary_lengths,
                &output,
                bytes.len() as i64,
                modified_millis(&output),
            ),
        };
        let name = index_feature_file::default_output(&output);
        std::fs::write(&name, index).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {name}: {error}"))
        })?;
    }
    Ok(None)
}

/// `ApplyBQSR.doWork`: every read that survives the filters, recalibrated and written back out.
///
/// It is `PrintReads` with a transformer between the traversal and the writer, and the transformer
/// is where nine arguments of this tool's own finally reach something. The recalibration table is
/// read WHOLE before the traversal starts, because the transformer needs its covariates before it
/// can judge a base.
pub fn apply_bqsr(parser: &Parser) -> Outcome {
    use gatk_engine::bqsr_transformer::ApplyBqsrArguments;

    let resolved_filters = resolve_read_filters(parser, "ApplyBQSR")?;

    let inputs = arguments(parser, "input");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --input is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let recal = argument(parser, "bqsr-recal-file").ok_or_else(|| {
        Thrown::command_line(
            "Argument bqsr-recal-file was missing: Argument 'bqsr-recal-file' is required",
        )
    })?;
    let recal_text = std::fs::read_to_string(&recal)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&recal, false)))?;

    let path = std::path::Path::new(&input);
    let bytes = std::fs::read(path)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&input, false)))?;
    let compressed = gatk_tools::read_walker_refusal::is_block_compressed(&bytes);
    let is_binary = if compressed {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .map(|inflated| inflated.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC))
            .unwrap_or(false)
    } else {
        bytes.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC)
    };
    let named_index = read_index(parser, 1, is_binary, compressed)?;
    let index = if is_binary {
        named_index.or_else(|| htsjdk_bam::sam_files::find_index(path))
    } else {
        None
    };
    // An index that is not one is refused by its MAGIC, and the refusal WAITS: `initializeReads`
    // opens the file, `initializeIntervals` runs next, and a query that does not resolve is
    // refused before anything asks the index what it is.
    let (source, bad_index) = match &index {
        Some(index) => match ReadsDataSource::open(path, index) {
            Ok(source) => (source, None),
            Err(_) => (
                ReadsDataSource::open_unindexed(path)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?,
                Some(
                    index
                        .file_name()
                        .map(|name| name.to_string_lossy().into_owned())
                        .unwrap_or_default(),
                ),
            ),
        },
        None => (
            ReadsDataSource::open_unindexed(path)
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
            None,
        ),
    };

    let header = source.header().clone();
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            validate_against_master(master, "reads", &header.sequences)?;
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
        }
    }
    let best = master
        .clone()
        .or_else(|| reference.clone())
        .unwrap_or_else(|| header.clone());
    let intervals = interval_arguments(parser, &best)?
        .map(|parameters| parameters.intervals)
        .unwrap_or_default();
    if !intervals.is_empty() && index.is_none() {
        return Err(Thrown::user(
            "Traversal by intervals was requested but some input files are not indexed.",
        ));
    }
    if let Some(name) = bad_index {
        return Err(Thrown::non_user(
            gatk_tools::read_walker_refusal::SAM_FORMAT,
            format!("Unknown BAM index file type: {name}"),
        ));
    }
    let filter = read_filter(parser, &resolved_filters, &header)?;

    let bqsr = ApplyBqsrArguments {
        preserve_qscores_less_than: number_or(parser, "preserve-qscores-less-than", 6),
        quantization_levels: number_or(parser, "quantize-quals", 0),
        static_quantization_quals: arguments(parser, "static-quantized-quals")
            .iter()
            .filter_map(|value| value.parse().ok())
            .collect(),
        round_down: flag(parser, "round-down-quantized"),
        emit_original_quals: flag(parser, "emit-original-quals"),
        use_original_base_qualities: flag(parser, "use-original-qualities"),
        global_qscore_prior: scalar(parser, "global-qscore-prior")
            .and_then(|value| value.parse().ok())
            .unwrap_or(-1.0),
        allow_missing_read_groups: flag(parser, "allow-missing-read-group"),
    };

    let command_line = crate::command_line::expanded("ApplyBQSR", parser);
    let options = gatk_tools::print_reads::Options {
        intervals,
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let (bam, bai) = gatk_tools::apply_bqsr::apply_bqsr_with(
        &source,
        &recal_text,
        &bqsr,
        &options,
        &filter,
        level,
        deflater,
    )
    .map_err(|error| Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message()),
    })?;

    std::fs::write(&output, &bam).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    if let Some(bai) = bai {
        let companion = std::path::Path::new(&output).with_extension("bai");
        std::fs::write(&companion, bai).map_err(|error| {
            Thrown::non_user(
                PORT_FAILURE,
                format!("could not write {}: {error}", companion.display()),
            )
        })?;
    }
    if flag(parser, "create-output-bam-md5") {
        let digest = format!("{output}.md5");
        std::fs::write(&digest, gatk_tools::gather_bam_files::md5_file(&bam)).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}"))
        })?;
    }
    Ok(None)
}

/// An integer argument with the tool's own default where it was not given.
fn number_or(parser: &Parser, long_name: &str, default: i32) -> i32 {
    scalar(parser, long_name)
        .and_then(|text| text.parse().ok())
        .unwrap_or(default)
}

/// `SelectVariants.doWork`: the records the arguments select, written as a VCF.
///
/// The tool is a variant walker like `CountVariants`, so the whole startup is shared -- and shared
/// rather than copied because the order of the refusals is what a covering-array row measures.
/// What follows the startup is the pipeline the five `select-variants-*` suites measure, in the
/// reference's own order: the queue is drained as far as the record about to be read, the record
/// is filtered, subset, filtered again, no-called and dropped from, and joins the queue rather
/// than the file.
///
/// # What this refuses rather than approximates
///
/// Six argument groups reach behaviour no static in this repository reproduces: a pedigree and
/// its Mendelian violations, the two random fractions, the GenomicsDB-only decoding, the
/// concordance tracks, `--variant-output-filtering` and `--fully-decode`. Each is refused when it
/// is SET, which is a port limitation with the tool's name on it rather than a silent difference:
/// a run that ignored `--select-random-fraction 0.5` would answer a question it was not asked.
pub fn select_variants(parser: &Parser) -> Outcome {
    let VariantStart {
        input,
        text,
        codec: _,
        intervals,
    } = variant_walker_startup(parser, "SelectVariants")?;

    select_variants_limits(parser)?;

    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;

    // `createSampleNameInclusionList(vcfHeaders)`, over the driving variants' own samples.
    let sample_arguments = gatk_tools::select_variants::SampleArguments {
        sample_names: arguments(parser, "sample-name"),
        sample_expressions: arguments(parser, "sample-expressions"),
        exclude_sample_names: arguments(parser, "exclude-sample-name"),
        exclude_sample_expressions: arguments(parser, "exclude-sample-expressions"),
        allow_nonoverlapping_command_line_samples: flag(
            parser,
            "allow-nonoverlapping-command-line-samples",
        ),
    };
    let selection = gatk_tools::select_variants::create_sample_name_inclusion_list(
        &file.header.samples,
        &sample_arguments,
    )
    .map_err(|refusal| Thrown {
        failure: Failure::User,
        exception: refusal.java_class(),
        message: Some(refusal.message()),
    })?;

    let subset_arguments = gatk_tools::select_variants::SubsetArguments {
        remove_unused_alternates: flag(parser, "remove-unused-alternates"),
        preserve_alleles: flag(parser, "preserve-alleles"),
        keep_original_chr_counts: flag(parser, "keep-original-ac"),
        keep_original_depth: flag(parser, "keep-original-dp"),
    };
    let filter_arguments = select_variants_filters(parser)?;
    let output_arguments = gatk_tools::select_variants::OutputArguments {
        set_filtered_genotypes_to_no_call: flag(parser, "set-filtered-gt-to-nocall"),
        info_annotations_to_drop: arguments(parser, "drop-info-annotation"),
        genotype_annotations_to_drop: arguments(parser, "drop-genotype-annotation"),
    };

    // `--sites-only-vcf-output` empties the sample columns, and it does so on the HEADER as well
    // as on every record, which is why it is read before the header is built.
    let sites_only = flag(parser, "sites-only-vcf-output");
    let header = gatk_tools::select_variants_header::output_header(
        &file.header,
        &gatk_tools::select_variants_header::HeaderArguments {
            keep_original_chr_counts: subset_arguments.keep_original_chr_counts,
            keep_original_depth: subset_arguments.keep_original_depth,
            info_annotations_to_drop: output_arguments.info_annotations_to_drop.clone(),
            genotype_annotations_to_drop: output_arguments.genotype_annotations_to_drop.clone(),
            add_output_vcf_command_line: flag(parser, "add-output-vcf-command-line"),
            tool_command_line: command_line_header_line(parser),
            samples: if sites_only {
                Vec::new()
            } else {
                selection.samples.clone()
            },
        },
    );

    // The traversal, which is the intervals' if there are any. A feature file with no index is
    // refused here rather than earlier, exactly as `CountVariants`' is.
    let located: Vec<LocatedRecord> = file
        .records
        .iter()
        .enumerate()
        .map(|(index, record)| LocatedRecord {
            index,
            contig: record.contig.clone(),
            start: record.start as i32,
            stop: record.stop as i32,
        })
        .collect();
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!(
                "Input {input} must support random access to enable traversal by intervals. \
                 If it's a file, please index it using the bundled tool IndexFeatureFile"
            )),
        });
    }

    let mut pending: gatk_tools::select_variants::PendingWriter<
        htsjdk_vcf::variant::VariantContext,
    > = gatk_tools::select_variants::PendingWriter::new();
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for located in gatk_engine::variant_source::traverse(&located, intervals.as_deref()) {
        let original = &file.records[located.index];
        // `apply` drains BEFORE it looks at the record, which is what lets a record trimmed onto a
        // later start be written first.
        for (_, vc) in pending.drain_before(&original.contig, original.start as i32) {
            written.push(vc);
        }

        let bridged = crate::variant_bridge::to_engine(original);
        if !gatk_tools::select_variants::keeps_before_subset(
            &bridged.record,
            &bridged.filter_record,
            &filter_arguments,
            &selection,
        )
        .map_err(select_error)?
        {
            continue;
        }

        let subset = gatk_tools::select_variants::subset_record(
            &bridged.record,
            &selection,
            &subset_arguments,
        )
        .map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.GATKException",
            message: Some(error.message()),
        })?;
        // The second round of JEXL sees the record the subset produced, not the one that was read.
        let after = crate::variant_bridge::to_engine(&crate::variant_bridge::from_engine(
            original, &subset,
        ));
        if !gatk_tools::select_variants::keeps_after_subset(
            &subset,
            &after.filter_record,
            &filter_arguments,
        )
        .map_err(select_error)?
        {
            continue;
        }

        let mut record = subset;
        if output_arguments.set_filtered_genotypes_to_no_call {
            gatk_tools::select_variants::set_filtered_genotypes_to_no_call(&mut record);
        }
        gatk_tools::select_variants::drop_annotations(&mut record, &output_arguments);
        // The file's own record is built HERE and carried through the queue: the queue reorders,
        // and nothing downstream could pair a reordered record with the one it was decoded from.
        let vc = crate::variant_bridge::from_engine(original, &record);
        pending.add(record, vc);
    }
    for (_, vc) in pending.drain() {
        written.push(vc);
    }

    if sites_only {
        for record in &mut written {
            record.genotypes.clear();
        }
    }

    let text = htsjdk_vcf::vcf_file::write_vcf(&header, &written).map_err(|error| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException",
        message: Some(format!("{error:?}")),
    })?;
    write_variant_output(parser, &output, &text)?;
    Ok(None)
}

/// A decoded record's position, which is all the traversal needs to select it.
struct LocatedRecord {
    index: usize,
    contig: String,
    start: i32,
    stop: i32,
}

impl gatk_engine::variant_source::Located for LocatedRecord {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn stop(&self) -> i32 {
        self.stop
    }
}

/// A `SelectError` as the reference throws it.
fn select_error(error: gatk_tools::select_variants::SelectError) -> Thrown {
    Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message()),
    }
}

/// The arguments that decide which records survive, read off the command line.
fn select_variants_filters(
    parser: &Parser,
) -> Result<gatk_tools::select_variants::FilterArguments, Thrown> {
    use gatk_tools::select_variants::{AlleleRestriction, VariantType};

    fn types(parser: &Parser, name: &str) -> Result<Vec<VariantType>, Thrown> {
        arguments(parser, name)
            .iter()
            .map(|value| match value.as_str() {
                "NO_VARIATION" => Ok(VariantType::NoVariation),
                "SNP" => Ok(VariantType::Snp),
                "MNP" => Ok(VariantType::Mnp),
                "INDEL" => Ok(VariantType::Indel),
                "SYMBOLIC" => Ok(VariantType::Symbolic),
                "MIXED" => Ok(VariantType::Mixed),
                other => Err(Thrown::command_line(format!(
                    "'{other}' is not a valid value for {name}."
                ))),
            })
            .collect()
    }

    let restriction = match argument(parser, "restrict-alleles-to").as_deref() {
        Some("BIALLELIC") => AlleleRestriction::Biallelic,
        Some("MULTIALLELIC") => AlleleRestriction::Multiallelic,
        _ => AlleleRestriction::All,
    };
    Ok(gatk_tools::select_variants::FilterArguments {
        types_to_include: types(parser, "select-type-to-include")?,
        types_to_exclude: types(parser, "select-type-to-exclude")?,
        allele_restriction: restriction,
        max_indel_size: number_or(parser, "max-indel-size", i32::MAX),
        min_indel_size: number_or(parser, "min-indel-size", 0),
        keep_ids: arguments(parser, "keep-ids"),
        exclude_ids: arguments(parser, "exclude-ids"),
        exclude_filtered: flag(parser, "exclude-filtered"),
        exclude_non_variants: flag(parser, "exclude-non-variants"),
        max_filtered_genotypes: number_or(parser, "max-filtered-genotypes", i32::MAX),
        min_filtered_genotypes: number_or(parser, "min-filtered-genotypes", 0),
        max_fraction_filtered_genotypes: fraction(parser, "max-fraction-filtered-genotypes", 1.0),
        min_fraction_filtered_genotypes: fraction(parser, "min-fraction-filtered-genotypes", 0.0),
        max_nocall_number: number_or(parser, "max-nocall-number", i32::MAX),
        max_nocall_fraction: fraction(parser, "max-nocall-fraction", 1.0),
        select_expressions: arguments(parser, "select"),
        select_genotype_expressions: arguments(parser, "select-genotype-expressions"),
        invert_select: flag(parser, "invertSelect"),
        apply_jexl_filters_first: flag(parser, "apply-jexl-filters-first"),
    })
}

/// A `double` argument, or the declared default when it was not given.
fn fraction(parser: &Parser, long_name: &str, default: f64) -> f64 {
    scalar(parser, long_name)
        .and_then(|text| text.parse().ok())
        .unwrap_or(default)
}

/// The six argument groups `SelectVariants` has and this port does not.
///
/// Each is refused when it is SET rather than ignored: a run that quietly dropped
/// `--select-random-fraction 0.5` would answer a question it was not asked, and a refusal that
/// names the port is the honest form of a gap (`gatk_rs::PortLimitation`).
fn select_variants_limits(parser: &Parser) -> Result<(), Thrown> {
    let mut refused: Vec<&str> = Vec::new();
    if argument(parser, "pedigree").is_some() {
        refused.push("--pedigree");
    }
    for flagged in [
        "mendelian-violation",
        "invert-mendelian-violation",
        "call-genotypes",
    ] {
        if flag(parser, flagged) {
            refused.push(match flagged {
                "mendelian-violation" => "--mendelian-violation",
                "invert-mendelian-violation" => "--invert-mendelian-violation",
                _ => "--call-genotypes",
            });
        }
    }
    if fraction(parser, "select-random-fraction", 1.0) != 1.0 {
        refused.push("--select-random-fraction");
    }
    if fraction(parser, "remove-fraction-genotypes", 0.0) != 0.0 {
        refused.push("--remove-fraction-genotypes");
    }
    if argument(parser, "concordance").is_some() {
        refused.push("--concordance");
    }
    if argument(parser, "discordance").is_some() {
        refused.push("--discordance");
    }
    if argument(parser, "variant-output-filtering").is_some() {
        refused.push("--variant-output-filtering");
    }
    if refused.is_empty() {
        return Ok(());
    }
    Err(Thrown::non_user(
        PORT_LIMITATION,
        format!(
            "SelectVariants in this port does not implement {}: a pedigree's Mendelian \
             violations, the two random fractions, the concordance tracks, the genotype caller \
             and the output filtering mode each reach behaviour no measured static reproduces, \
             and answering without them would be a different answer rather than a refusal",
            refused.join(", ")
        ),
    ))
}

/// `##GATKCommandLine=<ID=...,CommandLine="...",Version=...,Date=...>`, or nothing.
///
/// The four fields are in the reference's own order, and the value of the third is the toolkit
/// version this port claims. The fourth is the run's own wall-clock time, which is why the header
/// construction takes this as an input rather than building it: a golden of a file carrying it
/// would move on every run, and the `select-variants-header` suite elides it for that reason.
fn command_line_header_line(parser: &Parser) -> Option<htsjdk_vcf::header::HeaderLine> {
    if !flag(parser, "add-output-vcf-command-line") {
        return None;
    }
    Some(htsjdk_vcf::header::HeaderLine::Structured {
        key: "GATKCommandLine".to_string(),
        fields: vec![
            ("ID".to_string(), "SelectVariants".to_string()),
            (
                "CommandLine".to_string(),
                crate::command_line::expanded("SelectVariants", parser),
            ),
            ("Version".to_string(), crate::TOOLKIT_VERSION.to_string()),
            ("Date".to_string(), display_date_time()),
        ],
    })
}

/// `Utils.getDateTimeForDisplay(ZonedDateTime.now())`, which is
/// `DateTimeFormatter.ofLocalizedDateTime(FormatStyle.LONG)` under the US locale the reference
/// pins: `September 3, 2026 at 1:38:54 AM UTC`, measured from the pinned container.
///
/// The ZONE here is UTC where the reference uses the machine's own, and the difference is
/// deliberate rather than overlooked: the field holds the instant the run happened, so no two runs
/// agree on it and no golden can compare it. The pinned container runs UTC, which is where every
/// measurement of this port is made. Reproducing a local zone's abbreviation would need a tz
/// database for a field nothing checks.
fn display_date_time() -> String {
    const MONTHS: [&str; 12] = [
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ];
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|since| since.as_secs() as i64)
        .unwrap_or(0);
    let days = seconds.div_euclid(86_400);
    let time_of_day = seconds.rem_euclid(86_400);
    // `days` since 1970-01-01 to a civil date, by Howard Hinnant's `civil_from_days`.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let day_of_era = z - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let year = year_of_era + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let shifted_month = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * shifted_month + 2) / 5 + 1;
    let month = if shifted_month < 10 {
        shifted_month + 3
    } else {
        shifted_month - 9
    };
    let year = if month <= 2 { year + 1 } else { year };

    let hour24 = time_of_day / 3600;
    let minute = (time_of_day % 3600) / 60;
    let second = time_of_day % 60;
    // A twelve-hour clock, where midnight is 12 AM and noon is 12 PM.
    let hour = match hour24 % 12 {
        0 => 12,
        other => other,
    };
    let meridiem = if hour24 < 12 { "AM" } else { "PM" };
    format!(
        "{} {}, {} at {}:{:02}:{:02} {} UTC",
        MONTHS[(month - 1) as usize],
        day,
        year,
        hour,
        minute,
        second,
        meridiem
    )
}

/// The VCF a variant-writing tool leaves behind: the text, block compressed where the name says so,
/// with the index the arguments ask for beside it.
fn write_variant_output(parser: &Parser, output: &str, text: &str) -> Result<(), Thrown> {
    let block_compressed = output.ends_with(".gz") || output.ends_with(".bgz");
    let bytes = if block_compressed {
        let (level, deflater) = output_compression(parser);
        let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
        std::io::Write::write_all(&mut writer, text.as_bytes())
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
        writer
            .into_inner()
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?
    } else {
        text.as_bytes().to_vec()
    };
    std::fs::write(output, &bytes).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;

    if !flag(parser, "create-output-variant-index") {
        return Ok(());
    }
    let dictionary: Vec<(String, i32)> = text
        .lines()
        .take_while(|line| line.starts_with('#'))
        .filter_map(|line| {
            let body = line.strip_prefix("##contig=<")?.trim_end_matches('>');
            let name = body.split(',').find_map(|f| f.strip_prefix("ID="))?;
            let length = body
                .split(',')
                .find_map(|f| f.strip_prefix("length="))
                .and_then(|value| value.parse().ok())
                .unwrap_or(0);
            Some((name.to_string(), length))
        })
        .collect();
    let mut source = index_feature_file::Source::new(output);
    source.timestamp = modified_millis(output);
    let index = match index_feature_file::index_kind(output) {
        index_feature_file::IndexKind::Tabix => {
            let (level, deflater) = output_compression(parser);
            gatk_tools::index_feature_file::build_tabix(&bytes, &source, output, deflater, level)
                .map_err(|refusal| Thrown {
                    failure: Failure::User,
                    exception: refusal.java_class(),
                    message: Some(refusal.message()),
                })?
        }
        _ => on_the_fly_index(
            text,
            &dictionary,
            output,
            bytes.len() as i64,
            source.timestamp,
        ),
    };
    let companion = index_feature_file::default_output(output);
    std::fs::write(&companion, index).map_err(|error| {
        Thrown::non_user(
            PORT_FAILURE,
            format!("could not write {companion}: {error}"),
        )
    })
}
