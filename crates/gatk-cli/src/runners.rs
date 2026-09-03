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
    resolve_read_filters_in(parser, tool).map_err(|error| Thrown {
        // Every one of them is a `CommandLineException`, which is status ONE.
        failure: Failure::CommandLine,
        exception: error.class,
        message: Some(error.message),
    })
}

/// The same resolution as the parser's own, which is where it FIRST runs.
///
/// The descriptor validates while the command line is parsed, so by the time a runner asks, the
/// answer is either already a refusal or cannot become one. It is asked twice rather than cached
/// because the resolution is a pure function of four arguments and the tool's defaults.
pub(crate) fn resolve_read_filters_in(
    parser: &Parser,
    tool: &str,
) -> Result<Vec<gatk_tools::filter_resolution::ResolvedFilter>, gatk_barclay::Error> {
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
    .map_err(|error| gatk_barclay::Error {
        class: error.java_class(),
        message: error.message(),
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

/// Everything `GATKTool.onStartup` does for a READ walker, up to the record the traversal reads.
///
/// Shared rather than copied, because the ORDER is what a covering-array row over any of these
/// tools measures: `loadMasterSequenceDictionary`, `initializeReference`, `initializeReads`,
/// `initializeIntervals`, `validateSequenceDictionaries`, the traversal bounds, and only then the
/// record parse. A port that asked the reader first answered a later question first (#69), and
/// three tools of one archetype copying that order three times is three chances to get it wrong.
///
/// The filter is NOT built here. `read_filter` borrows the header, and a struct that owned both
/// would be self-referential; the caller builds it in one line from what this returns.
struct ReadWalkerStart {
    source: ReadsDataSource,
    header: SamHeader,
    intervals: Vec<gatk_engine::interval::SimpleInterval>,
    filters: Vec<gatk_tools::filter_resolution::ResolvedFilter>,
}

fn read_walker_startup(parser: &Parser, tool: &str) -> Result<ReadWalkerStart, Thrown> {
    let resolved_filters = resolve_read_filters(parser, tool)?;
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
    Ok(ReadWalkerStart {
        source,
        header,
        intervals,
        filters: resolved_filters,
    })
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
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CountReads")?;

    let filter = read_filter(parser, &filters, &header)?;
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
    let _ = resolve_read_filters(parser, "CountVariants")?;
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
    // The failure follows the CLASS rather than being assumed: the transformer throws a
    // `GATKException`, which leaves exit 3 and is printed as a stack trace rather than as
    // `A USER ERROR has occurred`.
    .map_err(|error| Thrown {
        failure: if error.is_user() {
            Failure::User
        } else {
            Failure::Other
        },
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

/// `CountBases.doWork`: the same traversal `CountReads` runs, summing lengths instead of records.
///
/// The whole startup is shared, which is the point of doing an archetype rather than a tool: these
/// two declare the same seventy arguments as `CountReads` and differ in one line of `apply` and in
/// what is printed.
pub fn count_bases(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CountBases")?;

    let filter = read_filter(parser, &filters, &header)?;
    let count = gatk_tools::count_reads::count_bases(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;

    if let Some(output) = argument(parser, "output") {
        // `print`, not `println`: the file is the number's digits and nothing else.
        std::fs::write(&output, gatk_tools::count_reads::output(count))
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    }
    Ok(Some(count.to_string()))
}

/// `FlagStat.doWork`: thirteen counters and their percentages, over the same traversal.
///
/// The counters need the record's contig AND its mate's, because two of them ask whether the mate
/// is on a different one, so the traversal carries the header's names rather than the indices.
pub fn flag_stat(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "FlagStat")?;

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

    if let Some(output) = argument(parser, "output") {
        std::fs::write(&output, &text)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    }
    // `onTraversalSuccess` returns the report itself, which `handleResult` prints.
    Ok(Some(text))
}

/// The contig a reference index names, or nothing where the index is `-1`.
fn contig_name(header: &SamHeader, index: i32) -> Option<&str> {
    usize::try_from(index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.as_str())
}

/// `CountBasesInReference.doWork`: every base of the traversal, counted by its byte.
///
/// A `ReferenceWalker` is the first archetype here whose traversal is the FASTA rather than a file
/// of records, so almost none of the read walker's startup applies: there are no reads to open, no
/// index to find and no dictionaries to compare. What is left is the reference itself and the
/// intervals over it, and the intervals are the tool's own -- `reference_walker::traverse` resolves
/// them against the FASTA's dictionary, which is the only dictionary a run of this tool has.
pub fn count_bases_in_reference(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "CountBasesInReference")?;
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut source =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    // `getBestAvailableSequenceDictionary` prefers a `--sequence-dictionary` over the reference's
    // own, and the intervals resolve against THAT. A run naming both therefore resolves `-L`
    // against the master and then queries the FASTA, which is how a contig the master declares and
    // the FASTA does not reaches the query rather than the parser (measured on rows 4, 6 and 9 of
    // this tool's array; rows 5, 7 and 8 are the same rule the other way round).
    let master = master_dictionary(parser)?;
    let own = gatk_tools::reference_walker::dictionary(&source);
    let best = master.unwrap_or(own);
    let intervals = match interval_arguments(parser, &best)? {
        Some(parameters) => parameters.intervals,
        // `getTraversalIntervals` with no interval argument at all: one interval per contig of the
        // dictionary, covering all of it.
        None => best
            .sequences
            .iter()
            .map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
                    .expect("a contig length is at least one")
            })
            .collect(),
    };
    let counts = gatk_tools::count_bases_in_reference::run_over(&mut source, &intervals)
        .map_err(reference_traversal_error)?;
    let report = counts.report();

    if let Some(output) = argument(parser, "output") {
        // `print`, not `println`: the rows already carry their own newlines.
        std::fs::write(&output, &report)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    }
    Ok(Some(report))
}

/// A reference walker's traversal failure, as the reference throws it.
///
/// The query's own refusal is a `UserException` naming the contig, which is what a run whose
/// master dictionary declares a contig the FASTA does not ends with.
fn reference_traversal_error(error: gatk_tools::reference_walker::TraversalError) -> Thrown {
    match error {
        gatk_tools::reference_walker::TraversalError::Reference(
            gatk_engine::reference::ReferenceError::UnknownContig(contig),
        ) => Thrown::user(format!(
            "Given reference file does not have data at the requested contig({contig})!"
        )),
        other => Thrown::user(format!("{other:?}")),
    }
}

/// `SplitIntervals.onTraversalStart`, which is the whole tool: `traverse()` is empty.
///
/// The third archetype here, and the first that writes a DIRECTORY. Its dictionary is the best
/// available one -- a `--sequence-dictionary`, else the reference, else the reads or the variants
/// -- and with no `-L` at all the intervals are every contig of it long enough to pass
/// `--min-contig-size`.
///
/// The file each shard is written to is the reference's own format, measured in the container: the
/// header is `@HD VN:1.6` and the `@SQ` lines of the dictionary, and each interval is the five
/// columns `IntervalListWriter` emits.
pub fn split_intervals(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "SplitIntervals")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?.map(|reference| reference.sequences);
    let best = match master.map(|header| header.sequences).or(reference) {
        Some(sequences) => sequences,
        None => {
            return Err(Thrown::user(
                "Reference sequence file or sequence dictionary required for this tool",
            ))
        }
    };
    let header = SamHeader {
        sequences: best.clone(),
        ..SamHeader::default()
    };
    let sequences: Vec<(String, i32)> = best
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();

    // `hasUserSuppliedIntervals()`: with none, the tool builds its own from the dictionary and
    // filters them by `--min-contig-size`, which is what `split` does when it is handed nothing.
    let given = interval_arguments(parser, &header)?.map(|parameters| {
        parameters
            .intervals
            .iter()
            .map(|interval| {
                htsjdk_bam::interval::Interval::new(&interval.contig, interval.start, interval.end)
            })
            .collect::<Vec<_>>()
    });

    let arguments = gatk_tools::split_intervals::Arguments {
        scatter_count: number_or(parser, "scatter-count", 1),
        min_contig_size: number_or(parser, "min-contig-size", 0),
        subdivision_mode: match scalar(parser, "subdivision-mode").as_deref() {
            Some("BALANCING_WITHOUT_INTERVAL_SUBDIVISION") => {
                gatk_engine::interval_list_scatter::ScatterMode::BalancingWithoutIntervalSubdivision
            }
            Some("BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW") => {
                gatk_engine::interval_list_scatter::ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow
            }
            Some("INTERVAL_COUNT") => gatk_engine::interval_list_scatter::ScatterMode::IntervalCount,
            Some("INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER") => {
                gatk_engine::interval_list_scatter::ScatterMode::IntervalCountWithDistributedRemainder
            }
            _ => gatk_engine::interval_list_scatter::ScatterMode::IntervalSubdivision,
        },
        prefix: argument(parser, "interval-file-prefix")
            .unwrap_or_else(|| gatk_tools::split_intervals::DEFAULT_PREFIX.to_string()),
        extension: argument(parser, "extension")
            .unwrap_or_else(|| gatk_tools::split_intervals::DEFAULT_EXTENSION.to_string()),
        num_digits: number_or(
            parser,
            "interval-file-num-digits",
            gatk_tools::split_intervals::DEFAULT_NUMBER_OF_DIGITS,
        ),
        dont_mix_contigs: flag(parser, "dont-mix-contigs"),
    };

    let shards = gatk_tools::split_intervals::split(given.as_deref(), &sequences, &arguments)
        .map_err(|error| Thrown {
            failure: Failure::User,
            exception: error.java_class(),
            message: Some(error.message()),
        })?;

    // `outputDir.mkdir()`, whose failure is a `RuntimeIOException` naming the absolute path.
    let directory = std::path::Path::new(&output);
    if !directory.exists() {
        std::fs::create_dir(directory).map_err(|_| {
            Thrown::non_user(
                "htsjdk.samtools.util.RuntimeIOException",
                format!("Unable to create directory: {}", directory.display()),
            )
        })?;
    }
    // A shard's header carries `SO:coordinate` for four of the five modes and not for the fifth:
    // `preprocessIntervalList` is `sorted()` everywhere but `INTERVAL_SUBDIVISION`, and `sorted()`
    // stamps the order on the copy it returns where `uniqued()` clones the original header and
    // stamps nothing. The port models that; the runner has to ask.
    // ...and `--dont-mix-contigs` takes it away again, whatever the mode: the regrouping builds
    // each shard as `new IntervalList(sequenceDictionary)` and `addall`, a FRESH list whose header
    // carries no sort order, so the stamp `sorted()` had put on the scatterer's output is gone.
    // Measured on rows 10, 14, 20, 21 and 23 of this tool's array, which are exactly the rows that
    // pair the flag with one of the four stamping modes.
    let sort_order =
        if arguments.subdivision_mode.stamps_sort_order() && !arguments.dont_mix_contigs {
            "\tSO:coordinate"
        } else {
            ""
        };
    for (name, list) in &shards {
        let mut text = format!("@HD\tVN:1.6{sort_order}\n");
        for sequence in &best {
            text.push_str(&format!(
                "@SQ\tSN:{}\tLN:{}\n",
                sequence.name, sequence.length
            ));
        }
        for interval in &list.intervals {
            text.push_str(&interval.to_file_line());
            text.push('\n');
        }
        let path = directory.join(name);
        std::fs::write(&path, text).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("{}: {error}", path.display()))
        })?;
    }
    Ok(None)
}

/// `PreprocessIntervals.onTraversalStart`, which like `SplitIntervals` is the whole tool.
///
/// The second interval utility, and it shares that one's shape: the best available dictionary, the
/// intervals over it, one file out. What is its own is the binning -- `--bin-length` chops the
/// padded intervals into fixed pieces -- and the filter that drops a bin whose every base is an N,
/// which is why this one needs the reference's BASES and not only its dictionary.
pub fn preprocess_intervals(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "PreprocessIntervals")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut source =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    // `validateIntervalArgumentCollection`, in its own order: this tool bins and pads by its own
    // two arguments, so it refuses every standard interval argument that would modify the input
    // intervals before it. Each is an `IllegalArgumentException` from `Utils.validateArg`, which
    // is a NON-user failure and leaves exit 3.
    for (wrong, message) in [
        (
            scalar(parser, "interval-set-rule").as_deref() == Some("INTERSECTION"),
            "Interval set rule must be set to UNION.",
        ),
        (
            number(parser, "interval-exclusion-padding") != 0,
            "Interval exclusion padding must be set to 0.",
        ),
        (
            number(parser, "interval-padding") != 0,
            "Interval padding must be set to 0.",
        ),
        (
            scalar(parser, "interval-merging-rule").as_deref() != Some("OVERLAPPING_ONLY"),
            "Interval merging rule must be set to OVERLAPPING_ONLY.",
        ),
    ] {
        if wrong {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                message,
            ));
        }
    }

    // The REFERENCE's dictionary, not the best available one: a `--sequence-dictionary` does not
    // reach this tool's output, and only the reference's `.dict` carries the attributes its `@SQ`
    // lines write through. A run naming both wrote no `M5` where the reference writes one, and a
    // dictionary built from the `.fai` drops it too.
    let from_dict = reference_dictionary(parser)?;
    let own = gatk_tools::reference_walker::dictionary(&source);
    let best = from_dict.unwrap_or(own);
    let sequences: Vec<gatk_tools::preprocess_intervals::Sequence> = best
        .sequences
        .iter()
        .map(|sequence| gatk_tools::preprocess_intervals::Sequence {
            name: sequence.name.clone(),
            length: sequence.length,
            md5: sequence.attributes.get("M5").map(str::to_string),
            uri: sequence.attributes.get("UR").map(str::to_string),
        })
        .collect();

    // This tool's own `Interval`, which is three fields and no strand: an interval list's other
    // two columns are written by the writer rather than carried by the interval.
    let given = interval_arguments(parser, &best)?.map(|parameters| {
        parameters
            .intervals
            .iter()
            .map(|interval| gatk_tools::filter_intervals::Interval {
                contig: interval.contig.clone(),
                start: interval.start,
                end: interval.end,
            })
            .collect::<Vec<_>>()
    });

    // One whole contig at a time, which is what the N filter asks for, and each is read once.
    let mut contigs: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    for sequence in &best.sequences {
        let bases = source
            .query(&sequence.name, 1, sequence.length)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        contigs.insert(sequence.name.clone(), bases);
    }

    let text = gatk_tools::preprocess_intervals::preprocess(
        given.as_deref(),
        &sequences,
        number_or(parser, "bin-length", 1000),
        number_or(parser, "padding", 250),
        |contig| contigs.get(contig).cloned().unwrap_or_default(),
    )
    .map_err(|error| Thrown {
        failure: Failure::CommandLine,
        exception: error.java_class(),
        message: Some(error.message()),
    })?;

    std::fs::write(&output, &text)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

/// `Pileup.apply`, once per locus of the traversal.
///
/// The fourth archetype: a LOCUS walker, whose unit is not a record or a base but the pileup of
/// every read covering one position. The read walker's startup still applies -- the reads are
/// opened, the dictionaries compared and the intervals resolved exactly as they are for
/// `CountReads` -- and what follows it is a different traversal entirely.
///
/// `--metadata` is refused rather than ignored: the features it annotates a locus with are a
/// second data source this runner does not open, and a silent empty column would be a different
/// answer rather than a refusal.
pub fn pileup(parser: &Parser) -> Outcome {
    if !arguments(parser, "metadata").is_empty() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "Pileup's --metadata is a feature source this port does not open yet, and a run that \
             ignored it would answer without the column it asks for. This message is the port's \
             own and not GATK's.",
        ));
    }
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "Pileup")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    // `hasReference()` decides the reference base: without one every locus reports `N`.
    let mut reference = match argument(parser, "reference") {
        Some(path) => Some(
            gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
        ),
        None => None,
    };
    // The contigs whole, so a locus can be answered without a query per position.
    let mut bases: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    if let Some(source) = reference.as_mut() {
        for (name, length) in source.sequences().to_vec() {
            let contig = source
                .query(&name, 1, length as i32)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            bases.insert(name, contig);
        }
    }

    // `MissingContigInSequenceDictionary`: the locus walker checks each interval's contig against
    // the REFERENCE's dictionary rather than the best available one, so a run whose master
    // declares a contig the FASTA does not is refused here rather than answering `N` for every
    // base of it. Measured on row 8 of this tool's array.
    if let Some(source) = reference.as_ref() {
        let known = gatk_tools::reference_walker::dictionary(source);
        for interval in &intervals {
            if !known
                .sequences
                .iter()
                .any(|sequence| sequence.name == interval.contig)
            {
                return Err(Thrown::user(format!(
                    "Contig {} not present in the sequence dictionary {}\n",
                    interval.contig,
                    gatk_tools::sequence_dictionary::pretty_print(&known.sequences)
                )));
            }
        }
    }

    let filter = read_filter(parser, &filters, &header)?;
    // The records the traversal would hand `apply`, unfiltered: the locus walker applies the
    // filter itself, and it does so BEFORE the loci are built, so a filtered read is absent from
    // the pileup rather than present and ignored.
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let applied = gatk_tools::locus_walker::traverse(
        &records,
        &header,
        None,
        if intervals.is_empty() {
            None
        } else {
            Some(&intervals)
        },
        gatk_tools::locus_walker::Options {
            max_depth_per_sample: number_or(parser, "max-depth-per-sample", 0),
            ..gatk_tools::locus_walker::Options::default()
        },
        &filter,
    )
    .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let output_insert_length = flag(parser, "output-insert-length");
    let show_verbose = flag(parser, "show-verbose");
    let mut text = String::new();
    for one in &applied {
        let base = bases
            .get(&one.context.contig)
            .and_then(|contig| contig.get((one.context.position - 1) as usize))
            .map(|base| *base as char)
            .unwrap_or('N');
        text.push_str(&gatk_tools::pileup::line(
            &one.context.pileup,
            base,
            &[],
            output_insert_length,
            show_verbose,
        ));
    }

    std::fs::write(&output, &text)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}
