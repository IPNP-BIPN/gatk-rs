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
        // With no output the report is written to `System.out` as the blocks are read, and a
        // refusal comes after what was printed.
        None => print!("{report}"),
    }
    match refusal {
        Some(refusal) => Err(Thrown::user(refusal.message())),
        // `doWork` returns 0, which `handleResult` prints.
        None => Ok(Some("0".to_string())),
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
    let mut agrees_with_header: Option<bool> = None;
    let mut parameterized: Vec<(gatk_readfilter::Parameterized, bool)> = Vec::new();
    for filter in resolved {
        let name = filter.name.as_str();
        if name == "WellformedReadFilter" {
            wellformed = Some(filter.negated);
        } else if name == "AlignmentAgreesWithHeaderReadFilter" {
            agrees_with_header = Some(filter.negated);
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
        } else if name == "ReadLengthReadFilter" {
            // Mutect2's chain carries this one, and its two bounds are the TOOL's defaults rather
            // than the filter's: thirty and `Integer.MAX_VALUE` where the filter alone would take
            // one and the same maximum. The declaration supplies both, so the fallbacks here are
            // only reached by a tool that declares neither.
            parameterized.push((
                gatk_readfilter::Parameterized::ReadLength {
                    min: number_or(parser, "min-read-length", 1),
                    max: number_or(parser, "max-read-length", i32::MAX),
                },
                filter.negated,
            ));
        } else if name == "MateDistantReadFilter" {
            // `PrintDistantMates`' own default, and the one argument that decides what "distant"
            // means. The filter was ported; only this branch was missing, so a row that kept the
            // tool's defaults refused instead of running.
            parameterized.push((
                gatk_readfilter::Parameterized::MateDistant {
                    threshold: number_or(parser, "mate-too-distant-length", 1000),
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
        if let Some(negated) = agrees_with_header {
            if gatk_readfilter::with_header::alignment_agrees_with_header(read, header) == negated {
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

/// `File.getAbsolutePath`: the working directory in front of a relative name, and nothing else.
///
/// Java does not normalise here, and the empty path is the working directory ITSELF rather than a
/// trailing separator on it, which is what makes `--annotated-intervals=` print `/work`.
fn java_absolute_path(path: &str) -> String {
    if path.starts_with('/') {
        return path.to_string();
    }
    let working = std::env::current_dir()
        .map(|directory| directory.display().to_string())
        .unwrap_or_default();
    if path.is_empty() {
        working
    } else {
        format!("{working}/{path}")
    }
}

/// `IOUtils.canReadFile`, which every copy-number input passes through before it is opened.
///
/// Three messages of one exception, and the port can tell the first two apart: a path with no file
/// at all, and a path that is a directory or a device. The unreadable-permissions case is the third
/// and is not reached here, because a file this port cannot open reads as one it can stat.
fn can_read_file(path: &str) -> Result<(), Thrown> {
    let reason = match std::fs::metadata(path) {
        Err(_) => "The input file does not exist.",
        Ok(metadata) if !metadata.is_file() => "The input file is not a regular file",
        Ok(_) => return Ok(()),
    };
    Err(Thrown {
        failure: Failure::User,
        exception: gatk_tools::read_walker_refusal::COULD_NOT_READ,
        message: Some(format!(
            "Couldn't read file {}. Error was: {reason}",
            java_absolute_path(path)
        )),
    })
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

/// Whether the tool's traversal calls `setTraversalBounds`, which only a walker's does.
///
/// The list is the tools here that extend `GATKTool` and override `traverse()` rather than
/// inheriting a walker's. For them `-L` is still parsed, still validated against the best available
/// dictionary and still refused when it names an unknown contig; what it does not do is bound the
/// reads, which is why an unindexed input reaches them.
fn sets_traversal_bounds(tool: &str) -> bool {
    !matches!(
        tool,
        "SplitIntervals"
            | "PreprocessIntervals"
            | "AnnotateIntervals"
            | "PrintReadsHeader"
            | "GetSampleName"
            | "TransferReadTags"
            | "PostProcessReadsForRSEM"
            | "CalibrateDragstrModel"
    )
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

    // `setTraversalBounds`, which the traversal calls before it reads anything -- if the traversal
    // is a WALKER's. `ReadWalker.traverse` and `LocusWalker.traverse` make that call; a `GATKTool`
    // that overrides `traverse()` never does, so `-L` does not bound its reads source and an
    // unindexed input is not refused at all. Measured on ten rows each of `TransferReadTags` and
    // `PostProcessReadsForRSEM`, whose corpus is query-name sorted and therefore has no index to
    // find: the reference traversed the whole file and the port refused every row.
    if intervals_given && index.is_none() && sets_traversal_bounds(tool) {
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
/// `IndexUtils.createSequenceDictionaryFromFeatureIndex`: the Tribble index's contig names, each
/// at `UNKNOWN_SEQUENCE_LENGTH`, which is zero. `None` without an index or with an empty one.
fn index_dictionary(path: &str) -> Option<SamHeader> {
    let index = index_feature_file::default_output(path);
    if !index.ends_with(".idx") {
        return None;
    }
    let bytes = std::fs::read(&index).ok()?;
    let parsed = htsjdk_tribble::index::TribbleIndex::read(&bytes).ok()?;
    let names = parsed.sequence_names();
    if names.is_empty() {
        return None;
    }
    Some(SamHeader {
        sequences: names
            .into_iter()
            .map(|name| htsjdk_bam::header::SequenceRecord::new(name, 0))
            .collect(),
        ..SamHeader::default()
    })
}

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
/// What a variant walker's startup produces, up to the record the traversal reads.
///
/// Shared rather than copied, for the reason [`read_walker_startup`] is: the ORDER is most of what
/// a covering-array row over any of these tools measures. `--read-index` is counted against the
/// READS inputs a variant walker still opens, the dictionaries are validated master-first and then
/// reference-first, and `-L` resolves against the DRIVING VARIANTS' dictionary rather than the
/// master when the VCF carries `##contig` lines of its own.
struct VariantWalkerStart {
    /// `--variant`, which is a scalar on these tools.
    input: String,
    /// The file's text, decompressed if it was block compressed.
    text: String,
    /// The codec the file's name resolved to.
    codec: gatk_tools::feature_codec::Codec,
    /// `-L`, resolved against the best available dictionary, or `None` when none was given.
    intervals: Option<Vec<gatk_engine::interval::SimpleInterval>>,
}

fn variant_walker_startup(parser: &Parser, tool: &str) -> Result<VariantWalkerStart, Thrown> {
    // A variant walker applies no read filter and still VALIDATES the ones a command line names:
    // the descriptor belongs to the command line rather than to the traversal, so `--read-filter`
    // and its three companions are refused here exactly as they are on a read walker.
    let _ = resolve_read_filters(parser, tool)?;
    // `--variant` is a SCALAR on this tool, where a read walker's `--input` is a collection: the
    // declaration says `collection: false`, and reading it as a list finds nothing at all.
    let input = argument(parser, "variant").ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;
    variant_walker_startup_over(parser, input)
}

/// [`variant_walker_startup`] over a driving file already chosen, which is what a
/// `MultiVariantWalker` with one `--variant` reaches: its collection holds the path, and the rest
/// of the startup is a variant walker's.
fn variant_walker_startup_over(
    parser: &Parser,
    input: String,
) -> Result<VariantWalkerStart, Thrown> {
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
    //
    // A VCF with no `##contig` line has its dictionary SYNTHESIZED from its Tribble index: the
    // index's contig names at an unknown length. That one is used only when no other source has a
    // dictionary (`VariantWalkerBase.getBestAvailableSequenceDictionary`), and an unknown length
    // is not checked against an interval's stop, so `-L` over such a VCF is accepted.
    let best = if header.sequences.is_empty() {
        match master.clone().or_else(|| reference.clone()).or_else(|| {
            reads_dictionaries.first().map(|sequences| SamHeader {
                sequences: sequences.clone(),
                ..SamHeader::default()
            })
        }) {
            Some(other) => other,
            None => index_dictionary(&input).unwrap_or_else(|| header.clone()),
        }
    } else {
        header.clone()
    };
    let intervals = interval_arguments(parser, &best)?.map(|parameters| parameters.intervals);

    Ok(VariantWalkerStart {
        input,
        text,
        codec,
        intervals,
    })
}

pub fn count_variants(parser: &Parser) -> Outcome {
    let VariantWalkerStart {
        input,
        text,
        codec,
        intervals,
        ..
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
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("0".to_string()))
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
    // `onTraversalSuccess` PRINTS the counts itself and returns 0, which `handleResult` prints
    // after them.
    print!("{report}");
    Ok(Some("0".to_string()))
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

    validate_copy_number_intervals(parser)?;

    // The REFERENCE's dictionary, not the best available one: a `--sequence-dictionary` does not
    // reach this tool's output, and only the reference's `.dict` carries the attributes its `@SQ`
    // lines write through. A run naming both wrote no `M5` where the reference writes one, and a
    // dictionary built from the `.fai` drops it too.
    let from_dict = reference_dictionary(parser)?;
    let own = gatk_tools::reference_walker::dictionary(&source);
    let written = from_dict.clone().unwrap_or(own.clone());
    // The INTERVALS resolve against a different dictionary from the one the file is written with:
    // `getBestAvailableSequenceDictionary` prefers a `--sequence-dictionary`, so a `-L` naming a
    // contig the master does not declare is refused even where the reference has it, while the
    // `@SQ` lines still come from the reference. One tool, two dictionaries.
    let best = master_dictionary(parser)?.or(from_dict).unwrap_or(own);
    let sequences: Vec<gatk_tools::preprocess_intervals::Sequence> = written
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

/// What a locus traversal refused with.
///
/// `--max-depth-per-sample` above zero used to be refused here as the port's own limitation. It is
/// not any more: the reservoir and leveling downsamplers are wired into the read-state managers, so
/// the argument thins a pileup the way the reference thins it (#1102).
fn locus_traversal_error(error: gatk_tools::locus_walker::LocusWalkerError) -> Thrown {
    Thrown::user(format!("{error:?}"))
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

    let filter = read_filter(parser, &filters, &header)?;
    // The records the traversal would hand `apply`, unfiltered: the locus walker applies the
    // filter itself, and it does so BEFORE the loci are built, so a filtered read is absent from
    // the pileup rather than present and ignored.
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;

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
    .map_err(locus_traversal_error)?;

    // `MissingContigInSequenceDictionary`, which is raised when a LOCUS asks the reference for its
    // base and not before. Three things follow from that, and all three are measured:
    //
    //   - it is the REFERENCE's dictionary that is consulted rather than the best available one;
    //   - the reads are read first, so a BAM handed another file's index answers htsjdk's failure
    //     (row 8 of `CheckPileup`'s array);
    //   - and a traversal that visits NO locus never asks, so a row whose filters keep no read at
    //     all writes an empty file rather than refusing. Measured on row 4 of this tool's array,
    //     where `--inverted-read-filter PrimaryLineReadFilter` keeps nothing and the reference
    //     wrote an empty pileup over a reference whose only contig is not the reads'.
    if let Some(source) = reference.as_ref() {
        let known = gatk_tools::reference_walker::dictionary(source);
        if let Some(unknown) = applied.iter().find(|one| {
            !known
                .sequences
                .iter()
                .any(|sequence| sequence.name == one.context.contig)
        }) {
            return Err(Thrown::user(format!(
                "Contig {} not present in the sequence dictionary {}\n",
                unknown.context.contig,
                gatk_tools::sequence_dictionary::pretty_print(&known.sequences)
            )));
        }
    }

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

/// `CheckPileup.apply`, once per locus, against the samtools mpileup file `--pileup` names.
///
/// `Pileup`'s traversal with a second data source beside it, and three things the other locus
/// walker does not have:
///
///   - the report is written to a FILE ONLY IF `--output` names one, and to stdout otherwise,
///     because `outFile` is optional and `new PrintStream(System.out)` is what a null becomes;
///   - a failing run still leaves the report behind: the line that explains the disagreement is
///     printed BEFORE the exception is thrown, and `closeTool` flushes the stream on the way out;
///   - and `onTraversalSuccess` returns the counters rather than writing them, so they are the
///     tool's RESULT and a refused run never reports them at all.
///
/// The truth file is read whole rather than queried per locus. The reference reaches it through a
/// `FeatureInput`, which is why the file has to be indexed for the reference to run at all; what
/// the index buys there is the query, and a port that holds every feature answers the same
/// question from memory.
pub fn check_pileup(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CheckPileup")?;
    let truth_path = argument(parser, "pileup").ok_or_else(|| {
        Thrown::command_line("Argument pileup was missing: Argument 'pileup' is required")
    })?;
    // `requiresReference()` is true, so the declaration itself is required and the parser refuses
    // a run without one before this is reached.
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    // The locus walker checks each interval's contig against the REFERENCE's dictionary, which is
    // the same check `Pileup` makes and for the same reason: a master dictionary declaring a
    // contig the FASTA does not is refused here rather than answered with `N`.

    let text = std::fs::read_to_string(&truth_path)
        .map_err(|error| Thrown::user(format!("{truth_path}: {error}")))?;
    let mut truth: Vec<gatk_engine::sam_pileup::SamPileupFeature> = Vec::new();
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let feature = gatk_engine::sam_pileup::decode(line).map_err(|error| Thrown {
            failure: Failure::User,
            exception: error.java_class(),
            message: Some(error.message()),
        })?;
        truth.push(feature);
    }

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;
    let applied = gatk_tools::locus_walker::traverse(
        &records,
        &header,
        Some(&mut reference),
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
    .map_err(locus_traversal_error)?;

    // `MissingContigInSequenceDictionary`, raised when a LOCUS asks the reference for its base: the
    // reads are read first, and a traversal that visits no locus never asks at all. Both halves are
    // measured, on row 8 of this tool's array and on row 4 of `Pileup`'s.
    let known = gatk_tools::reference_walker::dictionary(&reference);
    if let Some(unknown) = applied.iter().find(|one| {
        !known
            .sequences
            .iter()
            .any(|sequence| sequence.name == one.context.contig)
    }) {
        return Err(Thrown::user(format!(
            "Contig {} not present in the sequence dictionary {}\n",
            unknown.context.contig,
            gatk_tools::sequence_dictionary::pretty_print(&known.sequences)
        )));
    }

    let arguments = gatk_tools::check_pileup::CheckPileupArguments {
        ignore_overlaps: flag(parser, "ignore-overlaps"),
        continue_after_error: flag(parser, "continue-after-error"),
    };
    let output = argument(parser, "output");
    let mut report = String::new();
    let mut loci = 0i64;
    let mut bases = 0i64;

    for one in &applied {
        let mut context = one.reference.clone();
        let base = context
            .base(&mut reference)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        // `featureContext.getValues(mpileup)`, which is every feature OVERLAPPING the locus and of
        // which the tool takes the first. A samtools pileup feature is one base wide, so the
        // overlap is an equality.
        let feature = truth.iter().find(|feature| {
            feature.contig == one.context.contig && feature.position == one.context.position
        });
        let (line, error) =
            gatk_tools::check_pileup::apply(&one.context, base, feature, &arguments);
        if let Some(line) = line {
            report.push_str(&line);
        }
        if let Some(error) = error {
            if !arguments.continue_after_error {
                write_report(&output, &report)?;
                return Err(Thrown::user(format!("Bad input: {}", error.message())));
            }
        }
        loci += 1;
        bases += one.context.pileup.size() as i64;
    }

    write_report(&output, &report)?;
    Ok(Some(gatk_tools::check_pileup::summary(loci, bases)))
}

/// The report, to the file `--output` names or to stdout when it names none.
fn write_report(output: &Option<String>, report: &str) -> Result<(), Thrown> {
    match output {
        Some(path) => std::fs::write(path, report)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}"))),
        None => {
            print!("{report}");
            Ok(())
        }
    }
}

/// `FastaReferenceMaker.apply`, which appends one base of the reference per locus.
///
/// A second REFERENCE walker, and its startup is `CountBasesInReference`'s to the line: the
/// intervals resolve against the best available dictionary, which a `--sequence-dictionary`
/// outranks the reference in, and the traversal then queries the FASTA with them.
///
/// What is new is the OUTPUT, which is three files rather than one. `FastaReferenceWriterBuilder`
/// makes the `.fai` and the `.dict` beside the FASTA unless told otherwise, and htsjdk names them:
/// the index is the output plus `.fai`, and the dictionary replaces the FASTA's extension.
pub fn fasta_reference_maker(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "FastaReferenceMaker")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut source =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let master = master_dictionary(parser)?;
    let own = gatk_tools::reference_walker::dictionary(&source);
    let best = master.unwrap_or(own);
    let intervals = match interval_arguments(parser, &best)? {
        Some(parameters) => parameters.intervals,
        None => best
            .sequences
            .iter()
            .map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
                    .expect("a contig length is at least one")
            })
            .collect(),
    };

    let width = number_or(
        parser,
        "line-width",
        gatk_tools::fasta_reference_maker::DEFAULT_LINE_WIDTH as i32,
    )
    .max(0) as usize;
    // A refused traversal still leaves three files: the writer is built in `onTraversalStart` and
    // closed in `closeTool`, so a reference that has no data at a requested contig ends with an
    // empty FASTA, an empty `.fai` and a dictionary of nothing but its `@HD` line. Measured on rows
    // 8 and 18 of this tool's array, where the port wrote nothing at all.
    let refused = |error| -> Thrown {
        match gatk_tools::fasta_reference_maker::empty_outputs(width) {
            Ok(empty) => {
                let _ = write_outputs(&output, &empty);
                fasta_maker_error(error)
            }
            Err(failure) => fasta_maker_error(failure),
        }
    };
    let outputs = gatk_tools::fasta_reference_maker::run_over(&mut source, &intervals, width)
        .map_err(refused)?;

    write_outputs(&output, &outputs)?;
    Ok(None)
}

/// The FASTA and the two files htsjdk writes beside it, named the way htsjdk names them.
fn write_outputs(
    output: &str,
    outputs: &htsjdk_bam::fasta_writer::FastaOutputs,
) -> Result<(), Thrown> {
    write_file(output, &outputs.fasta)?;
    write_file(&format!("{output}.fai"), outputs.index.as_bytes())?;
    write_file(&dictionary_path(output), outputs.dictionary.as_bytes())
}

/// `ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence`: the FASTA's extension
/// replaced by `.dict`, and a name with no extension of that kind simply gains one.
fn dictionary_path(output: &str) -> String {
    for extension in [".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"] {
        if let Some(stem) = output.strip_suffix(extension) {
            return format!("{stem}.dict");
        }
    }
    format!("{output}.dict")
}

fn write_file(path: &str, bytes: &[u8]) -> Result<(), Thrown> {
    std::fs::write(path, bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))
}

/// What `FastaReferenceMaker` refused with, told apart by whose refusal it is.
fn fasta_maker_error(error: gatk_tools::fasta_reference_maker::MakerError) -> Thrown {
    match error {
        gatk_tools::fasta_reference_maker::MakerError::Traversal(traversal) => {
            reference_traversal_error(traversal)
        }
        // The writer's own refusals are `IllegalArgumentException`s, which the non-user handler
        // prints the class of: a `--line-width` of zero is refused before the reference is read.
        gatk_tools::fasta_reference_maker::MakerError::Writer(writer) => {
            Thrown::non_user("java.lang.IllegalArgumentException", format!("{writer:?}"))
        }
    }
}

/// `CheckReferenceCompatibility.traverse`: one BAM or one VCF checked against several references.
///
/// The tool's required argument is not the reference. What it checks is the INPUT, and
/// `initializeSequenceDictionaryForInput` refuses a command line naming both a BAM and a VCF, a
/// second reads input, or neither -- the last of those AFTER the dictionaries have been read,
/// which is why an unreadable reference is refused before an empty command line is.
///
/// One property of the input decides the whole algorithm: with an MD5 on every sequence the
/// comparison is `CompareReferences`' table, and without one it is name and length alone. A VCF
/// never reaches the first path, because `VCFContigHeaderLine.getSAMSequenceRecord` drops `M5`.
pub fn check_reference_compatibility(parser: &Parser) -> Outcome {
    use gatk_tools::check_reference_compatibility as check;
    use gatk_tools::compare_references as compare;

    let _ = resolve_read_filters(parser, "CheckReferenceCompatibility")?;
    let reads = arguments(parser, "input");
    let vcf = argument(parser, "variant");
    let references = arguments(parser, "references-to-compare");
    if references.is_empty() {
        return Err(Thrown::command_line(
            "Argument references-to-compare was missing: Argument 'references-to-compare' is required",
        ));
    }

    // `GATKTool.onStartup` opens the inputs and resolves `-L` against the best available
    // dictionary BEFORE `onTraversalStart` makes this tool's own checks, and the order is
    // observable: a `--sequence-dictionary` outranks everything, so an interval this tool would
    // never traverse is refused before the input pair is even looked at.
    let query_records = match reads.first() {
        Some(path) => reads_dictionary(parser, path, reads.len())?,
        None => Vec::new(),
    };
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    let best = master.or(reference).unwrap_or_else(|| SamHeader {
        sequences: query_records.clone(),
        ..SamHeader::default()
    });
    let _ = interval_arguments(parser, &best)?;

    // The first two refusals come before anything is read; the third comes after.
    if !reads.is_empty() && vcf.is_some() {
        return Err(input_refusal(check::InputError::BothBamAndVcf));
    }
    if reads.len() > 1 {
        return Err(input_refusal(check::InputError::ManyReadInputs));
    }

    let (query_name, query) = if let Some(path) = reads.first() {
        (file_name(path), query_records)
    } else if let Some(path) = &vcf {
        let text = feature_text(path)?;
        // A VCF's dictionary is its `##contig` lines, and the record htsjdk builds from one
        // carries no `M5` whatever the line says.
        (file_name(path), vcf_dictionary(&text).sequences)
    } else {
        (String::new(), Vec::new())
    };

    let mut dictionaries = Vec::new();
    for path in &references {
        let dictionary = std::path::Path::new(path).with_extension("dict");
        let text = std::fs::read_to_string(&dictionary)
            .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(path, false)))?;
        dictionaries.push((
            file_name(path),
            htsjdk_bam::reader::parse_header_text(&text).sequences,
        ));
    }

    if reads.is_empty() && vcf.is_none() {
        return Err(input_refusal(check::InputError::NoInput));
    }

    let md5_of = |records: &[htsjdk_bam::header::SequenceRecord]| -> Vec<Option<String>> {
        records
            .iter()
            .map(|record| record.attributes.get("M5").map(str::to_string))
            .collect()
    };
    let records: Vec<check::Record> = if check::md5s_present(&md5_of(&query)) {
        // `new ReferenceSequenceTable(dictionaries)` forces USE_DICT, so a REFERENCE with no `M5`
        // refuses the run even though the input has one on every sequence.
        let as_reference = |name: &String, records: &[htsjdk_bam::header::SequenceRecord]| {
            compare::Reference {
                column: name.clone(),
                sequences: records
                    .iter()
                    .map(|record| compare::Sequence {
                        name: record.name.clone(),
                        length: record.length as i64,
                        md5: record.attributes.get("M5").map(str::to_string),
                        // Never read: the mode is USE_DICT and every sequence has its `M5`.
                        calculated_md5: String::new(),
                    })
                    .collect(),
            }
        };
        let mut all = vec![as_reference(&query_name, &query)];
        for (name, records) in &dictionaries {
            all.push(as_reference(name, records));
        }
        // `TableError::message()` already carries `UserException$BadInput`'s own prefix.
        let table = compare::build(&all, compare::Md5Mode::UseDict).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            message: Some(error.message()),
        })?;
        // `compareAgainstKeyReference`: the pairs the key is the first half of, which are the ones
        // `compare_all` generates first because the key is index zero.
        let pairs = compare::compare_all(&table, &all).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            message: Some(error.message()),
        })?;
        pairs
            .iter()
            .take(dictionaries.len())
            .enumerate()
            .map(|(index, pair)| {
                check::evaluate_with_md5(pair, &missing_sequences(&query, &dictionaries[index].1))
            })
            .collect()
    } else {
        dictionaries
            .iter()
            // `!entry.getValue().equals(queryDictionary)`: a reference whose dictionary IS the
            // input's produces no row at all.
            .filter(|(_, records)| records != &query)
            .map(|(name, records)| {
                let status = match gatk_tools::sequence_dictionary::compare(records, &query, false)
                {
                    gatk_tools::sequence_dictionary::Compatibility::Identical => {
                        check::DictionaryCompatibility::Identical
                    }
                    gatk_tools::sequence_dictionary::Compatibility::Superset => {
                        check::DictionaryCompatibility::Superset
                    }
                    other => check::DictionaryCompatibility::Other(other.name()),
                };
                check::evaluate_without_md5(
                    name,
                    &query_name,
                    status,
                    &missing_sequences(&query, records),
                )
            })
            .collect()
    };

    let rendered = check::write_table(&query_name, &records);
    match argument(parser, "output") {
        Some(path) => write_file(&path, rendered.as_bytes())?,
        None => print!("{rendered}"),
    }
    Ok(None)
}

/// A feature file's text, decompressed when it is block compressed.
///
/// The codec is chosen by the file's NAME, which is `FeatureManager`'s rule, and a name no codec
/// claims is the refusal `IndexFeatureFile` gives.
fn feature_text(path: &str) -> Result<String, Thrown> {
    if gatk_tools::feature_codec::codec_for(path).is_none() {
        return Err(Thrown::user(
            index_feature_file::Refusal::NoSuitableCodecs {
                path: path.to_string(),
            }
            .message(),
        ));
    }
    let bytes = std::fs::read(path).map_err(|_| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: path.to_string(),
            }
            .message(),
        )
    })?;
    if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .ok()
            .and_then(|bytes| String::from_utf8(bytes).ok())
            .ok_or_else(|| {
                Thrown::non_user(
                    gatk_tools::read_walker_refusal::SAM_FORMAT,
                    format!("{path} is not a block compressed file"),
                )
            })
    } else {
        Ok(String::from_utf8_lossy(&bytes).into_owned())
    }
}

/// `getMissingSequencesIfSubset`: the reference's sequence names the input does not carry, in the
/// reference's own order.
fn missing_sequences(
    query: &[htsjdk_bam::header::SequenceRecord],
    reference: &[htsjdk_bam::header::SequenceRecord],
) -> Vec<String> {
    reference
        .iter()
        .filter(|record| !query.iter().any(|other| other.name == record.name))
        .map(|record| record.name.clone())
        .collect()
}

/// The file's name, which is what every message in this tool prints rather than the path.
fn file_name(path: &str) -> String {
    std::path::Path::new(path)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string())
}

/// `initializeSequenceDictionaryForInput`'s refusals, which are `UserException$BadInput` and whose
/// messages the port already carries with the prefix.
fn input_refusal(error: gatk_tools::check_reference_compatibility::InputError) -> Thrown {
    Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(error.message()),
    }
}

/// `CalculateMixingFractions.apply`: one bucket per sample, filled at singleton het SNPs.
///
/// A variant walker that ALSO opens the reads, which is why its startup is both: the driving
/// variants decide the traversal and the reads answer at each site. The counting itself is the
/// port's, including the two things that make this tool's table what it is -- a bucket nothing was
/// added to has an alt fraction of `0/0`, which is NaN, and the normalizer is the SUM of every
/// bucket's fraction, so one uncounted sample makes every row NaN; and the rows come out in a
/// `HashMap`'s iteration order, which is neither the header's nor alphabetical.
pub fn calculate_mixing_fractions(parser: &Parser) -> Outcome {
    use gatk_tools::calculate_mixing_fractions as mixing;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "CalculateMixingFractions")?;
    // `-L` is the traversal and the driving variants are queried by interval, so an input with no
    // index is refused before a record is read. Measured on two rows of this tool's array, where
    // the reference refused and the port answered with a table of NaN.
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;

    // The reads, read once. The reference asks a `ReadsContext` at every site, which is the same
    // set of records for a corpus this size; what it is not is the same WORK, and that difference
    // is a speed question rather than an answer question.
    let reads_paths = arguments(parser, "input");
    // A variant walker's reads go through the READ FILTERS like any walker's, and the tool counts
    // what survives them. Measured on two rows of this tool's array, where an inverted
    // `PrimaryLineReadFilter` left the reference with no read to count and every fraction NaN
    // while the port counted the lot.
    let resolved = resolve_read_filters(parser, "CalculateMixingFractions")?;
    // The records with the contig each one is on, because a read is filtered to the variant's
    // contig by its reference INDEX and the index is the reads header's, not the VCF's.
    let mut records: Vec<(String, htsjdk_bam::record::BamRecord)> = Vec::new();
    for path in &reads_paths {
        let source =
            gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let sequences = source.header().sequences.clone();
        let header = source.header().clone();
        let filter = read_filter(parser, &resolved, &header)?;
        for read in gatk_tools::read_walker::traverse(&source, &[], &filter)
            .map_err(|error| Thrown::user(format!("{error:?}")))?
        {
            let contig = sequences
                .get(read.reference_index as usize)
                .map(|sequence| sequence.name.clone())
                .unwrap_or_default();
            records.push((contig, read));
        }
    }

    let mut counts: std::collections::HashMap<String, mixing::AltAndTotalReadCounts> =
        std::collections::HashMap::new();
    let spans = gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref());
    for variant in &file.records {
        // `-L` bounds the traversal, so a site outside it is never applied and its reads never
        // counted -- which is one of the ways every row of the table becomes NaN.
        if let Some(spans) = spans {
            let inside = spans.iter().any(|interval| {
                interval.contig == variant.contig
                    && variant.start as i32 >= interval.start
                    && variant.start as i32 <= interval.end
            });
            if !inside {
                continue;
            }
        }
        if !mixing::is_biallelic_singleton_het_snp(variant) {
            continue;
        }
        let Some(sample) = mixing::variant_sample(variant) else {
            continue;
        };
        let alternate = variant.alleles.iter().find(|allele| !allele.is_reference());
        let Some(alt_base) = alternate
            .map(|allele| allele.base_string())
            .and_then(|bases| bases.as_bytes().first().copied())
        else {
            continue;
        };
        let at_site: Vec<htsjdk_bam::record::BamRecord> = records
            .iter()
            .filter(|(contig, _)| contig == &variant.contig)
            .map(|(_, read)| read.clone())
            .collect();
        let site = mixing::site_counts(&at_site, variant.start as i32, alt_base);
        let bucket = counts.entry(sample).or_default();
        bucket.alt += site.alt;
        bucket.total += site.total;
    }

    let rows = mixing::mixing_fractions(&file.header.samples, &counts)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
    write_file(&output, mixing::table(&rows).as_bytes())?;
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `AnnotateVcfWithExpectedAlleleFraction.apply`: one Float INFO field per record.
///
/// The arithmetic is the port's, and its two surprises are the point of the tool: the weights come
/// from the record's genotypes in the VCF's COLUMN order while the fractions come from
/// `getSampleNamesInOrder()`, which is SORTED, and the two arrays are multiplied element by
/// element; and the default tool header lines are added to the set AFTER the header was built from
/// it, so this tool's output carries no `##source=` and no `##GATKCommandLine` where its sibling
/// `AnnotateVcfWithBamDepth` carries both.
pub fn annotate_vcf_with_expected_allele_fraction(parser: &Parser) -> Outcome {
    use gatk_tools::annotate_vcf_with_expected_allele_fraction as expected;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "AnnotateVcfWithExpectedAlleleFraction")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let fractions_path = argument(parser, "mixing-fractions").ok_or_else(|| {
        Thrown::command_line(
            "Argument mixing-fractions was missing: Argument 'mixing-fractions' is required",
        )
    })?;
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }

    // `MixingFraction.readMixingFractions`, which is a `TableReader`: comment lines first, then a
    // column line, and a column the reader asks for by NAME. A table that is not this tool's is
    // refused by the reader rather than misread, and the refusal is an `IllegalArgumentException`
    // naming the column -- status three, not a user error. Measured on the row that hands this
    // tool `GetPileupSummaries`' table, where the port read two of its columns as a sample and a
    // fraction and refused for a reason of its own.
    let table_text = std::fs::read_to_string(&fractions_path)
        .map_err(|error| Thrown::user(format!("{fractions_path}: {error}")))?;
    let mut rows = table_text
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty());
    let columns: Vec<&str> = rows.next().unwrap_or_default().split('\t').collect();
    let column = |name: &str| -> Result<usize, Thrown> {
        columns
            .iter()
            .position(|found| *found == name)
            .ok_or_else(|| {
                Thrown::non_user(
                    "java.lang.IllegalArgumentException",
                    format!("there is no such column: {name}"),
                )
            })
    };
    let sample_column = column("SAMPLE")?;
    let fraction_column = column("MIXING_FRACTION")?;
    let mut table: Vec<(String, f64)> = Vec::new();
    for line in rows {
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(sample), Some(value)) = (fields.get(sample_column), fields.get(fraction_column))
        else {
            continue;
        };
        table.push((sample.to_string(), value.trim().parse().unwrap_or(f64::NAN)));
    }

    let mut file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let fractions =
        expected::fractions_in_sample_order(&table, &file.header.samples).map_err(|error| {
            Thrown {
                failure: Failure::User,
                exception: error.class(),
                message: Some(error.message()),
            }
        })?;

    for record in &mut file.records {
        let value = expected::annotation(record, &fractions);
        record.attributes.push((
            expected::AF_EXP.to_string(),
            htsjdk_vcf::variant::Value::Str(value),
        ));
    }

    // The `AF_EXP` line, and nothing else: the tool's default lines never reach the file.
    let declares_af_exp = file.header.lines.iter().any(|line| {
        matches!(line, htsjdk_vcf::header::HeaderLine::Compound { key, id, .. }
            if key == "INFO" && id == expected::AF_EXP)
    });
    if !declares_af_exp {
        file.header
            .lines
            .push(htsjdk_vcf::header::HeaderLine::Compound {
                key: "INFO".to_string(),
                id: expected::AF_EXP.to_string(),
                number: htsjdk_vcf::header::Cardinality::Fixed(1),
                line_type: htsjdk_vcf::header::LineType::Float,
                description: "expected allele fraction in pooled bam".to_string(),
                extra: Vec::new(),
            });
    }

    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = file
        .records
        .iter()
        .filter(|record| keep(record))
        .cloned()
        .collect();
    // `--sites-only-vcf-output` is the writer's, not the tool's: the samples leave the header and
    // the genotypes leave every record. Measured on eleven rows of this tool's array, where the
    // reference wrote a file with no FORMAT column and the port wrote the genotypes.
    apply_sites_only(parser, &mut file.header, &mut written);
    let rendered =
        htsjdk_vcf::vcf_file::write_vcf(&file.header, &written).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("{error:?}")),
        })?;
    write_variant_output(parser, &output, &rendered)?;
    Ok(None)
}

/// `AnnotateVcfWithBamDepth.apply`: one Integer INFO field per record, counted off the reads.
///
/// The sibling of [`annotate_vcf_with_expected_allele_fraction`], and the difference between them
/// is two statements in the other order: this tool adds its default header lines to the set BEFORE
/// the header is built from it, so `##source=` and `##GATKCommandLine` reach the file where the
/// other tool's do not.
///
/// The count is the port's, and CONTAINMENT is the condition: a read counts when it holds the
/// record's whole span, is not a duplicate, is mapped, passes vendor quality and is longer than
/// one base.
pub fn annotate_vcf_with_bam_depth(parser: &Parser) -> Outcome {
    use gatk_tools::annotate_vcf_with_bam_depth as depth;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "AnnotateVcfWithBamDepth")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }

    let resolved = resolve_read_filters(parser, "AnnotateVcfWithBamDepth")?;
    let mut reads: Vec<(String, htsjdk_bam::record::BamRecord)> = Vec::new();
    for path in arguments(parser, "input") {
        let source =
            gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let sequences = source.header().sequences.clone();
        let header = source.header().clone();
        let filter = read_filter(parser, &resolved, &header)?;
        for read in gatk_tools::read_walker::traverse(&source, &[], &filter)
            .map_err(|error| Thrown::user(format!("{error:?}")))?
        {
            let contig = sequences
                .get(read.reference_index as usize)
                .map(|sequence| sequence.name.clone())
                .unwrap_or_default();
            reads.push((contig, read));
        }
    }

    let mut file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let spans = gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref());
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for record in &file.records {
        if let Some(spans) = spans {
            let inside = spans.iter().any(|interval| {
                interval.contig == record.contig
                    && record.stop as i32 >= interval.start
                    && record.start as i32 <= interval.end
            });
            if !inside {
                continue;
            }
        }
        let at_site: Vec<depth::Read> = reads
            .iter()
            .map(|(contig, read)| depth::Read {
                contig,
                start: read.alignment_start,
                // `getEnd()`: the start plus the cigar's reference length, less one.
                end: read.alignment_start + read.cigar.reference_length() as i32 - 1,
                flags: read.flags,
            })
            .collect();
        written.push(depth::annotate(record, depth::bam_depth(&at_site, record)));
    }

    let declares = file.header.lines.iter().any(|line| {
        matches!(line, htsjdk_vcf::header::HeaderLine::Compound { key, id, .. }
            if key == "INFO" && id == depth::BAM_DEPTH)
    });
    if !declares {
        file.header
            .lines
            .push(htsjdk_vcf::header::HeaderLine::Compound {
                key: "INFO".to_string(),
                id: depth::BAM_DEPTH.to_string(),
                number: htsjdk_vcf::header::Cardinality::Fixed(1),
                line_type: htsjdk_vcf::header::LineType::Integer,
                description: "pooled bam depth".to_string(),
                extra: Vec::new(),
            });
    }

    let keep = variant_output_filter(parser, intervals.as_deref())?;
    written.retain(|record| keep(record));
    apply_sites_only(parser, &mut file.header, &mut written);
    let rendered =
        htsjdk_vcf::vcf_file::write_vcf(&file.header, &written).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("{error:?}")),
        })?;
    write_variant_output(parser, &output, &rendered)?;
    Ok(None)
}

/// `CountFalsePositives.doWork`: two counters, a denominator and a six-column table.
///
/// The first tool here whose answer depends on the SIZE of the interval argument rather than on
/// which records it selects: the denominator is the merged intervals' bases, so two overlapping
/// `-L` arguments contribute their union once. `requiresIntervals()` is true, so a run without
/// `-L` is refused by the parser and never reaches this.
///
/// The counting is the port's, and the `snp` column is not what its name says: `isIndel()` is
/// `getType() == INDEL` and nothing looser, so an MNP, a symbolic allele and a record with no
/// alternate all land in the other bucket.
pub fn count_false_positives(parser: &Parser) -> Outcome {
    use gatk_tools::count_false_positives as counting;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "CountFalsePositives")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let resolved = intervals.unwrap_or_default();
    let selected: Vec<htsjdk_vcf::variant::VariantContext> = file
        .records
        .iter()
        .filter(|record| {
            resolved.iter().any(|interval| {
                interval.contig == record.contig
                    && record.stop as i32 >= interval.start
                    && record.start as i32 <= interval.end
            })
        })
        .cloned()
        .collect();

    let counts = counting::count(&selected);
    let territory = counting::target_territory(&resolved);
    let id = counting::id_from_path(&input);
    write_file(&output, counting::table(&id, counts, territory).as_bytes())?;
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `EvaluateInfoFieldConcordance.apply`: one INFO key of the eval file against one of the truth.
///
/// The first runner here with a SECOND variant input beside the driving one, walked in lockstep by
/// [`gatk_engine::concordance_walker`]. Only true positives are looked at, a record whose key is
/// absent is counted but contributes nothing, and the "absolute" difference is
/// `Math.sqrt(delta * delta)` where `Math.abs` was meant -- all three are the port's, measured by
/// its own suite.
pub fn evaluate_info_field_concordance(parser: &Parser) -> Outcome {
    use gatk_tools::evaluate_info_field_concordance as concordance;

    let _ = resolve_read_filters(parser, "EvaluateInfoFieldConcordance")?;
    let summary = argument(parser, "summary").ok_or_else(|| {
        Thrown::command_line("Argument summary was missing: Argument 'summary' is required")
    })?;
    let eval_path = argument(parser, "evaluation").ok_or_else(|| {
        Thrown::command_line("Argument evaluation was missing: Argument 'evaluation' is required")
    })?;
    let truth_path = argument(parser, "truth").ok_or_else(|| {
        Thrown::command_line("Argument truth was missing: Argument 'truth' is required")
    })?;
    let eval_key = argument(parser, "eval-info-key").ok_or_else(|| {
        Thrown::command_line(
            "Argument eval-info-key was missing: Argument 'eval-info-key' is required",
        )
    })?;
    let truth_key = argument(parser, "truth-info-key").ok_or_else(|| {
        Thrown::command_line(
            "Argument truth-info-key was missing: Argument 'truth-info-key' is required",
        )
    })?;

    let read = |path: &str| -> Result<htsjdk_vcf::reader::VcfFile, Thrown> {
        let text = std::fs::read_to_string(path)
            .map_err(|error| Thrown::user(format!("{path}: {error}")))?;
        htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: failure.error.class(),
            message: Some(failure.error.message()),
        })
    };
    let eval_text = std::fs::read_to_string(&eval_path)
        .map_err(|error| Thrown::user(format!("{eval_path}: {error}")))?;
    let eval_file = read(&eval_path)?;
    let truth_file = read(&truth_path)?;

    // `onTraversalStart`'s two checks, which read the HEADERS and not the records: a key no header
    // declares is refused before the walk.
    let declares = |file: &htsjdk_vcf::reader::VcfFile, key: &str| {
        file.header.lines.iter().any(|line| {
            matches!(line, htsjdk_vcf::header::HeaderLine::Compound { key: kind, id, .. }
                if kind == "INFO" && id == key)
        })
    };
    concordance::check_keys(
        declares(&eval_file, &eval_key),
        &eval_key,
        &eval_path,
        declares(&truth_file, &truth_key),
        &truth_key,
        &truth_path,
    )
    .map_err(|error| Thrown {
        failure: Failure::User,
        exception: error.class(),
        message: Some(error.message()),
    })?;

    let dictionary: Vec<String> = eval_file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Contig { fields, .. } => fields
                .iter()
                .find(|(key, _)| key == "ID")
                .map(|(_, value)| value.clone()),
            _ => None,
        })
        .collect();

    // `-L` bounds the traversal, so a record outside the window is not walked at all and its
    // difference never reaches the mean. Measured on eight rows of this tool's array, where every
    // window gave the reference a different mean and the port the same one.
    let intervals = interval_arguments(parser, &vcf_dictionary(&eval_text))?
        .map(|parameters| parameters.intervals);
    let inside = |record: &htsjdk_vcf::variant::VariantContext| -> bool {
        match &intervals {
            None => true,
            Some(list) => list.iter().any(|interval| {
                interval.contig == record.contig
                    && record.stop as i32 >= interval.start
                    && record.start as i32 <= interval.end
            }),
        }
    };
    let truth_kept: Vec<htsjdk_vcf::variant::VariantContext> = truth_file
        .records
        .iter()
        .filter(|record| inside(record))
        .cloned()
        .collect();
    let eval_kept: Vec<htsjdk_vcf::variant::VariantContext> = eval_file
        .records
        .iter()
        .filter(|record| inside(record))
        .cloned()
        .collect();
    let as_records = |records: &[htsjdk_vcf::variant::VariantContext]| -> Vec<ConcordanceVariant> {
        records
            .iter()
            .map(|record| ConcordanceVariant {
                contig: record.contig.clone(),
                start: record.start as i32,
                filtered: record.is_filtered(),
            })
            .collect()
    };
    let truth_records = as_records(&truth_kept);
    let eval_records = as_records(&eval_kept);
    // `areVariantsAtSameLocusConcordant` is allele equality on this tool: two records at one locus
    // are concordant when their alternates match.
    let steps = gatk_engine::concordance_walker::concordance(
        &truth_records,
        &eval_records,
        &dictionary,
        |_, _| true,
    );

    let mut totals = concordance::Concordance::default();
    for step in &steps {
        if step.state != gatk_engine::concordance_walker::ConcordanceState::TruePositive {
            continue;
        }
        let (Some(truth_index), Some(eval_index)) = (step.truth, step.eval) else {
            continue;
        };
        let truth = &truth_kept[truth_index];
        let eval = &eval_kept[eval_index];
        // `isSNP()` then `isIndel()`, both `getType() ==` and nothing looser, which is the type
        // `RemoveNearbyIndels` measured.
        let eval_type = match gatk_tools::remove_nearby_indels::variant_type(eval) {
            gatk_tools::remove_nearby_indels::VariantType::Snp => concordance::EvalType::Snp,
            gatk_tools::remove_nearby_indels::VariantType::Indel => concordance::EvalType::Indel,
            _ => concordance::EvalType::Other,
        };
        totals.add(
            eval_type,
            info_as_double(eval, &eval_key),
            info_as_double(truth, &truth_key),
        );
    }

    write_file(&summary, totals.table(&eval_key, &truth_key).as_bytes())?;
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// As much of a record as the concordance iterator looks at.
struct ConcordanceVariant {
    contig: String,
    start: i32,
    filtered: bool,
}

impl gatk_engine::concordance_walker::ConcordanceRecord for ConcordanceVariant {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn is_filtered(&self) -> bool {
        self.filtered
    }
}

/// `vc.getAttributeAsDouble(key, 0)`, which is absent rather than zero here: a record with no such
/// key is counted and contributes nothing.
fn info_as_double(record: &htsjdk_vcf::variant::VariantContext, key: &str) -> Option<f64> {
    record
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .and_then(|(_, value)| value.format())
        .and_then(|text| text.parse().ok())
}

/// `CallCopyRatioSegments.doWork`: segments in, called segments out, and a second file beside them.
///
/// No walker, no reference, no reads: eighteen arguments and a table transformed by arithmetic.
/// The statistics are computed TWICE over the copy-neutral segments -- once over all of them to
/// find the outliers, once over what is left to call against -- and the port keeps both passes
/// because a port filtering with the recomputed pair would drop a different set.
///
/// The output is two files: `-O` itself, and the legacy `.igv.seg` beside it, whose path is `-O`
/// with ONE extension removed and `.igv.seg` appended, and whose columns are IGV's rather than the
/// table's: the call comes before the mean.
pub fn call_copy_ratio_segments(parser: &Parser) -> Outcome {
    use gatk_tools::call_copy_ratio_segments as calling;

    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let bound = |name: &str, default: f64| -> f64 {
        scalar(parser, name)
            .and_then(|text| text.parse().ok())
            .unwrap_or(default)
    };
    let lower = bound("neutral-segment-copy-ratio-lower-bound", 0.9);
    let upper = bound("neutral-segment-copy-ratio-upper-bound", 1.1);
    let outlier = bound("outlier-neutral-segment-copy-ratio-z-score-threshold", 2.0);
    let calling_threshold = bound("calling-copy-ratio-z-score-threshold", 2.0);

    let text = std::fs::read_to_string(&input)
        .map_err(|error| Thrown::user(format!("{input}: {error}")))?;
    // The file is a SAM header, a column line and the rows. The header travels into the output
    // unchanged, which is what carries the sample name and the dictionary through.
    let mut header = String::new();
    let mut rows = Vec::new();
    let mut columns_seen = false;
    let mut sample = String::new();
    for line in text.lines() {
        if line.starts_with('@') {
            header.push_str(line);
            header.push('\n');
            if let Some(rest) = line.strip_prefix("@RG\t") {
                for field in rest.split('\t') {
                    if let Some(name) = field.strip_prefix("SM:") {
                        sample = name.to_string();
                    }
                }
            }
            continue;
        }
        if !columns_seen {
            columns_seen = true;
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        rows.push(calling::CopyRatioSegment {
            contig: fields[0].to_string(),
            start: fields[1].parse().unwrap_or_default(),
            end: fields[2].parse().unwrap_or_default(),
            num_points: fields[3].parse().unwrap_or_default(),
            mean_log2_copy_ratio: fields[4].trim().parse().unwrap_or(f64::NAN),
        });
    }

    let calls = calling::make_calls(&rows, lower, upper, outlier, calling_threshold)
        .map_err(|error| Thrown::non_user(error.java_class(), error.message().to_string()))?;

    write_file(
        &output,
        calling::write_called(&header, &rows, &calls).as_bytes(),
    )?;
    write_file(
        &legacy_segments_path(&output),
        calling::write_legacy(&sample, &rows, &calls).as_bytes(),
    )?;
    Ok(None)
}

/// `FilenameUtils.removeExtension(path) + ".igv.seg"`: ONE extension removed, whatever it was, so
/// a `-O` with none simply gains the suffix.
fn legacy_segments_path(output: &str) -> String {
    let (directory, name) = match output.rfind('/') {
        Some(slash) => (&output[..=slash], &output[slash + 1..]),
        None => ("", output),
    };
    let stem = match name.rfind('.') {
        Some(dot) if dot > 0 => &name[..dot],
        _ => name,
    };
    format!("{directory}{stem}.igv.seg")
}

/// `VariantFiltration.apply`: JEXL over the record, JEXL over each genotype, and a mask beside
/// them.
///
/// The largest namespace declared here, eighty-nine arguments, and the filtering itself is the
/// port's: the FILTER column is sorted where a filter was applied, a genotype's FT is `PASS` when
/// the record has an FT column at all, and the cluster window looks at the SNPs around a record
/// rather than at the record.
pub fn variant_filtration(parser: &Parser) -> Outcome {
    use gatk_tools::variant_filtration as filtration;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "VariantFiltration")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }

    let mut file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;

    // `filterExpressions` and `filterNames` are read in step: the nth expression carries the nth
    // name, and a name list shorter than the expression list is the parser's refusal rather than
    // this tool's.
    let compile = |expressions: Vec<String>,
                   names: Vec<String>|
     -> Result<Vec<filtration::MatchExp>, Thrown> {
        expressions
            .iter()
            .enumerate()
            .map(|(index, text)| {
                let name = names.get(index).cloned().unwrap_or_default();
                filtration::MatchExp::new(&name, text).map_err(|error| {
                    Thrown::user(format!("Invalid JEXL expression detected: {error:?}"))
                })
            })
            .collect()
    };
    let site = compile(
        arguments(parser, "filter-expression"),
        arguments(parser, "filter-name"),
    )?;
    let genotype = compile(
        arguments(parser, "genotype-filter-expression"),
        arguments(parser, "genotype-filter-name"),
    )?;

    // `--mask` is a feature file read whole: what the tool asks of it is where its features are,
    // and a BED and a VCF answer that the same way.

    let flag_of = |name: &str| flag(parser, name);
    let tool_arguments = filtration::Arguments {
        cluster_size: number_or(parser, "cluster-size", 3),
        cluster_window: number_or(parser, "cluster-window-size", 0),
        mask_name: argument(parser, "mask-name").unwrap_or_else(|| "Mask".to_string()),
        filter_records_not_in_mask: flag_of("filter-not-in-mask"),
        invert_filter_expression: flag_of("invert-filter-expression"),
        invert_genotype_filter_expression: flag_of("invert-genotype-filter-expression"),
        missing_values_evaluate_as_failing: flag_of("missing-values-evaluate-as-failing"),
        invalidate_previous_filters: flag_of("invalidate-previous-filters"),
        set_filtered_genotypes_to_no_call: flag_of("set-filtered-genotype-to-no-call"),
        mask_extension: number_or(parser, "mask-extension", 0),
    };

    let spans = gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref());
    let kept: Vec<htsjdk_vcf::variant::VariantContext> = file
        .records
        .iter()
        .filter(|record| match spans {
            None => true,
            Some(list) => list.iter().any(|interval| {
                interval.contig == record.contig
                    && record.stop as i32 >= interval.start
                    && record.start as i32 <= interval.end
            }),
        })
        .cloned()
        .collect();

    // `FeatureDataSource` asks the mask for its index at the FIRST query, so a traversal that
    // reaches no record never asks. The refusal below is conditioned on that.
    let queries_the_mask = !kept.is_empty();
    let mask_path = argument(parser, "mask");
    let mask: Vec<(String, i32, i32)> = match &mask_path {
        None => Vec::new(),
        Some(path) => {
            // The mask is QUERIED by interval, and `FeatureDataSource` asks for the index at the
            // first query rather than at startup: a bounded traversal that reaches no record never
            // queries and never refuses. Measured on rows of this tool's array where `-L` names a
            // window the file has no record in, and the reference read the unindexed BED happily.
            if !has_feature_index(path) && queries_the_mask {
                return Err(Thrown::user(format!(
                    "Input {path} must support random access to enable queries by interval. If \
                     it's a file, please index it using the bundled tool IndexFeatureFile"
                )));
            }
            let text = std::fs::read_to_string(path)
                .map_err(|error| Thrown::user(format!("{path}: {error}")))?;
            if path.ends_with(".bed") {
                text.lines()
                    .filter_map(|line| {
                        htsjdk_tribble::bed::decode(line, htsjdk_tribble::bed::StartOffset::One)
                            .ok()
                            .flatten()
                    })
                    .map(|feature| (feature.contig, feature.start, feature.end))
                    .collect()
            } else {
                htsjdk_vcf::reader::read_vcf(&text)
                    .map(|masked| {
                        masked
                            .records
                            .iter()
                            .map(|record| {
                                (
                                    record.contig.clone(),
                                    record.start as i32,
                                    record.stop as i32,
                                )
                            })
                            .collect()
                    })
                    .unwrap_or_default()
            }
        }
    };

    let records: Vec<filtration::Record> = kept.iter().map(filtration_record).collect();
    let filtered = filtration::filter_records(&records, &site, &genotype, &mask, &tool_arguments);

    let mut written = kept.clone();
    for (index, record) in written.iter_mut().enumerate() {
        let rendered = filtration::rendered_filters(&filtered[index]);
        record.filters = match rendered.as_str() {
            "." => None,
            "PASS" => Some(Vec::new()),
            names => Some(names.split(';').map(str::to_string).collect()),
        };
        for sample in 0..record.genotypes.len() {
            if let Some(ft) = filtration::rendered_genotype_filter(&filtered[index], sample) {
                record.genotypes[sample].filters = Some(ft);
            }
            if filtered[index]
                .no_called
                .get(sample)
                .copied()
                .unwrap_or(false)
            {
                record.genotypes[sample].alleles = Vec::new();
            }
        }
    }

    // The header gains one `##FILTER` line per filter this run can apply, in the reference's own
    // order: the clustered-SNP name first when the window is open, then the site expressions, then
    // the genotype ones, then the mask. An inverted expression describes itself as `Inverse of:`.
    let mut filter_lines: Vec<(String, String)> = Vec::new();
    if tool_arguments.cluster_window > 0 {
        filter_lines.push((
            "SnpCluster".to_string(),
            "SNPs found in clusters".to_string(),
        ));
    }
    let described = |text: &str, inverted: bool| -> String {
        if inverted {
            format!("Inverse of: {text}")
        } else {
            text.to_string()
        }
    };
    for (index, expression) in arguments(parser, "filter-expression").iter().enumerate() {
        if let Some(name) = arguments(parser, "filter-name").get(index) {
            filter_lines.push((
                name.clone(),
                described(expression, tool_arguments.invert_filter_expression),
            ));
        }
    }
    // `possiblyInvertFilterExpression` reads `invertFilterExpression` and nothing else, and BOTH
    // loops call it: a genotype filter's DESCRIPTION is inverted by the site flag even though its
    // decision is inverted by the genotype one. Measured on this tool's array.
    for (index, expression) in arguments(parser, "genotype-filter-expression")
        .iter()
        .enumerate()
    {
        if let Some(name) = arguments(parser, "genotype-filter-name").get(index) {
            filter_lines.push((
                name.clone(),
                described(expression, tool_arguments.invert_filter_expression),
            ));
        }
    }
    if mask_path.is_some() {
        let description = argument(parser, "mask-description").unwrap_or_else(|| {
            if tool_arguments.filter_records_not_in_mask {
                "Doesn't overlap a user-input mask".to_string()
            } else {
                "Overlaps a user-input mask".to_string()
            }
        });
        filter_lines.push((tool_arguments.mask_name.clone(), description));
    }
    // `VCFStandardHeaderLines.getFormatLine(FT)` when any genotype expression was given, the
    // chromosome counts when filtered genotypes become no-calls, and the allele-specific status
    // line when that mode is on. All three are the ENGINE's lines rather than the tool's filters.
    let mut compound: Vec<(
        &str,
        &str,
        htsjdk_vcf::header::Cardinality,
        htsjdk_vcf::header::LineType,
        &str,
    )> = Vec::new();
    if !arguments(parser, "genotype-filter-expression").is_empty() {
        compound.push((
            "FORMAT",
            "FT",
            htsjdk_vcf::header::Cardinality::Unbounded,
            htsjdk_vcf::header::LineType::String,
            "Genotype-level filter",
        ));
    }
    if tool_arguments.set_filtered_genotypes_to_no_call {
        compound.push((
            "INFO",
            "AC",
            htsjdk_vcf::header::Cardinality::A,
            htsjdk_vcf::header::LineType::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        ));
        compound.push((
            "INFO",
            "AF",
            htsjdk_vcf::header::Cardinality::A,
            htsjdk_vcf::header::LineType::Float,
            "Allele Frequency, for each ALT allele, in the same order as listed",
        ));
        compound.push((
            "INFO",
            "AN",
            htsjdk_vcf::header::Cardinality::Fixed(1),
            htsjdk_vcf::header::LineType::Integer,
            "Total number of alleles in called genotypes",
        ));
    }
    if flag(parser, "apply-allele-specific-filters") {
        compound.push((
            "INFO",
            "AS_FilterStatus",
            htsjdk_vcf::header::Cardinality::A,
            htsjdk_vcf::header::LineType::String,
            "Filter status for each allele, as assessed by ApplyVQSR. Note that the VCF filter \
             field will reflect the most lenient/sensitive status across all alleles.",
        ));
    }
    // The header is a `HashSet` of LINES: an identical line collapses and a line with the same ID
    // and a different description does NOT, so a file already declaring `AF` ends up with two
    // `##INFO=<ID=AF,...>` lines. Measured on a row of this tool's array.
    for (kind, id, number, line_type, description) in compound {
        let declared = file.header.lines.iter().any(|line| {
            matches!(line, htsjdk_vcf::header::HeaderLine::Compound { key, id: found, description: text, .. }
                if key == kind && found == id && text == description)
        });
        if !declared {
            file.header
                .lines
                .push(htsjdk_vcf::header::HeaderLine::Compound {
                    key: kind.to_string(),
                    id: id.to_string(),
                    number,
                    line_type,
                    description: description.to_string(),
                    extra: Vec::new(),
                });
        }
    }

    for (name, description) in filter_lines {
        let declared = file.header.lines.iter().any(
            |line| matches!(line, htsjdk_vcf::header::HeaderLine::Filter { id, .. } if id == &name),
        );
        if !declared {
            file.header
                .lines
                .push(htsjdk_vcf::header::HeaderLine::Filter {
                    id: name,
                    description,
                });
        }
    }

    let keep = variant_output_filter(parser, intervals.as_deref())?;
    written.retain(|record| keep(record));
    apply_sites_only(parser, &mut file.header, &mut written);
    let rendered =
        htsjdk_vcf::vcf_file::write_vcf(&file.header, &written).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("{error:?}")),
        })?;
    write_variant_output(parser, &output, &rendered)?;
    Ok(None)
}

/// A record as the filtering reads it: the INFO map, the genotype fields and the FILTER column.
fn filtration_record(
    record: &htsjdk_vcf::variant::VariantContext,
) -> gatk_tools::variant_filtration::Record {
    use gatk_tools::variant_filtration as filtration;
    // `VariantJEXLContext`'s fixed names come FIRST and the INFO attributes after them, because an
    // expression says `QUAL > 50` and no INFO field is called QUAL. A context of attributes alone
    // refuses the expression as an unknown variable, which reads as "the filter did not match":
    // measured on rows of this tool's array, where the reference applied the site filter and the
    // port applied only the mask.
    // The values are the OBJECTS the reference's map hands JEXL rather than their text: `POS` is an
    // `Integer`, `QUAL` a `Double` and `N_ALLELES` an `Integer`, which is what makes `QUAL > 50` an
    // answer rather than a `NumberFormatException` (#1142). An INFO attribute is the `String` the
    // codec decoded it to, and `FILTER` is text on both of its branches.
    use gatk_engine::jexl::Value as JexlValue;
    let mut info: filtration::Context = filtration::Context::new();
    info.insert("CHROM".to_string(), JexlValue::Str(record.contig.clone()));
    info.insert("POS".to_string(), JexlValue::Int(record.start as i32));
    info.insert(
        "QUAL".to_string(),
        JexlValue::Double(-10.0 * record.log10_p_error),
    );
    info.insert(
        "N_ALLELES".to_string(),
        JexlValue::Int(record.alleles.len() as i32),
    );
    let filtered = record
        .filters
        .as_ref()
        .is_some_and(|filters| !filters.is_empty());
    info.insert(
        "FILTER".to_string(),
        JexlValue::Str(if filtered { "1" } else { "0" }.to_string()),
    );
    for (key, value) in &record.attributes {
        if let Some(text) = value.format() {
            info.insert(key.clone(), JexlValue::Str(text));
        }
    }
    for filter in record.filters.iter().flatten() {
        info.entry(filter.clone())
            .or_insert_with(|| JexlValue::Str("1".to_string()));
    }
    let genotypes = record
        .genotypes
        .iter()
        .map(|genotype| filtration::GenotypeFields {
            fields: genotype
                .extended
                .iter()
                .filter_map(|(key, value)| {
                    value
                        .format()
                        .map(|text| (key.clone(), JexlValue::Str(text)))
                })
                .collect(),
            filters: genotype
                .filters
                .as_ref()
                .filter(|text| *text != "PASS" && *text != ".")
                .map(|text| text.split(';').map(str::to_string).collect())
                .unwrap_or_default(),
        })
        .collect();
    filtration::Record {
        contig: record.contig.clone(),
        start: record.start as i32,
        stop: record.stop as i32,
        is_snp: gatk_tools::remove_nearby_indels::variant_type(record)
            == gatk_tools::remove_nearby_indels::VariantType::Snp,
        filters: record.filters.clone(),
        info,
        genotypes,
    }
}

/// `CollectAllelicCounts.apply`: the reference and the alternate base counted at every locus.
///
/// The first copy-number runner here that reads READS. Three things are the port's and each is a
/// decision the name does not show: the alternate count is `total - reference` rather than the
/// alternate BASE's count, the minimum base quality is the collector's own threshold rather than a
/// read filter, and a locus whose reference base is not one of `ACGT` produces no row at all
/// instead of an empty one.
pub fn collect_allelic_counts(parser: &Parser) -> Outcome {
    use gatk_tools::collect_allelic_counts as allelic;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CollectAllelicCounts")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    // `ReadUtils.getSamplesFromHeader(...)`: the SM of the first read group, which is the column
    // name the table carries.
    let sample = header
        .read_groups
        .iter()
        .find_map(|group| group.attributes.get("SM").map(str::to_string))
        .unwrap_or_default();

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;
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
            // `emitEmptyLoci()` is TRUE on this tool, which is what makes the table a row per
            // POSITION of the interval rather than a row per pileup: a locus no read reaches is
            // written with zero counts and an `N` alternate. Measured on this tool's array, where
            // the reference wrote 5,505 rows over `chr1:501-6000` and the port wrote seventy-five.
            emit_empty_loci: true,
            ..gatk_tools::locus_walker::Options::default()
        },
        &filter,
    )
    .map_err(locus_traversal_error)?;

    // The reference's own dictionary, which is what its refusal prints.
    let reference_sequences = gatk_tools::reference_walker::dictionary(&reference).sequences;
    let minimum_base_quality = number_or(parser, "minimum-base-quality", 20).max(0) as u8;
    let mut counts = Vec::new();
    for one in &applied {
        // `referenceContext.getBase()`: the one base under the locus, upper-cased and with its
        // IUPAC codes flattened like every other reference query in this engine.
        let bases = reference
            .query(
                &one.context.contig,
                one.context.position,
                one.context.position,
            )
            .map_err(|error| match error {
                // `MissingContigInSequenceDictionary`, raised by the reference query itself, and
                // the dictionary it prints is the REFERENCE's. Measured on a row of this tool's
                // array where the reads are on `chr1` and `--reference` carries `chrOther`.
                gatk_engine::reference::ReferenceError::UnknownContig(contig) => {
                    Thrown::user(format!(
                        "Contig {contig} not present in the sequence dictionary {}\n",
                        gatk_tools::sequence_dictionary::pretty_print(&reference_sequences)
                    ))
                }
                other => Thrown::user(format!("{other:?}")),
            })?;
        let Some(base) = bases.first().copied() else {
            continue;
        };
        if let Some(count) = allelic::collect_at_locus(
            base,
            &one.context.pileup,
            &one.context.contig,
            one.context.position,
            minimum_base_quality,
        ) {
            counts.push(count);
        }
    }

    let sequences: Vec<(String, i32)> = header
        .sequences
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();
    write_file(
        &output,
        allelic::write(&sequences, &sample, &counts).as_bytes(),
    )?;
    Ok(None)
}

/// One copy-number table as its reader sees it: the metadata, the column names, and the rows.
///
/// The header is kept rather than skipped because it IS the collection's sequence dictionary, and
/// the columns are kept by name because the reader asks for them by name.
struct Table {
    header: SamHeader,
    columns: Vec<String>,
    body: Vec<Vec<String>>,
}

/// `FilterIntervals.doWork`: three filters over one shared mask, and the order is the arithmetic.
///
/// The second CHAIN here: the annotated intervals are what `AnnotateIntervals` writes and the
/// counts are what `CollectReadCounts` writes, so the corpus carries both produced by the
/// reference itself. The filtering is the port's, including the last thing it does and the easiest
/// to miss: a contig left with a single surviving interval loses it, so a run that filters down to
/// exactly one ends with none and is then refused for having nothing left.
pub fn filter_intervals(parser: &Parser) -> Outcome {
    use gatk_tools::filter_intervals as filtering;

    let _ = resolve_read_filters(parser, "FilterIntervals")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // `CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection`, which this tool makes
    // like every other copy-number tool: it bins and pads by its own arguments and refuses the
    // standard interval ones that would modify its input first. Measured on every row of this
    // tool's array, where the reference refused at exit three and the port filtered.
    validate_copy_number_intervals(parser)?;
    let annotated = argument(parser, "annotated-intervals");
    let counts_paths = arguments(parser, "input");

    // A table of this shape is a SAM header, a column line and the rows; the columns are named and
    // the reader asks for them by name, which is what `CollectAllelicCounts`' reader does too.
    //
    // `AbstractRecordCollection`'s constructor calls `IOUtils.canReadFile` before it opens
    // anything, so a path that is not a file is refused by the CHECK rather than by the read, and
    // the message names the file's absolute path. Measured on a row of this tool's array that
    // passed `--annotated-intervals=` with no value, where the reference printed the working
    // directory: "Couldn't read file /work. Error was: The input file does not exist."
    let read_table = |path: &str| -> Result<Table, Thrown> {
        can_read_file(path)?;
        let text = std::fs::read_to_string(path)
            .map_err(|error| Thrown::user(format!("{path}: {error}")))?;
        // The `@SQ` block at the head of the table is the collection's METADATA, which is the
        // dictionary `-L` resolves against and the dictionary the output interval list carries.
        // It is the file's own and not the reference's: this tool takes no `--reference`, and a
        // run that resolved `-L` against an empty dictionary refused every interval as an unknown
        // contig.
        let header = htsjdk_bam::reader::parse_header_text(
            &text
                .lines()
                .take_while(|line| line.starts_with('@'))
                .map(|line| format!("{line}\n"))
                .collect::<String>(),
        );
        let mut rows = text
            .lines()
            .filter(|line| !line.starts_with('@') && !line.trim().is_empty());
        let columns: Vec<String> = rows
            .next()
            .unwrap_or_default()
            .split('\t')
            .map(str::to_string)
            .collect();
        let body = rows
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect();
        Ok(Table {
            header,
            columns,
            body,
        })
    };

    // `validateArguments` runs before a single input is opened, and its order is observable: the
    // missing-inputs refusal is a `UserException` at two, and the duplicate check under it is a
    // `Utils.validateArg` at three.
    if annotated.is_none() && counts_paths.is_empty() {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(filtering::FilterError::NoInputs.message().to_string()),
        });
    }
    let mut seen = std::collections::HashSet::new();
    if !counts_paths.iter().all(|path| seen.insert(path.clone())) {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "List of input read-count files cannot contain duplicates.",
        ));
    }

    let mut intervals: Vec<filtering::Interval> = Vec::new();
    let mut annotations: Vec<(String, Vec<f64>)> = Vec::new();
    // `metadata`, which is the FIRST input's dictionary: the annotated intervals when they were
    // given, and the first counts file otherwise.
    let mut metadata: Option<SamHeader> = None;
    if let Some(path) = &annotated {
        let Table {
            header,
            columns,
            body: rows,
        } = read_table(path)?;
        metadata = Some(header);
        for row in &rows {
            intervals.push(filtering::Interval {
                contig: row[0].clone(),
                start: row[1].parse().unwrap_or_default(),
                end: row[2].parse().unwrap_or_default(),
            });
        }
        for (index, name) in columns.iter().enumerate().skip(3) {
            annotations.push((
                name.clone(),
                rows.iter()
                    .map(|row| {
                        row.get(index)
                            .and_then(|v| v.parse().ok())
                            .unwrap_or(f64::NAN)
                    })
                    .collect(),
            ));
        }
    }

    // Each counts file keeps its OWN intervals, because the matrix is built by looking each
    // interval up per file rather than by position: a file whose rows are a different set is
    // refused, and a file whose rows are a superset is subset down.
    let mut counts_tables: Vec<(String, Vec<filtering::Interval>, Vec<f64>)> = Vec::new();
    for path in &counts_paths {
        let Table {
            header, body: rows, ..
        } = read_table(path)?;
        if metadata.is_none() {
            metadata = Some(header);
        }
        let file_intervals: Vec<filtering::Interval> = rows
            .iter()
            .map(|row| filtering::Interval {
                contig: row[0].clone(),
                start: row[1].parse().unwrap_or_default(),
                end: row[2].parse().unwrap_or_default(),
            })
            .collect();
        if intervals.is_empty() && annotated.is_none() {
            intervals = file_intervals.clone();
        }
        counts_tables.push((
            path.clone(),
            file_intervals,
            rows.iter()
                .map(|row| row.get(3).and_then(|v| v.parse().ok()).unwrap_or(f64::NAN))
                .collect(),
        ));
    }

    // The three refusals are not the same KIND: the empty intersection is an
    // `IllegalArgumentException` from `Utils.validateArg`, which leaves exit three, and the other
    // two are `UserException`s at two. The port's own `java_class` says which is which.
    let refuse = |error: filtering::FilterError| {
        let class = error.java_class();
        let message = error.message().to_string();
        if class == "java.lang.IllegalArgumentException" {
            Thrown::non_user(class, message)
        } else {
            Thrown {
                failure: Failure::User,
                exception: class,
                message: Some(message),
            }
        }
    };
    if intervals.is_empty() {
        return Err(refuse(filtering::FilterError::EmptyIntersection));
    }

    let dictionary = metadata.unwrap_or_default();
    // The dictionary goes into the output list whole: `IntervalList` is built from the metadata's
    // own `SAMSequenceDictionary`, so an `M5` the annotated intervals carried is written out with
    // it. Measured on this tool's array, where the only differing byte was that field.
    let sequences: Vec<gatk_tools::preprocess_intervals::Sequence> = dictionary
        .sequences
        .iter()
        .map(|sequence| gatk_tools::preprocess_intervals::Sequence {
            name: sequence.name.clone(),
            length: sequence.length,
            md5: sequence.attributes.get("M5").map(str::to_string),
            uri: sequence.attributes.get("UR").map(str::to_string),
        })
        .collect();

    // `ListUtils.intersection`, which is a LIST intersection and not a genomic one: an input bin
    // survives only when a requested interval EQUALS it, so `-L chr1:1-7000` over four
    // thousand-base bins intersects to nothing and the run is refused. Measured on every row of
    // this tool's array, where the port read the window as containment and kept all four.
    //
    // With both kinds of input the intersection is taken TWICE, the annotated intervals first and
    // the first counts file second, so a bin either side is missing is gone before any filter runs.
    let requested = interval_arguments(parser, &dictionary)?.map(|parameters| parameters.intervals);
    let kept: Vec<usize> = (0..intervals.len())
        .filter(|&index| {
            let interval = &intervals[index];
            let asked = requested.as_ref().is_none_or(|windows| {
                windows.iter().any(|window| {
                    window.contig == interval.contig
                        && window.start == interval.start
                        && window.end == interval.end
                })
            });
            let counted = annotated.is_none()
                || counts_tables
                    .first()
                    .is_none_or(|(_, file, _)| file.contains(interval));
            asked && counted
        })
        .collect();
    if kept.is_empty() {
        return Err(refuse(filtering::FilterError::EmptyIntersection));
    }
    intervals = kept.iter().map(|&index| intervals[index].clone()).collect();
    for (_, values) in annotations.iter_mut() {
        *values = kept.iter().map(|&index| values[index]).collect();
    }

    // `constructReadCountMatrix`: one row per file, each pulled out by INTERVAL and not by row
    // number, and a file that does not carry every intersected interval is refused at three.
    let mut counts: Vec<Vec<f64>> = Vec::new();
    for (path, file_intervals, values) in &counts_tables {
        let row: Vec<f64> = intervals
            .iter()
            .filter_map(|interval| {
                file_intervals
                    .iter()
                    .position(|candidate| candidate == interval)
                    .map(|index| values[index])
            })
            .collect();
        if row.len() != intervals.len() {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!(
                    "Intervals for read-count file {path} do not contain all specified intervals."
                ),
            ));
        }
        counts.push(row);
    }

    let mut mask = vec![false; intervals.len()];
    let bound = |name: &str, default: f64| -> f64 {
        scalar(parser, name)
            .and_then(|text| text.parse().ok())
            .unwrap_or(default)
    };
    // The annotations run first and in the header's own order, because the mask is SHARED: an
    // interval a GC filter failed is not in the population the count percentiles are taken over.
    for (name, values) in &annotations {
        let (minimum, maximum) = match name.as_str() {
            "GC_CONTENT" => (
                bound("minimum-gc-content", 0.1),
                bound("maximum-gc-content", 0.9),
            ),
            "MAPPABILITY" => (
                bound("minimum-mappability", 0.9),
                bound("maximum-mappability", 1.0),
            ),
            "SEGMENTAL_DUPLICATION_CONTENT" => (
                bound("minimum-segmental-duplication-content", 0.0),
                bound("maximum-segmental-duplication-content", 0.5),
            ),
            _ => continue,
        };
        filtering::update_mask_by_annotation(&mut mask, values, minimum, maximum)
            .map_err(refuse)?;
    }
    if !counts.is_empty() {
        filtering::update_mask_by_low_counts(
            &mut mask,
            &counts,
            number_or(parser, "low-count-filter-count-threshold", 5),
            bound("low-count-filter-percentage-of-samples", 90.0),
        )
        .map_err(refuse)?;
        filtering::update_mask_by_extreme_counts(
            &mut mask,
            &counts,
            bound("extreme-count-filter-minimum-percentile", 1.0),
            bound("extreme-count-filter-maximum-percentile", 99.0),
            bound("extreme-count-filter-percentage-of-samples", 90.0),
        )
        .map_err(refuse)?;
    }
    let contigs: Vec<String> = intervals
        .iter()
        .map(|interval| interval.contig.clone())
        .collect();
    filtering::update_mask_by_solitary_intervals(&mut mask, &contigs).map_err(refuse)?;

    let kept: Vec<filtering::Interval> = intervals
        .iter()
        .enumerate()
        .filter(|(index, _)| !mask[*index])
        .map(|(_, interval)| interval.clone())
        .collect();
    write_file(&output, filtering::write(&sequences, &kept).as_bytes())?;
    Ok(None)
}

/// `LeftAlignAndTrimVariants.apply`: each record split, trimmed, and walked left.
///
/// The tool three engine bricks were built for, and the runner's share of it is small: what is
/// its own is the WINDOW, and the window is a function of the record before it. `lastVariant` is
/// the record as WRITTEN, which has already moved left, so aligning one record frees the next; and
/// a record skipped for being too long still becomes the bound the next one is measured against.
/// Both live in [`gatk_tools::left_align_and_trim_variants`], where a golden measures them.
///
/// The header is the input's, merged with itself, plus the chromosome-count lines when the run
/// splits: `addChromosomeCountsToHeader` is called on that branch alone, so a run that only
/// trims writes the lines the input had.
pub fn left_align_and_trim_variants(parser: &Parser) -> Outcome {
    use gatk_tools::left_align_and_trim_variants as align;

    let VariantWalkerStart {
        input,
        text,
        codec: _,
        intervals,
    } = variant_walker_startup(parser, "LeftAlignAndTrimVariants")?;

    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;

    let arguments = align::Arguments {
        dont_trim_alleles: flag(parser, "dont-trim-alleles"),
        split_multiallelics: flag(parser, "split-multi-allelics"),
        max_indel_size: number_or(parser, "max-indel-length", align::DEFAULT_MAX_INDEL_SIZE),
        max_leading_bases: number_or(
            parser,
            "max-leading-bases",
            align::DEFAULT_MAX_LEADING_BASES,
        ),
    };
    let keep_original_counts = flag(parser, "keep-original-ac");
    // A sharded output is a WRITER, not a transformation: the reference writes
    // `out.shard_00000.vcf.gz` and one file per shard, which is a feature of its own rather than
    // this tool's. Refusing says so; writing one file would be a different answer.
    if number_or(parser, "max-variants-per-shard", 0) > 0 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--max-variants-per-shard writes one file per shard, which this port does not write \
             yet. This message is the port's own and not GATK's.",
        ));
    }

    // `createVCFHeaderLineList`: the input's own lines merged with themselves, the tool's default
    // lines, and the count lines the splitting adds.
    let (mut lines, _warnings) = htsjdk_vcf::merge::smart_merge_headers(
        &[htsjdk_vcf::merge::Source {
            header: &file.header,
            version: None,
        }],
        true,
    )
    .unwrap_or_else(|_| (file.header.lines.clone(), Vec::new()));
    if flag(parser, "add-output-vcf-command-line") {
        lines.push(htsjdk_vcf::header::HeaderLine::Unstructured {
            key: "source".to_string(),
            value: "LeftAlignAndTrimVariants".to_string(),
        });
        if let Some(command_line) = command_line_header_line(parser, "LeftAlignAndTrimVariants") {
            lines.push(command_line);
        }
    }
    if keep_original_counts {
        // `GATKVCFHeaderLines.getInfoLine` for the three keys the subsetter can set, which the
        // header declares whether or not a record turns out to carry them.
        for (id, number, line_type, description) in [
            (
                "AC_Orig",
                htsjdk_vcf::header::Cardinality::A,
                htsjdk_vcf::header::LineType::Integer,
                "Original AC",
            ),
            (
                "AF_Orig",
                htsjdk_vcf::header::Cardinality::A,
                htsjdk_vcf::header::LineType::Float,
                "Original AF",
            ),
            (
                "AN_Orig",
                htsjdk_vcf::header::Cardinality::Fixed(1),
                htsjdk_vcf::header::LineType::Integer,
                "Original AN",
            ),
        ] {
            lines.push(htsjdk_vcf::header::HeaderLine::Compound {
                key: "INFO".to_string(),
                id: id.to_string(),
                number,
                line_type,
                description: description.to_string(),
                extra: Vec::new(),
            });
        }
    }
    if arguments.split_multiallelics {
        for key in ["AC", "AF", "AN"] {
            lines.retain(|line| {
                !matches!(
                    line,
                    htsjdk_vcf::header::HeaderLine::Compound { key: k, id, .. }
                        if k == "INFO" && id == key
                )
            });
            if let Some(standard) = htsjdk_vcf::standard_header_lines::standard_info_line(key) {
                lines.push(standard);
            }
        }
    }
    // `VcfUtils.getSortedSampleSet`, which is the input's samples in their natural order.
    //
    // `--sites-only-vcf-output` empties the set here as well as on every record: the writer is
    // built with no samples at all, so the output has no FORMAT column rather than a no-call in
    // one. Measured on eight rows of this tool's array, where the port wrote `GT ./.`.
    let sites_only = flag(parser, "sites-only-vcf-output");
    let mut samples = if sites_only {
        Vec::new()
    } else {
        file.header.samples.clone()
    };
    samples.sort();
    let header =
        update_header_contig_lines(parser, htsjdk_vcf::header::VcfHeader { lines, samples })?;

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

    // The whole contig, once per contig: the alignment indexes the reference by one-based position
    // and walks LEFT from the record, so a slice around the record would have to be re-based.
    let mut contigs: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    let mut last: Option<gatk_engine::variant_context_utils::Variant> = None;
    for located in gatk_engine::variant_source::traverse(&located, intervals.as_deref()) {
        let original = &file.records[located.index];
        let bridged = crate::variant_bridge::to_engine(original);

        // The counts the subsetter would carry over are only reachable on a record that HAS them;
        // the corpus has none, so the flag is the header's business alone. A record that carries
        // one and is split would need the allele subsetting this port does not have.
        if keep_original_counts
            && arguments.split_multiallelics
            && bridged
                .record
                .variant
                .attributes
                .iter()
                .any(|(key, _)| key == "AC" || key == "AF" || key == "AN")
        {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                "--keep-original-ac over a record that carries AC asks the subsetter for the \
                 counts it had, which this port does not carry yet. This message is the port's \
                 own and not GATK's.",
            ));
        }

        let pieces = if arguments.split_multiallelics {
            if original
                .genotypes
                .iter()
                .any(|genotype| genotype.extended.iter().any(|(key, _)| key == "AF"))
            {
                return Err(Thrown::non_user(
                    PORT_LIMITATION,
                    "A record whose genotypes carry AF is split by the SOMATIC splitter, which \
                     this port does not carry yet. This message is the port's own and not GATK's.",
                ));
            }
            gatk_engine::variant_context_utils::split_variant_context_to_biallelics(
                &bridged.record.variant,
                false,
            )
            .map_err(split_refusal)?
        } else {
            vec![bridged.record.variant.clone()]
        };

        for piece in pieces {
            if align::largest_indel_length(&piece) > arguments.max_indel_size {
                // Written untouched, and it is still the record the next one is measured against.
                written.push(crate::variant_bridge::from_engine(
                    original,
                    &gatk_tools::select_variants::Record {
                        variant: piece.clone(),
                        samples: file.header.samples.clone(),
                    },
                ));
                last = Some(piece);
                continue;
            }
            let distance = match &last {
                Some(previous) if previous.contig == piece.contig => piece.start - previous.stop,
                _ => i32::MAX,
            };
            let window = arguments.max_leading_bases.min(distance - 1);
            // `leftAlignAndTrim` returns before it reads a base when the record is not an indel or
            // the window is empty, so the reference is opened only where the alignment needs it.
            // Measured on a row whose `--reference` carries another contig entirely and whose
            // records are all SNVs or too long to align: the reference wrote the file.
            let needs_reference = piece.is_indel() && window > 0;
            if needs_reference && !contigs.contains_key(&piece.contig) {
                let length = reference
                    .sequences()
                    .iter()
                    .find(|(name, _)| *name == piece.contig)
                    .map(|(_, length)| *length as i32)
                    .unwrap_or(0);
                // `ReferenceContext` over a contig the FASTA does not carry is the walker's own
                // refusal, and it is the one `reference_traversal_error` renders. Measured on a
                // row of this tool's array where `--reference` is the corpus's other contig.
                if length == 0 {
                    return Err(Thrown::user(format!(
                        "Given reference file does not have data at the requested contig({})!",
                        piece.contig
                    )));
                }
                // ZERO-based: `leftAlignAndTrim` slices `[start - 1 .. stop]`, so the vector is
                // the contig's bases with position one at index zero. A one-based vector with a
                // pad in front shifted every comparison by a base, and the records moved one base
                // left instead of the thousand the window allowed.
                let mut bases: Vec<u8> = Vec::new();
                if length > 0 {
                    // The engine's own query, which upper-cases and flattens IUPAC exactly as
                    // `ReferenceDataSource.of(path)` does for a `ReferenceContext`: the alignment
                    // compares these bases against the record's alleles.
                    bases.extend(
                        reference
                            .query(&piece.contig, 1, length)
                            .map_err(|error| match error {
                                gatk_engine::reference::ReferenceError::UnknownContig(contig) => {
                                    Thrown::user(format!(
                                        "Contig {contig} not present in the sequence dictionary {}\n",
                                        gatk_tools::sequence_dictionary::pretty_print(
                                            &gatk_tools::reference_walker::dictionary(&reference)
                                                .sequences
                                        )
                                    ))
                                }
                                other => Thrown::user(format!("{other:?}")),
                            })?,
                    );
                }
                contigs.insert(piece.contig.clone(), bases);
            }
            let empty: Vec<u8> = Vec::new();
            let bases = contigs.get(&piece.contig).unwrap_or(&empty);
            let aligned = gatk_engine::variant_context_utils::left_align_and_trim(
                &piece,
                bases,
                window,
                !arguments.dont_trim_alleles,
            )
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
            written.push(crate::variant_bridge::from_engine(
                original,
                &gatk_tools::select_variants::Record {
                    variant: aligned.clone(),
                    samples: file.header.samples.clone(),
                },
            ));
            last = Some(aligned);
        }
    }

    let keep = variant_output_filter(parser, intervals.as_deref())?;
    written.retain(|record| keep(record));
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
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// What the biallelic splitter refuses: the trimming underneath, or the genotype combinatorics.
///
/// Both are `GATKException`s in the reference rather than user errors, which is what makes them
/// status three.
fn split_refusal(error: gatk_engine::variant_context_utils::SplitError) -> Thrown {
    Thrown::non_user(
        "org.broadinstitute.hellbender.exceptions.GATKException",
        error.message(),
    )
}

/// `GatherBQSRReports.doWork`: several recalibration reports summed into one.
///
/// A `CommandLineProgram` rather than a `GATKTool`, so there is no reference, no intervals and no
/// read filters to resolve: the whole run is read the files, gather, write. What the gather does is
/// [`gatk_tools::gather_bqsr_reports`], where a golden measures it; the runner's share is the
/// order the files are read in, which is the order `--input` was given.
pub fn gather_bqsr_reports(parser: &Parser) -> Outcome {
    use gatk_tools::gather_bqsr_reports as gather;

    let inputs = arguments(parser, "input");
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let mut texts = Vec::new();
    for path in &inputs {
        texts.push(std::fs::read_to_string(path).map_err(|error| Thrown {
            failure: Failure::User,
            exception: gatk_tools::read_walker_refusal::COULD_NOT_READ,
            message: Some(format!(
                "Couldn't read file {}. Error was: {error}",
                java_absolute_path(path)
            )),
        })?);
    }
    let borrowed: Vec<&str> = texts.iter().map(String::as_str).collect();
    // `GATKException` is not a user error, which is what puts it at status three; the empty input
    // list is an `IllegalArgumentException`, which is three as well.
    let gathered = gather::gather(&borrowed)
        .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?;
    write_file(&output, gathered.as_bytes())?;
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("0".to_string()))
}

/// One record of either side, as the concordance iterator reads it.
struct ConcordanceLocus {
    index: usize,
    contig: String,
    start: i32,
    filtered: bool,
}

impl gatk_engine::concordance_walker::ConcordanceRecord for ConcordanceLocus {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn is_filtered(&self) -> bool {
        self.filtered
    }
}

/// `Concordance.apply`: a truth callset and an evaluation one walked together.
///
/// The first `AbstractConcordanceWalker` with a runner here, and what makes it one is that it
/// drives TWO feature inputs at once: the iterator in [`gatk_engine::concordance_walker`] steps
/// them by the dictionary's contig order, and every step is labelled with one of five states. The
/// table those states become is [`gatk_tools::concordance`], where a golden measures it; the
/// runner reads the files, applies the truth side's filter, and writes what the labels ask for.
///
/// The truth side has a filter of its own (`makeTruthVariantFilter`: not filtered, not symbolic)
/// and the eval side has none, which is why a filtered eval record is a STATE rather than a
/// dropped record.
pub fn concordance(parser: &Parser) -> Outcome {
    use gatk_tools::concordance as conc;

    let _ = resolve_read_filters(parser, "Concordance")?;
    let truth_path = argument(parser, "truth").ok_or_else(|| {
        Thrown::command_line("Argument truth was missing: Argument 'truth' is required")
    })?;
    // `--evaluation` is the long name and `-eval` the short one; the parser knows only the long
    // one, and a runner that asked for the short one read nothing at all.
    let eval_path = argument(parser, "evaluation").ok_or_else(|| {
        Thrown::command_line("Argument evaluation was missing: Argument 'evaluation' is required")
    })?;
    let summary_path = argument(parser, "summary").ok_or_else(|| {
        Thrown::command_line("Argument summary was missing: Argument 'summary' is required")
    })?;

    let read_vcf = |path: &str| -> Result<(String, htsjdk_vcf::reader::VcfFile), Thrown> {
        let bytes = std::fs::read(path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: path.to_string(),
                }
                .message(),
            )
        })?;
        let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
            htsjdk_bgzf::read::decompress_all(&bytes)
                .ok()
                .and_then(|bytes| String::from_utf8(bytes).ok())
                .ok_or_else(|| {
                    Thrown::non_user(
                        gatk_tools::read_walker_refusal::SAM_FORMAT,
                        format!("{path} is not a block compressed file"),
                    )
                })?
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })?;
        Ok((text, file))
    };
    let (truth_text, truth_file) = read_vcf(&truth_path)?;
    let (_, eval_file) = read_vcf(&eval_path)?;

    // `validateSequenceDictionaries` still runs, and it is the ENGINE's: the master against the
    // reads, the reference and the features, then the reference against the reads. Only the best
    // available dictionary is this tool's own. Measured on a row whose `--sequence-dictionary`
    // shares no contig with the BAM the command line named: the reference refused and the port
    // walked.
    let reads_inputs = arguments(parser, "input").len();
    let reads_dictionaries: Vec<Vec<htsjdk_bam::header::SequenceRecord>> =
        arguments(parser, "input")
            .iter()
            .map(|reads| reads_dictionary(parser, reads, reads_inputs))
            .collect::<Result<_, _>>()?;
    let features = vcf_dictionary(&truth_text);
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            for reads in &reads_dictionaries {
                validate_against_master(master, "reads", reads)?;
            }
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
            validate_against_master(master, "features", &features.sequences)?;
        }
        let compare = |left_name: &str,
                       left: &[htsjdk_bam::header::SequenceRecord],
                       right_name: &str,
                       right: &[htsjdk_bam::header::SequenceRecord]|
         -> Result<(), Thrown> {
            gatk_tools::sequence_dictionary::validate(
                left_name, left, right_name, right, false, false,
            )
            .map_err(|refusal| Thrown {
                failure: Failure::User,
                exception: refusal.java_class(),
                message: Some(refusal.message()),
            })
        };
        if let Some(reference) = &reference {
            for reads in &reads_dictionaries {
                compare("reference", &reference.sequences, "reads", reads)?;
            }
            compare(
                "reference",
                &reference.sequences,
                "features",
                &features.sequences,
            )?;
        }
        for reads in &reads_dictionaries {
            compare("reads", reads, "features", &features.sequences)?;
        }
    }

    // `AbstractConcordanceWalker.getBestAvailableSequenceDictionary` is FINAL and returns the
    // TRUTH file's dictionary. Not the master, not the reference: a run whose
    // `--sequence-dictionary` carries another contig entirely still resolves `-L chr1` against the
    // truth file, and the comparator orders by the contig's index in that same dictionary.
    // Measured on three rows of this tool's array, where the port refused the interval the
    // reference walked.
    let sequences: Vec<htsjdk_bam::header::SequenceRecord> =
        vcf_dictionary(&truth_text).sequences.clone();

    // The comparator orders by the contig's INDEX in that dictionary, so only the names matter to
    // it; the lengths matter to `-L`, which validates against them.
    let dictionary: Vec<String> = sequences
        .iter()
        .map(|sequence| sequence.name.clone())
        .collect();

    // `-L` bounds BOTH sides: the base class hands each file the same traversal parameters, so a
    // window that ends before a call drops that call from the comparison entirely. Measured on a
    // row whose window is chr1:1-6000 and whose eval file calls at 6001: the reference counted one
    // false positive fewer.
    let dictionary_header = SamHeader {
        sequences: sequences.clone(),
        ..SamHeader::default()
    };
    let intervals =
        interval_arguments(parser, &dictionary_header)?.map(|parameters| parameters.intervals);
    for (path, file) in [(&truth_path, &truth_file), (&eval_path, &eval_file)] {
        let _ = file;
        if intervals.is_some() && !has_feature_index(path) {
            return Err(Thrown {
                failure: Failure::User,
                exception: "org.broadinstitute.hellbender.exceptions.UserException",
                message: Some(format!(
                    "Input {path} must support random access to enable traversal by intervals. \
                     If it's a file, please index it using the bundled tool IndexFeatureFile"
                )),
            });
        }
    }
    let in_traversal = |contig: &str, start: i32, stop: i32| -> bool {
        match &intervals {
            None => true,
            Some(windows) => windows.iter().any(|window| {
                window.contig == contig && window.start <= stop && start <= window.end
            }),
        }
    };

    let is_symbolic_or_sv = |record: &htsjdk_vcf::variant::VariantContext| {
        record.alleles[1..]
            .iter()
            .any(|allele| allele.is_symbolic() || allele.display_string().starts_with('<'))
    };
    // `makeTruthVariantFilter`, which the base class applies before the iterator sees a record.
    let truth: Vec<ConcordanceLocus> = truth_file
        .records
        .iter()
        .enumerate()
        .filter(|(_, record)| {
            in_traversal(&record.contig, record.start as i32, record.stop as i32)
                && conc::truth_variant_filter(record.is_filtered(), is_symbolic_or_sv(record))
        })
        .map(|(index, record)| ConcordanceLocus {
            index,
            contig: record.contig.clone(),
            start: record.start as i32,
            filtered: record.is_filtered(),
        })
        .collect();
    let eval: Vec<ConcordanceLocus> = eval_file
        .records
        .iter()
        .enumerate()
        .filter(|(_, record)| in_traversal(&record.contig, record.start as i32, record.stop as i32))
        .map(|(index, record)| ConcordanceLocus {
            index,
            contig: record.contig.clone(),
            start: record.start as i32,
            filtered: record.is_filtered(),
        })
        .collect();

    let alleles = |record: &htsjdk_vcf::variant::VariantContext| -> (String, Vec<String>) {
        (
            record
                .alleles
                .first()
                .map(|allele| allele.display_string())
                .unwrap_or_default(),
            record.alleles[1..]
                .iter()
                .map(|allele| allele.display_string())
                .collect(),
        )
    };
    let steps =
        gatk_engine::concordance_walker::concordance(&truth, &eval, &dictionary, |left, right| {
            let (truth_reference, truth_alternates) = alleles(&truth_file.records[left.index]);
            let (eval_reference, eval_alternates) = alleles(&eval_file.records[right.index]);
            conc::variants_at_same_locus_are_concordant(
                &truth_reference,
                &truth_alternates,
                &eval_reference,
                &eval_alternates,
            )
        });

    // One record per FILTER line of the EVAL header, whatever any record carries.
    //
    // A `##FILTER` line carries an ID and a description and no Number or Type, so the header model
    // holds it as a STRUCTURED line rather than a compound one: reading only the compound ones
    // found no filters at all, and the first filtered record then dereferenced a record that was
    // never created.
    let declared: Vec<String> = eval_file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Compound { key, id, .. } if key == "FILTER" => {
                Some(id.clone())
            }
            htsjdk_vcf::header::HeaderLine::Structured { key, fields } if key == "FILTER" => fields
                .iter()
                .find(|(name, _)| name == "ID")
                .map(|(_, value)| value.clone()),
            _ => None,
        })
        .collect();
    let mut analysis = conc::FilterAnalysis::new(&declared);
    let filter_analysis_path = argument(parser, "filter-analysis");

    let mut summary = conc::Summary::default();
    let mut annotated: std::collections::HashMap<
        &'static str,
        Vec<htsjdk_vcf::variant::VariantContext>,
    > = std::collections::HashMap::new();
    for step in &steps {
        // `getTruthIfPresentElseEval().isSNP()`: the truth record decides the bucket when there is
        // one, so a false positive is bucketed by the EVAL record and nothing else is.
        let representative = match (step.truth, step.eval) {
            (Some(index), _) => &truth_file.records[truth[index].index],
            (None, Some(index)) => &eval_file.records[eval[index].index],
            (None, None) => continue,
        };
        let is_snp = gatk_tools::remove_nearby_indels::variant_type(representative)
            == gatk_tools::remove_nearby_indels::VariantType::Snp;
        summary.add(step.state, is_snp);

        if let Some(index) = step.eval {
            let record = &eval_file.records[eval[index].index];
            analysis
                .apply(
                    step.state,
                    &record.filters.clone().unwrap_or_default(),
                    filter_analysis_path.is_some(),
                )
                .map_err(|error| Thrown::non_user(error.class(), error.message()))?;
        }

        for (file, side) in conc::writes(step.state) {
            let name = match file {
                conc::AnnotatedVcf::TruePositivesAndFalseNegatives => {
                    "true-positives-and-false-negatives"
                }
                conc::AnnotatedVcf::TruePositivesAndFalsePositives => {
                    "true-positives-and-false-positives"
                }
                conc::AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives => {
                    "filtered-true-negatives-and-false-negatives"
                }
            };
            if argument(parser, name).is_none() {
                continue;
            }
            let record = match side {
                conc::Side::Truth => step
                    .truth
                    .map(|index| &truth_file.records[truth[index].index]),
                conc::Side::Eval => step.eval.map(|index| &eval_file.records[eval[index].index]),
            };
            let Some(record) = record else {
                continue;
            };
            // `annotateWithConcordanceState`: the state's abbreviation under `STATUS`, on a copy.
            let mut copy = record.clone();
            copy.attributes
                .retain(|(key, _)| key != conc::TRUTH_STATUS_VCF_ATTRIBUTE);
            copy.attributes.push((
                conc::TRUTH_STATUS_VCF_ATTRIBUTE.to_string(),
                htsjdk_vcf::variant::Value::Str(step.state.abbreviation().to_string()),
            ));
            annotated.entry(name).or_default().push(copy);
        }
    }

    write_file(&summary_path, summary.table().as_bytes())?;
    if let Some(path) = &filter_analysis_path {
        let table = analysis
            .table()
            .map_err(|error| Thrown::non_user(error.class(), error.message()))?;
        write_file(path, table.as_bytes())?;
    }
    for (name, side) in [
        (
            "true-positives-and-false-negatives",
            conc::AnnotatedVcf::TruePositivesAndFalseNegatives,
        ),
        (
            "true-positives-and-false-positives",
            conc::AnnotatedVcf::TruePositivesAndFalsePositives,
        ),
        (
            "filtered-true-negatives-and-false-negatives",
            conc::AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives,
        ),
    ] {
        let Some(path) = argument(parser, name) else {
            continue;
        };
        // The header is the SIDE's own, with the STATUS line and the tool's default lines added.
        let source = match side.header() {
            conc::Side::Truth => &truth_file.header,
            conc::Side::Eval => &eval_file.header,
        };
        let mut header = source.clone();
        header.lines.push(htsjdk_vcf::header::HeaderLine::Compound {
            key: "INFO".to_string(),
            id: conc::TRUTH_STATUS_VCF_ATTRIBUTE.to_string(),
            number: htsjdk_vcf::header::Cardinality::Fixed(1),
            line_type: htsjdk_vcf::header::LineType::String,
            description: "Truth status: TP/FP/FN for true positive/false positive/false negative."
                .to_string(),
            extra: Vec::new(),
        });
        if flag(parser, "add-output-vcf-command-line") {
            header
                .lines
                .push(htsjdk_vcf::header::HeaderLine::Unstructured {
                    key: "source".to_string(),
                    value: "Concordance".to_string(),
                });
            if let Some(command_line) = command_line_header_line(parser, "Concordance") {
                header.lines.push(command_line);
            }
        }
        let mut records = annotated.remove(name).unwrap_or_default();
        // `--sites-only-vcf-output` empties the sample set on the header and the genotypes on every
        // record, the same way it does on every other tool here that writes a VCF.
        if flag(parser, "sites-only-vcf-output") {
            header.samples.clear();
            for record in &mut records {
                record.genotypes.clear();
            }
        }
        let text = htsjdk_vcf::vcf_file::write_vcf(&header, &records).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("{error:?}")),
        })?;
        write_variant_output(parser, &path, &text)?;
    }
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `DepthOfCoverage.apply`: every base of every interval, counted per sample.
///
/// The first tool here that writes a FAMILY of files, and which of them appear is a function of the
/// four `omit` arguments alone: no omission depends on whether the data would have filled the file.
/// This port writes ONE of the seven, the per-base table, and refuses a command line that asks for
/// any of the other six: their quantiles come from the partitioned data store, which is not ported,
/// and a file of plausible numbers would be worse than a refusal.
///
/// What IS ported is what a row of that table holds: every base of the interval is a row whether a
/// read reaches it or not, the partition is the SAMPLE rather than the read group, and
/// `--min-base-quality` filters a BASE rather than a read.
pub fn depth_of_coverage(parser: &Parser) -> Outcome {
    use gatk_tools::depth_of_coverage as coverage;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "DepthOfCoverage")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    // The two arguments whose values decide the FORMAT of what is written, held at their defaults.
    if scalar(parser, "output-format").as_deref().unwrap_or("CSV") != "CSV" {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--output-format TABLE writes the per-base file as a padded table, which this port \
             does not write yet. This message is the port's own and not GATK's.",
        ));
    }
    let partitions = arguments(parser, "partition-type");
    if !partitions.is_empty() && partitions != vec!["sample".to_string()] {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--partition-type beyond `sample` partitions the counts by read group, library or \
             centre as well, which this port does not carry yet. This message is the port's own \
             and not GATK's.",
        ));
    }
    if argument(parser, "calculate-coverage-over-genes").is_some() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--calculate-coverage-over-genes writes a per-gene table this port does not write \
             yet. This message is the port's own and not GATK's.",
        ));
    }

    let omissions = coverage::Omissions {
        locus_table: flag(parser, "omit-locus-table"),
        depth_output_at_each_base: flag(parser, "omit-depth-output-at-each-base"),
        per_sample_statistics: flag(parser, "omit-per-sample-statistics"),
        interval_statistics: flag(parser, "omit-interval-statistics"),
    };
    // The per-base file's suffix is the EMPTY one; every other suffix is a table this port does
    // not compute.
    let unported: Vec<&str> = coverage::written_suffixes(&omissions)
        .into_iter()
        .filter(|suffix| !suffix.is_empty())
        .collect();
    if !unported.is_empty() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "This run asks for {}, whose quantiles come from a partitioned data store this \
                 port does not carry yet. Pass --omit-locus-table, --omit-interval-statistics and \
                 --omit-per-sample-statistics. This message is the port's own and not GATK's.",
                unported.join(", ")
            ),
        ));
    }

    // `calculateCoverageHistogramBinEndpoints`, which runs at startup and refuses the three
    // binning arguments together rather than one at a time.
    let (start, stop, bins) = (
        number_or(parser, "start", 1),
        number_or(parser, "stop", 500),
        number_or(parser, "nBins", 499),
    );
    if bins > stop - start || start < 1 {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            message: Some(
                "Bad input: the start must be at least 1 and the number of bins may not exceed \
                 stop - start"
                    .to_string(),
            ),
        });
    }

    let minimum = number_or(parser, "min-base-quality", 0) as i64;
    let maximum = number_or(parser, "max-base-quality", 127) as i64;
    for (name, value) in [("min-base-quality", minimum), ("max-base-quality", maximum)] {
        // Every one of the three is the PARSER's refusal, which is a `CommandLineException` and
        // therefore status one, whichever side of the byte range the value fell out on.
        coverage::check_base_quality(name, value)
            .map_err(|error| Thrown::command_line(error.message()))?;
    }

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(reads_traversal_error)?;

    // `ReadUtils.getSamplesFromHeader`: the distinct SM values, in the header's own order, which is
    // the order the sample columns are written in.
    let mut samples: Vec<String> = Vec::new();
    for group in &header.read_groups {
        if let Some(sample) = group.attributes.get("SM") {
            if !samples.iter().any(|seen| seen == sample) {
                samples.push(sample.to_string());
            }
        }
    }
    let sample_of = |record: &htsjdk_bam::record::BamRecord| -> String {
        gatk_engine::read_pileup::sample_name(record, &header).unwrap_or_default()
    };
    let reads: Vec<coverage::Read> = records
        .iter()
        .filter_map(|record| {
            contig_name(&header, record.reference_index).map(|contig| coverage::Read {
                name: record.read_name.clone(),
                sample: sample_of(record),
                contig: contig.to_string(),
                start: record.alignment_start,
                bases: record.read_bases.clone(),
                base_qualities: record
                    .base_qualities
                    .iter()
                    .map(|quality| i32::from(*quality))
                    .collect(),
            })
        })
        .collect();

    // The traversal still runs when no file will be written, and two refusals live inside it.
    //
    // The REFERENCE is asked for the bases under each interval, so a `--reference` that does not
    // carry the interval's contig is refused here rather than at startup: measured on a row whose
    // dictionary validation is turned off and whose reference is the corpus's other contig.
    //
    // And it is asked only when `--print-base-counts` is on: the base-count field is the one thing
    // in this table that needs a reference base, so a run without it writes its rows whatever the
    // reference carries and a run with it is refused. Both rows are in this tool's array.
    if flag(parser, "print-base-counts") {
        if let Some(path) = argument(parser, "reference") {
            let mut reference =
                gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&path))
                    .map_err(|error| Thrown::user(format!("{error:?}")))?;
            let sequences = gatk_tools::reference_walker::dictionary(&reference).sequences;
            for interval in &intervals {
                if let Err(gatk_engine::reference::ReferenceError::UnknownContig(contig)) =
                    reference.query(&interval.contig, interval.start, interval.start)
                {
                    return Err(Thrown::user(format!(
                        "Contig {contig} not present in the sequence dictionary {}\n",
                        gatk_tools::sequence_dictionary::pretty_print(&sequences)
                    )));
                }
            }
        }
    }
    // And `CoverageUtils.getBaseCounts` throws on the two fragment modes before it counts a base:
    // the feature is disabled in the reference itself (gatk#6491), so the refusal is the answer.
    // It comes AFTER the binning check and after the reference query, which is what two rows of
    // this tool's array measure.
    if matches!(
        scalar(parser, "count-type").as_deref(),
        Some("COUNT_FRAGMENTS") | Some("COUNT_FRAGMENTS_REQUIRE_SAME_BASE")
    ) {
        return Err(Thrown::non_user(
            "java.lang.UnsupportedOperationException",
            "Fragment based counting is currently unsupported",
        ));
    }
    // `--include-deletions` adds a sixth column, `D`, to every base-count field and counts a
    // deletion into it. The ported counter has the five bases and no deletion slot, so a run that
    // asks for it would write a table one column short of the reference's.
    if flag(parser, "include-deletions") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--include-deletions adds a D column to the base counts, which this port does not \
             count yet. This message is the port's own and not GATK's.",
        ));
    }

    // `--omit-depth-output-at-each-base` removes the per-base file, which is the only file this
    // port writes: the run then writes nothing at all, and so does the reference once the other
    // three omissions have taken their files away.
    if omissions.depth_output_at_each_base {
        return Ok(Some("success".to_string()));
    }

    let print_base_counts = flag(parser, "print-base-counts");
    let mut text = coverage::per_locus_header(&samples, print_base_counts);
    text.push('\n');
    for interval in &intervals {
        for locus in coverage::per_locus(
            &reads,
            &samples,
            &interval.contig,
            interval.start,
            interval.end,
            minimum as i32,
            maximum as i32,
        ) {
            text.push_str(&coverage::per_locus_row(&locus, print_base_counts));
            text.push('\n');
        }
    }
    write_file(&output, text.as_bytes())?;
    // `onTraversalSuccess` returns the word, lower case, which `handleResult` prints.
    Ok(Some("success".to_string()))
}

/// One record as the posteriors read it: the alleles, the counts in INFO, and the genotypes.
fn posterior_record(
    vc: &htsjdk_vcf::variant::VariantContext,
) -> gatk_tools::calculate_genotype_posteriors::Record {
    use gatk_tools::calculate_genotype_posteriors as posteriors;

    let mut attributes = std::collections::BTreeMap::new();
    for key in ["AC", "AN", "MLEAC"] {
        if let Some((_, value)) = vc.attributes.iter().find(|(name, _)| name == key) {
            if let Some(text) = value.format() {
                attributes.insert(key.to_string(), text);
            }
        }
    }
    posteriors::Record {
        id: vc.id.clone(),
        start: vc.start as i32,
        alleles: vc
            .alleles
            .iter()
            .map(|allele| posteriors::Allele {
                bases: allele.display_string(),
                is_ref: allele.is_reference(),
            })
            .collect(),
        attributes,
        genotypes: vc
            .genotypes
            .iter()
            .map(|genotype| posteriors::Genotype {
                sample: genotype.sample_name.clone(),
                alleles: genotype
                    .alleles
                    .iter()
                    .filter_map(|allele| {
                        vc.alleles.iter().position(|candidate| candidate == allele)
                    })
                    .collect(),
                depth: genotype.dp,
                likelihoods: genotype.pl.clone(),
                posteriors: genotype
                    .extended
                    .iter()
                    .find(|(key, _)| key == "PP")
                    .and_then(|(_, value)| value.format())
                    .map(|text| {
                        text.split(',')
                            .filter_map(|piece| piece.trim().parse().ok())
                            .collect()
                    }),
            })
            .collect(),
    }
}

/// `CalculateGenotypePosteriors.apply`: likelihoods turned into posteriors under a prior built
/// from allele counts.
///
/// Where the counts come from is the whole of the tool, and three arguments move it: a site no
/// supporting callset carries falls back to the INPUT's own samples, but only when there are ten of
/// them or `--num-reference-samples-if-no-call` was given, and `--ignore-input-samples` takes even
/// that away. A site with no counts at all gets a FLAT prior, which is `PG` all zeros and `PP`
/// equal to `PL`.
///
/// `calculateChromosomeCounts` runs on every record BEFORE the priors, so the counts the posteriors
/// read are the recomputed ones rather than whatever the file carried.
///
/// The family half is a second algorithm (`FamilyLikelihoods`) that this port does not carry, so a
/// command line with a `--pedigree` is refused rather than answered without it.
pub fn calculate_genotype_posteriors(parser: &Parser) -> Outcome {
    use gatk_tools::calculate_genotype_posteriors as posteriors;

    let VariantWalkerStart {
        input,
        text,
        codec: _,
        intervals,
    } = variant_walker_startup(parser, "CalculateGenotypePosteriors")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if argument(parser, "pedigree").is_some() && !flag(parser, "skip-family-priors") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--pedigree asks for the family priors, a second algorithm this port does not carry \
             yet. Pass --skip-family-priors. This message is the port's own and not GATK's.",
        ));
    }

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;

    // Every supporting callset, read whole: `featureContext.getValues` asks each of them for the
    // records at the driving record's locus, and only those with the SAME START are used.
    let mut supporting: Vec<posteriors::Record> = Vec::new();
    for path in arguments(parser, "supporting-callsets") {
        let bytes = std::fs::read(&path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile { path: path.clone() }.message(),
            )
        })?;
        let support_text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
            htsjdk_bgzf::read::decompress_all(&bytes)
                .ok()
                .and_then(|bytes| String::from_utf8(bytes).ok())
                .unwrap_or_default()
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        let support = htsjdk_vcf::reader::read_vcf(&support_text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })?;
        supporting.extend(support.records.iter().map(posterior_record));
    }

    let number = |name: &str, default: f64| -> f64 {
        scalar(parser, name)
            .and_then(|text| text.parse().ok())
            .unwrap_or(default)
    };
    let options = posteriors::Options {
        snp_prior_dirichlet: number("global-prior-snp", 0.001),
        indel_prior_dirichlet: number("global-prior-indel", 0.001),
        use_input_samples_allele_counts: !flag(parser, "discovered-allele-count-priors-off"),
        use_mleac: !flag(parser, "default-to-allele-count"),
        ignore_input_samples_for_missing_resources: flag(parser, "ignore-input-samples"),
        use_flat_priors_for_indels: flag(parser, "use-flat-priors-for-indels"),
    };
    let missing_reference_samples = number_or(parser, "num-reference-samples-if-no-call", 0);
    let skip_population_priors = flag(parser, "skip-population-priors");

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

    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for located in gatk_engine::variant_source::traverse(&located, intervals.as_deref()) {
        let original = &file.records[located.index];
        // `calculateChromosomeCounts(builder, false)` on every record, whatever the priors do next.
        let mut bridged = crate::variant_bridge::to_engine(original);
        gatk_tools::select_variants::calculate_chromosome_counts(&mut bridged.record.variant);
        let mut vc = crate::variant_bridge::from_engine(original, &bridged.record);
        if skip_population_priors {
            written.push(vc);
            continue;
        }

        let record = posterior_record(&vc);
        let matching: Vec<posteriors::Record> = supporting
            .iter()
            .filter(|resource| resource.start == record.start)
            .cloned()
            .collect();
        let answer = posteriors::calculate_posterior_probs(
            &record,
            &matching,
            if matching.is_empty() {
                missing_reference_samples
            } else {
                0
            },
            &options,
        );

        if let Some(prior) = &answer.prior {
            vc.attributes.retain(|(key, _)| key != "PG");
            vc.attributes.push((
                "PG".to_string(),
                htsjdk_vcf::variant::Value::Str(
                    prior
                        .iter()
                        .map(|value| value.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                ),
            ));
        }
        for called in &answer.genotypes {
            let Some(genotype) = vc
                .genotypes
                .iter_mut()
                .find(|genotype| genotype.sample_name == called.sample)
            else {
                continue;
            };
            genotype.alleles = called
                .alleles
                .iter()
                .filter_map(|index| original.alleles.get(*index).cloned())
                .collect();
            if let Some(gq) = called.gq {
                genotype.gq = Some(gq);
            }
            genotype.extended.retain(|(key, _)| key != "PP");
            if let Some(values) = &called.posteriors {
                genotype.extended.push((
                    "PP".to_string(),
                    htsjdk_vcf::variant::Value::Str(
                        values
                            .iter()
                            .map(|value| value.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    ),
                ));
            }
        }
        written.push(vc);
    }

    // The two header lines the tool adds, plus the three the RECOMPUTATION needs: every record goes
    // through `calculateChromosomeCounts`, so `AC`, `AF` and `AN` are written whether the input
    // declared them or not. Measured on every accepted row of this tool's array, where the port
    // wrote a count into a header that did not declare it.
    let mut header = file.header.clone();
    for key in ["AC", "AF", "AN"] {
        header.lines.retain(|line| {
            !matches!(
                line,
                htsjdk_vcf::header::HeaderLine::Compound { key: k, id, .. }
                    if k == "INFO" && id == key
            )
        });
        if let Some(standard) = htsjdk_vcf::standard_header_lines::standard_info_line(key) {
            header.lines.push(standard);
        }
    }
    for (key, id, number, line_type, description) in [
        (
            "INFO",
            "PG",
            htsjdk_vcf::header::Cardinality::G,
            htsjdk_vcf::header::LineType::Integer,
            "Genotype Likelihood Prior",
        ),
        (
            "FORMAT",
            "PP",
            htsjdk_vcf::header::Cardinality::G,
            htsjdk_vcf::header::LineType::Integer,
            "Phred-scaled Posterior Genotype Probabilities",
        ),
    ] {
        header.lines.push(htsjdk_vcf::header::HeaderLine::Compound {
            key: key.to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: description.to_string(),
            extra: Vec::new(),
        });
    }
    if flag(parser, "add-output-vcf-command-line") {
        header
            .lines
            .push(htsjdk_vcf::header::HeaderLine::Unstructured {
                key: "source".to_string(),
                value: "CalculateGenotypePosteriors".to_string(),
            });
        if let Some(command_line) = command_line_header_line(parser, "CalculateGenotypePosteriors")
        {
            header.lines.push(command_line);
        }
    }
    // And NOT `VcfUtils.updateHeaderContigLines`: this tool writes the input's own contig lines
    // whatever `--reference` says, so a port that rebuilt them added an `assembly=` field and a
    // `##reference` line the reference never writes. Measured on every accepted row here.

    let keep = variant_output_filter(parser, intervals.as_deref())?;
    written.retain(|record| keep(record));
    if flag(parser, "sites-only-vcf-output") {
        header.samples.clear();
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

/// `DenoiseReadCounts.doWork`, in the branch with no panel: the counts standardised.
///
/// The next link of the copy-number chain, and the shortest: read the counts `CollectReadCounts`
/// wrote, take them to fractional coverage, divide by the sample median, take the log, and subtract
/// the median of THAT. Two files come out and they are the same file, because a run with no panel
/// has nothing to denoise with: `writeResult` is handed the standardised values twice.
///
/// The panel itself is HDF5 and the GC correction reads the annotated intervals' GC column; neither
/// is written here, so both arguments are refused rather than ignored.
pub fn denoise_read_counts(parser: &Parser) -> Outcome {
    use gatk_tools::denoise_read_counts as denoise;

    if argument(parser, "count-panel-of-normals").is_some() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--count-panel-of-normals is an HDF5 panel this port does not read yet. This message \
             is the port's own and not GATK's.",
        ));
    }
    if argument(parser, "annotated-intervals").is_some() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--annotated-intervals turns on the GC-bias correction, which this port does not \
             carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let standardized = argument(parser, "standardized-copy-ratios").ok_or_else(|| {
        Thrown::command_line(
            "Argument standardized-copy-ratios was missing: Argument 'standardized-copy-ratios' \
             is required",
        )
    })?;
    let denoised = argument(parser, "denoised-copy-ratios").ok_or_else(|| {
        Thrown::command_line(
            "Argument denoised-copy-ratios was missing: Argument 'denoised-copy-ratios' is required",
        )
    })?;
    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;

    can_read_file(&input)?;
    let text = std::fs::read_to_string(&input)
        .map_err(|error| Thrown::user(format!("{input}: {error}")))?;
    // The `@SQ` block is the metadata, the `@RG`'s `SM` the sample, and the rows the counts.
    let header = htsjdk_bam::reader::parse_header_text(
        &text
            .lines()
            .take_while(|line| line.starts_with('@'))
            .map(|line| format!("{line}\n"))
            .collect::<String>(),
    );
    let sequences: Vec<(String, i32)> = header
        .sequences
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();
    let sample = header
        .read_groups
        .iter()
        .find_map(|group| group.attributes.get("SM").map(str::to_string))
        .unwrap_or_default();

    let mut intervals: Vec<(String, i32, i32)> = Vec::new();
    let mut counts: Vec<f64> = Vec::new();
    for line in text
        .lines()
        .filter(|line| !line.starts_with('@') && !line.trim().is_empty())
        .skip(1)
    {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 4 {
            continue;
        }
        intervals.push((
            columns[0].to_string(),
            columns[1].parse().unwrap_or_default(),
            columns[2].parse().unwrap_or_default(),
        ));
        counts.push(columns[3].parse().unwrap_or(f64::NAN));
    }

    let values = denoise::standardize(&counts).map_err(|error| Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message().to_string()),
    })?;
    let table = denoise::write(&sequences, &sample, &intervals, &values);
    // The same bytes twice: with no panel the denoised result IS the standardised one.
    write_file(&standardized, table.as_bytes())?;
    write_file(&denoised, table.as_bytes())?;
    Ok(None)
}

/// `ValidateBasicSomaticShortMutations.apply`: a discovery call asked of a second pair of BAMs.
///
/// The second tool here that drives TWO samples of reads at once, and unlike `GetNormalArtifactData`
/// it names them: `--val-case-sample-name` and `--val-control-sample-name` pick the pileups out of
/// whatever `--input` holds, and `--discovery-sample-name` picks the genotype out of the VCF.
///
/// Three files can come out and each is optional but the first: the validation table, the annotated
/// VCF and the summary. What decides a judgment is
/// [`gatk_tools::validate_basic_somatic_short_mutations`], where a golden measures it.
pub fn validate_basic_somatic_short_mutations(parser: &Parser) -> Outcome {
    use gatk_tools::validate_basic_somatic_short_mutations as validate;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "ValidateBasicSomaticShortMutations")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }

    let options = validate::Arguments {
        discovery_sample: argument(parser, "discovery-sample-name").unwrap_or_default(),
        validation_case_name: argument(parser, "val-case-sample-name").unwrap_or_default(),
        validation_control_name: argument(parser, "val-control-sample-name").unwrap_or_default(),
        min_power: scalar(parser, "min-power")
            .and_then(|text| text.parse().ok())
            .unwrap_or(validate::DEFAULT_MIN_POWER),
        max_validation_normal_count: number_or(
            parser,
            "max-validation-normal-count",
            validate::DEFAULT_MAX_VALIDATION_NORMAL_COUNT,
        ),
        min_bq_cutoff: number_or(
            parser,
            "min-base-quality-cutoff",
            validate::DEFAULT_MIN_BQ_CUTOFF,
        ),
    };

    // The reads of every `--input`, kept with the header they came from: the sample a read belongs
    // to is its read group's, and the two names above pick which pileup a read lands in.
    let resolved = resolve_read_filters(parser, "ValidateBasicSomaticShortMutations")?;
    let mut reads: Vec<(String, htsjdk_bam::record::BamRecord)> = Vec::new();
    for path in arguments(parser, "input") {
        let source =
            gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let header = source.header().clone();
        let filter = read_filter(parser, &resolved, &header)?;
        for read in gatk_tools::read_walker::traverse(&source, &[], &filter)
            .map_err(|error| Thrown::user(format!("{error:?}")))?
        {
            let sample = gatk_engine::read_pileup::sample_name(&read, &header).unwrap_or_default();
            reads.push((sample, read));
        }
    }

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;
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

    let mut table: Vec<gatk_engine::basic_somatic_short_mutation_validator::BasicValidationResult> =
        Vec::new();
    let mut summary = gatk_tools::concordance::Summary::default();
    let mut annotated: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for located in gatk_engine::variant_source::traverse(&located, intervals.as_deref()) {
        let record = &file.records[located.index];
        let Some(genotype) = record
            .genotypes
            .iter()
            .find(|genotype| genotype.sample_name == options.discovery_sample)
        else {
            continue;
        };
        let allele =
            |allele: &htsjdk_vcf::allele::Allele| gatk_engine::variant_context_utils::Allele {
                bases: allele.display_string().into_bytes(),
                is_reference: allele.is_reference(),
            };
        let validation_genotype =
            gatk_engine::basic_somatic_short_mutation_validator::ValidationGenotype {
                alleles: genotype.alleles.iter().map(allele).collect(),
                ad: genotype.ad.clone(),
                filters: genotype.filters.clone(),
            };
        // The pileup of one named sample at the record's START, which is the locus the validator
        // counts at.
        let pileup_of = |wanted: &str| -> Option<gatk_engine::read_pileup::ReadPileup<'_>> {
            let elements: Vec<gatk_engine::pileup::PileupElement> = reads
                .iter()
                .filter(|(sample, _)| sample == wanted)
                .filter_map(|(_, read)| {
                    let offset = record.start as i32 - read.alignment_start;
                    gatk_engine::pileup::PileupElement::for_read_and_offset(read, offset)
                })
                .collect();
            if elements.is_empty() {
                None
            } else {
                Some(gatk_engine::read_pileup::ReadPileup::new(
                    &record.contig,
                    record.start as i32,
                    elements,
                ))
            }
        };
        let case = pileup_of(&options.validation_case_name);
        let control = pileup_of(&options.validation_control_name);

        let applied = validate::apply(
            &record.contig,
            record.start as i32,
            record.stop as i32,
            &allele(&record.alleles[0]),
            &record.alleles[1..].iter().map(allele).collect::<Vec<_>>(),
            &record.filters.clone().unwrap_or_default(),
            &validation_genotype,
            case.as_ref(),
            control.as_ref(),
            &options,
        )
        .map_err(|error| Thrown {
            failure: match error {
                validate::ToolError::NullResult => Failure::Other,
                _ => Failure::User,
            },
            exception: error.java_class(),
            message: Some(error.message()),
        })?;
        let Some(applied) = applied else {
            continue;
        };
        if let Some(result) = &applied.result {
            table.push(result.clone());
        }
        validate::count_towards_summary(
            &mut summary,
            &applied,
            gatk_tools::remove_nearby_indels::variant_type(record)
                == gatk_tools::remove_nearby_indels::VariantType::Snp,
        );
        // The annotated record is the input's with the judgment in INFO; the genotypes travel with
        // it, because the writer's header carries the input's samples.
        let mut copy = record.clone();
        copy.attributes.retain(|(key, _)| {
            key != validate::JUDGMENT_KEY
                && key != validate::POWER_KEY
                && key != validate::VALIDATION_AD_KEY
        });
        copy.attributes.push((
            validate::JUDGMENT_KEY.to_string(),
            htsjdk_vcf::variant::Value::Str(applied.judgment.name().to_string()),
        ));
        if let (Some(power), Some((reference_count, alternate_count))) =
            (applied.power, applied.validation_ad)
        {
            // `VCFEncoder`'s own rendering, which is three decimals for a double in this range:
            // the port wrote the full `0.41258741258741227` where the reference wrote `0.413`.
            copy.attributes.push((
                validate::POWER_KEY.to_string(),
                htsjdk_vcf::variant::Value::Str(htsjdk_vcf::variant::format_vcf_double(power)),
            ));
            copy.attributes.push((
                validate::VALIDATION_AD_KEY.to_string(),
                htsjdk_vcf::variant::Value::Str(format!("{reference_count},{alternate_count}")),
            ));
        }
        annotated.push(copy);
    }

    write_file(
        &output,
        gatk_engine::basic_somatic_short_mutation_validator::write_table(&table).as_bytes(),
    )?;
    if let Some(path) = argument(parser, "summary") {
        write_file(&path, summary.table().as_bytes())?;
    }
    if let Some(path) = argument(parser, "annotated-vcf") {
        // `new VCFHeader(headerLines, inputHeader.getGenotypeSamples())`: the input's lines, the
        // three this tool declares, the tool's default lines, and the input's OWN samples, so the
        // annotated file keeps its genotype columns.
        let mut header = file.header.clone();
        for (id, number, line_type, description) in [
            (
                validate::JUDGMENT_KEY,
                htsjdk_vcf::header::Cardinality::Fixed(1),
                htsjdk_vcf::header::LineType::String,
                "Validation judgment: validated, unvalidated, or skipped.",
            ),
            (
                validate::POWER_KEY,
                htsjdk_vcf::header::Cardinality::Fixed(1),
                htsjdk_vcf::header::LineType::Float,
                "Power to validate variant in validation bam.",
            ),
            (
                validate::VALIDATION_AD_KEY,
                htsjdk_vcf::header::Cardinality::A,
                htsjdk_vcf::header::LineType::Integer,
                "Ref and alt allele count in validation bam.",
            ),
        ] {
            header.lines.push(htsjdk_vcf::header::HeaderLine::Compound {
                key: "INFO".to_string(),
                id: id.to_string(),
                number,
                line_type,
                description: description.to_string(),
                extra: Vec::new(),
            });
        }
        if flag(parser, "add-output-vcf-command-line") {
            header
                .lines
                .push(htsjdk_vcf::header::HeaderLine::Unstructured {
                    key: "source".to_string(),
                    value: "ValidateBasicSomaticShortMutations".to_string(),
                });
            if let Some(command_line) =
                command_line_header_line(parser, "ValidateBasicSomaticShortMutations")
            {
                header.lines.push(command_line);
            }
        }
        // `--sites-only-vcf-output` empties the sample set on the header and the genotypes on the
        // records, the same way it does on every other writer here.
        let mut records = annotated.clone();
        if flag(parser, "sites-only-vcf-output") {
            header.samples.clear();
            for record in &mut records {
                record.genotypes.clear();
            }
        }
        let text = htsjdk_vcf::vcf_file::write_vcf(&header, &records).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("{error:?}")),
        })?;
        write_variant_output(parser, &path, &text)?;
    }
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `PrintReadCounts.apply`: a depth-evidence or counts file rewritten for the CNV callers.
///
/// A `FeatureWalker`, and the two feature types it accepts disagree about what a header is: an
/// `.rd.txt` carries one line of column names and NO dictionary, so the run needs
/// `--sequence-dictionary` and refuses without one; a `.counts.tsv` carries a whole SAM header and
/// never consults it. Which files come out, and what is in them when the run does not finish, is
/// [`gatk_tools::print_read_counts`], where a golden measures it.
pub fn print_read_counts(parser: &Parser) -> Outcome {
    use gatk_tools::print_read_counts as counts;

    let _ = resolve_read_filters(parser, "PrintReadCounts")?;
    let input = argument(parser, "input-counts").ok_or_else(|| {
        Thrown::command_line(
            "Argument input-counts was missing: Argument 'input-counts' is required",
        )
    })?;
    let prefix = argument(parser, "output-prefix").unwrap_or_default();
    // `FeatureManager` asks every codec whether it can decode the path, and the two that matter
    // here answer by EXTENSION: `.counts.tsv` for the simple counts and `.rd.txt` for the depth
    // evidence. A file called `counts.tsv` has neither, so it is refused for having no codec
    // rather than read as the table it holds. Measured on four rows of this tool's array.
    if !(input.ends_with(".counts.tsv") || input.ends_with(".rd.txt")) {
        return Err(Thrown::user(format!(
            "Cannot read file://{} because no suitable codecs found",
            java_absolute_path(&input)
        )));
    }
    let text = std::fs::read_to_string(&input).map_err(|_| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: input.clone(),
            }
            .message(),
        )
    })?;

    // The two headers this tool accepts, told apart the way the codecs tell them apart: a SAM
    // header begins with `@`, a depth header with its column names.
    let parsed = if text.starts_with('@') {
        let header = htsjdk_bam::reader::parse_header_text(
            &text
                .lines()
                .take_while(|line| line.starts_with('@'))
                .map(|line| format!("{line}\n"))
                .collect::<String>(),
        );
        let records = text
            .lines()
            .filter(|line| !line.starts_with('@') && !line.trim().is_empty())
            .skip(1)
            .filter_map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                (columns.len() >= 4).then(|| counts::SimpleCount {
                    contig: columns[0].to_string(),
                    start: columns[1].parse().unwrap_or_default(),
                    end: columns[2].parse().unwrap_or_default(),
                    count: columns[3].parse().unwrap_or_default(),
                })
            })
            .collect();
        counts::Input::Counts(counts::CountsFile { header, records })
    } else {
        let mut lines = text.lines().filter(|line| !line.trim().is_empty());
        let header = lines.next().unwrap_or_default();
        counts::Input::Depth(counts::DepthFile {
            samples: counts::depth_header_samples(header),
            records: lines.filter_map(counts::decode_depth).collect(),
        })
    };

    let dictionary = master_dictionary(parser)?.map(|header| header.sequences);
    let intervals: Vec<counts::Interval> = arguments(parser, "intervals")
        .iter()
        .filter_map(|query| {
            let (contig, range) = query.split_once(':')?;
            let (start, end) = range.split_once('-')?;
            Some(counts::Interval {
                contig: contig.to_string(),
                start: start.parse().ok()?,
                end: end.parse().ok()?,
            })
        })
        .collect();

    let run = counts::run(
        &parsed,
        dictionary.as_deref(),
        &prefix,
        argument(parser, "output-file-list").as_deref(),
        &intervals,
    );
    // A refused run still WROTE: the files the tool had opened before it threw are on disk, half a
    // header and all, which is what the golden measures and why they are written before the
    // refusal is raised.
    for (path, contents) in run.disk.files() {
        write_file(&path, contents.as_bytes())?;
    }
    if let Some(error) = &run.error {
        // The class is borrowed from the error rather than a constant, so it is matched to one of
        // the three the module can name: a `Thrown` holds a `&'static str`.
        let class = match error.java_class() {
            "org.broadinstitute.hellbender.exceptions.UserException" => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            "java.lang.ArrayIndexOutOfBoundsException" => {
                "java.lang.ArrayIndexOutOfBoundsException"
            }
            _ => "java.lang.IllegalArgumentException",
        };
        return Err(Thrown {
            failure: if class == "org.broadinstitute.hellbender.exceptions.UserException" {
                Failure::User
            } else {
                Failure::Other
            },
            exception: class,
            message: Some(error.message()),
        });
    }
    Ok(None)
}

/// One read as the four SV evidence writers read it.
fn sv_read(
    record: &htsjdk_bam::record::BamRecord,
    header: &SamHeader,
) -> gatk_tools::collect_sv_evidence::Read {
    let contig_of = |index: i32| -> Option<String> {
        header
            .sequences
            .get(usize::try_from(index).ok()?)
            .map(|sequence| sequence.name.clone())
    };
    gatk_tools::collect_sv_evidence::Read {
        name: record.read_name.clone(),
        contig_index: usize::try_from(record.reference_index).unwrap_or(usize::MAX),
        contig: contig_of(record.reference_index).unwrap_or_default(),
        start: record.alignment_start,
        mapping_quality: i32::from(record.mapping_quality),
        cigar: record
            .cigar
            .elements
            .iter()
            .map(|element| (element.op.to_char() as char, element.length as i32))
            .collect(),
        paired: record.flags & 0x1 != 0,
        properly_paired: record.flags & 0x2 != 0,
        mate_unmapped: record.flags & 0x8 != 0,
        mate_contig_index: usize::try_from(record.mate_reference_index).ok(),
        mate_contig: contig_of(record.mate_reference_index),
        mate_start: Some(record.mate_alignment_start),
        reverse_strand: record.flags & 0x10 != 0,
        mate_reverse_strand: record.flags & 0x20 != 0,
        supplementary: record.flags & 0x800 != 0,
        secondary: record.flags & 0x100 != 0,
        duplicate: record.flags & 0x400 != 0,
        unmapped: record.flags & 0x4 != 0,
        bases: record.read_bases.clone(),
        base_qualities: record
            .base_qualities
            .iter()
            .map(|quality| i32::from(*quality))
            .collect(),
    }
}

/// `CollectSVEvidence.apply`: one BAM walked once for four kinds of evidence.
///
/// Every tool with a runner here so far wrote one kind of thing; this one writes up to four, and
/// each has its own rule for which read it will look at. What a read contributes is
/// [`gatk_tools::collect_sv_evidence`], where a golden measures it; the runner opens the BAM, the
/// sites and the intervals, and checks the file NAMES, which the reference does before it reads a
/// record: a writer refuses a name it could not read back.
pub fn collect_sv_evidence(parser: &Parser) -> Outcome {
    use gatk_tools::collect_sv_evidence as evidence;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CollectSVEvidence")?;
    let _ = argument(parser, "sample-name").ok_or_else(|| {
        Thrown::command_line("Argument sample-name was missing: Argument 'sample-name' is required")
    })?;

    let pair_file = argument(parser, "pe-file");
    let split_file = argument(parser, "sr-file");
    let site_file = argument(parser, "sd-file");
    let depth_file = argument(parser, "depth-evidence-file");
    if pair_file.is_none() && split_file.is_none() && site_file.is_none() && depth_file.is_none() {
        return Err(Thrown::user(evidence::NO_OUTPUT_MESSAGE));
    }
    // Each writer tests the name it was given BEFORE the traversal, and refuses one it could not
    // read back: which names it can is the codecs' own `canDecode`.
    let outputs = [
        ("pe", &pair_file),
        ("sr", &split_file),
        ("sd", &site_file),
        ("rd", &depth_file),
    ];
    for (kind, path) in outputs {
        if let Some(name) = path {
            if evidence::encoding(kind, name).is_none() {
                return Err(Thrown::user(evidence::bad_name_message(kind, name)));
            }
        }
    }
    // A name the codecs accept can still select a writer this port does not carry: GATK writes a
    // block-compressed name through BGZF at `--compression-level` with a tabix index beside it,
    // and a `.bci` name through the binary container. Only plain text is written below, so the
    // other two refuse here, after every name was tested and before anything is written.
    for (kind, path) in outputs {
        if let Some(name) = path {
            if evidence::encoding(kind, name)
                != Some(evidence::Encoding::Text {
                    block_compressed: false,
                })
            {
                return Err(Thrown::non_user(
                    PORT_LIMITATION,
                    format!(
                        "{name} asks for block-compressed or .bci evidence, which this port does \
                         not write yet. This message is the port's own and not GATK's."
                    ),
                ));
            }
        }
    }

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(reads_traversal_error)?;
    let reads: Vec<evidence::Read> = records
        .iter()
        .map(|record| sv_read(record, &header))
        .collect();

    if let Some(path) = &pair_file {
        let mut text = String::new();
        for pair in evidence::discordant_pairs(&reads) {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                pair.contig,
                pair.position,
                if pair.strand { "+" } else { "-" },
                pair.mate_contig,
                pair.mate_position,
                if pair.mate_strand { "+" } else { "-" }
            ));
        }
        write_file(path, text.as_bytes())?;
    }
    if let Some(path) = &split_file {
        let mut text = String::new();
        for split in evidence::split_reads(&reads) {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                split.contig,
                split.position,
                match split.side {
                    evidence::Side::Left => "left",
                    evidence::Side::Right => "right",
                    evidence::Side::Middle => "middle",
                },
                split.count
            ));
        }
        write_file(path, text.as_bytes())?;
    }
    if let Some(path) = &site_file {
        let sites = match argument(parser, "site-depth-locs-vcf") {
            Some(vcf) => {
                let text = std::fs::read_to_string(&vcf).map_err(|_| {
                    Thrown::user(
                        index_feature_file::Refusal::CouldNotReadInputFile { path: vcf.clone() }
                            .message(),
                    )
                })?;
                let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
                    failure: Failure::User,
                    exception: "htsjdk.tribble.TribbleException",
                    message: Some(failure.error.message()),
                })?;
                file.records
                    .iter()
                    .map(|record| evidence::Site {
                        contig: record.contig.clone(),
                        position: record.start as i32,
                        reference: record.alleles[0].display_string(),
                        alternates: record.alleles[1..]
                            .iter()
                            .map(|allele| allele.display_string())
                            .collect(),
                    })
                    .collect()
            }
            None => Vec::new(),
        };
        // Each row carries the SAMPLE between the position and the four counts: the file is read
        // back per sample, so the name travels with every record rather than in a header.
        let sample = argument(parser, "sample-name").unwrap_or_default();
        let mut text = String::new();
        for depth in evidence::site_depths(
            &reads,
            &sites,
            number_or(
                parser,
                "site-depth-min-mapq",
                evidence::DEFAULT_SITE_DEPTH_MIN_MAPQ,
            ),
            number_or(
                parser,
                "site-depth-min-baseq",
                evidence::DEFAULT_SITE_DEPTH_MIN_BASEQ,
            ),
        ) {
            text.push_str(&format!(
                "{}\t{}\t{sample}\t{}\t{}\t{}\t{}\n",
                depth.contig,
                depth.position,
                depth.counts[0],
                depth.counts[1],
                depth.counts[2],
                depth.counts[3]
            ));
        }
        write_file(path, text.as_bytes())?;
    }
    if let Some(path) = &depth_file {
        let windows: Vec<(String, i32, i32)> = match argument(parser, "depth-evidence-intervals") {
            Some(list) => {
                let text = std::fs::read_to_string(&list).map_err(|_| {
                    Thrown::user(
                        index_feature_file::Refusal::CouldNotReadInputFile { path: list.clone() }
                            .message(),
                    )
                })?;
                let parsed: Vec<(String, i32, i32)> = text
                    .lines()
                    .filter(|line| !line.starts_with('@') && !line.trim().is_empty())
                    .filter_map(|line| {
                        let columns: Vec<&str> = line.split('\t').collect();
                        Some((
                            columns.first()?.to_string(),
                            columns.get(1)?.parse().ok()?,
                            columns.get(2)?.parse().ok()?,
                        ))
                    })
                    .collect();
                if parsed.is_empty() {
                    return Err(Thrown::user(evidence::empty_intervals_message(&list)));
                }
                parsed
            }
            None => Vec::new(),
        };
        // The depth file opens with a column line naming the sample, which is the one evidence
        // file of the four that has a header at all.
        let mut text = format!(
            "#Chr\tStart\tEnd\t{}\n",
            argument(parser, "sample-name").unwrap_or_default()
        );
        for depth in evidence::depth_evidence(
            &reads,
            &windows,
            number_or(
                parser,
                "depth-evidence-min-mapq",
                evidence::DEFAULT_DEPTH_EVIDENCE_MIN_MAPQ,
            ),
        ) {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                depth.contig, depth.start, depth.end, depth.count
            ));
        }
        write_file(path, text.as_bytes())?;
    }
    Ok(None)
}

/// The name `VariantContext.Type` prints, which the two GVCF messages quote.
fn variant_context_type_name(kind: gatk_tools::remove_nearby_indels::VariantType) -> &'static str {
    use gatk_tools::remove_nearby_indels::VariantType;
    match kind {
        VariantType::NoVariation => "NO_VARIATION",
        VariantType::Snp => "SNP",
        VariantType::Mnp => "MNP",
        VariantType::Indel => "INDEL",
        VariantType::Symbolic => "SYMBOLIC",
        VariantType::Mixed => "MIXED",
    }
}

/// One record as `ValidateVariants` reads it, which is less of it than a writer needs.
///
/// The attributes are SORTED, because the only place they are read is a message that prints them
/// through a `TreeMap`, and the filters are the applied ones: a record that passed carries an empty
/// list here whether its column said `PASS` or `.`.
fn validation_record(
    vc: &htsjdk_vcf::variant::VariantContext,
) -> gatk_tools::validate_variants::Record {
    let integers = |key: &str| -> Vec<i32> {
        vc.attributes
            .iter()
            .find(|(name, _)| name == key)
            .and_then(|(_, value)| value.format())
            .map(|text| {
                text.split(',')
                    .filter_map(|piece| piece.trim().parse().ok())
                    .collect()
            })
            .unwrap_or_default()
    };
    let mut attributes: Vec<(String, String)> = vc
        .attributes
        .iter()
        .map(|(key, value)| (key.clone(), value.format().unwrap_or_default()))
        .collect();
    attributes.sort_by(|left, right| left.0.cmp(&right.0));

    gatk_tools::validate_variants::Record {
        contig: vc.contig.clone(),
        start: vc.start as i32,
        reference: vc
            .alleles
            .first()
            .map(|allele| allele.display_string())
            .unwrap_or_default(),
        alternates: vc.alleles[1..]
            .iter()
            .map(|allele| allele.display_string())
            .collect(),
        filters: vc.filters.clone().unwrap_or_default(),
        allele_counts: integers("AC"),
        allele_number: integers("AN").first().copied(),
        genotypes: vc
            .genotypes
            .iter()
            .map(|genotype| {
                genotype
                    .alleles
                    .iter()
                    .map(|allele| {
                        if allele.is_no_call() {
                            None
                        } else {
                            vc.alleles.iter().position(|candidate| candidate == allele)
                        }
                    })
                    .collect()
            })
            .collect(),
        qual: if vc.has_log10_p_error() {
            Some(vc.phred_scaled_qual())
        } else {
            None
        },
        variant_type: variant_context_type_name(gatk_tools::remove_nearby_indels::variant_type(vc))
            .to_string(),
        attributes,
    }
}

/// `ValidateVariants.apply` and `onTraversalSuccess`: a VCF checked, and NOTHING written.
///
/// The only tool with a runner here whose whole output is its refusal. Three things it decides are
/// the tool rather than the checks, and all three are measured:
///
///   * `--validation-type-to-exclude` is the argument that turns checks ON. With nothing excluded
///     the run is `ALL`, which is what the inputs allow; exclude anything at all and the concrete
///     set is built and the exclusions removed from it, so `REF` comes back and a run with no
///     `--reference` is refused before a record is read;
///   * `--validate-GVCF` excludes the allele check on its own account, which is what sends a plain
///     GVCF run down that same branch;
///   * `--warn-on-errors` turns every refusal into a log line, so the run SUCCEEDS and writes
///     nothing, which is indistinguishable from a valid file except in the log.
///
/// The GVCF coverage check runs at the end over the whole traversal, and what it subtracts from is
/// the interval argument when there is one and the whole dictionary when there is not.
pub fn validate_variants(parser: &Parser) -> Outcome {
    use gatk_tools::validate_variants as validate;

    let VariantWalkerStart {
        input,
        text,
        codec: _,
        intervals,
    } = variant_walker_startup(parser, "ValidateVariants")?;

    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;

    let mut types_to_exclude = Vec::new();
    for name in arguments(parser, "validation-type-to-exclude") {
        match validate::ValidationType::parse(&name) {
            Some(kind) => types_to_exclude.push(kind),
            // The parser refuses an unknown constant before the tool runs, so this is unreachable
            // from a command line; a value that got here is the port's own failure to say so.
            None => {
                return Err(Thrown::non_user(
                    PORT_FAILURE,
                    format!("{name} is not a ValidationType this port knows."),
                ))
            }
        }
    }

    let reference_path = argument(parser, "reference");
    let validate_gvcf = flag(parser, "validate-GVCF");
    let arguments = validate::Arguments {
        types_to_exclude,
        do_not_validate_filtered_records: flag(parser, "do-not-validate-filtered-records"),
        warn_on_errors: flag(parser, "warn-on-errors"),
        validate_gvcf,
        has_reference: reference_path.is_some(),
        has_dbsnp: argument(parser, "dbsnp").is_some(),
    };

    // `throwOrWarn`: with `--warn-on-errors` the message is logged and the traversal carries on, so
    // every refusal below goes through this and a warning run ends at zero.
    let refuse = |error: validate::ValidationError| -> Option<Thrown> {
        if arguments.warn_on_errors {
            return None;
        }
        Some(Thrown {
            failure: Failure::User,
            exception: error.java_class(),
            message: Some(error.message()),
        })
    };

    let types = match validate::types_to_apply(&arguments) {
        Ok(types) => types,
        Err(error) => match refuse(error) {
            Some(thrown) => return Err(thrown),
            // `calculateValidationTypesToApply` is called from `onTraversalStart`, where the warn
            // path swallows nothing: the refusal is raised before `throwOrWarn` exists. Measured
            // by the golden this module carries.
            None => {
                return Err(Thrown::user(
                    validate::ValidationError::MissingReference.message(),
                ))
            }
        },
    };

    let mut reference = match &reference_path {
        Some(path) => Some(
            gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
        ),
        None => None,
    };

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

    let mut order = validate::OrderCheck::new();
    let mut covered: Vec<gatk_engine::interval::SimpleInterval> = Vec::new();
    // The overlap watch, which is `--fail-gvcf-on-overlap`'s whole input. `previous` is the
    // MERGED interval rather than the record's own, because that is what the reference compares
    // against, and `overlapping` is overwritten on every hit: the message calls it the first
    // overlapping interval and the assignment makes it the last.
    let mut previous: Option<(gatk_engine::interval::SimpleInterval, bool)> = None;
    let mut overlapping: Option<gatk_engine::interval::SimpleInterval> = None;
    for located in gatk_engine::variant_source::traverse(&located, intervals.as_deref()) {
        let vc = &file.records[located.index];
        let record = validation_record(vc);

        if validate_gvcf {
            if let Err(error) = order.check(&record) {
                if let Some(thrown) = refuse(error) {
                    return Err(thrown);
                }
            }
            // The blocks a GVCF writes are adjacent rather than overlapping, so the reference
            // merges an adjacent pair (margin ONE) into one interval before adding it, and counts
            // an actually overlapping pair (margin ZERO) as an overlap. A reference block is a
            // record whose only alternate is `<NON_REF>`, and only a pair with one of those in it
            // is an overlap worth reporting.
            let this_is_reference =
                record.alternates.len() == 1 && record.alternates[0] == "<NON_REF>";
            if let Some(interval) = gatk_engine::interval::SimpleInterval::new(
                &vc.contig,
                vc.start as i32,
                vc.stop as i32,
            ) {
                let overlaps = |margin: i32| -> bool {
                    match &previous {
                        Some((before, _)) => {
                            before.contig == interval.contig
                                && before.start <= interval.end + margin
                                && interval.start <= before.end + margin
                        }
                        None => false,
                    }
                };
                if overlaps(0)
                    && (previous.as_ref().is_some_and(|(_, was)| *was) || this_is_reference)
                {
                    overlapping = Some(interval.clone());
                }
                let merged = if overlaps(1) {
                    let (before, _) = previous.as_ref().expect("a previous interval to overlap");
                    gatk_engine::interval::SimpleInterval::new(
                        &interval.contig,
                        before.start,
                        before.end.max(interval.end),
                    )
                    .unwrap_or_else(|| interval.clone())
                } else {
                    interval.clone()
                };
                covered.push(merged.clone());
                previous = Some((merged, this_is_reference));
            }
        }

        // `getRefBasesAtPosition(reader, contig, start, refLength)`: the reference under the whole
        // reference allele, which is what the REF check compares against.
        let observed = match &mut reference {
            Some(source) => {
                let length = record.reference.len() as i32;
                let bases = source
                    .query(&record.contig, record.start, record.start + length - 1)
                    .map_err(|error| match error {
                        // A reference that does not carry the record's contig is a refusal and not
                        // a skipped check, and the dictionary it prints is the REFERENCE's.
                        // Measured on a row of this tool's array where `--reference` is the
                        // corpus's other contig and dictionary validation is turned off.
                        gatk_engine::reference::ReferenceError::UnknownContig(contig) => {
                            Thrown::user(format!(
                                "Contig {contig} not present in the sequence dictionary {}\n",
                                gatk_tools::sequence_dictionary::pretty_print(
                                    &gatk_tools::reference_walker::dictionary(source).sequences
                                )
                            ))
                        }
                        other => Thrown::user(format!("{other:?}")),
                    })?;
                Some(String::from_utf8_lossy(&bases).into_owned())
            }
            None => None,
        };
        if let Err(error) =
            validate::validate_record(&record, &input, &types, observed.as_deref(), &arguments)
        {
            if let Some(thrown) = refuse(error) {
                return Err(thrown);
            }
        }
    }

    if validate_gvcf {
        // The whole region is the interval argument when there is one and the dictionary when there
        // is not, and what is subtracted from it is the merged span of every record traversed.
        let dictionary = match master_dictionary(parser)?.or(reference_dictionary(parser)?) {
            Some(header) => header,
            None => vcf_dictionary(&text),
        };
        let wanted = match intervals.as_deref() {
            Some(given) if !given.is_empty() => given.to_vec(),
            _ => gatk_engine::interval_args::whole_reference(&dictionary),
        };
        let uncovered =
            gatk_engine::interval_args::subtract_regions(&wanted, &covered, &dictionary);
        let loci: i64 = uncovered
            .iter()
            .map(|interval| i64::from(interval.end - interval.start + 1))
            .sum();
        if loci > 0 {
            let first = &uncovered[0];
            // `GenomeLoc.toString`, which collapses a one-base interval to a single position.
            let rendered = if first.start == first.end {
                format!("{}:{}", first.contig, first.start)
            } else {
                format!("{}:{}-{}", first.contig, first.start, first.end)
            };
            if let Some(thrown) = refuse(validate::ValidationError::NotCovering {
                loci,
                first_gap: rendered,
            }) {
                return Err(thrown);
            }
        }
        // And the overlap, which is checked AFTER the coverage: a GVCF that both overlaps and
        // leaves a gap is refused for the gap.
        if flag(parser, "fail-gvcf-on-overlap") {
            if let Some(interval) = &overlapping {
                let message = format!(
                    "This GVCF contained overlapping reference blocks.  The first overlapping \
                     interval is {}:{}-{}",
                    interval.contig, interval.start, interval.end
                );
                if !arguments.warn_on_errors {
                    return Err(Thrown {
                        failure: Failure::User,
                        exception:
                            "org.broadinstitute.hellbender.exceptions.UserException$ValidationFailure",
                        message: Some(message),
                    });
                }
            }
        }
    }
    Ok(None)
}

/// `GetNormalArtifactData.apply`: one locus, two pileups of the same reads, and a seeded draw.
///
/// The first Mutect tool with a runner here, and the RNG is what makes it one: the whole traversal
/// shares one `java.util.Random(47382911)`, drawn on once per CANDIDATE locus, so which loci reach
/// the table depends on how many came before. The draw also happens BEFORE the last rejection
/// rule, which is why a locus rejected for too much tumour support still consumes a number.
///
/// The split is by SAMPLE and not by file: `--normal-sample` names the samples whose reads are the
/// normal, and everything else in the same pileup is the tumour. One input carrying two samples is
/// therefore the shape this port reads, and the tool's own arithmetic lives in
/// [`gatk_tools::get_normal_artifact_data`], where a golden already measures it.
pub fn get_normal_artifact_data(parser: &Parser) -> Outcome {
    use gatk_tools::get_normal_artifact_data as artifact;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "GetNormalArtifactData")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let normal_samples = arguments(parser, "normal-sample");
    let error_probability = scalar(parser, "error-prob")
        .and_then(|text| text.parse().ok())
        .unwrap_or(artifact::DEFAULT_ERROR_PROBABILITY);

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;
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
    .map_err(locus_traversal_error)?;

    let reference_sequences = gatk_tools::reference_walker::dictionary(&reference).sequences;
    let mut generator = gatk_engine::java_random::JavaRandom::gatk();
    let mut rows = Vec::new();
    for one in &applied {
        // `ReadUtils.getSampleName` per element, against the list: a read with no read group has
        // no sample, which is in neither list and therefore counts as tumour, exactly as the
        // reference's `contains(null)` on a list of names does.
        let is_normal = |element: &gatk_engine::pileup::PileupElement| {
            gatk_engine::read_pileup::sample_name(element.read, &header)
                .is_some_and(|sample| normal_samples.contains(&sample))
        };
        let normal = one.context.pileup.filtered(is_normal);
        let tumor = one.context.pileup.filtered(|element| !is_normal(element));

        let bases = reference
            .query(
                &one.context.contig,
                one.context.position,
                one.context.position,
            )
            .map_err(|error| match error {
                gatk_engine::reference::ReferenceError::UnknownContig(contig) => {
                    Thrown::user(format!(
                        "Contig {contig} not present in the sequence dictionary {}\n",
                        gatk_tools::sequence_dictionary::pretty_print(&reference_sequences)
                    ))
                }
                other => Thrown::user(format!("{other:?}")),
            })?;
        let Some(base) = bases.first().copied() else {
            continue;
        };
        if let artifact::Outcome::Kept(record) =
            artifact::apply(&normal, &tumor, base, error_probability, || {
                generator.next_double()
            })
        {
            rows.push(*record);
        }
    }

    write_file(&output, artifact::write(&rows).as_bytes())?;
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `CollectF1R2Counts`, a locus walker whose one output is an ARCHIVE written by `closeTool`.
///
/// The collector exists from `onTraversalStart`, and with it one alt table writer per sample; the
/// histograms are only written by `onTraversalSuccess`. `closeTool` runs in `doWork`'s `finally`,
/// so a traversal that is refused part way still leaves a `.tar.gz` behind, holding the alt tables
/// as far as they got and no histogram at all. A refusal in the startup, before the collector
/// exists, leaves nothing.
///
/// The members are named `./<url-encoded sample><extension>`, in `File.listFiles()` order, which is
/// the file system's and not reproducible; the port writes them sorted, and the covering array
/// compares the archive member by member rather than byte by byte.
pub fn collect_f1r2_counts(parser: &Parser) -> Outcome {
    use gatk_tools::collect_f1r2_counts as f1r2;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CollectF1R2Counts")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let args = f1r2::Args {
        min_median_map_qual: number_or(parser, "f1r2-median-mq", 50),
        min_base_quality: number_or(parser, "f1r2-min-bq", 20),
        max_depth: number_or(parser, "f1r2-max-depth", f1r2::DEFAULT_MAX_DEPTH),
    };

    // `ReadUtils.getSamplesFromHeader`: the read groups' samples, once each.
    let mut samples: Vec<String> = Vec::new();
    for group in &header.read_groups {
        if let Some(sample) = group.attributes.get("SM") {
            if !samples.iter().any(|known| known == sample) {
                samples.push(sample.to_string());
            }
        }
    }
    let mut collector = f1r2::Collector::new(args, &samples);
    let archive = |collector: &f1r2::Collector, histograms: bool| -> Result<(), Thrown> {
        let mut files = collector.files(histograms);
        files.sort_by(|left, right| left.0.cmp(&right.0));
        let entries: Vec<(String, Vec<u8>)> = files
            .into_iter()
            .map(|(name, text)| (format!("./{name}"), text.into_bytes()))
            .collect();
        write_file(&output, &tar_gz(&entries))
    };

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;
    let applied = match gatk_tools::locus_walker::traverse(
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
    ) {
        Ok(applied) => applied,
        Err(error) => {
            archive(&collector, false)?;
            return Err(locus_traversal_error(error));
        }
    };

    let sequences = gatk_tools::reference_walker::dictionary(&reference).sequences;
    let one_sample = samples.len() == 1;
    for one in &applied {
        let Some(length) = sequences
            .iter()
            .find(|sequence| sequence.name == one.context.contig)
            .map(|sequence| sequence.length)
        else {
            archive(&collector, false)?;
            return Err(Thrown::user(format!(
                "Contig {} not present in the sequence dictionary {}\n",
                one.context.contig,
                gatk_tools::sequence_dictionary::pretty_print(&sequences)
            )));
        };
        let mut elements = Vec::new();
        for element in &one.context.pileup.elements {
            let sample = gatk_engine::read_pileup::sample_name(element.read, &header);
            // `splitBySample(header, null)`: a read with no sample is refused, but only on the
            // path that splits, which is every sample count but one.
            if sample.is_none() && !one_sample {
                archive(&collector, false)?;
                return Err(Thrown::user(format!(
                    "SAM/BAM/CRAM file (unknown) is malformed: Read {} is missing the read group \
                     (RG) tag, which is required by the GATK.  Please use \
                     http://gatkforums.broadinstitute.org/discussion/59/\
                     companion-utilities-replacereadgroups to fix this problem",
                    element.read.read_name
                )));
            }
            let flags = element.read.flags;
            elements.push(f1r2::Element {
                sample: sample.unwrap_or_default(),
                base: element.base(),
                qual: element.qual(),
                reverse_strand: flags & 0x10 != 0,
                first_of_pair: flags & 0x1 != 0 && flags & 0x40 != 0,
                mapping_quality: element.mapping_qual() as i32,
                deletion: element.is_deletion(),
                after_insertion: element.is_after_insertion(),
                before_deletion_start: element.is_before_deletion_start(),
            });
        }
        // `getKmerAround(position, 1)`: the window widened within the contig, and no k-mer at all
        // when the contig's end cut it short.
        let position = one.context.position;
        let start = (position - 1).max(1);
        let end = (position + 1).min(length);
        let kmer = if end - start < 2 * f1r2::REF_CONTEXT_PADDING as i32 {
            None
        } else {
            let bases = reference
                .query(&one.context.contig, start, end)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            Some(String::from_utf8_lossy(&bases).to_ascii_uppercase())
        };
        collector.process(&elements, kmer.as_deref());
    }

    archive(&collector, true)?;
    // What `onTraversalSuccess` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `CreateSomaticPanelOfNormals`, a variant walker over a multi-sample VCF.
///
/// The site rules, the germline test and the beta fit are [`gatk_tools::create_somatic_panel_of_normals`]'s;
/// what the runner adds is the header, which is built fresh rather than copied:
///
///   - the tool's own two INFO lines, one `##normal_sample` line per sample of the input, and the
///     default tool lines, all in a `HashSet` the writer sorts;
///   - and the contig lines of the INPUT's dictionary, set after the writer exists, so an input
///     with no contig line at all is a `NullPointerException` that leaves an empty file behind.
///
/// Each record written is a site with the input's alleles, no ID, no QUAL, no FILTER and no
/// genotype, holding `FRACTION` and `BETA` and nothing else. The germline resource is queried by
/// overlap and only its FIRST record is read, its `AF` values summed.
pub fn create_somatic_panel_of_normals(parser: &Parser) -> Outcome {
    use gatk_tools::create_somatic_panel_of_normals as pon;
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    use htsjdk_vcf::variant::Value;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "CreateSomaticPanelOfNormals")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let min_sample_count = number_or(
        parser,
        "min-sample-count",
        pon::DEFAULT_MIN_SAMPLE_COUNT as i32,
    );
    let max_germline_probability = scalar(parser, "max-germline-probability")
        .and_then(|text| text.parse::<f64>().ok())
        .unwrap_or(pon::DEFAULT_MAX_GERMLINE_PROBABILITY);
    let germline = match argument(parser, "germline-resource") {
        Some(path) => Some(feature_variants(&path)?.0),
        None => None,
    };

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let samples = file.header.samples.clone();

    let mut header = htsjdk_vcf::header::VcfHeader::new();
    header.lines.extend(default_tool_vcf_header_lines(
        parser,
        "CreateSomaticPanelOfNormals",
    ));
    header.lines.push(HeaderLine::Compound {
        key: "INFO".to_string(),
        id: "FRACTION".to_string(),
        number: Cardinality::Fixed(1),
        line_type: LineType::Float,
        description: "Fraction of samples exhibiting artifact".to_string(),
        extra: Vec::new(),
    });
    header.lines.push(HeaderLine::Compound {
        key: "INFO".to_string(),
        id: "BETA".to_string(),
        number: Cardinality::Fixed(2),
        line_type: LineType::Float,
        description: "Beta distribution parameters to fit artifact allele fractions".to_string(),
        extra: Vec::new(),
    });
    for sample in &samples {
        let line = HeaderLine::Unstructured {
            key: "normal_sample".to_string(),
            value: sample.clone(),
        };
        if !header.lines.contains(&line) {
            header.lines.push(line);
        }
    }
    let contigs: Vec<Vec<(String, String)>> = file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Contig { fields, .. } => Some(fields.clone()),
            _ => None,
        })
        .collect();
    if contigs.is_empty() {
        // `setSequenceDictionary(null)`, after `createVCFWriter` opened the file and before the
        // header was written: `closeTool` closes a writer that wrote nothing.
        write_variant_output(parser, &output, "")?;
        return Err(Thrown::non_user(
            "java.lang.NullPointerException",
            "Cannot invoke \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \
             \"dictionary\" is null",
        ));
    }
    for (index, fields) in contigs.iter().enumerate() {
        let field = |key: &str| {
            fields
                .iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value.clone())
        };
        // `getSAMSequenceRecord` then `new VCFContigHeaderLine(record, record.getAssembly())`: the
        // ID, the length (zero when the line had none) and the assembly survive, nothing else.
        let mut rebuilt = vec![
            ("ID".to_string(), field("ID").unwrap_or_default()),
            (
                "length".to_string(),
                field("length")
                    .and_then(|length| length.parse::<i32>().ok())
                    .unwrap_or(0)
                    .to_string(),
            ),
        ];
        if let Some(assembly) = field("assembly") {
            rebuilt.push(("assembly".to_string(), assembly));
        }
        header.lines.push(HeaderLine::Contig {
            index: index as i32,
            fields: rebuilt,
        });
    }

    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut records = written.to_vec();
        apply_sites_only(parser, &mut header, &mut records);
        let out = write_vcf_honouring_lenient(parser, &header, &records)?;
        write_variant_output(parser, &output, &out)
    };

    let germline_frequency = |record: &htsjdk_vcf::variant::VariantContext| -> f64 {
        let Some(resource) = germline.as_ref() else {
            return 0.0;
        };
        let Some(first) = resource.iter().find(|candidate| {
            candidate.contig == record.contig
                && candidate.start <= record.stop
                && record.start <= candidate.stop
        }) else {
            return 0.0;
        };
        // `getAttributeAsDoubleList(vc, AF, 0.0)`: a missing value reads as the default.
        let text_of = |value: &Value| value.format().unwrap_or_default();
        let values: Vec<String> = match first.attributes.iter().find(|(key, _)| key == "AF") {
            None => Vec::new(),
            Some((_, Value::List(items))) => items.iter().map(text_of).collect(),
            Some((_, value)) => text_of(value).split(',').map(str::to_string).collect(),
        };
        values
            .iter()
            .map(|value| {
                if value == "." {
                    0.0
                } else {
                    value.parse::<f64>().unwrap_or(0.0)
                }
            })
            .sum()
    };

    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for record in kept {
        record.genotypes.decode();
        let site = pon::Site {
            contig: record.contig.clone(),
            position: record.start as i32,
            reference: record.reference().display_string(),
            alternates: record
                .alternate_alleles()
                .iter()
                .map(|allele| allele.display_string())
                .collect(),
            genotypes: record
                .genotypes
                .iter()
                .map(|genotype| pon::Genotype {
                    sample: genotype.sample_name.clone(),
                    allele_depths: genotype.ad.clone(),
                })
                .collect(),
        };
        let entries = pon::build_panel(
            std::slice::from_ref(&site),
            samples.len(),
            |_| germline_frequency(record),
            min_sample_count.max(0) as usize,
            max_germline_probability,
        );
        let Some(entry) = entries.into_iter().next() else {
            continue;
        };
        let mut out = htsjdk_vcf::variant::VariantContext::new(
            &record.contig,
            record.start,
            record.alleles.clone(),
        );
        out.stop = record.stop;
        out.attributes = vec![
            ("FRACTION".to_string(), Value::Double(entry.fraction)),
            (
                "BETA".to_string(),
                Value::List(vec![
                    Value::Double(entry.beta.alpha),
                    Value::Double(entry.beta.beta),
                ]),
            ),
        ];
        written.push(out);
    }

    finish(&written)?;
    Ok(Some("SUCCESS".to_string()))
}

/// `SplitCRAM`, a `CommandLineProgram` that cuts a CRAM at container boundaries.
///
/// Where the cuts fall is [`gatk_tools::split_cram::plan`]'s. What the runner adds is the bytes of
/// each shard, which are htsjdk's writers rather than a copy of the input's first bytes:
///
///   - the file definition, written back from the one read;
///   - the SAM header container, REBUILT: the header text is parsed and encoded again, the block is
///     GZIP at `Defaults.COMPRESSION_LEVEL`, which `GATKConfig` sets to two, and the container
///     header is `makeSAMFileHeaderContainer`'s, unmapped with one block and no landmark;
///   - every data container as it was read, which `Container.write` reproduces for a file htsjdk
///     wrote, its blocks keeping their compressed bytes;
///   - and the version 3 EOF container.
///
/// Nothing is printed: `doWork` returns null.
pub fn split_cram(parser: &Parser) -> Outcome {
    use gatk_tools::split_cram as split;
    use htsjdk_cram::varint::{write_unsigned_itf8, write_unsigned_ltf8};

    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let template =
        argument(parser, "output").unwrap_or_else(|| split::DEFAULT_TEMPLATE.to_string());
    let shard_records = scalar(parser, "shard-records")
        .and_then(|text| text.parse::<i64>().ok())
        .unwrap_or(split::DEFAULT_SHARD_RECORDS);
    let max_output_count = number_or(parser, "shard-max-output-count", 0);

    // `onStartup`, before the input is opened.
    if !split::accepts_template(&template) {
        let error = split::SplitError::TemplateMissingFormatter {
            template: template.clone(),
        };
        // `SplitError::java_class` borrows the error; the class is a constant.
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            error.message(),
        ));
    }

    let cram = std::fs::read(&input).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file {}. Error was: It doesn't exist.",
            std::path::Path::new(&input).display()
        ))
    })?;
    if cram.len() < 4 || &cram[..4] != b"CRAM" {
        return Err(Thrown::non_user(
            "java.lang.RuntimeException",
            "Input does not have a valid CRAM header.",
        ));
    }
    let walk = htsjdk_cram::file::read_file(&cram)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, error.message()))?;
    let major = walk.definition.major;

    // `CramIO.readSAMFileHeader`: a little-endian length, then that many bytes of text.
    let text_length = walk
        .sam_header
        .get(..4)
        .map(|bytes| i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]).max(0) as usize)
        .unwrap_or(0);
    let text_end = (4 + text_length).min(walk.sam_header.len());
    let text = String::from_utf8_lossy(walk.sam_header.get(4..text_end).unwrap_or(&[]));
    let header = htsjdk_bam::reader::parse_header_text(&text);
    let encoded = header.encode_replacing_version();

    // `samHeaderToByteArray`, then `createGZIPFileHeaderBlock` and the block's own `write`.
    let mut raw = (encoded.len() as i32).to_le_bytes().to_vec();
    raw.extend_from_slice(encoded.as_bytes());
    let compressed = java_gzip(&raw, 2);
    let mut block = vec![1u8, 0u8];
    block.extend_from_slice(&write_unsigned_itf8(0).0);
    block.extend_from_slice(&write_unsigned_itf8(compressed.len() as i32).0);
    block.extend_from_slice(&write_unsigned_itf8(raw.len() as i32).0);
    block.extend_from_slice(&compressed);
    if major >= 3 {
        let crc = htsjdk_cram::compression_header::crc32(&block);
        block.extend_from_slice(&crc.to_le_bytes());
    }
    // `makeSAMFileHeaderContainer(blockSize)`, written by `ContainerHeader.write`.
    let mut container = (block.len() as i32).to_le_bytes().to_vec();
    container.extend_from_slice(&write_unsigned_itf8(-1).0);
    container.extend_from_slice(&write_unsigned_itf8(0).0);
    container.extend_from_slice(&write_unsigned_itf8(0).0);
    container.extend_from_slice(&write_unsigned_itf8(0).0);
    container.extend_from_slice(&write_unsigned_ltf8(0).0);
    container.extend_from_slice(&write_unsigned_ltf8(0).0);
    container.extend_from_slice(&write_unsigned_itf8(1).0);
    container.extend_from_slice(&write_unsigned_itf8(0).0);
    if major >= 3 {
        let crc = htsjdk_cram::compression_header::crc32(&container);
        container.extend_from_slice(&crc.to_le_bytes());
    }
    let mut preamble = walk.definition.write();
    preamble.extend_from_slice(&container);
    preamble.extend_from_slice(&block);

    // The data containers, as bytes, stopping at the EOF container the iterator does not return.
    let containers: Vec<(&[u8], i32)> = walk
        .containers
        .iter()
        .filter(|one| !one.header.is_eof())
        .map(|one| {
            let end =
                one.offset + one.header.byte_length + one.header.blocks_byte_size.max(0) as usize;
            (
                &cram[one.offset..end.min(cram.len())],
                one.header.record_count,
            )
        })
        .collect();
    let counts: Vec<i32> = containers.iter().map(|(_, count)| *count).collect();
    let shards = split::plan(&counts, shard_records, max_output_count, &template)
        .map_err(|error| Thrown::non_user("java.lang.IllegalArgumentException", error.message()))?;

    let mut next = 0usize;
    for shard in &shards {
        let mut out = preamble.clone();
        for _ in &shard.containers {
            out.extend_from_slice(containers[next].0);
            next += 1;
        }
        out.extend_from_slice(&CRAM_V3_EOF);
        write_file(&shard.name, &out)?;
    }
    Ok(None)
}

/// `CramIO.ZERO_F_EOF_MARKER`, the version 3 EOF container.
const CRAM_V3_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// `java.util.zip.GZIPOutputStream` at `level`: the fixed ten-byte header with no name and no
/// time, whose operating system byte is 255 (unknown, measured on the oracle's JDK rather than the
/// zero older JDKs wrote), the JDK's zlib deflate, then the CRC and the length.
fn java_gzip(data: &[u8], level: u32) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    let mut compressor = flate2::Compress::new(flate2::Compression::new(level), false);
    let mut deflated = Vec::with_capacity(data.len() + 64);
    loop {
        let status = compressor
            .compress_vec(
                &data[compressor.total_in() as usize..],
                &mut deflated,
                flate2::FlushCompress::Finish,
            )
            .expect("deflating into a vector does not fail");
        if status == flate2::Status::StreamEnd {
            break;
        }
        deflated.reserve(deflated.capacity().max(64));
    }
    out.extend_from_slice(&deflated);
    let mut crc = flate2::Crc::new();
    crc.update(data);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

/// `FilterMutectCalls`, a `MultiplePassVariantWalker` over Mutect2's calls.
///
/// The four passes, the filters and the two outputs' content are
/// [`gatk_tools::filter_mutect_calls::run`]'s. What the runner adds is the translation both ways:
///
///   - each `VariantContext` becomes the `Record` the filters read, its genotypes told apart by the
///     header's `##normal_sample` lines, every other sample being a tumour;
///   - the header is rebuilt in `onTraversalStart` and written as soon as the writer exists, BEFORE
///     the stats table is looked for, so a missing table is a `CouldNotReadInputFile` that leaves a
///     VCF of nothing but its header;
///   - and each output record is the input's with its FILTER replaced, `AS_FilterStatus` set and the
///     phred-scaled posteriors the filters annotate, its genotypes written from the file.
///
/// The three inputs that bring a second model in, `--contamination-table`, `--tumor-segmentation`
/// and `--ob-priors`, are refused rather than ignored.
pub fn filter_mutect_calls(parser: &Parser) -> Outcome {
    use gatk_engine::accumulate_data::AccumulationAllele;
    use gatk_engine::allele_filter::GenotypeData;
    use gatk_engine::filtering_engine::{EngineArguments, Record};
    use gatk_engine::mutect_filter_list::FilterArguments;
    use gatk_engine::somatic_clustering_model::AlternateAllele;
    use gatk_engine::threshold_calculator::Strategy;
    use gatk_tools::filter_mutect_calls as fmc;
    use htsjdk_vcf::header::HeaderLine;
    use htsjdk_vcf::variant::Value;

    for held in ["contamination-table", "tumor-segmentation", "ob-priors"] {
        if !arguments(parser, held).is_empty() {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "FilterMutectCalls' --{held} is a second model this port does not read yet, \
                     and a run that ignored it would filter differently. This message is the \
                     port's own and not GATK's."
                ),
            ));
        }
    }

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "FilterMutectCalls")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let double = |name: &str, default: f64| -> f64 {
        scalar(parser, name)
            .and_then(|text| text.parse::<f64>().ok())
            .unwrap_or(default)
    };

    let defaults = EngineArguments::default();
    let list_defaults = FilterArguments::default();
    let mut list = FilterArguments {
        mitochondria: flag(parser, "mitochondria-mode"),
        microbial: flag(parser, "microbial-mode"),
        min_median_mapping_quality: number_or(
            parser,
            "min-median-mapping-quality",
            list_defaults.min_median_mapping_quality,
        ),
        log_snv_prior: double("log-snv-prior", list_defaults.log_snv_prior),
        log_indel_prior: double("log-indel-prior", list_defaults.log_indel_prior),
        read_orientation_priors: false,
    };
    let min_median_mapping_quality = list.min_median_mapping_quality();
    let engine = EngineArguments {
        list,
        min_median_base_quality: number_or(
            parser,
            "min-median-base-quality",
            defaults.min_median_base_quality,
        ),
        min_median_mapping_quality,
        long_indel_length: number_or(parser, "long-indel-length", defaults.long_indel_length),
        unique_alt_read_count: number_or(
            parser,
            "unique-alt-read-count",
            defaults.unique_alt_read_count,
        ),
        contamination_estimate: double("contamination-estimate", defaults.contamination_estimate),
        min_reads_on_each_strand: number_or(
            parser,
            "min-reads-per-strand",
            defaults.min_reads_on_each_strand,
        ),
        min_median_read_position: number_or(
            parser,
            "min-median-read-position",
            defaults.min_median_read_position,
        ),
        min_af: double("min-allele-fraction", defaults.min_af),
        normal_pileup_p_value_threshold: double(
            "normal-p-value-threshold",
            defaults.normal_pileup_p_value_threshold,
        ),
        n_ratio: double("max-n-ratio", defaults.n_ratio),
        max_events_in_region: number_or(
            parser,
            "max-events-in-region",
            defaults.max_events_in_region,
        ),
        max_events_in_haplotype: number_or(
            parser,
            "max-events-in-haplotype",
            defaults.max_events_in_haplotype,
        ),
        num_alt_alleles_threshold: number_or(
            parser,
            "max-alt-allele-count",
            defaults.num_alt_alleles_threshold as i32,
        )
        .max(0) as usize,
        max_median_fragment_length_difference: number_or(
            parser,
            "max-median-fragment-length-difference",
            defaults.max_median_fragment_length_difference,
        ),
        min_slippage_length: number_or(parser, "min-slippage-length", defaults.min_slippage_length),
        slippage_rate: double("pcr-slippage-rate", defaults.slippage_rate),
        max_distance_to_filtered_call_on_same_haplotype: number_or(
            parser,
            "distance-on-haplotype",
            defaults.max_distance_to_filtered_call_on_same_haplotype,
        ),
    };
    let tool_defaults = fmc::ToolArguments::default();
    let strategy = match scalar(parser, "threshold-strategy").as_deref() {
        Some("CONSTANT") => Strategy::Constant,
        Some("FALSE_DISCOVERY_RATE") => Strategy::FalseDiscoveryRate,
        _ => Strategy::OptimalFScore,
    };

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;

    // `onTraversalStart`'s header: the input's lines less `filtering_status`, a new
    // `filtering_status`, the STRQ and AS_FilterStatus INFO lines, every Mutect FILTER line and the
    // default tool lines, in one set the writer sorts.
    let mut header = file.header.clone();
    header.lines.retain(
        |line| !matches!(line, HeaderLine::Unstructured { key, .. } if key == "filtering_status"),
    );
    let mut added: Vec<HeaderLine> = vec![HeaderLine::Unstructured {
        key: "filtering_status".to_string(),
        value: "These calls have been filtered by FilterMutectCalls to label false positives \
                with a list of failed filters and true positives with PASS."
            .to_string(),
    }];
    let parsed_line = |text: &str| -> Option<HeaderLine> {
        htsjdk_vcf::reader::read_vcf(&format!(
            "##fileformat=VCFv4.2\n##{text}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        ))
        .ok()
        .and_then(|parsed| {
            parsed.header.lines.into_iter().find(
                |line| !matches!(line, HeaderLine::Unstructured { key, .. } if key == "fileformat"),
            )
        })
    };
    for (_, line) in gatk_engine::mutect_filter_list::INFO_LINES {
        added.extend(parsed_line(line));
    }
    for name in gatk_engine::mutect_filter_list::MUTECT_FILTER_NAMES {
        if let Some(line) = gatk_engine::mutect_filter_list::filter_line(name) {
            added.extend(parsed_line(&line));
        }
    }
    added.extend(default_tool_vcf_header_lines(parser, "FilterMutectCalls"));
    for line in added {
        if !header.lines.contains(&line) {
            header.lines.push(line);
        }
    }
    let write_records = |records: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut records = records.to_vec();
        apply_sites_only(parser, &mut header, &mut records);
        let out = write_vcf_honouring_lenient(parser, &header, &records)?;
        write_variant_output(parser, &output, &out)
    };

    // `new File(statsTable == null ? drivingVariantFile + ".stats" : statsTable)`, looked for
    // after the header is written.
    let stats_path = argument(parser, "stats").unwrap_or_else(|| format!("{input}.stats"));
    let Ok(stats) = std::fs::read_to_string(&stats_path) else {
        write_records(&[])?;
        let missing = fmc::MissingStatsTable { path: stats_path };
        return Err(Thrown {
            failure: Failure::User,
            exception: missing.class(),
            message: Some(missing.message()),
        });
    };
    // `MutectStats.readFromFile`, of which the clustering model reads `callable` alone.
    let callable_sites = stats.lines().skip(1).find_map(|line| {
        let mut fields = line.split('\t');
        match (fields.next(), fields.next()) {
            (Some("callable"), Some(value)) => value.trim().parse::<f64>().ok(),
            _ => None,
        }
    });

    let normal_samples: Vec<String> = file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Unstructured { key, value } if key == "normal_sample" => {
                Some(value.clone())
            }
            _ => None,
        })
        .collect();

    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let text_of = |value: &Value| value.format().unwrap_or_default();
    let values_of = |value: &Value| -> Vec<String> {
        match value {
            Value::List(items) => items.iter().map(text_of).collect(),
            other => text_of(other).split(',').map(str::to_string).collect(),
        }
    };
    let info = |record: &htsjdk_vcf::variant::VariantContext, key: &str| -> Option<Vec<String>> {
        record
            .attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| values_of(value))
    };
    let doubles = |record: &htsjdk_vcf::variant::VariantContext, key: &str| -> Option<Vec<f64>> {
        info(record, key).map(|values| {
            values
                .iter()
                .map(|value| value.trim().parse::<f64>().unwrap_or(0.0))
                .collect()
        })
    };
    let ints = |record: &htsjdk_vcf::variant::VariantContext, key: &str| -> Option<Vec<i32>> {
        info(record, key).map(|values| {
            values
                .iter()
                .map(|value| value.trim().parse::<i32>().unwrap_or(0))
                .collect()
        })
    };

    let mut records: Vec<Record> = Vec::new();
    for variant in &kept {
        variant.genotypes.decode();
        let reference = variant.reference();
        let reference_length = reference.len() as i32;
        let alternates: Vec<AccumulationAllele> = variant
            .alternate_alleles()
            .iter()
            .map(|allele| AccumulationAllele {
                allele: AlternateAllele {
                    length: if allele.is_symbolic() {
                        0
                    } else {
                        allele.len() as i32
                    },
                    symbolic: allele.is_symbolic(),
                },
                non_ref: allele.display_string() == "<NON_REF>",
            })
            .collect();
        // `getIndelLengths`, which answers only for an INDEL or MIXED record.
        let kinds: Vec<&str> = variant
            .alternate_alleles()
            .iter()
            .map(|allele| {
                if allele.is_symbolic() {
                    "SYMBOLIC"
                } else if allele.len() == reference.len() {
                    if allele.len() == 1 {
                        "SNP"
                    } else {
                        "MNP"
                    }
                } else {
                    "INDEL"
                }
            })
            .collect();
        let mixed = kinds.windows(2).any(|pair| pair[0] != pair[1]);
        let indel_lengths = if mixed || kinds.first() == Some(&"INDEL") {
            Some(
                alternates
                    .iter()
                    .map(|alternate| alternate.allele.length - reference_length)
                    .collect(),
            )
        } else {
            None
        };
        let mut genotypes = Vec::new();
        let mut allele_fractions = Vec::new();
        let mut phasing = Vec::new();
        for genotype in variant.genotypes.iter() {
            genotypes.push(GenotypeData {
                tumor: !normal_samples.contains(&genotype.sample_name),
                allele_depths: genotype.ad.clone().unwrap_or_default(),
                values: Vec::new(),
            });
            allele_fractions.push(
                genotype
                    .get("AF")
                    .map(|value| {
                        values_of(value)
                            .iter()
                            .map(|text| text.trim().parse::<f64>().unwrap_or(0.0))
                            .collect()
                    })
                    .unwrap_or_default(),
            );
            let field = |key: &str| {
                genotype
                    .get(key)
                    .map(text_of)
                    .filter(|text| !text.is_empty() && text != ".")
            };
            phasing.push((field("PGT"), field("PID")));
        }
        let single_int = |key: &str| ints(variant, key).and_then(|values| values.first().copied());
        records.push(Record {
            start: variant.start as i32,
            reference_length,
            alternates,
            genotypes,
            allele_fractions,
            phasing,
            tumor_log_10_odds: doubles(variant, "TLOD"),
            normal_artifact_log_10_odds: doubles(variant, "NALOD"),
            normal_log_10_odds: doubles(variant, "NLOD"),
            population_af: doubles(variant, "POPAF"),
            median_base_quality: ints(variant, "MBQ"),
            median_mapping_quality: ints(variant, "MMQ"),
            median_fragment_length: ints(variant, "MFRL"),
            median_read_position: ints(variant, "MPOS"),
            unique_alt_read_count: ints(variant, "AS_UNIQ_ALT_READ_COUNT"),
            strand_bias_table: info(variant, "AS_SB_TABLE").map(|values| values.join(",")),
            n_count: single_int("NCount"),
            event_count_in_region: single_int("ECNT"),
            event_count_in_haplotype: single_int("ECNTH"),
            repeats_per_allele: info(variant, "RPA"),
            repeat_unit: info(variant, "RU").map(|values| values.join(",")),
            in_panel_of_normals: variant.attributes.iter().any(|(key, _)| key == "PON"),
            indel_lengths,
        });
    }

    let arguments = fmc::ToolArguments {
        engine,
        strategy,
        initial_posterior_threshold: double(
            "initial-threshold",
            tool_defaults.initial_posterior_threshold,
        ),
        max_false_discovery_rate: double(
            "false-discovery-rate",
            tool_defaults.max_false_discovery_rate,
        ),
        f_score_beta: double("f-score-beta", tool_defaults.f_score_beta),
        callable_sites,
        log_artifact_prior: double("log-artifact-prior", tool_defaults.log_artifact_prior),
    };
    let result = match fmc::run(&records, &arguments) {
        Ok(result) => result,
        Err(error) => {
            write_records(&[])?;
            return Err(Thrown::non_user(error.class, error.message));
        }
    };

    let mut written = Vec::new();
    for (variant, applied) in kept.iter().zip(&result.records) {
        let mut out = (*variant).clone();
        let mut filters = applied.filters.clone();
        filters.sort();
        filters.dedup();
        out.filters = Some(filters);
        out.attributes.retain(|(key, _)| key != "AS_FilterStatus");
        out.attributes.push((
            "AS_FilterStatus".to_string(),
            Value::Str(applied.as_filter_status.clone()),
        ));
        for (key, quality) in &applied.annotations {
            out.attributes.retain(|(name, _)| name != key);
            out.attributes
                .push((key.clone(), Value::Int(i64::from(*quality))));
        }
        written.push(out);
    }
    write_records(&written)?;
    let stats_output = argument(parser, "filtering-stats")
        .unwrap_or_else(|| format!("{output}.filteringStats.tsv"));
    write_file(&stats_output, result.filtering_stats.as_bytes())?;
    Ok(None)
}

/// `ComposeSTRTableFile`, a `GATKTool` that scans the reference for tandem repeats and writes a zip.
///
/// The scan and the decimation are [`gatk_tools::compose_str_table`]'s. What the runner adds is the
/// zip `STRTableFileBuilder.store` writes: `reference.dict` (the best available dictionary, under an
/// `@HD` line of its own), `decimation.txt`, `sites.bin` and `sites.idx`, `sites.txt` when asked for,
/// and `summary.txt`, whose annotations carry the tool's command line.
///
/// The intervals are merged again under `ALL`, whatever `--interval-merging-rule` said, and are
/// traversed contig by contig in the dictionary's order. The members are written stored rather
/// than deflated and in name order; the covering array compares a zip member by member.
pub fn compose_str_table_file(parser: &Parser) -> Outcome {
    use gatk_tools::compose_str_table as str_table;

    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let spec = scalar(parser, "decimation").unwrap_or_else(|| "DEFAULT".to_string());
    let table = if spec.eq_ignore_ascii_case("NONE") {
        str_table::DecimationTable::none()
    } else if spec.eq_ignore_ascii_case("DEFAULT") {
        str_table::DecimationTable::default_table()
    } else {
        let text = std::fs::read_to_string(&spec)
            .map_err(|_| Thrown::user(format!("Couldn't read file {spec}. Error was: {spec}")))?;
        str_table::DecimationTable::parse(&text, &spec).map_err(|error| match error {
            str_table::DecimationError::BadInput(message) => bad_input(message),
        })?
    };
    let settings = str_table::Settings {
        max_period: number_or(parser, "max-period", 8).max(1) as usize,
        max_repeat: number_or(parser, "max-repeats", 20).max(1) as usize,
    };

    let _ = resolve_read_filters(parser, "ComposeSTRTableFile")?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    // `getBestAvailableSequenceDictionary`: the master one, else the reference's own `.dict`.
    let dictionary = match master_dictionary(parser)? {
        Some(master) => master,
        None => reference_dictionary(parser)?.unwrap_or_default(),
    };
    let intervals = interval_arguments(parser, &dictionary)?.map(|parameters| parameters.intervals);

    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();
    for sequence in &dictionary.sequences {
        let bases = reference
            .query(&sequence.name, 1, sequence.length)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        contigs.push((sequence.name.clone(), bases));
    }
    // `sortAndMergeIntervals(..., IntervalMergingRule.ALL)`: overlapping and adjacent intervals
    // become one, per contig.
    let mut chosen: Vec<(String, i64, i64)> = Vec::new();
    if let Some(intervals) = &intervals {
        let mut sorted: Vec<(usize, i64, i64, String)> = intervals
            .iter()
            .filter_map(|interval| {
                dictionary
                    .sequences
                    .iter()
                    .position(|sequence| sequence.name == interval.contig)
                    .map(|index| {
                        (
                            index,
                            i64::from(interval.start),
                            i64::from(interval.end),
                            interval.contig.clone(),
                        )
                    })
            })
            .collect();
        sorted.sort();
        for (_, start, end, contig) in sorted {
            match chosen.last_mut() {
                Some((last_contig, _, last_end))
                    if *last_contig == contig && start <= *last_end + 1 =>
                {
                    *last_end = (*last_end).max(end);
                }
                _ => chosen.push((contig, start, end)),
            }
        }
    }
    let scan = if intervals.is_some() && chosen.is_empty() {
        str_table::Scan::default()
    } else {
        str_table::scan(&contigs, &chosen, settings, &table)
    };

    let names: Vec<String> = dictionary
        .sequences
        .iter()
        .map(|sequence| sequence.name.clone())
        .collect();
    let dictionary_text = htsjdk_bam::header::SamHeader {
        sequences: dictionary.sequences.clone(),
        ..htsjdk_bam::header::SamHeader::default()
    }
    .encode_replacing_version();
    let (sites, index) = str_table::sites_binary(&scan.emitted);
    let annotations = vec![(
        "commandLine".to_string(),
        crate::command_line::expanded("ComposeSTRTableFile", parser),
    )];
    let summary = str_table::summary(
        &scan,
        settings,
        &annotations,
        &table,
        gatk_barclay::java_double_to_string,
    );
    let mut members: Vec<(String, Vec<u8>)> = vec![
        ("decimation.txt".to_string(), table.print().into_bytes()),
        ("reference.dict".to_string(), dictionary_text.into_bytes()),
        ("sites.bin".to_string(), sites),
        ("sites.idx".to_string(), index),
        ("summary.txt".to_string(), summary.into_bytes()),
    ];
    if flag(parser, "generate-sites-text-output") {
        members.push((
            "sites.txt".to_string(),
            str_table::sites_text(&scan.emitted, &names).into_bytes(),
        ));
    }
    members.sort_by(|left, right| left.0.cmp(&right.0));
    write_file(&output, &stored_zip(&members))?;
    Ok(None)
}

/// A zip of `members`, each stored rather than deflated: local headers, the central directory and
/// its end record, with no time on any entry.
fn stored_zip(members: &[(String, Vec<u8>)]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut central = Vec::new();
    for (name, data) in members {
        let offset = out.len() as u32;
        let mut crc = flate2::Crc::new();
        crc.update(data);
        let crc = crc.sum();
        let mut header = Vec::new();
        header.extend_from_slice(&20u16.to_le_bytes()); // version needed
        header.extend_from_slice(&0u16.to_le_bytes()); // flags
        header.extend_from_slice(&0u16.to_le_bytes()); // stored
        header.extend_from_slice(&0u16.to_le_bytes()); // time
        header.extend_from_slice(&0x21u16.to_le_bytes()); // date: 1980-01-01
        header.extend_from_slice(&crc.to_le_bytes());
        header.extend_from_slice(&(data.len() as u32).to_le_bytes());
        header.extend_from_slice(&(data.len() as u32).to_le_bytes());
        header.extend_from_slice(&(name.len() as u16).to_le_bytes());
        header.extend_from_slice(&0u16.to_le_bytes()); // extra
        out.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        out.extend_from_slice(&header);
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(data);
        central.extend_from_slice(&0x0201_4b50u32.to_le_bytes());
        central.extend_from_slice(&20u16.to_le_bytes()); // version made by
        central.extend_from_slice(&header);
        central.extend_from_slice(&0u16.to_le_bytes()); // comment
        central.extend_from_slice(&0u16.to_le_bytes()); // disk
        central.extend_from_slice(&0u16.to_le_bytes()); // internal attributes
        central.extend_from_slice(&0u32.to_le_bytes()); // external attributes
        central.extend_from_slice(&offset.to_le_bytes());
        central.extend_from_slice(name.as_bytes());
    }
    let central_offset = out.len() as u32;
    out.extend_from_slice(&central);
    out.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&(members.len() as u16).to_le_bytes());
    out.extend_from_slice(&(members.len() as u16).to_le_bytes());
    out.extend_from_slice(&(central.len() as u32).to_le_bytes());
    out.extend_from_slice(&central_offset.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out
}

/// `CalibrateDragstrModel`, a `GATKTool` that piles the reads up over the STR table's sites and
/// estimates the DRAGstr parameters from them.
///
/// The estimation and the file layout are [`gatk_tools::calibrate_dragstr_model`]'s. The runner
/// adds the traversal the tool does itself, interval by interval: the sites whose START falls in
/// the interval, in the table's order, and for each one the reads overlapping it that are neither
/// unmapped, secondary nor QC-failed, with `XQ` read as the mapping quality. Then the downsampling by
/// the table's decimation bits, the qualifying filter, and either the estimate or Illumina's
/// defaults, which is decided by a minimum count per period and repeat length.
///
/// `--parallel` and `--threads` are refused: the parallel collection merges shards in an order the
/// sequential one does not, and that order reaches both the estimate and the sites file.
pub fn calibrate_dragstr_model(parser: &Parser) -> Outcome {
    use gatk_tools::calibrate_dragstr_model as cdm;

    if flag(parser, "parallel") || number_or(parser, "threads", 0) > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "CalibrateDragstrModel's parallel collection is not ported. This message is the \
             port's own and not GATK's.",
        ));
    }
    let bad_value = |name: &str, message: String| {
        Thrown::command_line(format!("Argument {name} has a bad value: {message}"))
    };
    let sequence = |name: &str, default: &str| -> Result<(Vec<f64>, bool), Thrown> {
        let text = scalar(parser, name).unwrap_or_else(|| default.to_string());
        let set = scalar(parser, name).is_some();
        cdm::double_sequence(&text)
            .map(|values| (values, set))
            .map_err(|message| bad_value(name, message))
    };
    let (gp, gp_set) = sequence("gp-values", "10:1.0:50")?;
    let (api, api_set) = sequence("api-values", "0:1.0:40")?;
    let (gop, _) = sequence("gop-values", "10:.25:50")?;
    let parameters = cdm::HyperParameters {
        phred_gp_values: gp,
        phred_api_values: api,
        phred_gop_values: gop,
        het_to_hom_ratio: scalar(parser, "het-to-hom-ratio")
            .and_then(|text| text.parse().ok())
            .unwrap_or(2.0),
        min_loci_count: number_or(parser, "min-loci-count", 50).max(0) as usize,
        api_mono_threshold: f64::from(number_or(parser, "api-mono-threshold", 3)),
        max_period: number_or(parser, "max-period", 8).max(1) as usize,
        max_repeat_length: number_or(parser, "max-repeats", 20).max(1) as usize,
    };
    // `DragstrHyperParameters.validate`, an else-if chain: a GP sequence given on the command line
    // is checked and nothing after it is.
    let not_phred = |values: &[f64]| values.iter().find(|d| !d.is_finite() || **d < 0.0).copied();
    if gp_set {
        if let Some(d) = not_phred(&parameters.phred_gp_values) {
            return Err(bad_value(
                "gp-values",
                format!(
                    "Not a valid Phred value: {}",
                    gatk_barclay::java_double_to_string(d)
                ),
            ));
        }
    } else if api_set {
        if let Some(d) = not_phred(&parameters.phred_api_values) {
            return Err(bad_value(
                "api-values",
                format!(
                    "Not a valid Phred value: {}",
                    gatk_barclay::java_double_to_string(d)
                ),
            ));
        }
    } else if !parameters.het_to_hom_ratio.is_finite() || parameters.het_to_hom_ratio <= 0.0 {
        return Err(bad_value(
            "het-to-hom-ratio",
            format!(
                "must be finite and greater than 0 but found {}",
                gatk_barclay::java_double_to_string(parameters.het_to_hom_ratio)
            ),
        ));
    } else if number_or(parser, "min-loci-count", 50) < 1 {
        return Err(bad_value(
            "min-loci-count",
            format!(
                "must be greater than 0 but found {}",
                number_or(parser, "min-loci-count", 50)
            ),
        ));
    }
    let min_depth = number_or(parser, "minimum-depth", 10);
    let padding = i64::from(number_or(parser, "pileup-padding", 5));
    let min_mq = number_or(parser, "sampling-min-mq", 20);
    let downsample_size = number_or(parser, "down-sample-size", 4096).max(0) as usize;
    let force = flag(parser, "force-estimation");
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let table_path = argument(parser, "str-table-path").ok_or_else(|| {
        Thrown::command_line(
            "Argument str-table-path was missing: Argument 'str-table-path' is required",
        )
    })?;
    let sites_output = scalar(parser, "debug-sites-output");

    let ReadWalkerStart {
        source,
        header,
        intervals,
        ..
    } = read_walker_startup(parser, "CalibrateDragstrModel")?;
    // `getBestAvailableSequenceDictionary`: the master, else the reference, else the reads.
    let dictionary: Vec<htsjdk_bam::header::SequenceRecord> = match master_dictionary(parser)? {
        Some(master) => master.sequences,
        None => match reference_dictionary(parser)? {
            Some(reference) => reference.sequences,
            None => header.sequences.clone(),
        },
    };
    let read_groups: Vec<String> = header
        .read_groups
        .iter()
        .map(|group| group.id.clone())
        .collect();
    let mut samples: Vec<Option<String>> = Vec::new();
    for group in &header.read_groups {
        let sample = group.attributes.get("SM").map(str::to_string);
        if !samples.contains(&sample) {
            samples.push(sample);
        }
    }
    if samples.len() > 1 {
        let names: Vec<String> = samples
            .iter()
            .map(|s| s.clone().unwrap_or_else(|| "null".to_string()))
            .collect();
        return Err(Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            format!(
                "the input alignment(s) have more than one sample: {}",
                names.join(", ")
            ),
        ));
    }
    let sample = samples.first().cloned().flatten();

    // The sites file is opened before the table is.
    let mut sites_lines: Vec<String> = Vec::new();
    if let Some(path) = &sites_output {
        write_file(path, b"")?;
    }

    let zip = std::fs::read(&table_path).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file {table_path}. Error was: {table_path}"
        ))
    })?;
    let members = read_zip(&zip).map_err(|message| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            message,
        )
    })?;
    let member = |name: &str| -> Result<&Vec<u8>, Thrown> {
        members.get(name).ok_or_else(|| {
            Thrown::non_user(
                "org.broadinstitute.hellbender.exceptions.GATKException",
                format!("str-table-file {table_path} is missing {name}"),
            )
        })
    };
    let table_dictionary =
        htsjdk_bam::reader::parse_header_text(&String::from_utf8_lossy(member("reference.dict")?))
            .sequences;
    let decimation = gatk_tools::compose_str_table::DecimationTable::parse(
        &String::from_utf8_lossy(member("decimation.txt")?),
        "decimation.txt",
    )
    .map_err(|error| match error {
        gatk_tools::compose_str_table::DecimationError::BadInput(message) => bad_input(message),
    })?;
    let sites = cdm::read_sites(member("sites.bin")?);

    match gatk_tools::sequence_dictionary::compare(&dictionary, &table_dictionary, false) {
        gatk_tools::sequence_dictionary::Compatibility::Identical
        | gatk_tools::sequence_dictionary::Compatibility::Superset
        | gatk_tools::sequence_dictionary::Compatibility::NonCanonicalHumanOrder
        | gatk_tools::sequence_dictionary::Compatibility::OutOfOrder => {}
        other => {
            return Err(Thrown::non_user(
                "org.broadinstitute.hellbender.exceptions.GATKException",
                format!(
                    "the reference and str-table sequence dictionary are incompatible: {}",
                    other.name()
                ),
            ));
        }
    }

    // `getTraversalIntervals`: the ones given, else every contig of the best dictionary.
    let traversal: Vec<gatk_engine::interval::SimpleInterval> = if intervals.is_empty() {
        dictionary
            .iter()
            .filter_map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
            })
            .collect()
    } else {
        intervals
    };
    let mut all = cdm::Stratified::new(parameters.max_period, parameters.max_repeat_length);
    for interval in &traversal {
        let Some(contig_index) = table_dictionary
            .iter()
            .position(|s| s.name == interval.contig)
        else {
            continue;
        };
        let contig_length = dictionary
            .iter()
            .find(|s| s.name == interval.contig)
            .map_or(0, |s| i64::from(s.length));
        let records =
            gatk_tools::read_walker::traverse(&source, std::slice::from_ref(interval), &|_| true)
                .map_err(reads_traversal_error)?;
        let reads: Vec<cdm::PileRead> = records
            .iter()
            .filter(|read| {
                read.flags & (0x4 | 0x100 | 0x200) == 0
                    && read.alignment_start <= read.alignment_end()
            })
            .map(|read| {
                let xq = match read.tags.get(htsjdk_bam::tag::Tag::new(b"XQ")) {
                    Some(htsjdk_bam::tag::TagValue::Int(value)) => Some(*value as i32),
                    _ => None,
                };
                cdm::PileRead {
                    start: i64::from(read.alignment_start),
                    end: i64::from(read.alignment_end()),
                    mapping_quality: xq.unwrap_or(i32::from(read.mapping_quality)),
                    supplementary: read.flags & 0x800 != 0,
                    cigar: read
                        .cigar
                        .elements
                        .iter()
                        .map(|element| (element.op.to_char(), i64::from(element.length)))
                        .collect(),
                }
            })
            .collect();
        for site in sites.iter().filter(|site| {
            site.contig == contig_index as i32
                && site.start >= i64::from(interval.start)
                && site.start <= i64::from(interval.end)
        }) {
            let overlapping: Vec<&cdm::PileRead> = reads
                .iter()
                .filter(|read| read.start <= site.end() && read.end >= site.start)
                .collect();
            let case = cdm::collect(site, &overlapping, padding, contig_length);
            all.add(case).map_err(|(index, length)| {
                Thrown::non_user(
                    "java.lang.ArrayIndexOutOfBoundsException",
                    format!("Index {index} out of bounds for length {length}"),
                )
            })?;
        }
    }

    // `downSample`, combination by combination.
    let mut kept = cdm::Stratified::new(parameters.max_period, parameters.max_repeat_length);
    for period in 1..=parameters.max_period {
        for repeats in 1..=parameters.max_repeat_length {
            let bit = decimation.decimation_bit(period, repeats).max(0) as usize;
            let cell = all.cells[period - 1][repeats - 1].clone();
            let survivors = cdm::downsample_cell(&cell, bit, downsample_size, &mut sites_lines)
                .map_err(|(index, length)| {
                    Thrown::non_user(
                        "java.lang.ArrayIndexOutOfBoundsException",
                        format!("Index {index} out of bounds for length {length}"),
                    )
                })?;
            for case in survivors {
                let _ = kept.add(case);
            }
        }
    }
    let final_sites = kept.qualifying(min_depth, min_mq, 0);
    if sites_output.is_some() {
        for row in &kept.cells {
            for cell in row {
                for case in cell {
                    let fate = if case.qualifies(min_depth, min_mq, 0) {
                        "used"
                    } else {
                        "skipped"
                    };
                    sites_lines.push(case.line(fate));
                }
            }
        }
    }

    // `isThereEnoughCases` answers true under `--force-estimation` whatever the counts, so the file
    // says `estimated` and `estimatedByForce` is never written.
    let enough = force
        || cdm::enough_cases(
            &final_sites,
            parameters.max_period,
            parameters.max_repeat_length,
        );
    let using_defaults = !enough && !force;
    let annotations = vec![
        (
            "sample",
            sample.unwrap_or_else(|| "<unspecified>".to_string()),
        ),
        (
            "readGroups",
            if read_groups.is_empty() {
                "<unspecified>".to_string()
            } else {
                read_groups.join(", ")
            },
        ),
        (
            "estimatedOrDefaults",
            if using_defaults {
                "defaults"
            } else if enough {
                "estimated"
            } else {
                "estimatedByForce"
            }
            .to_string(),
        ),
        (
            "commandLine",
            crate::command_line::expanded("CalibrateDragstrModel", parser),
        ),
    ];
    let text = if using_defaults {
        cdm::params_file(&annotations, 20, &cdm::default_rows())
    } else {
        let rows = cdm::estimate(&parameters, &final_sites.as_cases());
        cdm::params_file(&annotations, parameters.max_repeat_length, &rows)
    };
    if let Some(path) = &sites_output {
        let mut body = sites_lines.join("\n");
        if !sites_lines.is_empty() {
            body.push('\n');
        }
        write_file(path, body.as_bytes())?;
    }
    write_file(&output, text.as_bytes())?;
    Ok(None)
}

/// The members of a zip, by name, from its central directory: stored or deflated.
fn read_zip(bytes: &[u8]) -> Result<std::collections::HashMap<String, Vec<u8>>, String> {
    use std::io::Read;
    let u16_at = |at: usize| -> Option<usize> {
        bytes
            .get(at..at + 2)
            .map(|b| u16::from_le_bytes([b[0], b[1]]) as usize)
    };
    let u32_at = |at: usize| -> Option<usize> {
        bytes
            .get(at..at + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
    };
    let end = (0..bytes.len().saturating_sub(21))
        .rev()
        .find(|&at| u32_at(at) == Some(0x0605_4b50))
        .ok_or_else(|| "not a zip file".to_string())?;
    let count = u16_at(end + 10).ok_or("truncated zip")?;
    let mut at = u32_at(end + 16).ok_or("truncated zip")?;
    let mut members = std::collections::HashMap::new();
    for _ in 0..count {
        if u32_at(at) != Some(0x0201_4b50) {
            return Err("a broken central directory".to_string());
        }
        let method = u16_at(at + 10).ok_or("truncated zip")?;
        let compressed = u32_at(at + 20).ok_or("truncated zip")?;
        let name_length = u16_at(at + 28).ok_or("truncated zip")?;
        let extra_length = u16_at(at + 30).ok_or("truncated zip")?;
        let comment_length = u16_at(at + 32).ok_or("truncated zip")?;
        let local = u32_at(at + 42).ok_or("truncated zip")?;
        let name = String::from_utf8_lossy(
            bytes
                .get(at + 46..at + 46 + name_length)
                .ok_or("truncated zip")?,
        )
        .into_owned();
        at += 46 + name_length + extra_length + comment_length;
        let local_name = u16_at(local + 26).ok_or("truncated zip")?;
        let local_extra = u16_at(local + 28).ok_or("truncated zip")?;
        let data_start = local + 30 + local_name + local_extra;
        let data = bytes
            .get(data_start..data_start + compressed)
            .ok_or("truncated zip")?;
        let content = match method {
            0 => data.to_vec(),
            8 => {
                let mut out = Vec::new();
                flate2::read::DeflateDecoder::new(data)
                    .read_to_end(&mut out)
                    .map_err(|error| error.to_string())?;
                out
            }
            other => return Err(format!("zip method {other} is not supported")),
        };
        members.insert(name, content);
    }
    Ok(members)
}

/// `LearnReadOrientationModel`, a `CommandLineProgram` over one or more `CollectF1R2Counts` archives.
///
/// The EM is [`gatk_tools::learn_read_orientation_model::learn_prior`]'s. The runner is the rest of
/// `doWork`: every archive's reference and depth-one histograms summed per sample, its alt tables
/// gathered per sample in input order, each canonical context combined with its reverse
/// complement, and one prior table per sample written into a `.tar.gz`. A context is skipped, and
/// keeps its flat prior with no examples, when its reference histogram is empty or it has no alt
/// site of either strand.
pub fn learn_read_orientation_model(parser: &Parser) -> Outcome {
    use gatk_tools::learn_read_orientation_model as lrom;
    use std::collections::{BTreeMap, HashMap};

    let inputs = arguments(parser, "input");
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let absolute_output = java_absolute_path(&output);
    if !absolute_output.ends_with(".tar.gz") {
        return Err(Thrown::user(format!(
            "Couldn't write file {absolute_output} because Output file must end in .tar.gz"
        )));
    }
    let threshold = scalar(parser, "convergence-threshold")
        .and_then(|text| text.parse::<f64>().ok())
        .unwrap_or(1e-4);
    let max_iterations = number_or(parser, "num-em-iterations", 20);
    let max_depth = number_or(parser, "max-depth", 200);

    // Every archive, extracted: the histograms by sample and the alt tables in input order.
    type Histograms = LabelledHistograms;
    let mut ref_by_sample: HashMap<String, Vec<Histograms>> = HashMap::new();
    let mut alt_by_sample: HashMap<String, Vec<Histograms>> = HashMap::new();
    let mut records_by_sample: HashMap<String, Vec<lrom::AltSite>> = HashMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    for input in &inputs {
        let bytes = std::fs::read(input)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{input}: {error}")))?;
        let tar = {
            use std::io::Read;
            let mut out = Vec::new();
            flate2::read::MultiGzDecoder::new(bytes.as_slice())
                .read_to_end(&mut out)
                .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{input}: {error}")))?;
            out
        };
        let mut members = tar_members(&tar);
        members.sort_by(|left, right| left.0.cmp(&right.0));
        for (name, data) in members {
            let text = String::from_utf8_lossy(&data).into_owned();
            if name.ends_with(".ref_histogram") || name.ends_with(".alt_histogram") {
                let (sample, histograms) = parse_metrics_histograms(&text);
                let target = if name.ends_with(".ref_histogram") {
                    &mut ref_by_sample
                } else {
                    &mut alt_by_sample
                };
                target.entry(sample).or_default().push(histograms);
            } else if name.ends_with(".alt_table") {
                let (sample, records) = parse_alt_table(&text);
                if !sample_order.contains(&sample) {
                    sample_order.push(sample.clone());
                }
                records_by_sample.entry(sample).or_default().extend(records);
            }
        }
    }
    // `sumHistogramsFromFiles`: the first file's histograms, the others added bin by bin.
    let sum_files = |files: &[Histograms]| -> Histograms {
        let mut sum = files.first().cloned().unwrap_or_default();
        for other in files.iter().skip(1) {
            for (label, bins) in other {
                if let Some((_, target)) = sum.iter_mut().find(|(name, _)| name == label) {
                    for (depth, count) in bins {
                        *target.entry(*depth).or_insert(0.0) += count;
                    }
                }
            }
        }
        sum
    };

    let canonical = lrom::canonical_kmers();
    let order = gatk_tools::collect_f1r2_counts::ref_context_order();
    let mut files: Vec<(String, Vec<u8>)> = Vec::new();
    for sample in &sample_order {
        let refs = sum_files(ref_by_sample.get(sample).map_or(&[][..], Vec::as_slice));
        let alts = sum_files(alt_by_sample.get(sample).map_or(&[][..], Vec::as_slice));
        let records = &records_by_sample[sample];
        let mut priors: HashMap<String, lrom::Prior> = HashMap::new();
        for context in &canonical {
            let revcomp = lrom::reverse_complement(context);
            let ref_bins = |label: &str| -> BTreeMap<i32, f64> {
                refs.iter()
                    .find(|(name, _)| name == label)
                    .map(|(_, bins)| bins.clone())
                    .unwrap_or_default()
            };
            // `combineRefHistogramWithRC`: the context's bins with the reverse complement's.
            let mut combined: BTreeMap<i32, f64> =
                (1..=max_depth).map(|depth| (depth, 0.0)).collect();
            let forward = ref_bins(context);
            let backward = ref_bins(&revcomp);
            for (depth, count) in &forward {
                *combined.entry(*depth).or_insert(0.0) +=
                    count + backward.get(depth).unwrap_or(&0.0);
            }
            // `combineAltDepthOneHistogramWithRC`, alt bases in order, F1R2 before F2R1.
            let mut combined_alts = Vec::new();
            let middle = context.as_bytes()[1];
            for alt in [lrom::Base::A, lrom::Base::C, lrom::Base::G, lrom::Base::T] {
                if alt.name().as_bytes()[0] == middle {
                    continue;
                }
                for f1r2 in [true, false] {
                    let label = format!(
                        "{context}_{}_{}",
                        alt.name(),
                        if f1r2 { "F1R2" } else { "F2R1" }
                    );
                    let other = format!(
                        "{revcomp}_{}_{}",
                        alt.complement().name(),
                        if f1r2 { "F2R1" } else { "F1R2" }
                    );
                    let find = |label: &str| {
                        alts.iter()
                            .find(|(name, _)| name == label)
                            .map(|(_, bins)| bins.clone())
                            .unwrap_or_default()
                    };
                    let (a, b) = (find(&label), find(&other));
                    let mut counts: BTreeMap<i32, f64> =
                        (1..=max_depth).map(|depth| (depth, 0.0)).collect();
                    for (depth, count) in &a {
                        *counts.entry(*depth).or_insert(0.0) +=
                            count + b.get(depth).unwrap_or(&0.0);
                    }
                    combined_alts.push(lrom::AltHistogram { alt, f1r2, counts });
                }
            }
            // `mergeDesignMatrices`: the context's own records, then the reverse complement's
            // turned around.
            let mut design: Vec<lrom::AltSite> = records
                .iter()
                .filter(|r| &r.context == context)
                .cloned()
                .collect();
            design.extend(
                records
                    .iter()
                    .filter(|r| r.context == revcomp)
                    .map(lrom::AltSite::reverse_complement),
            );
            let ref_sum: f64 = combined.values().sum();
            if ref_sum == 0.0 || design.is_empty() {
                continue;
            }
            let prior = lrom::learn_prior(
                context,
                &combined,
                &combined_alts,
                &design,
                threshold,
                max_iterations,
                max_depth,
            );
            priors.insert(revcomp.clone(), prior.reverse_complement());
            priors.insert(context.clone(), prior);
        }
        let columns = [
            "context",
            "rev_comp",
            "f1r2_a",
            "f1r2_c",
            "f1r2_g",
            "f1r2_t",
            "f2r1_a",
            "f2r1_c",
            "f2r1_g",
            "f2r1_t",
            "hom_ref",
            "germline_het",
            "somatic_het",
            "hom_var",
            "num_examples",
            "num_alt_examples",
        ];
        let rows: Vec<Vec<String>> = order
            .iter()
            .map(|kmer| {
                let prior = priors.get(kmer).cloned().unwrap_or_else(|| {
                    let middle = lrom::Base::from_byte(kmer.as_bytes()[1]).expect("a k-mer");
                    lrom::Prior {
                        context: kmer.clone(),
                        pi: lrom::flat_prior(middle),
                        examples: 0,
                        alt_examples: 0,
                    }
                });
                let mut row = vec![
                    prior.context.clone(),
                    lrom::reverse_complement(&prior.context),
                ];
                row.extend(
                    prior
                        .pi
                        .iter()
                        .map(|value| gatk_engine::tsv_table::java_double_to_string(*value)),
                );
                row.push(prior.examples.to_string());
                row.push(prior.alt_examples.to_string());
                row
            })
            .collect();
        let text = gatk_engine::tsv_table::write_table(&columns, &rows, &[("SAMPLE", sample)]);
        files.push((
            format!(
                "./{}.orientation_priors",
                gatk_tools::get_sample_name::url_encode_utf8(sample)
            ),
            text.into_bytes(),
        ));
    }
    files.sort_by(|left, right| left.0.cmp(&right.0));
    write_file(&output, &tar_gz(&files))?;
    Ok(Some("SUCCESS".to_string()))
}

/// Every regular entry of an uncompressed tar stream, by name, GNU long names followed.
fn tar_members(tar: &[u8]) -> Vec<(String, Vec<u8>)> {
    let field = |block: &[u8], from: usize, to: usize| -> String {
        let raw = &block[from..to];
        let end = raw.iter().position(|byte| *byte == 0).unwrap_or(raw.len());
        String::from_utf8_lossy(&raw[..end]).into_owned()
    };
    let mut out = Vec::new();
    let mut offset = 0;
    let mut long_name: Option<String> = None;
    while offset + 512 <= tar.len() {
        let block = &tar[offset..offset + 512];
        if block.iter().all(|byte| *byte == 0) {
            break;
        }
        let size = usize::from_str_radix(field(block, 124, 136).trim(), 8).unwrap_or(0);
        let data_start = offset + 512;
        let data = &tar[data_start..(data_start + size).min(tar.len())];
        match block[156] {
            b'L' => long_name = Some(field(data, 0, data.len())),
            b'0' | 0 => {
                let name = long_name.take().unwrap_or_else(|| field(block, 0, 100));
                out.push((name, data.to_vec()));
            }
            _ => {}
        }
        offset = data_start + size.div_ceil(512) * 512;
    }
    out
}

/// Labelled histograms: each one's value label and its counts by bin.
type LabelledHistograms = Vec<(String, std::collections::BTreeMap<i32, f64>)>;

/// A Picard metrics file's first header, which is the sample, and its histograms by label.
fn parse_metrics_histograms(text: &str) -> (String, LabelledHistograms) {
    let lines: Vec<&str> = text.lines().collect();
    let mut sample = String::new();
    for (index, line) in lines.iter().enumerate() {
        if line.starts_with("## htsjdk.samtools.metrics.StringHeader") {
            sample = lines
                .get(index + 1)
                .and_then(|value| value.strip_prefix("# "))
                .unwrap_or("")
                .to_string();
            break;
        }
    }
    let mut histograms = Vec::new();
    if let Some(start) = lines
        .iter()
        .position(|line| line.starts_with("## HISTOGRAM"))
    {
        if let Some(header) = lines.get(start + 1) {
            let labels: Vec<&str> = header.split('\t').skip(1).collect();
            histograms = labels
                .iter()
                .map(|label| (label.to_string(), std::collections::BTreeMap::new()))
                .collect();
            for line in lines.iter().skip(start + 2) {
                if line.trim().is_empty() {
                    break;
                }
                let mut cells = line.split('\t');
                let Some(depth) = cells.next().and_then(|d| d.parse::<i32>().ok()) else {
                    continue;
                };
                for (index, cell) in cells.enumerate() {
                    if let (Some((_, bins)), Ok(value)) =
                        (histograms.get_mut(index), cell.parse::<f64>())
                    {
                        bins.insert(depth, value);
                    }
                }
            }
        }
    }
    (sample, histograms)
}

/// An alt table's sample and rows.
fn parse_alt_table(
    text: &str,
) -> (
    String,
    Vec<gatk_tools::learn_read_orientation_model::AltSite>,
) {
    let mut sample = String::new();
    let mut rows = Vec::new();
    let mut seen_header = false;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("#<METADATA>SAMPLE=") {
            sample = rest.to_string();
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if !seen_header {
            seen_header = true;
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() < 7 {
            continue;
        }
        let number = |index: usize| cells[index].parse::<i32>().unwrap_or(0);
        if let Some(alt) = gatk_tools::learn_read_orientation_model::Base::parse(cells[6]) {
            rows.push(gatk_tools::learn_read_orientation_model::AltSite {
                context: cells[0].to_string(),
                ref_count: number(1),
                alt_count: number(2),
                ref_f1r2: number(3),
                alt_f1r2: number(4),
                alt,
            });
        }
    }
    (sample, rows)
}

/// `GeneExpressionEvaluation`, a read walker counting fragments over a gff3's features.
///
/// The counting is [`gatk_tools::gene_expression_evaluation::count`]'s. The runner is
/// `onTraversalStart` and the traversal: the grouping features of the gff3 overlapping each
/// traversal interval, each with the intervals of its descendants of the overlap types, keyed by
/// the feature shrunk to its label attribute so that one found twice is one feature; then the
/// reads through the tool's own ten default filters, the mapping quality one being the instance
/// EQUAL multi-mapping drops to zero.
pub fn gene_expression_evaluation(parser: &Parser) -> Outcome {
    use gatk_tools::gene_expression_evaluation as gee;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "GeneExpressionEvaluation")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let gff_path = argument(parser, "gff-file").ok_or_else(|| {
        Thrown::command_line("Argument gff-file was missing: Argument 'gff-file' is required")
    })?;
    let label = match scalar(parser, "feature-label-key").as_deref() {
        Some("ID") => gee::FeatureLabel::Id,
        _ => gee::FeatureLabel::Name,
    };
    let multi_overlap_method = match scalar(parser, "multi-overlap-method").as_deref() {
        Some("EQUAL") => gee::MultiOverlapMethod::Equal,
        _ => gee::MultiOverlapMethod::Proportional,
    };
    let multi_map_method = match scalar(parser, "multi-map-method").as_deref() {
        Some("EQUAL") => gee::MultiMapMethod::Equal,
        _ => gee::MultiMapMethod::Ignore,
    };
    let read_strands = match scalar(parser, "read-strands").as_deref() {
        Some("FORWARD_FORWARD") => gee::ReadStrands::ForwardForward,
        Some("REVERSE_FORWARD") => gee::ReadStrands::ReverseForward,
        Some("REVERSE_REVERSE") => gee::ReadStrands::ReverseReverse,
        _ => gee::ReadStrands::ForwardReverse,
    };
    let grouping: Vec<String> = {
        let given = arguments(parser, "grouping-type");
        if given.is_empty() {
            vec!["gene".to_string()]
        } else {
            given
        }
    };
    let overlap: Vec<String> = {
        let given = arguments(parser, "overlap-type");
        if given.is_empty() {
            vec!["exon".to_string()]
        } else {
            given
        }
    };
    let minimum = number_or(parser, "minimum-mapping-quality", 10);
    let settings = gee::Settings {
        multi_overlap_method,
        multi_map_method,
        read_strands,
        unspliced: flag(parser, "unspliced"),
        feature_label: label,
        minimum_mapping_quality: minimum,
        filter_mapping_quality: false,
    };
    let effective_minimum = settings.effective_minimum_mapping_quality();

    // `onTraversalStart`: the sample, one across every read group.
    let mut sample: Option<String> = None;
    for group in &header.read_groups {
        let this = group.attributes.get("SM").map(str::to_string);
        match &sample {
            None => sample = this,
            Some(first) => {
                if this.as_deref() != Some(first.as_str()) {
                    return Err(Thrown::non_user(
                        "org.broadinstitute.hellbender.exceptions.GATKException",
                        "Cannot run GeneExpressionEvaluation on multi-sample bam.",
                    ));
                }
            }
        }
    }
    let dictionary: Vec<htsjdk_bam::header::SequenceRecord> = match master_dictionary(parser)? {
        Some(master) => master.sequences,
        None => match reference_dictionary(parser)? {
            Some(reference) => reference.sequences,
            None => header.sequences.clone(),
        },
    };
    let all_intervals: Vec<gatk_engine::interval::SimpleInterval> = if intervals.is_empty() {
        dictionary
            .iter()
            .filter_map(|s| gatk_engine::interval::SimpleInterval::new(&s.name, 1, s.length))
            .collect()
    } else {
        intervals.clone()
    };

    let text = {
        let bytes = std::fs::read(&gff_path).map_err(|_| {
            Thrown::user(format!(
                "Couldn't read file {gff_path}. Error was: {gff_path}"
            ))
        })?;
        if bytes.len() > 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
            gunzip(&bytes).map_err(|error| Thrown::non_user(PORT_FAILURE, error.to_string()))?
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        }
    };
    let gff = gee::parse_gff3(&text).map_err(|message| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(message),
    })?;
    let mut features: Vec<gee::GroupingFeature> = Vec::new();
    for interval in &all_intervals {
        for (index, feature) in gff.iter().enumerate() {
            let base = &feature.base;
            if base.contig != interval.contig
                || base.start > interval.end
                || base.end < interval.start
                || !grouping.contains(&base.kind)
            {
                continue;
            }
            let overlaps: Vec<gee::Interval> = gee::descendants(&gff, index)
                .into_iter()
                .filter(|child| overlap.contains(&gff[*child].base.kind))
                .map(|child| gee::Interval {
                    contig: gff[child].base.contig.clone(),
                    start: gff[child].base.start,
                    end: gff[child].base.end,
                })
                .collect();
            // `shrinkBaseData`: only the label attribute survives.
            let mut shrunk = base.clone();
            shrunk.attributes.retain(|key, _| key == label.key());
            if label.value(&shrunk).is_none() {
                return Err(Thrown::user(format!(
                    "no geneid field {} found in feature at {}:{}-{}",
                    if label == gee::FeatureLabel::Id {
                        "ID"
                    } else {
                        "NAME"
                    },
                    shrunk.contig,
                    shrunk.start,
                    shrunk.end
                )));
            }
            match features.iter_mut().find(|known| known.base == shrunk) {
                Some(known) => {
                    for interval in overlaps {
                        if !known.overlaps.contains(&interval) {
                            known.overlaps.push(interval);
                        }
                    }
                }
                None => features.push(gee::GroupingFeature {
                    base: shrunk,
                    overlaps,
                }),
            }
        }
    }

    // The tool's own filters, the mapping quality one with EQUAL's zero.
    let mut plain_filters = filters.clone();
    let mapping_quality = plain_filters
        .iter()
        .position(|filter| filter.name == "MappingQualityReadFilter")
        .map(|index| plain_filters.remove(index));
    let maximum =
        scalar(parser, "maximum-mapping-quality").and_then(|text| text.parse::<i32>().ok());
    let filter = read_filter(parser, &plain_filters, &header)?;
    let keep = |read: &BamRecord| -> bool {
        if !filter(read) {
            return false;
        }
        match &mapping_quality {
            None => true,
            Some(resolved) => {
                let mq = i32::from(read.mapping_quality);
                let passes = mq >= effective_minimum && maximum.is_none_or(|max| mq <= max);
                passes != resolved.negated
            }
        }
    };
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &keep)
        .map_err(reads_traversal_error)?;
    let name_of = |index: i32| -> Option<String> {
        header.sequences.get(index as usize).map(|s| s.name.clone())
    };
    let int_tag = |read: &BamRecord, tag: &[u8; 2]| -> Option<i32> {
        match read.tags.get(htsjdk_bam::tag::Tag::new(tag)) {
            Some(htsjdk_bam::tag::TagValue::Int(value)) => Some(*value as i32),
            _ => None,
        }
    };
    let blocks_of =
        |cigar: &htsjdk_bam::cigar::Cigar, start: i32, contig: &str| -> Vec<gee::Interval> {
            htsjdk_bam::alignment_block::alignment_blocks(cigar, start)
                .into_iter()
                .map(|block| gee::Interval {
                    contig: contig.to_string(),
                    start: block.reference_start,
                    end: block.reference_start + block.length - 1,
                })
                .collect()
        };
    let reads: Vec<gee::Read> = records
        .iter()
        .map(|read| {
            let contig = name_of(read.reference_index).unwrap_or_default();
            let mate_contig = name_of(read.mate_reference_index);
            let mate_blocks = match read.tags.get(htsjdk_bam::tag::Tag::new(b"MC")) {
                Some(htsjdk_bam::tag::TagValue::Str(text)) => {
                    htsjdk_bam::text_parse::parse_cigar(text).ok().map(|cigar| {
                        blocks_of(
                            &cigar,
                            read.mate_alignment_start,
                            mate_contig.as_deref().unwrap_or(""),
                        )
                    })
                }
                _ => None,
            };
            gee::Read {
                name: read.read_name.clone(),
                contig: contig.clone(),
                start: read.alignment_start,
                blocks: blocks_of(&read.cigar, read.alignment_start, &contig),
                end: read.alignment_end(),
                reverse: read.flags & 0x10 != 0,
                paired: read.flags & 0x1 != 0,
                proper_pair: read.flags & 0x2 != 0,
                first_of_pair: read.flags & 0x1 != 0 && read.flags & 0x40 != 0,
                mate_unmapped: read.flags & 0x8 != 0,
                mate_contig,
                mate_start: (read.mate_alignment_start > 0).then_some(read.mate_alignment_start),
                mate_blocks,
                mate_reverse: read.flags & 0x20 != 0,
                mate_quality: int_tag(read, b"MQ"),
                mapping_quality: i32::from(read.mapping_quality),
                hits: int_tag(read, b"NH"),
                fragment_length: read.inferred_insert_size,
            }
        })
        .collect();
    let coverages = gee::count(&features, &reads, &settings)
        .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?;
    let inputs = arguments(parser, "input");
    let text = gee::write_counts(
        &features,
        &coverages,
        sample.as_deref().unwrap_or("null"),
        label,
        &inputs,
        &gff_path,
    );
    write_file(&output, text.as_bytes())?;
    Ok(None)
}

/// `FastaAlternateReferenceMaker.apply`, which is the maker's with a VCF applied at every locus.
///
/// The startup is `FastaReferenceMaker`'s to the line, because it IS that class's: `-L` resolves
/// against the best available dictionary, which a `--sequence-dictionary` outranks the reference
/// in, and the traversal then queries the FASTA. What this tool adds is read before the traversal
/// and in the reference's own order:
///
///   - `super.onTraversalStart()` builds the writer, so a `--line-width` the writer refuses is
///     refused before either check below and leaves no files at all;
///   - `--snp-mask-priority` without `--snp-mask` is a `CommandLineException`, thrown by the tool
///     rather than by the parser;
///   - `--use-iupac-sample` is checked against the DRIVING variants' header samples, so a sample
///     that appears in no record still passes.
///
/// Both of those refusals happen after the writer exists, which is why they leave the same three
/// empty files a refused traversal does.
///
/// The two feature inputs are read whole rather than queried. htsjdk's `FeatureDataSource` queries
/// them by interval and therefore wants an index; nothing in this tool's array reaches an
/// unindexed one, so the refusal that would produce is not modelled here rather than guessed at.
pub fn fasta_alternate_reference_maker(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "FastaAlternateReferenceMaker")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    // `--variant` is a scalar `FeatureInput` here, as it is on every variant walker.
    let variant = argument(parser, "variant").ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;

    let (variants, samples) = feature_variants(&variant)?;
    let mask = match argument(parser, "snp-mask") {
        Some(path) => Some(feature_variants(&path)?.0),
        None => None,
    };

    let mut source =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let master = master_dictionary(parser)?;
    let own = gatk_tools::reference_walker::dictionary(&source);
    let best = master.unwrap_or(own);
    let intervals = match interval_arguments(parser, &best)? {
        Some(parameters) => parameters.intervals,
        None => best
            .sequences
            .iter()
            .map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
                    .expect("a contig length is at least one")
            })
            .collect(),
    };

    let width = number_or(
        parser,
        "line-width",
        gatk_tools::fasta_reference_maker::DEFAULT_LINE_WIDTH as i32,
    )
    .max(0) as usize;

    let arguments = gatk_tools::fasta_alternate_reference_maker::AlternateArguments {
        mask: mask.as_deref(),
        mask_priority: flag(parser, "snp-mask-priority"),
        iupac_sample: argument(parser, "use-iupac-sample"),
    };

    // The writer is built in `onTraversalStart` and closed in `closeTool`, so every refusal after
    // it exists still leaves a FASTA, a `.fai` and a dictionary of nothing but its `@HD` line.
    let refused = |error| -> Thrown {
        match gatk_tools::fasta_reference_maker::empty_outputs(width) {
            Ok(empty) => {
                let _ = write_outputs(&output, &empty);
                alternate_maker_error(error)
            }
            // The width itself is what the writer refused, and that happens before any file is
            // opened: `FastaReferenceWriterBuilder.build` checks it first for exactly that reason.
            Err(failure) => alternate_maker_error(
                gatk_tools::fasta_alternate_reference_maker::AlternateError::Maker(failure),
            ),
        }
    };
    let outputs = gatk_tools::fasta_alternate_reference_maker::run_over(
        &mut source,
        &intervals,
        width,
        &variants,
        &arguments,
        &samples,
    )
    .map_err(refused)?;

    write_outputs(&output, &outputs)?;
    Ok(None)
}

/// A `FeatureInput<VariantContext>`: the file's records and the samples its header declares.
///
/// The codec is chosen by the file's NAME, which is `FeatureManager`'s rule, and a name no codec
/// claims is the same refusal `IndexFeatureFile` gives.
fn feature_variants(
    path: &str,
) -> Result<(Vec<htsjdk_vcf::variant::VariantContext>, Vec<String>), Thrown> {
    if gatk_tools::feature_codec::codec_for(path).is_none() {
        return Err(Thrown::user(
            index_feature_file::Refusal::NoSuitableCodecs {
                path: path.to_string(),
            }
            .message(),
        ));
    }
    let bytes = std::fs::read(path).map_err(|_| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: path.to_string(),
            }
            .message(),
        )
    })?;
    let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
        htsjdk_bgzf::read::decompress_all(&bytes)
            .ok()
            .and_then(|bytes| String::from_utf8(bytes).ok())
            .ok_or_else(|| {
                Thrown::non_user(
                    gatk_tools::read_walker_refusal::SAM_FORMAT,
                    format!("{path} is not a block compressed file"),
                )
            })?
    } else {
        String::from_utf8_lossy(&bytes).into_owned()
    };
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;
    Ok((file.records, file.header.samples))
}

/// What `FastaAlternateReferenceMaker` refused with, told apart by whose refusal it is.
fn alternate_maker_error(
    error: gatk_tools::fasta_alternate_reference_maker::AlternateError,
) -> Thrown {
    match error {
        // The tool's own two checks, one of which Barclay's exception class names even though the
        // parser accepted the command line.
        gatk_tools::fasta_alternate_reference_maker::AlternateError::Argument(argument) => Thrown {
            failure: Failure::User,
            exception: argument.java_class(),
            message: Some(argument.message()),
        },
        gatk_tools::fasta_alternate_reference_maker::AlternateError::Maker(maker) => {
            fasta_maker_error(maker)
        }
    }
}

/// `CompareReferences.traverse`, which is a `GATKTool` that overrides the traversal entirely.
///
/// The engine's reference is the FIRST column and `--references-to-compare` the rest, in a
/// `LinkedHashMap` keyed by path: the same path twice is ONE entry, which is how a reference
/// compared with itself produces no pair at all and the tool then walks off the empty list.
/// That `IndexOutOfBoundsException` is the reference's answer and it is reproduced here.
///
/// The output is split: the table goes to `--output`, or to stdout when there is none, and the
/// analysis always goes to stdout behind a line of asterisks.
pub fn compare_references(parser: &Parser) -> Outcome {
    use gatk_tools::compare_references as compare;

    let _ = resolve_read_filters(parser, "CompareReferences")?;
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let others = arguments(parser, "references-to-compare");
    if others.is_empty() {
        return Err(Thrown::command_line(
            "Argument references-to-compare was missing: Argument 'references-to-compare' is required",
        ));
    }

    // `GATKTool.onStartup` resolves the intervals against the best available dictionary before any
    // tool code runs, and this tool then ignores them: `traverse()` is its own. What survives is
    // the refusal, so an interval on a contig the reference does not carry is refused here exactly
    // as it is on a tool that traverses.
    let dictionary = reference_dictionary(parser)?.unwrap_or_else(SamHeader::default);
    let master = master_dictionary(parser)?;
    let best = master.unwrap_or(dictionary);
    let _ = interval_arguments(parser, &best)?;

    // The enums arrive as their constant names, and `scalar` is what reads one: `argument` reads
    // the `Tagged` and `Str` values a path or a string argument holds, and an enum is neither.
    let mode = match scalar(parser, "md5-calculation-mode").as_deref() {
        Some("USE_DICT") => compare::Md5Mode::UseDict,
        Some("ALWAYS_RECALCULATE") => compare::Md5Mode::AlwaysRecalculate,
        _ => compare::Md5Mode::RecalculateIfMissing,
    };
    let base_comparison =
        scalar(parser, "base-comparison").unwrap_or_else(|| "NO_BASE_COMPARISON".to_string());

    // `LinkedHashMap.put`: the engine's reference first, then the others, and a path already in
    // the map keeps its position and adds no column.
    let mut paths: Vec<String> = vec![reference];
    for path in others {
        if !paths.contains(&path) {
            paths.push(path);
        }
    }

    // `onTraversalStart`, in its own order: the directory before the count.
    if base_comparison != "NO_BASE_COMPARISON" {
        let directory = argument(parser, "base-comparison-output");
        match &directory {
            None => {
                return Err(could_not_create_output_file(
                    "null",
                    &format!(
                        "Output directory not provided but required in -base-comparison {base_comparison} mode."
                    ),
                ))
            }
            Some(path) if !std::path::Path::new(path).exists() => {
                return Err(could_not_create_output_file(
                    path,
                    "Output directory non-existent.",
                ))
            }
            Some(_) => {}
        }
        if paths.len() != 2 {
            return Err(bad_input(
                "Base comparison modes can only be run on 2 references.".to_string(),
            ));
        }
    }

    let references = compare_references_read(&paths, mode)?;
    let table = compare::build(&references, mode).map_err(compare_table_error)?;

    let rendered = compare::write_table(&table);
    let mut stdout = String::new();
    match argument(parser, "output") {
        Some(path) => write_file(&path, rendered.as_bytes())?,
        None => stdout.push_str(&rendered),
    }
    if flag(parser, "display-sequences-by-name") {
        stdout.push_str(&compare::write_by_sequence_name(
            &table,
            &references,
            flag(parser, "display-only-differing-sequences"),
        ));
    }

    let pairs = compare::compare_all(&table, &references).map_err(compare_table_error)?;
    stdout.push_str("*********************************************************\n");
    for pair in &pairs {
        stdout.push_str(&pair.rendered());
        stdout.push('\n');
    }
    // Everything above is printed as it is produced, so a refusal below it keeps what came before,
    // and `onTraversalSuccess` returns null: there is no `Tool returned:` line on this tool.
    print!("{stdout}");
    // `referencePairs.get(0)` is read before the switch on the base-comparison mode, so a run with
    // one reference dies here having already printed everything above.
    if pairs.is_empty() {
        return Err(Thrown::non_user(
            "java.lang.IndexOutOfBoundsException",
            "Index 0 out of bounds for length 0",
        ));
    }
    // FULL_ALIGNMENT runs mummer, which no row of this tool's array reaches: the mode's two
    // non-default constants are held out of the fixtures for that reason.
    Ok(None)
}

/// Each reference as the table reads it: the file's NAME as the column, and its `.dict`.
///
/// The MD5 is recalculated only where the mode asks for it, which is what keeps `USE_DICT` from
/// reading a base at all. The recalculation is NOT the tool's own reference query: it is
/// `ReferenceUtils.calculateMD5`, which opens the FASTA with `preserveCase` and `preserveIUPAC`
/// both on and then upper-cases each byte, so a soft-masked base counts and an ambiguity code is
/// not flattened to `N`.
fn compare_references_read(
    paths: &[String],
    mode: gatk_tools::compare_references::Md5Mode,
) -> Result<Vec<gatk_tools::compare_references::Reference>, Thrown> {
    use gatk_tools::compare_references as compare;

    let mut references = Vec::new();
    for path in paths {
        let dictionary = std::path::Path::new(path).with_extension("dict");
        let text = std::fs::read_to_string(&dictionary)
            .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(path, false)))?;
        let header = htsjdk_bam::reader::parse_header_text(&text);
        let column = std::path::Path::new(path)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.clone());

        let mut fasta: Option<htsjdk_bam::fasta_index::IndexedFasta<std::fs::File>> = None;
        let mut sequences = Vec::new();
        for record in &header.sequences {
            let md5 = record.attributes.get("M5").map(str::to_string);
            let recalculate = match mode {
                compare::Md5Mode::AlwaysRecalculate => true,
                compare::Md5Mode::RecalculateIfMissing => md5.as_deref().is_none_or(str::is_empty),
                compare::Md5Mode::UseDict => false,
            };
            let calculated = if recalculate {
                if fasta.is_none() {
                    fasta = Some(
                        htsjdk_bam::fasta_index::IndexedFasta::open(std::path::Path::new(path))
                            .map_err(|error| Thrown::user(error.message()))?,
                    );
                }
                let reader = fasta.as_mut().expect("the FASTA opened just above");
                let bases = reader
                    .query(&record.name, 1, record.length as i64)
                    .map_err(|error| Thrown::user(error.message()))?;
                compare::calculate_md5(&bases)
            } else {
                // Never read in this mode, and the reference has not opened the file either.
                String::new()
            };
            sequences.push(compare::Sequence {
                name: record.name.clone(),
                length: record.length as i64,
                md5,
                calculated_md5: calculated,
            });
        }
        references.push(compare::Reference { column, sequences });
    }
    Ok(references)
}

/// `UserException.CouldNotCreateOutputFile`, whose message names the file first.
fn could_not_create_output_file(path: &str, reason: &str) -> Thrown {
    Thrown {
        failure: Failure::User,
        exception:
            "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile",
        message: Some(format!("Could not create file {path}. {reason}")),
    }
}

/// What the table refused. Its `message()` already carries `UserException$BadInput`'s own
/// `Bad input: ` prefix, so this does not add a second one.
fn compare_table_error(error: gatk_tools::compare_references::TableError) -> Thrown {
    // The class is written out rather than taken from `java_class()`, which borrows the error:
    // both of that enum's variants are the same exception, and the port says so.
    Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(error.message()),
    }
}

/// `CollectReadCounts.apply`, which counts one read into the interval its START falls in.
///
/// A read walker whose traversal is `CountReads`', and three things around it that are the tool's
/// own:
///
///   - `requiresIntervals()` is true, so the interval argument is REQUIRED by its declaration and
///     a run without `-L` is refused by the parser rather than by the tool;
///   - `validateIntervalArgumentCollection` is the copy-number check `PreprocessIntervals` and
///     `AnnotateIntervals` make, which is why it lives in one place here;
///   - and the sample name is read off the header's read groups, one distinct `SM` or a refusal.
///
/// `--format HDF5` is the default and is refused: the port writes the TSV and HDF5 is a file format
/// rather than a spelling of it.
pub fn collect_read_counts(parser: &Parser) -> Outcome {
    if scalar(parser, "format").as_deref() != Some("TSV") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "CollectReadCounts writes HDF5 by default, which this port does not write; pass \
             --format TSV. This message is the port's own and not GATK's.",
        ));
    }
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CollectReadCounts")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if arguments(parser, "intervals").is_empty()
        && arguments(parser, "exclude-intervals").is_empty()
    {
        return Err(Thrown::command_line(
            "Argument intervals was missing: Argument 'intervals' is required",
        ));
    }
    validate_copy_number_intervals(parser)?;

    let sample = sample_name(&header)?;
    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &filter)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    // A read with no contig cannot fall in an interval, so it counts nowhere: `apply` reads
    // `read.getContig()`, which the wellformed filter has already guaranteed is a mapped one.
    let starts: Vec<(&str, i32)> = records
        .iter()
        .filter_map(|record| {
            contig_name(&header, record.reference_index)
                .map(|contig| (contig, record.alignment_start))
        })
        .collect();
    let counts = gatk_tools::collect_read_counts::count(&starts, &intervals);
    // The TSV's `@SQ` lines are the METADATA's dictionary, which is the reads' own header rather
    // than the best available one: a `--sequence-dictionary` resolves the intervals and is only
    // warned about here when it disagrees.
    let sequences: Vec<(String, i32)> = header
        .sequences
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();
    let text = gatk_tools::collect_read_counts::write(&sequences, &sample, &intervals, &counts);
    write_file(&output, text.as_bytes())?;
    Ok(None)
}

/// `MetadataUtils.readSampleName`: one distinct `SM` over the read groups, or an
/// `IllegalArgumentException` naming what it found instead.
fn sample_name(header: &SamHeader) -> Result<String, Thrown> {
    let refuse = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);
    if header.read_groups.is_empty() {
        return Err(refuse(
            "The input header does not contain any read groups.  Cannot determine a sample name."
                .to_string(),
        ));
    }
    let mut samples: Vec<String> = Vec::new();
    for group in &header.read_groups {
        if let Some(sample) = group.attributes.get("SM").map(str::to_string) {
            if !samples.contains(&sample) {
                samples.push(sample);
            }
        }
    }
    match samples.len() {
        0 => Err(refuse(
            "The input header does not contain a sample name.".to_string(),
        )),
        1 => Ok(samples.remove(0)),
        _ => Err(refuse(format!(
            "The input header contains more than one unique sample name: {}",
            samples.join(", ")
        ))),
    }
}

/// `GetSampleName.traverse`, which is the whole tool: it opens the reads to read their header and
/// never asks for a record.
///
/// The read walker's startup still runs in front of it -- the dictionaries are compared and the
/// intervals resolved -- because the arguments that ask for those are declared whether or not a
/// traversal uses them.
pub fn get_sample_name(parser: &Parser) -> Outcome {
    let ReadWalkerStart { source, .. } = read_walker_startup(parser, "GetSampleName")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let text =
        gatk_tools::get_sample_name::get_sample_name(&source, flag(parser, "use-url-encoding"))
            .map_err(|error| Thrown::user(format!("Bad input: {}", error.message())))?;
    std::fs::write(&output, &text)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

/// `PrintDistantMates.apply`, which writes the reads whose mate is far away or on another contig.
///
/// `PrintReads`' plumbing with a filter in front of it, and one alteration on the way out: a
/// printed read carries its original alignment in an `OA` tag and is written unmapped-adjacent,
/// which `do_distant_mate_alterations` is the port of.
pub fn print_distant_mates(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "PrintDistantMates")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let filter = read_filter(parser, &filters, &header)?;
    // The tool's own `doWork`, which is three steps rather than one: every read the filter chain
    // kept is ALTERED -- moved to its mate's position and unmapped -- and the result is RE-SORTED,
    // because an alteration that moves a read breaks the order it arrived in. The port had all of
    // that already; the runner was reimplementing two of the three and getting both wrong.
    let command_line = crate::command_line::expanded("PrintDistantMates", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let (bytes, bai) = gatk_tools::print_distant_mates::print_distant_mates_with(
        &source, &options, &filter, level, deflater,
    )
    .map_err(|error| Thrown::user(format!("{error:?}")))?;
    std::fs::write(&output, &bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    if let Some(bai) = bai {
        // The index REPLACES the output's extension, and the digest APPENDS to it.
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
        std::fs::write(&digest, gatk_tools::gather_bam_files::md5_file(&bytes)).map_err(
            |error| Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}")),
        )?;
    }
    Ok(None)
}

/// `CompareIntervalLists.doWork`: two interval lists, and whether they cover the same bases.
///
/// The first tool here that is a `CommandLineProgram` rather than a `GATKTool`, which the counts
/// say plainly: fifteen arguments declared and fifteen on the instance, where every walker has
/// seventy against thirty-eight. So there is no startup at all -- no reads, no plugin descriptor,
/// no interval collection -- and `-L` and `-L2` are this tool's own arguments rather than the
/// engine's.
pub fn compare_interval_lists(parser: &Parser) -> Outcome {
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let dictionary = reference_dictionary(parser)?.ok_or_else(|| {
        Thrown::user(gatk_tools::read_walker_refusal::cannot_read(
            &reference, false,
        ))
    })?;

    // Each file is parsed and `ALL`-merged against the reference's dictionary before the two are
    // compared, which is what the tool's own `getGenomeLocs` does to each.
    let mut lists = Vec::new();
    for name in ["L", "L2"] {
        // SCALARS, not collections: this tool declares `-L` and `-L2` itself, as one string each,
        // where the engine's own `--intervals` is a list. Reading them as a collection found
        // nothing at all and refused a command line the reference runs.
        let given = argument(parser, name).ok_or_else(|| {
            Thrown::command_line(format!(
                "Argument {name} was missing: Argument '{name}' is required"
            ))
        })?;
        let parameters = gatk_engine::interval_arguments::traversal_parameters(
            std::slice::from_ref(&given),
            &[],
            &dictionary,
            SetRule::Union,
            MergingRule::All,
            0,
            0,
        )
        .map_err(|refusal| Thrown {
            failure: Failure::User,
            exception: refusal.java_class(),
            message: Some(refusal.message()),
        })?;
        lists.push(parameters.intervals);
    }

    let comparison = gatk_tools::compare_interval_lists::equate_intervals(&lists[0], &lists[1]);
    match comparison {
        gatk_tools::compare_interval_lists::Comparison::Equal => {
            // `doWork` prints the verdict and returns 0.
            println!("Intervals are equal");
            Ok(Some("0".to_string()))
        }
        gatk_tools::compare_interval_lists::Comparison::Different(difference) => Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("Intervals are not equal: \n{difference}")),
        }),
        // `test.pop()` on an empty list, which leaves the engine as itself and carries NO message
        // at all. `None` rather than an empty string: the handler prints `<class>: <message>`, so
        // an empty one leaves a trailing colon the reference never writes.
        gatk_tools::compare_interval_lists::Comparison::TestExhausted => Err(Thrown {
            failure: Failure::Other,
            exception: "java.util.NoSuchElementException",
            message: None,
        }),
    }
}

/// `FixMisencodedBaseQualityReads.apply`: every quality with 31 taken off it.
///
/// A read walker that writes a BAM, which is `PrintDistantMates`' shape; what is its own is the
/// refusal, which fires on the FIRST quality below the offset and abandons the run rather than
/// writing what it had.
pub fn fix_misencoded_base_quality_reads(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "FixMisencodedBaseQualityReads")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("FixMisencodedBaseQualityReads", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let outcome =
        gatk_tools::fix_misencoded_base_quality_reads::fix_misencoded_base_quality_reads_with(
            &source, &options, &filter, level, deflater,
        )
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let (bytes, bai) = match outcome {
        Ok(written) => written,
        Err(refusal) => {
            return Err(Thrown {
                failure: Failure::User,
                exception: refusal.class(),
                message: Some(refusal.message()),
            })
        }
    };

    std::fs::write(&output, &bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
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
        std::fs::write(&digest, gatk_tools::gather_bam_files::md5_file(&bytes)).map_err(
            |error| Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}")),
        )?;
    }
    Ok(None)
}

/// `UnmarkDuplicates.apply`, which clears one flag bit and writes the read back.
///
/// `PrintReads`' plumbing with a mutation in the middle, and the same three companions on the way
/// out: the index when `--create-output-bam-index` asks for one, the MD5 when
/// `--create-output-bam-md5` does, and the `@PG` line the writer adds.
pub fn unmark_duplicates(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "UnmarkDuplicates")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("UnmarkDuplicates", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let (bytes, bai) = gatk_tools::unmark_duplicates::unmark_duplicates_with(
        &source, &options, &filter, level, deflater,
    )
    .map_err(reads_traversal_error)?;
    write_bam(parser, &output, &bytes, bai)
}

/// `RevertBaseQualityScores`, which is `UnmarkDuplicates`' plumbing with an abort in the middle.
///
/// The refusal is the tool's own `UserException` and it happens PART WAY: the reference writes as
/// it goes, so a run that hits a read with no `OQ` leaves whatever had been flushed behind and
/// exits non-zero. This port writes nothing in that case, and the difference is the one thing here
/// that is not the reference's: a partial BAM is not an answer a covering array can compare, and
/// the row is the refusal either way.
pub fn revert_base_quality_scores(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "RevertBaseQualityScores")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("RevertBaseQualityScores", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let run = gatk_tools::revert_base_quality_scores::revert_base_quality_scores_with(
        &source, &options, &filter, level, deflater,
    )
    .map_err(reads_traversal_error)?;
    match run {
        Ok((bytes, bai)) => write_bam(parser, &output, &bytes, bai),
        // Two exceptions, two handlers: the tool's own `UserException` is decorated and exits at
        // two, while `fastqToPhred`'s `IllegalArgumentException` is a bug rather than a refusal and
        // prints its class before the message at three.
        Err(refusal) => Err(
            if refusal.class() == gatk_tools::main_entry::USER_EXCEPTION {
                Thrown::user(refusal.message())
            } else {
                Thrown::non_user(refusal.class(), refusal.message())
            },
        ),
    }
}

/// `AddOriginalAlignmentTags`, the first of the archetype that writes TAGS rather than changing
/// the read.
///
/// Its refusal is htsjdk's rather than the tool's: `getMateReferenceName` on an unpaired read
/// throws `IllegalStateException`, so the handler prints the class in front of the message and the
/// run ends at status three rather than two.
pub fn add_original_alignment_tags(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "AddOriginalAlignmentTags")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("AddOriginalAlignmentTags", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let run = gatk_tools::add_original_alignment_tags::add_original_alignment_tags_with(
        &source, &options, &filter, level, deflater,
    )
    .map_err(reads_traversal_error)?;
    match run {
        Ok((bytes, bai)) => write_bam(parser, &output, &bytes, bai),
        Err(refusal) => Err(Thrown::non_user(refusal.class(), refusal.message())),
    }
}

/// `LeftAlignIndels`, the first read walker here whose REFERENCE is required.
///
/// The window each read is left-aligned in is the read's own span, which the walker builds as
/// `new ReferenceContext(reference, new SimpleInterval(read))`, so the reference is queried per
/// read and not per interval. `AlignmentUtils.leftAlignIndels` raises `IllegalArgumentException`
/// on a cigar it cannot align, which is a bug rather than a refusal: the handler prints the class
/// and the run ends at three.
pub fn left_align_indels(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "LeftAlignIndels")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // `requiresReference()` is true, so the parser refuses a run without one before this.
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("LeftAlignIndels", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let dictionary = gatk_tools::reference_walker::dictionary(&reference);
    let run = gatk_tools::left_align_indels::left_align_indels_with(
        &source,
        &mut reference,
        &options,
        &filter,
        level,
        deflater,
    )
    .map_err(|error| match error {
        // `MissingContigInSequenceDictionary`, raised by the per-read reference query: the
        // dictionary it prints is the REFERENCE's, which is why this is formatted here. Measured on
        // row 6 of this tool's array, where the reads are on `chr1` and `--reference` is the fasta
        // that carries `chrOther`; the port answered a `SAMFormatException` at three.
        gatk_engine::reads::ReadsError::ContigNotInDictionary(contig) => Thrown::user(format!(
            "Contig {contig} not present in the sequence dictionary {}\n",
            gatk_tools::sequence_dictionary::pretty_print(&dictionary.sequences)
        )),
        other => reads_traversal_error(other),
    })?;
    match run {
        Ok((bytes, bai)) => write_bam(parser, &output, &bytes, bai),
        Err(error) => Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            format!("{error:?}"),
        )),
    }
}

/// `DumpTabixIndex`, which is no walker at all: a `.tbi` in, its text out.
///
/// Three refusals, in the order the tool reaches them, and the first is the file's NAME.
///
///   - a path that does not end in `.tbi` is refused before anything is opened, whatever it
///     holds: `Expected a .tbi file as input.`;
///   - a `.tbi` that is not gzipped fails inside `java.util.zip`, and the tool CATCHES that and
///     raises its own `Trouble reading index.` -- the bare `java.util.zip.ZipException` in the
///     dump-tabix-index golden is the same failure reached through the tool's method rather than
///     through `Main`, which is a different door and a different handler;
///   - and a gzipped `.tbi` whose magic is not `TBI\1` is `Incorrect magic number for tabix
///     index`.
///
/// All three are `UserException` at status two.
pub fn dump_tabix_index(parser: &Parser) -> Outcome {
    let input = argument(parser, "tabix-index").ok_or_else(|| {
        Thrown::command_line("Argument tabix-index was missing: Argument 'tabix-index' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // The NAME is checked before anything is opened, so a file that is not called `.tbi` is
    // refused for its name and never reaches the gzip layer whatever it holds.
    if !input.ends_with(".tbi") {
        return Err(Thrown::user("Expected a .tbi file as input."));
    }
    let bytes = std::fs::read(&input)
        .map_err(|error| Thrown::user(format!("Couldn't read {input}: {error}")))?;
    let decompressed = htsjdk_bgzf::read::decompress_all(&bytes)
        .map_err(|_| Thrown::user("Trouble reading index."))?;
    let text = gatk_tools::dump_tabix_index::dump_tabix_index(&decompressed)
        .map_err(|error| Thrown::user(error.message()))?;
    std::fs::write(&output, text).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    Ok(None)
}

/// `ReadAnonymizer`, the first declared walker that REWRITES the bases it reads.
///
/// Every base a match consumes is replaced by the reference's, a base that already agreed keeps its
/// own quality and one that was replaced takes `--ref-base-quality`, and the cigar's `M` becomes
/// `=` unless `--use-simple-cigar` says otherwise. The window is the read's own span and is built
/// per read, so a reference that does not carry the read's contig is refused during the traversal.
pub fn read_anonymizer(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "ReadAnonymizer")?;
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
    let command_line = crate::command_line::expanded("ReadAnonymizer", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let arguments = gatk_tools::read_anonymizer::AnonymizerArguments {
        ref_base_quality: u8::try_from(number_or(
            parser,
            "ref-base-quality",
            i32::from(gatk_tools::read_anonymizer::DEFAULT_REF_BASE_QUALITY),
        ))
        .unwrap_or(gatk_tools::read_anonymizer::DEFAULT_REF_BASE_QUALITY),
        use_simple_cigar: flag(parser, "use-simple-cigar"),
    };
    let dictionary = gatk_tools::reference_walker::dictionary(&reference);
    let (bytes, bai) = gatk_tools::read_anonymizer::read_anonymizer_with(
        &source,
        &mut reference,
        &arguments,
        &options,
        &filter,
        output_compression(parser).0,
        output_compression(parser).1,
    )
    .map_err(|error| match error {
        gatk_engine::reads::ReadsError::ContigNotInDictionary(contig) => Thrown::user(format!(
            "Contig {contig} not present in the sequence dictionary {}\n",
            gatk_tools::sequence_dictionary::pretty_print(&dictionary.sequences)
        )),
        other => reads_traversal_error(other),
    })?;
    write_bam(parser, &output, &bytes, bai)
}

/// `PrintFileDiagnostics`, which chooses an analyzer by the input's EXTENSION and nothing else.
///
/// `.bai` is the only branch this port carries, and it is the textual index htsjdk writes; `.cram`
/// and `.crai` are the reference's other two. A name no analyzer claims is a `RuntimeException`
/// quoting the argument, and a `.bai` that does not read is htsjdk's `SAMException`, so the two
/// refusals exit at three and print their classes.
pub fn print_file_diagnostics(parser: &Parser) -> Outcome {
    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // `HTSAnalyzerFactory.getFileAnalyzer` prints the input as given before it chooses anything.
    println!("{input}");
    let analyzer = gatk_tools::print_file_diagnostics::analyzer_for(&input)
        .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?;
    let bytes = std::fs::read(&input)
        .map_err(|error| Thrown::user(format!("Couldn't read {input}: {error}")))?;
    let report = match analyzer {
        gatk_tools::print_file_diagnostics::Analyzer::Bai => {
            gatk_tools::print_file_diagnostics::bai_report(&bytes)
                .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?
        }
        // The CRAM analyzers read a container structure this port does not carry. A row that asks
        // for one says so rather than answering with a report that is not the reference's.
        other => {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "The {other:?} analyzer is a GATK feature that this port does not carry yet. \
                     This message is the port's own and not GATK's."
                ),
            ))
        }
    };
    std::fs::write(&output, report).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    // `BAIAnalyzer.doAnalysis` prints where it wrote, `analyze` then emits an empty line, and
    // `doWork` returns 0, which `handleResult` prints.
    println!("\nOutput written to {output}\n");
    println!();
    Ok(Some("0".to_string()))
}

/// `SplitReads`, whose `--output` names a DIRECTORY and whose file names it builds itself.
///
/// One file per key of the splitters that were asked for, named
/// `<input base name><key><input extension>`, and with no `--split-*` at all that is one file named
/// after the input. The key is built per read, so a read whose read group cannot answer a splitter
/// is the tool's refusal rather than a file named `unknown` -- except where the whole key is
/// `.unknown`, which is a file.
pub fn split_reads(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "SplitReads")?;
    let directory = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let input = arguments(parser, "input")
        .into_iter()
        .next()
        .unwrap_or_default();
    // `getOutputFileName`: the base name and the extension of the INPUT, which is what the key is
    // inserted between.
    let file_name = std::path::Path::new(&input)
        .file_name()
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_else(|| input.clone());
    let (base_name, extension) = match file_name.rfind('.') {
        Some(dot) => (file_name[..dot].to_string(), file_name[dot..].to_string()),
        None => (file_name.clone(), String::new()),
    };

    let mut splitters = Vec::new();
    if flag(parser, "split-sample") {
        splitters.push(gatk_tools::split_reads::Splitter::Sample);
    }
    if flag(parser, "split-read-group") {
        splitters.push(gatk_tools::split_reads::Splitter::ReadGroupId);
    }
    if flag(parser, "split-library-name") {
        splitters.push(gatk_tools::split_reads::Splitter::LibraryName);
    }

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("SplitReads", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let run = gatk_tools::split_reads::split_reads_with(
        &source, &options, &splitters, &base_name, &extension, &filter, level, deflater,
    )
    .map_err(reads_traversal_error)?;
    let files = match run {
        Ok(files) => files,
        Err(refusal) => return Err(Thrown::non_user(refusal.class(), refusal.message())),
    };
    for file in &files {
        let path = std::path::Path::new(&directory).join(&file.name);
        std::fs::write(&path, &file.bam).map_err(|error| {
            Thrown::non_user(
                PORT_FAILURE,
                format!("could not write {}: {error}", path.display()),
            )
        })?;
        if let Some(index) = &file.index {
            let companion = path.with_extension("bai");
            std::fs::write(&companion, index).map_err(|error| {
                Thrown::non_user(
                    PORT_FAILURE,
                    format!("could not write {}: {error}", companion.display()),
                )
            })?;
        }
        // The digest is written per FILE, like the index: a run that splits into six files and
        // asks for md5s leaves six of them, each APPENDED to its own name.
        if flag(parser, "create-output-bam-md5") {
            let digest = format!("{}.md5", path.display());
            std::fs::write(&digest, gatk_tools::gather_bam_files::md5_file(&file.bam)).map_err(
                |error| {
                    Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}"))
                },
            )?;
        }
    }
    Ok(None)
}

/// `ClipReads`, which writes a BAM and, when asked, a statistics file beside it.
///
/// `--cycles-to-trim` is kept unparsed until the tool runs, because the parse fails with a message
/// of its own: every malformed spelling raises the same `RuntimeException`.
pub fn clip_reads(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "ClipReads")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    // `-XF`, read here rather than in the tool: the reference reads a FASTA and the port takes the
    // records, so a file that cannot be read is this runner's refusal.
    let mut clip_sequence_file = Vec::new();
    if let Some(path) = argument(parser, "clip-sequences-file") {
        let text = std::fs::read_to_string(&path)
            .map_err(|error| Thrown::user(format!("Couldn't read {path}: {error}")))?;
        clip_sequence_file = gatk_tools::clip_reads::parse_clip_sequence_file(&text);
    }
    let arguments_for_clip = gatk_tools::clip_reads::ClipArguments {
        q_trimming_threshold: number_or(parser, "q-trimming-threshold", -1),
        cycles_to_clip: argument(parser, "cycles-to-trim"),
        clip_sequences: arguments(parser, "clip-sequence"),
        clip_sequence_file,
        clipping_representation: {
            use gatk_engine::clipping::ClippingRepresentation as Representation;
            // The parser refuses a constant this list does not carry, so the default arm is the
            // argument left unset rather than a spelling nothing recognises.
            match scalar(parser, "clip-representation").as_deref() {
                Some("WRITE_Q0S") => Representation::WriteQ0s,
                Some("WRITE_NS_Q0S") => Representation::WriteNsQ0s,
                Some("SOFTCLIP_BASES") => Representation::SoftclipBases,
                Some("HARDCLIP_BASES") => Representation::HardclipBases,
                Some("REVERT_SOFTCLIPPED_BASES") => Representation::RevertSoftclippedBases,
                _ => Representation::WriteNs,
            }
        },
        only_do_read: argument(parser, "read"),
        clip_adapter: flag(parser, "clip-adapter"),
        min_read_length: number_or(parser, "min-read-length-to-output", 0),
    };

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("ClipReads", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let run = gatk_tools::clip_reads::clip_reads_with(
        &source,
        &options,
        &arguments_for_clip,
        &filter,
        level,
        deflater,
    )
    .map_err(reads_traversal_error)?;
    let (bytes, bai, statistics) = match run {
        Ok(produced) => produced,
        // Both of these are `RuntimeException`s rather than the tool's own refusal, so the class is
        // printed and the run ends at three. The cycles message is the reference's own; a clip the
        // `ReadClipper` refuses is not reached by any row of this tool's array, and the corpus is
        // where that would show.
        Err(gatk_tools::clip_reads::ClipReadsError::BadlyFormattedCycles(argument)) => {
            return Err(Thrown::non_user(
                "java.lang.RuntimeException",
                format!("Badly formatted cyclesToClip argument: {argument}"),
            ))
        }
        Err(gatk_tools::clip_reads::ClipReadsError::Clip(error)) => {
            return Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!("{error:?}"),
            ))
        }
    };
    if let Some(path) = argument(parser, "output-statistics") {
        std::fs::write(&path, &statistics).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {path}: {error}"))
        })?;
    }
    write_bam(parser, &output, &bytes, bai)?;
    // `onTraversalSuccess` returns the clipping statistics, which `handleResult` prints.
    Ok(Some(statistics))
}

/// `SplitNCigarReads`, which splits a read at every `N` of its cigar.
///
/// The overhang-fixing manager holds a queue of reads and asks the reference for the bases around a
/// splice, so the runner hands it a query closure over the reference source rather than a slice: the
/// window is per SPLICE and not per read.
pub fn split_n_cigar_reads(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "SplitNCigarReads")?;
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
    let command_line = crate::command_line::expanded("SplitNCigarReads", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let arguments_for_split = gatk_tools::split_n_cigar_reads::SplitArguments {
        refactor_ndn_cigar_reads: flag(parser, "refactor-cigar-string"),
        skip_mq_transform: flag(parser, "skip-mapping-quality-transform"),
        process_secondary_alignments: flag(parser, "process-secondary-alignments"),
        overhang: gatk_engine::overhang_fixing_manager::OverhangArguments {
            max_records_in_memory: usize::try_from(number_or(
                parser,
                "max-reads-in-memory",
                150_000,
            ))
            .unwrap_or(150_000),
            max_mismatches_in_overhang: number_or(parser, "max-mismatches-in-overhang", 1),
            max_bases_in_overhang: number_or(parser, "max-bases-in-overhang", 40),
            do_not_fix_overhangs: flag(parser, "do-not-fix-overhangs"),
            process_secondary_reads: flag(parser, "process-secondary-alignments"),
        },
    };

    let dictionary = gatk_tools::reference_walker::dictionary(&reference);
    let mut query = |contig: &str, start: i32, end: i32| -> Result<Vec<u8>, String> {
        reference
            .query(contig, start, end)
            .map_err(|error| format!("{error:?}"))
    };
    let (level, deflater) = output_compression(parser);
    let produced = gatk_tools::split_n_cigar_reads::split_n_cigar_reads_with(
        &source,
        &arguments_for_split,
        &options,
        &filter,
        &mut query,
        level,
        deflater,
    );
    let (bytes, bai) = match produced {
        Ok(produced) => produced,
        Err(gatk_tools::split_n_cigar_reads::SplitToolError::Reads(error)) => {
            return Err(match error {
                gatk_engine::reads::ReadsError::ContigNotInDictionary(contig) => {
                    Thrown::user(format!(
                        "Contig {contig} not present in the sequence dictionary {}\n",
                        gatk_tools::sequence_dictionary::pretty_print(&dictionary.sequences)
                    ))
                }
                other => reads_traversal_error(other),
            })
        }
        // `splitReadBasedOnCigar` and the manager both raise GATK's own unchecked exception, so the
        // class is printed and the run ends at three.
        Err(gatk_tools::split_n_cigar_reads::SplitToolError::Split(error)) => {
            return Err(Thrown::non_user(
                "org.broadinstitute.hellbender.exceptions.GATKException",
                error.message(),
            ))
        }
    };
    write_bam(parser, &output, &bytes, bai)
}

/// `MethylationTypeCaller`, the first declared read walker that writes a VCF.
///
/// `--add-output-vcf-command-line` adds `##source` and `##GATKCommandLine`, and the second carries
/// the run's own DATE, which is why the suite's comparable runs turn it off. The line is built here
/// because the tool takes the default lines as a parameter for that reason.
pub fn methylation_type_caller(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "MethylationTypeCaller")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let mut default_lines = Vec::new();
    if flag(parser, "add-output-vcf-command-line") {
        let command_line = crate::command_line::expanded("MethylationTypeCaller", parser);
        default_lines.push(htsjdk_vcf::header::HeaderLine::Unstructured {
            key: "source".to_string(),
            value: "MethylationTypeCaller".to_string(),
        });
        default_lines.push(htsjdk_vcf::header::HeaderLine::Structured {
            key: "GATKCommandLine".to_string(),
            fields: vec![
                ("ID".to_string(), "MethylationTypeCaller".to_string()),
                ("CommandLine".to_string(), command_line),
                ("Version".to_string(), crate::TOOLKIT_VERSION.to_string()),
            ],
        });
    }

    let text = gatk_tools::methylation_type_caller::methylation_type_caller(
        &source,
        &mut reference,
        if intervals.is_empty() {
            None
        } else {
            Some(&intervals)
        },
        default_lines,
        flag(parser, "sites-only-vcf-output"),
        // A command line adds to the tool's default filters and can invert or disable them, and a
        // row that keeps no read writes a header and nothing else. Measured on row 13 of this
        // tool's array, where `--inverted-read-filter PrimaryLineReadFilter` keeps only the
        // non-primary reads and the corpus has none.
        &read_filter(parser, &filters, &header)?,
    )
    .map_err(|error| Thrown::user(error.message()))?;
    std::fs::write(&output, &text).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;

    // A variant output carries the same two companions a BAM does, under their own arguments: the
    // index the file's name implies, and the digest APPENDED to the whole name.
    if flag(parser, "create-output-variant-index") {
        // `getBestAvailableSequenceDictionary`, which the index's `DICT:` properties come from: a
        // `--sequence-dictionary` OUTRANKS the reference's own. Measured on row 6 of this tool's
        // array, where `--reference other.fasta` and `--sequence-dictionary matching.dict`
        // disagree and the reference indexed `chr1` where the port indexed `chrOther`.
        let master = master_dictionary(parser)?;
        let lengths: Vec<(String, i32)> = match &master {
            Some(header) => header
                .sequences
                .iter()
                .map(|sequence| (sequence.name.clone(), sequence.length))
                .collect(),
            None => gatk_tools::reference_walker::dictionary(&reference)
                .sequences
                .iter()
                .map(|sequence| (sequence.name.clone(), sequence.length))
                .collect(),
        };
        let index = on_the_fly_index(
            &text,
            &lengths,
            &output,
            text.len() as i64,
            modified_millis(&output),
        );
        let name = format!("{output}.idx");
        std::fs::write(&name, index).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {name}: {error}"))
        })?;
    }
    if flag(parser, "create-output-variant-md5") {
        let digest = format!("{output}.md5");
        std::fs::write(
            &digest,
            gatk_tools::gather_bam_files::md5_file(text.as_bytes()),
        )
        .map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}"))
        })?;
    }
    Ok(None)
}

/// `BaseRecalibrator`, whose output is a GATKReport and whose second input is a set of KNOWN SITES.
///
/// The sites are read whole rather than queried: the counting pass asks, per base, whether the
/// locus is known, and a port that holds them all answers that from memory. Both formats GATK
/// registers for this argument are read here, a BED and a VCF, and the golden's two runs over the
/// same sites in the two formats produce the same table.
///
/// The reference contig is read whole for the same reason `ReadAnonymizer`'s window is: the
/// counting pass compares the read's bases against it, and BAQ looks wider than the read.
pub fn base_recalibrator(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "BaseRecalibrator")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let dictionary = gatk_tools::reference_walker::dictionary(&reference);

    let sites_paths = arguments(parser, "known-sites");
    if sites_paths.is_empty() {
        return Err(Thrown::command_line(
            "Argument known-sites was missing: Argument 'known-sites' is required",
        ));
    }
    let mut known_sites = Vec::new();
    for path in &sites_paths {
        let text = std::fs::read_to_string(path)
            .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(path, false)))?;
        // The codec is chosen by the file's NAME, which is what `FeatureManager` does. Neither is
        // read as an interval argument: these are features, so they are not checked against the
        // reference's dictionary and a site on a contig the reference does not carry is simply a
        // site no read is at.
        if path.ends_with(".bed") {
            for line in text.lines() {
                // `BEDCodec` with `StartOffset.ONE`: the file is half-open and zero-based and the
                // locus it decodes to is neither.
                if let Ok(Some(feature)) =
                    htsjdk_tribble::bed::decode(line, htsjdk_tribble::bed::StartOffset::One)
                {
                    known_sites.push(gatk_engine::interval::SimpleInterval {
                        contig: feature.contig,
                        start: feature.start,
                        end: feature.end,
                    });
                }
            }
        } else {
            // The whole file through the VCF reader, header included: a record's INFO is decoded
            // against the declarations above it, and a reader given an empty header answers
            // nothing at all, which is a run with no known sites and no way to tell.
            let file = htsjdk_vcf::reader::read_vcf(&text)
                .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
            for record in &file.records {
                known_sites.push(gatk_engine::interval::SimpleInterval {
                    contig: record.contig.clone(),
                    start: record.start as i32,
                    end: record.stop as i32,
                });
            }
        }
    }

    // The contig the READS are on, whole: the counting pass compares a read's bases against the
    // reference, so the contig it needs is the reads' and not whichever the reference lists first.
    // Taking the reference's first contig instead sliced a thousand bases of `chrOther` with `chr1`
    // coordinates and panicked.
    //
    // The query is a CLOSURE because the order is observable: the reference is opened at startup
    // and read only once the traversal has produced reads, so a row that hands a BAM the wrong
    // index AND a reference without the reads' contig answers the read failure. Querying first
    // answered the contig instead.
    let wanted = header
        .sequences
        .first()
        .map(|sequence| sequence.name.clone())
        .unwrap_or_default();
    let mut bases = || -> Result<Vec<u8>, gatk_tools::base_recalibrator::BaseRecalibratorError> {
        let length = dictionary
            .sequences
            .iter()
            .find(|sequence| sequence.name == wanted)
            .map(|sequence| sequence.length)
            .ok_or_else(|| {
                gatk_tools::base_recalibrator::BaseRecalibratorError::Reads(
                    gatk_engine::reads::ReadsError::ContigNotInDictionary(wanted.clone()),
                )
            })?;
        reference.query(&wanted, 1, length).map_err(|error| {
            gatk_tools::base_recalibrator::BaseRecalibratorError::Reads(
                gatk_engine::reads::ReadsError::Malformed(format!("{error:?}")),
            )
        })
    };

    let arguments_for_engine = gatk_engine::base_recalibration_engine::EngineArguments {
        covariates: gatk_engine::covariates::RecalibrationArguments {
            mismatches_context_size: number_or(parser, "mismatches-context-size", 2),
            indels_context_size: number_or(parser, "indels-context-size", 3),
            maximum_cycle_value: number_or(parser, "maximum-cycle-value", 500),
            low_qual_tail: u8::try_from(number_or(parser, "low-quality-tail", 2)).unwrap_or(2),
        },
        enable_baq: flag(parser, "enable-baq"),
        compute_indel_bqsr_tables: flag(parser, "compute-indel-bqsr-tables"),
        preserve_qscores_less_than: number_or(parser, "preserve-qscores-less-than", 6),
        default_base_qualities: i8::try_from(number_or(parser, "default-base-qualities", -1))
            .unwrap_or(-1),
        use_original_base_qualities: flag(parser, "use-original-qualities"),
    };
    let filter = read_filter(parser, &filters, &header)?;
    let table = gatk_tools::base_recalibrator::base_recalibrator(
        &source,
        &mut bases,
        &known_sites,
        &arguments_for_engine,
        number_or(parser, "quantizing-levels", 16),
        &filter,
        &intervals,
    )
    .map_err(|error| match error {
        // A read failure is the reference's own exception, class and all: the EOF on a record body
        // is `RuntimeEOFException` at status three, not a user banner carrying a Rust debug string.
        gatk_tools::base_recalibrator::BaseRecalibratorError::Reads(
            gatk_engine::reads::ReadsError::ContigNotInDictionary(contig),
        ) => Thrown::user(format!(
            "Contig {contig} not present in the sequence dictionary {}\n",
            gatk_tools::sequence_dictionary::pretty_print(&dictionary.sequences)
        )),
        gatk_tools::base_recalibrator::BaseRecalibratorError::Reads(error) => {
            reads_traversal_error(error)
        }
        other => Thrown::user(other.message()),
    })?;
    std::fs::write(&output, table).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    // What `doWork` returns, which `handleResult` prints after `Tool returned:`.
    Ok(Some("SUCCESS".to_string()))
}

/// `GtfToBed`, whose BED is one-based because nothing converts the GTF's own coordinates.
///
/// The dictionary is REQUIRED and comes from `--sequence-dictionary`: the tool sorts its rows by
/// the contig's index in it, so a run without one is refused before a line is read.
pub fn gtf_to_bed(parser: &Parser) -> Outcome {
    let input = argument(parser, "gtf-path").ok_or_else(|| {
        Thrown::command_line("Argument gtf-path was missing: Argument 'gtf-path' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let text = std::fs::read_to_string(&input)
        .map_err(|_| Thrown::user(gatk_tools::read_walker_refusal::cannot_read(&input, false)))?;
    // `validateSequenceDictionaries` runs at STARTUP, before a line of the annotation is read: a
    // master dictionary and a `--reference` that disagree are refused there, and the two messages
    // are the dictionary comparison's own -- "No overlapping contigs found" for two that share
    // nothing and "Found contigs with the same name but different lengths" for two that do.
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let (Some(master), Some(path)) =
            (master_dictionary(parser)?, argument(parser, "reference"))
        {
            let mut reference =
                gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&path))
                    .map_err(|error| Thrown::user(format!("{error:?}")))?;
            let theirs = gatk_tools::reference_walker::dictionary(&reference);
            validate_against_master(&master, "reference", &theirs.sequences)?;
            let _ = &mut reference;
        }
    }
    let dictionary = match master_dictionary(parser)? {
        Some(header) => Some(
            header
                .sequences
                .iter()
                .map(|sequence| sequence.name.clone())
                .collect::<Vec<String>>(),
        ),
        None => None,
    };
    let features = gatk_tools::gtf_to_bed::parse_features(&text);
    let bed = gatk_tools::gtf_to_bed::run(
        &features,
        dictionary.as_deref(),
        flag(parser, "sort-by-transcript"),
        // `--use-basic-transcript`, singular. The plural spelling is not an argument at all, so a
        // runner that asked for it read false on every row and answered with the transcripts the
        // reference had dropped.
        flag(parser, "use-basic-transcript"),
    )
    .map_err(|error| match error {
        // The comparator's refusal is an `IllegalArgumentException` and not the tool's own: it
        // reaches the handler as a bug rather than a user error, so the class is printed and the
        // run ends at three.
        gatk_tools::gtf_to_bed::GtfError::UnknownContig { .. } => {
            Thrown::non_user("java.lang.IllegalArgumentException", error.message())
        }
        other => Thrown::user(other.message()),
    })?;
    std::fs::write(&output, bed).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    })?;
    Ok(None)
}

/// A traversal's refusal, told apart by whose exception it is.
///
/// A record that does not decode is htsjdk's `SAMFormatException` and no `UserException` at all,
/// so the handler prints its CLASS and the run ends at status three rather than two. Measured on
/// row 9 of `UnmarkDuplicates`' array, where a BAM is handed the OTHER BAM's index: the query
/// lands mid-file, the next record's length reads as zero, and the port answered a user banner
/// where the reference answered `htsjdk.samtools.SAMFormatException: Invalid record length: 0`.
fn reads_traversal_error(error: gatk_engine::reads::ReadsError) -> Thrown {
    match &error {
        gatk_engine::reads::ReadsError::Malformed(detail) => {
            let message = match detail.strip_prefix("InvalidRecordLength(") {
                Some(rest) => format!("Invalid record length: {}", rest.trim_end_matches(')')),
                None => detail.clone(),
            };
            Thrown::non_user(gatk_tools::read_walker_refusal::SAM_FORMAT, message)
        }
        // A record whose LENGTH decodes and whose body is not there is a different exception in a
        // different package: `BinaryCodec.readBytes` counts what it asked for and what it got, and
        // names the file it was reading. Measured on row 7 of `AddOriginalAlignmentTags`' array,
        // where reads.bam is handed pairs.bai: the seek lands mid-record, four bytes of a read
        // name read as a length of 1178218778, and 320 bytes are left in the file.
        gatk_engine::reads::ReadsError::PrematureEof {
            expected,
            received,
            file,
        } => Thrown::non_user(
            gatk_tools::read_walker_refusal::RUNTIME_EOF,
            format!(
                "Premature EOF. Expected {expected} but only received {received}; \
                 BinaryCodec in readmode; file: {file}"
            ),
        ),
        other => Thrown::user(format!("{other:?}")),
    }
}

/// A BAM and the two files that follow it, written the way every read walker here writes them.
fn write_bam(
    parser: &Parser,
    output: &str,
    bytes: &[u8],
    bai: Option<Vec<u8>>,
) -> Result<Option<String>, Thrown> {
    std::fs::write(output, bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    // htsjdk's writer builds an index only for a `SO:coordinate` file, whatever
    // `--create-output-bam-index` says: an unsorted or queryname BAM is written without one.
    let bai = bai.filter(|_| is_coordinate_sorted(bytes));
    if let Some(bai) = bai {
        let companion = std::path::Path::new(output).with_extension("bai");
        std::fs::write(&companion, bai).map_err(|error| {
            Thrown::non_user(
                PORT_FAILURE,
                format!("could not write {}: {error}", companion.display()),
            )
        })?;
    }
    if flag(parser, "create-output-bam-md5") {
        let digest = format!("{output}.md5");
        std::fs::write(&digest, gatk_tools::gather_bam_files::md5_file(bytes)).map_err(
            |error| Thrown::non_user(PORT_FAILURE, format!("could not write {digest}: {error}")),
        )?;
    }
    Ok(None)
}

/// `GetPileupSummaries.apply`, which writes one row per biallelic SNP the population VCF carries.
///
/// The first tool here whose traversal drives a FEATURE source: `-V` is not an input the tool reads
/// once but the thing each locus is looked up in, and the run is refused before the first locus if
/// that file's header declares no `AF`.
///
/// Three refusals in three different places, which is most of what a row of this tool's array
/// measures: the header check is `onTraversalStart`, the biallelic-SNP and frequency-range tests
/// are `apply`, and "no variant carried an AF at all" is `onTraversalSuccess` -- a run whose
/// records all lack the field is refused AFTER the table is written rather than before it.
pub fn get_pileup_summaries(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "GetPileupSummaries")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let variants = argument(parser, "variant").ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;

    // `read_vcf` takes the file's TEXT and not its path: the reader is the thing being compared,
    // so the filesystem stays out of it and the runner does the reading.
    let text = std::fs::read_to_string(&variants)
        .map_err(|error| Thrown::user(format!("{variants}: {error}")))?;
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    // `getHeaderForFeatures(variants).getInfoHeaderLines()`, which is the check the tool makes
    // before it opens its writer.
    let declares_allele_frequency = file.header.lines.iter().any(|line| {
        matches!(
            line,
            htsjdk_vcf::header::HeaderLine::Compound { key, id, .. } if key == "INFO" && id == "AF"
        )
    });

    // `ReadUtils.getSamplesFromHeader(...).stream().findFirst().get()`, which is the SM of the
    // first read group in the header's own order.
    let sample = header
        .read_groups
        .iter()
        .find_map(|group| group.attributes.get("SM").map(str::to_string))
        .unwrap_or_default();

    let filter = read_filter(parser, &filters, &header)?;
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
    .map_err(locus_traversal_error)?;

    // `featureContext.getValues(variants)` at each locus, which is every record OVERLAPPING it in
    // the file's own order; the tool reads the first and ignores the rest.
    let sites: Vec<(Vec<htsjdk_vcf::variant::VariantContext>, [i32; 4])> = applied
        .iter()
        .map(|one| {
            let overlapping: Vec<htsjdk_vcf::variant::VariantContext> = file
                .records
                .iter()
                .filter(|record| {
                    record.contig == one.context.contig
                        && record.start <= one.context.position as i64
                        && one.context.position as i64 <= record.stop
                })
                .cloned()
                .collect();
            (overlapping, one.context.pileup.base_counts())
        })
        .collect();

    let text = gatk_tools::get_pileup_summaries::run(
        &sites,
        declares_allele_frequency,
        &sample,
        fraction_or(parser, "minimum-population-allele-frequency", 0.01),
        fraction_or(parser, "maximum-population-allele-frequency", 0.2),
    );
    match text {
        Ok(text) => {
            std::fs::write(&output, &text)
                .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
            // `onTraversalSuccess` returns the word rather than a count, and the dispatcher prints
            // whatever a tool returns.
            Ok(Some("SUCCESS".to_string()))
        }
        // `message()` already carries the `Bad input: ` the exception adds, so the runner does not
        // add a second one: the reference prints the prefix once.
        Err(refusal) => Err(Thrown::user(refusal.message())),
    }
}

/// A double-valued argument, or the default the declaration carries.
fn fraction_or(parser: &Parser, long_name: &str, default: f64) -> f64 {
    scalar(parser, long_name)
        .and_then(|text| text.parse().ok())
        .unwrap_or(default)
}

/// `PrintReadsHeader.traverse`, which is the whole tool: it opens the reads for their header.
///
/// `GetSampleName`'s shape with a different answer, and the same startup in front of it: a
/// `GATKTool` declares the interval and dictionary arguments whether or not a traversal uses them,
/// so the refusals they carry are refusals here too.
pub fn print_reads_header(parser: &Parser) -> Outcome {
    let ReadWalkerStart { source, .. } = read_walker_startup(parser, "PrintReadsHeader")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let bytes = gatk_tools::print_reads_header::print_reads_header(&source);
    std::fs::write(&output, &bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

/// `CalculateContamination.doWork`, which reads the table `GetPileupSummaries` writes.
///
/// The first CHAIN this port can run: one tool's output is the other's input, and both ends are
/// compared against the reference. It is a `CommandLineProgram` and no `GATKTool`, so there is no
/// startup to speak of -- no reads, no dictionary, no intervals -- and the whole run is two tables
/// in and one or two out.
///
/// The sample name is the INPUT table's metadata and is what both outputs are written under; a
/// table with no `SAMPLE` metadata tag leaves it empty, which is what the reference does with it.
pub fn calculate_contamination(parser: &Parser) -> Outcome {
    let input = argument(parser, "input").ok_or_else(|| {
        Thrown::command_line("Argument input was missing: Argument 'input' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let (sample, sites) = read_pileup_table(&input)?;
    let matched = match argument(parser, "matched-normal") {
        Some(path) => Some(read_pileup_table(&path)?.1),
        None => None,
    };

    let segmentation_path = argument(parser, "tumor-segmentation");
    let answer = gatk_tools::calculate_contamination::run_from_command_line(
        &sites,
        matched.as_deref(),
        segmentation_path.is_some(),
        fraction_or(parser, "low-coverage-ratio-threshold", 0.1),
        fraction_or(parser, "high-coverage-ratio-threshold", 3.0),
    );

    // The segmentation is written FIRST, before the contamination is calculated at all, so a run
    // that dies in the model still leaves it.
    if let (Some(path), Some(records)) = (&segmentation_path, &answer.segmentation) {
        // The model's record carries the three interval fields flat; the writer's carries a
        // `SimpleInterval`. They are the same record either side of the table layer.
        let rows: Vec<gatk_engine::contamination_tables::MinorAlleleFractionRecord> = records
            .iter()
            .map(
                |record| gatk_engine::contamination_tables::MinorAlleleFractionRecord {
                    segment: gatk_engine::interval::SimpleInterval {
                        contig: record.contig.clone(),
                        start: record.start,
                        end: record.end,
                    },
                    minor_allele_fraction: record.minor_allele_fraction,
                },
            )
            .collect();
        let text = gatk_engine::contamination_tables::write_segments(
            sample.as_deref().unwrap_or_default(),
            &rows,
        );
        std::fs::write(path, &text)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
    }

    let text = gatk_engine::contamination_tables::write_contamination(&[
        gatk_engine::contamination_tables::ContaminationRecord {
            sample: sample.unwrap_or_default(),
            contamination: answer.contamination,
            error: answer.error,
        },
    ]);
    std::fs::write(&output, &text)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(Some("SUCCESS".to_string()))
}

/// `PileupSummary.readFromFile`, with the file read here rather than inside the port.
fn read_pileup_table(
    path: &str,
) -> Result<
    (
        Option<String>,
        Vec<gatk_engine::pileup_summary::PileupSummary>,
    ),
    Thrown,
> {
    let text =
        std::fs::read_to_string(path).map_err(|error| Thrown::user(format!("{path}: {error}")))?;
    gatk_engine::pileup_summary::read_from_file(&text, path).map_err(|error| Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message()),
    })
}

/// `CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection`, in its own order.
///
/// The copy-number tools bin and pad by their OWN arguments, so they refuse every standard
/// interval argument that would modify the input intervals before them. Each is an
/// `IllegalArgumentException` from `Utils.validateArg`, which is a non-user failure and leaves
/// exit 3, and `AnnotateIntervals` calls the same method `PreprocessIntervals` does.
fn validate_copy_number_intervals(parser: &Parser) -> Result<(), Thrown> {
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
    Ok(())
}
/// `AnnotateIntervals.onTraversalStart`, which like the other copy-number tools is the whole run.
///
/// A third interval utility, and the second to call
/// `CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection`: the four standard
/// interval arguments it refuses are the same four `PreprocessIntervals` refuses, for the same
/// reason.
///
/// `--mappability-track` and `--segmental-duplication-track` are refused rather than ignored: each
/// adds a COLUMN to the table from a feature file this runner does not open, and a table missing a
/// column its arguments asked for is a different answer rather than a refusal.
pub fn annotate_intervals(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "AnnotateIntervals")?;
    for track in ["mappability-track", "segmental-duplication-track"] {
        if argument(parser, track).is_some() {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "AnnotateIntervals' --{track} annotates each interval from a feature file this \
                     port does not open yet, and a table without the column it asks for would be a \
                     different answer. This message is the port's own and not GATK's."
                ),
            ));
        }
    }
    validate_copy_number_intervals(parser)?;

    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut source =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    // ONE dictionary here, unlike `PreprocessIntervals`: `getBestAvailableSequenceDictionary`
    // answers both questions, so a `--sequence-dictionary` resolves the intervals AND is written
    // into the `@SQ` lines. Writing the reference's own put an `M5` in ten rows where the
    // reference writes none, which is the master's dictionary showing through.
    let from_dict = reference_dictionary(parser)?;
    let own = gatk_tools::reference_walker::dictionary(&source);
    let best = master_dictionary(parser)?.or(from_dict).unwrap_or(own);
    let written = best.clone();
    let intervals = match interval_arguments(parser, &best)? {
        Some(parameters) => parameters.intervals,
        None => best
            .sequences
            .iter()
            .map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
                    .expect("a contig length is at least one")
            })
            .collect(),
    };

    let mut text = String::from("@HD\tVN:1.6\n");
    for sequence in &written.sequences {
        text.push_str(&format!(
            "@SQ\tSN:{}\tLN:{}",
            sequence.name, sequence.length
        ));
        for key in ["M5", "UR"] {
            if let Some(value) = sequence.attributes.get(key) {
                text.push_str(&format!("\t{key}:{value}"));
            }
        }
        text.push('\n');
    }
    text.push_str(&gatk_tools::annotate_intervals::columns(&["GC_CONTENT"]));
    text.push('\n');
    for interval in &intervals {
        let bases = source
            .query(&interval.contig, interval.start, interval.end)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        text.push_str(&gatk_tools::annotate_intervals::row(
            &interval.contig,
            interval.start,
            interval.end,
            &[gatk_tools::annotate_intervals::gc_content(&bases)],
        ));
        text.push('\n');
    }

    std::fs::write(&output, &text)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(None)
}

/// `CallableLoci`, a locus walker whose two outputs are a BED and a summary of six counts.
///
/// Four things here are the tool's own and none of them are decoration:
///
///   - `emitEmptyLoci()` and `includeNs()` are both true, so the traversal reports EVERY base of
///     its intervals, uncovered ones included, and a locus over an `N` is answered `REF_N` before
///     any depth is counted;
///   - without `-L` the loci come from the REFERENCE's dictionary and not the reads'
///     (`getTraversalIntervals` asks `hasReference()`), so a reference longer than the reads' header
///     claims still produces a line per base;
///   - the single-sample check runs in `onTraversalStart` BEFORE either stream is opened, so a
///     refused run leaves no output file at all, unlike the tools whose writer is built first;
///   - and its six default read filters do not include the walker's own, so `--disable-read-filter`
///     lists them and nothing else.
pub fn callable_loci(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "CallableLoci")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let summary_output = argument(parser, "summary").ok_or_else(|| {
        Thrown::command_line("Argument summary was missing: Argument 'summary' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let known = gatk_tools::reference_walker::dictionary(&reference);

    // `onTraversalStart`'s own check, on the DISTINCT samples in read-group order.
    let mut samples: Vec<String> = Vec::new();
    for group in &header.read_groups {
        if let Some(sample) = group.attributes.get("SM") {
            if !samples.iter().any(|seen| seen == sample) {
                samples.push(sample.to_string());
            }
        }
    }
    if samples.len() != 1 {
        return Err(bad_input(format!(
            "CallableLoci only works for a single sample.  Found {} samples ({}).",
            samples.len(),
            samples.join(", ")
        )));
    }

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;

    // The contigs whole: a state per base would otherwise be a reference query per base.
    let mut bases: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    for (name, length) in reference.sequences().to_vec() {
        let contig = reference
            .query(&name, 1, length as i32)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        bases.insert(name, contig);
    }
    // `getTraversalIntervals()`: the user's intervals, or every interval of the reference.
    let requested: Vec<gatk_engine::interval::SimpleInterval> = if intervals.is_empty() {
        known
            .sequences
            .iter()
            .map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
                    .expect("a contig length is at least one")
            })
            .collect()
    } else {
        intervals.clone()
    };
    let applied = gatk_tools::locus_walker::traverse(
        &records,
        &header,
        None,
        Some(&requested),
        gatk_tools::locus_walker::Options {
            include_deletions: true,
            include_ns: true,
            emit_empty_loci: true,
            max_depth_per_sample: number_or(parser, "max-depth-per-sample", 0),
        },
        &filter,
    )
    .map_err(locus_traversal_error)?;

    // `MissingContigInSequenceDictionary`, raised when a LOCUS asks the reference for its base: the
    // reads come first, and a traversal that visits no locus never asks. `Pileup`'s note carries the
    // rows that measure both halves.
    if let Some(unknown) = applied.iter().find(|one| {
        !known
            .sequences
            .iter()
            .any(|sequence| sequence.name == one.context.contig)
    }) {
        return Err(Thrown::user(format!(
            "Contig {} not present in the sequence dictionary {}\n",
            unknown.context.contig,
            gatk_tools::sequence_dictionary::pretty_print(&known.sequences)
        )));
    }

    let thresholds = gatk_tools::callable_loci::Arguments {
        max_low_mapq: number_or(parser, "max-low-mapq", 1),
        min_mapping_quality: number_or(parser, "min-mapping-quality", 10),
        min_base_quality: number_or(parser, "min-base-quality", 20),
        min_depth: number_or(parser, "min-depth", 4),
        // `maxDepth` is an `Integer` and not an `int`: unset is null, and null is not a threshold
        // rather than a threshold of zero. A run with `--max-depth 0` calls every covered locus
        // excessive, and a run without it calls none.
        max_depth: scalar(parser, "max-depth").and_then(|text| text.parse().ok()),
        min_depth_low_mapq: number_or(parser, "min-depth-for-low-mapq", 10),
        max_low_mapq_fraction: fraction_or(parser, "max-fraction-of-reads-with-low-mapq", 0.1),
    };
    // A locus with no reference base at all cannot be read from the map, and the reference would
    // have refused before the traversal; `N` is what the state machine answers for one.
    let loci: Vec<(String, i32, gatk_tools::callable_loci::State)> = applied
        .iter()
        .map(|one| {
            let base = bases
                .get(&one.context.contig)
                .and_then(|contig| contig.get((one.context.position - 1) as usize))
                .copied()
                .unwrap_or(b'N');
            let pileup: Vec<gatk_tools::callable_loci::Element> = one
                .context
                .pileup
                .elements
                .iter()
                .map(|element| gatk_tools::callable_loci::Element {
                    mapping_quality: i32::from(element.read.mapping_quality),
                    base_quality: element.qual() as i32,
                    is_deletion: element.is_deletion(),
                })
                .collect();
            (
                one.context.contig.clone(),
                one.context.position,
                gatk_tools::callable_loci::state_at(base, &pileup, &thresholds),
            )
        })
        .collect();

    let format = match scalar(parser, "format").as_deref() {
        Some("STATE_PER_BASE") => gatk_tools::callable_loci::OutputFormat::StatePerBase,
        _ => gatk_tools::callable_loci::OutputFormat::Bed,
    };
    let (bed, summary) = gatk_tools::callable_loci::write(&loci, format);
    write_file(&output, bed.as_bytes())?;
    write_file(&summary_output, summary.as_bytes())?;
    Ok(None)
}

/// `ShiftFasta`, which is a `GATKTool` with `traverse` overridden and therefore no walker at all.
///
/// Its four outputs are written from three arguments: `-O` takes the FASTA and the index and
/// dictionary htsjdk writes beside it, `--shift-back-output` the chain, and `--interval-file-name`
/// the base name of a PAIR of interval lists. Two of those are opened in `onTraversalStart`, so a
/// refused traversal still leaves them behind, empty: see `docs`'s note on a writer built at
/// startup.
pub fn shift_fasta(parser: &Parser) -> Outcome {
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let chain_output = argument(parser, "shift-back-output").ok_or_else(|| {
        Thrown::command_line(
            "Argument shift-back-output was missing: Argument 'shift-back-output' is required",
        )
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let interval_name = argument(parser, "interval-file-name");

    // A `GATKTool` loads the master dictionary, opens the reference, resolves the intervals and
    // validates the dictionaries whether or not its traversal uses any of them, so those refusals
    // are this tool's refusals too.
    let master = master_dictionary(parser)?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let theirs = gatk_tools::reference_walker::dictionary(&reference);
    // `initializeIntervals` runs BEFORE `validateSequenceDictionaries`, so an interval that the best
    // available dictionary does not carry is refused before the two dictionaries are compared at
    // all. Measured on row 7 of this tool's array, where `--sequence-dictionary other.dict` and
    // `--reference reference.fasta` share no contig AND `-L chr1:50000-60000` is off the master: the
    // reference answered `Badly formed genome unclippedLoc` and the port answered the mismatch.
    let best = master.clone().unwrap_or_else(|| theirs.clone());
    let _ = interval_arguments(parser, &best)?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            validate_against_master(master, "reference", &theirs.sequences)?;
        }
    }

    let offsets: Vec<i32> = arguments(parser, "shift-offset-list")
        .iter()
        .filter_map(|text| text.parse().ok())
        .collect();
    let width = number_or(parser, "line-width", 60).max(0) as usize;

    // The three writers `onTraversalStart` builds, which exist before the offset list is checked:
    // the FASTA's (with its `.fai` and `.dict`), the chain's, and the pair the interval base name
    // asks for. A refused traversal closes all of them, so the run leaves an empty FASTA, an empty
    // index, a dictionary of nothing but its `@HD` line, and three empty text files.
    let write_companions = |chain: &str, plain: &str, shifted: &str| -> Result<(), Thrown> {
        write_file(&chain_output, chain.as_bytes())?;
        if let Some(name) = &interval_name {
            write_file(&format!("{name}.intervals"), plain.as_bytes())?;
            write_file(&format!("{name}.shifted.intervals"), shifted.as_bytes())?;
        }
        Ok(())
    };
    let refused = |error: gatk_tools::shift_fasta::ShiftError| -> Thrown {
        match gatk_tools::fasta_reference_maker::empty_outputs(width) {
            Ok(empty) => {
                let _ = write_outputs(&output, &empty);
                let _ = write_companions("", "", "");
                shift_error(error)
            }
            Err(failure) => fasta_maker_error(failure),
        }
    };
    let outputs = gatk_tools::shift_fasta::run(&mut reference, &offsets, width).map_err(refused)?;

    write_outputs(&output, &outputs.reference)?;
    write_companions(
        &outputs.chain,
        &outputs.intervals,
        &outputs.shifted_intervals,
    )?;
    Ok(None)
}

/// What `ShiftFasta` refused with, told apart by whose refusal it is.
fn shift_error(error: gatk_tools::shift_fasta::ShiftError) -> Thrown {
    match &error {
        // The writer's refusals are htsjdk's exceptions and no `UserException` at all, so the handler
        // prints the class and the run ends at three. A `--shift-offset-list 0` shifts no contig,
        // which leaves the writer with nothing to close: `no sequences were added to the reference`.
        // Measured on rows 4 and 9 of this tool's array, where the port had answered at zero with
        // four empty files.
        gatk_tools::shift_fasta::ShiftError::Writer(_) => {
            Thrown::non_user(error.java_class(), error.message())
        }
        // The tool's own two, whose banner the handler prints at two.
        _ => Thrown {
            failure: Failure::User,
            exception: error.java_class(),
            message: Some(error.message()),
        },
    }
}

/// `TransferReadTags`, which is a `GATKTool` with `traverse()` overridden and TWO reads sources.
///
/// The second one is opened by the tool and not by the engine: `new ReadsPathDataSource(path)` with
/// no index, no interval bound and no filter, which is why the unmapped side is read whole here.
/// The engine's own source is reached through `directlyAccessEngineReadsDataSource().iterator()`,
/// so the filter chain the command line resolved is SELECTED and never consulted: a
/// `--read-filter` on this tool changes what `--disable-read-filter` lists and nothing else.
pub fn transfer_read_tags(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "TransferReadTags")?;
    // Resolved for its refusals, which are the parser's, and then not applied: see above.
    let _ = read_filter(parser, &filters, &header)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let unmapped_path = argument(parser, "unmapped-sam").ok_or_else(|| {
        Thrown::command_line(
            "Argument unmapped-sam was missing: Argument 'unmapped-sam' is required",
        )
    })?;
    let read_tags = arguments(parser, "read-tags");
    let unmapped =
        gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(&unmapped_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let command_line = crate::command_line::expanded("TransferReadTags", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let run = gatk_tools::transfer_read_tags::transfer_read_tags_with(
        &source, &unmapped, &read_tags, &options, level, deflater,
    )
    .map_err(reads_traversal_error)?;
    match run {
        Ok((bytes, bai)) => {
            write_bam(parser, &output, &bytes, bai)?;
            Ok(Some("SUCCESS".to_string()))
        }
        Err(refusal) => Err(transfer_error(refusal)),
    }
}

/// Which exception each of `TransferReadTags`' refusals arrives as.
///
/// Only one of the five is a `UserException`. `Utils.validate` throws `IllegalStateException` and
/// `Utils.nonNull` and `Utils.nonEmpty` throw `IllegalArgumentException`, so four of the five reach
/// the handler as a bug rather than a refusal: the class is printed and the run ends at three.
fn transfer_error(error: gatk_tools::transfer_read_tags::TransferError) -> Thrown {
    use gatk_tools::transfer_read_tags::TransferError;
    let message = error.message();
    match error {
        TransferError::UnmappedEmptyAndAlignedIsNot => Thrown::user(message),
        TransferError::AlignedNotQueryNameSorted | TransferError::NotInUnmapped { .. } => {
            Thrown::non_user("java.lang.IllegalStateException", message)
        }
        TransferError::NoReadTags | TransferError::AttributeEmpty { .. } => {
            Thrown::non_user("java.lang.IllegalArgumentException", message)
        }
    }
}

/// `PostProcessReadsForRSEM`, the second `GATKTool` here that groups the traversal itself.
///
/// Its default read filter is a SINGLETON that is not the walker's -- `NOT_SUPPLEMENTARY_ALIGNMENT`
/// alone, with no wellformed filter at all -- and `traverse()` calls `makeReadFilter()`, so the
/// chain a row builds is the chain that decides which reads reach the query-name groups.
pub fn post_process_reads_for_rsem(parser: &Parser) -> Outcome {
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "PostProcessReadsForRSEM")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("PostProcessReadsForRSEM", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    let run = gatk_tools::post_process_reads_for_rsem::post_process_reads_for_rsem_with(
        &source, &options, &filter, level, deflater,
    )
    .map_err(reads_traversal_error)?;
    match run {
        Ok((bytes, bai)) => {
            write_bam(parser, &output, &bytes, bai)?;
            Ok(Some("SUCCESS".to_string()))
        }
        // Four of the five refusals are the tool's own `UserException`s; the fifth is the JVM's
        // `NullPointerException`, raised from inside a guard that exists because the value may be
        // null and dereferences it anyway.
        Err(refusal) => Err(match refusal {
            gatk_tools::post_process_reads_for_rsem::RsemError::NullDereference(_) => {
                Thrown::non_user("java.lang.NullPointerException", refusal.message())
            }
            gatk_tools::post_process_reads_for_rsem::RsemError::PrimaryAlreadySet { .. } => {
                Thrown::non_user("java.lang.IllegalStateException", refusal.message())
            }
            other => Thrown::user(other.message()),
        }),
    }
}

/// `VariantsToTable`, the first `VariantWalker` here whose output is a TABLE.
///
/// The driving variants are the traversal and `-F`, `-GF`, `-ASF` and `-ASGF` decide the columns,
/// so the tool reads what a VCF DECLARES as much as what it holds: a field's `Number` says whether
/// its value is a list, and an `R` there is the only thing the allele-specific arguments treat
/// differently.
///
/// Three things `onTraversalStart` does that the table's shape depends on:
///
///   - with none of the four field arguments given, the columns become every mandatory field except
///     INFO, then every INFO id the header declares, then every FORMAT id with `GT` moved FIRST;
///   - the samples are a SORTED set, and asking for no genotype field at all empties it, which is
///     what keeps a genotype column out of a table nobody asked one for;
///   - and a file with no samples and no fields is refused rather than answered with an empty table.
pub fn variants_to_table(parser: &Parser) -> Outcome {
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "VariantsToTable")?;

    let declared = vcf_declarations(&text);
    let mut fields = arguments(parser, "fields");
    let mut genotype_fields = arguments(parser, "genotype-fields");
    let allele_specific_fields = arguments(parser, "asFieldsToTake");
    let mut allele_specific_genotype_fields = arguments(parser, "asGenotypeFieldsToTake");
    if fields.is_empty()
        && genotype_fields.is_empty()
        && allele_specific_fields.is_empty()
        && allele_specific_genotype_fields.is_empty()
    {
        // `VCFHeader.HEADER_FIELDS.values()` minus INFO, in the enum's own order.
        fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
            .iter()
            .map(|name| (*name).to_string())
            .collect();
        fields.extend(declared.info.clone());
        for id in &declared.format {
            if id == "GT" {
                genotype_fields.insert(0, id.clone());
            } else {
                genotype_fields.push(id.clone());
            }
        }
    }

    // `VcfUtils.getSortedSampleSet`, which is a `TreeSet` and therefore sorted.
    let mut samples: Vec<String> =
        if genotype_fields.is_empty() && allele_specific_genotype_fields.is_empty() {
            Vec::new()
        } else {
            let mut names = vcf_samples(&text);
            names.sort();
            names.dedup();
            names
        };
    if samples.is_empty()
        && !(genotype_fields.is_empty() && allele_specific_genotype_fields.is_empty())
    {
        genotype_fields.clear();
        allele_specific_genotype_fields.clear();
        if fields.is_empty() && allele_specific_fields.is_empty() {
            return Err(Thrown::user(
                "There are no samples and no fields - no output will be produced",
            ));
        }
    }
    if genotype_fields.is_empty() && allele_specific_genotype_fields.is_empty() {
        samples.clear();
    }

    let table_arguments = gatk_tools::variants_to_table::Arguments {
        fields,
        genotype_fields,
        allele_specific_fields,
        allele_specific_genotype_fields,
        split_multi_allelic: flag(parser, "split-multi-allelic"),
        show_filtered: flag(parser, "show-filtered"),
        moltenize: flag(parser, "moltenize"),
        error_if_missing_data: flag(parser, "error-if-missing-data"),
    };

    let mut out = String::new();
    for column in gatk_tools::variants_to_table::header(&samples, &table_arguments) {
        if !out.is_empty() && !out.ends_with('\n') {
            out.push('\t');
        }
        out.push_str(&column);
    }
    out.push('\n');

    let all = vcf_records(&text, &samples_in_file_order(&text), &declared);
    // `-L` is the traversal and not a filter the tool applies: the driving variants are queried by
    // interval, which is why an input with no random access is refused before a record is read.
    // Measured on rows 0, 4 and 12 of this tool's array, where the reference emitted the ONE record
    // the window holds and the port emitted all eight.
    let spans: Vec<Locus> = all
        .iter()
        .map(|(record, _)| Locus {
            contig: record.contig.clone(),
            start: record.start,
            stop: record.start + record.reference.len() as i32 - 1,
        })
        .collect();
    if gatk_engine::variant_source::intervals_for_traversal(intervals.as_deref()).is_some()
        && !has_feature_index(&input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.clone(),
            }
            .message(),
        ));
    }
    let kept: Vec<usize> = gatk_engine::variant_source::traverse(&spans, intervals.as_deref())
        .into_iter()
        .map(|span| {
            spans
                .iter()
                .position(|other| std::ptr::eq(other, span))
                .expect("a span of this list")
        })
        .collect();
    let records: Vec<&(gatk_tools::variants_to_table::Record, String)> =
        kept.into_iter().map(|index| &all[index]).collect();

    let mut emitted = 0_usize;
    for (record, raw) in records {
        // `showFiltered || vc.isNotFiltered()`: a record whose FILTER is neither `.` nor `PASS` is
        // skipped, and the counter that numbers the moltenized rows never sees it.
        if !table_arguments.show_filtered && !record.filters.is_empty() {
            continue;
        }
        emitted += 1;
        let rows = gatk_tools::variants_to_table::extract_fields(
            record,
            &samples,
            &table_arguments,
            &declared.per_allele,
        )
        .map_err(|missing| {
            // `String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc)`: the
            // source is the driving input's NAME, which for a `-V` with no logical name is
            // `Unknown`, and the third is the whole `VariantContext.toString()`.
            Thrown::user(format!(
                "Missing field {} in vc Unknown at {}",
                missing.field,
                variant_context_to_string(record, raw)
            ))
        })?;
        for row in rows {
            if table_arguments.moltenize {
                let mut index = 0;
                for field in &table_arguments.fields {
                    out.push_str(&format!("{emitted}\tsite\t{field}\t{}\n", row[index]));
                    index += 1;
                }
                for sample in &samples {
                    for field in &table_arguments.genotype_fields {
                        out.push_str(&format!(
                            "{emitted}\t{}\t{field}\t{}\n",
                            sample.replace(' ', "_"),
                            row[index]
                        ));
                        index += 1;
                    }
                }
            } else {
                out.push_str(&row.join("\t"));
                out.push('\n');
            }
        }
    }

    // `-O` is optional and a null one is `System.out`, which is the same shape `CheckPileup` has.
    write_report(&argument(parser, "output"), &out)?;
    Ok(None)
}

/// What a VCF's header DECLARES, which is what decides a value's shape.
struct VcfDeclarations {
    /// The INFO ids, in the order the header wrote them.
    info: Vec<String>,
    /// The FORMAT ids, in the same order.
    format: Vec<String>,
    /// Every id whose `Number` is `R`, which is what `-ASF` and `-ASGF` split.
    per_allele: gatk_tools::variants_to_table::CountTypes,
    /// Every id whose `Number` is not `1`, whose value is therefore a list.
    lists: std::collections::HashSet<String>,
}

fn vcf_declarations(text: &str) -> VcfDeclarations {
    let mut declared = VcfDeclarations {
        info: Vec::new(),
        format: Vec::new(),
        per_allele: gatk_tools::variants_to_table::CountTypes::new(),
        lists: std::collections::HashSet::new(),
    };
    for line in text.lines() {
        let (kind, rest) = if let Some(rest) = line.strip_prefix("##INFO=<") {
            ("INFO", rest)
        } else if let Some(rest) = line.strip_prefix("##FORMAT=<") {
            ("FORMAT", rest)
        } else {
            continue;
        };
        let field = |key: &str| -> Option<String> {
            rest.split(',')
                .find_map(|entry| entry.trim().strip_prefix(&format!("{key}=")))
                .map(|value| value.trim_end_matches('>').to_string())
        };
        let Some(id) = field("ID") else { continue };
        let number = field("Number").unwrap_or_default();
        if kind == "INFO" {
            declared.info.push(id.clone());
        } else {
            declared.format.push(id.clone());
        }
        declared.per_allele.insert(id.clone(), number == "R");
        if number != "1" {
            declared.lists.insert(id);
        }
    }
    declared
}

/// The sample names, in the order the `#CHROM` line writes them.
fn samples_in_file_order(text: &str) -> Vec<String> {
    text.lines()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| {
            line.split('\t')
                .skip(9)
                .map(|name| name.to_string())
                .collect()
        })
        .unwrap_or_default()
}

fn vcf_samples(text: &str) -> Vec<String> {
    samples_in_file_order(text)
}

/// The records, projected to what the table reads.
///
/// A value is a LIST where the file wrote commas and the header said its `Number` is not one, which
/// is the same rule the tool's own suite uses: the table prints a list differently from a string
/// that happens to hold a comma.
fn vcf_records(
    text: &str,
    samples: &[String],
    declared: &VcfDeclarations,
) -> Vec<(gatk_tools::variants_to_table::Record, String)> {
    let value = |text: &str, key: &str| -> gatk_tools::variants_to_table::Value {
        if declared.lists.contains(key) && text.contains(',') {
            gatk_tools::variants_to_table::Value::Many(
                text.split(',').map(|part| part.to_string()).collect(),
            )
        } else {
            gatk_tools::variants_to_table::Value::One(text.to_string())
        }
    };
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let info = field
                .get(7)
                .map(|column| {
                    column
                        .split(';')
                        .filter_map(|entry| entry.split_once('='))
                        .map(|(key, text)| (key.to_string(), value(text, key)))
                        .collect()
                })
                .unwrap_or_default();
            let keys: Vec<&str> = field
                .get(8)
                .map(|f| f.split(':').collect())
                .unwrap_or_default();
            let genotypes = (0..samples.len())
                .map(|index| {
                    field
                        .get(9 + index)
                        .map(|column| {
                            column
                                .split(':')
                                .enumerate()
                                .filter_map(|(at, text)| {
                                    keys.get(at)
                                        .map(|key| ((*key).to_string(), value(text, key)))
                                })
                                .collect()
                        })
                        .unwrap_or_default()
                })
                .collect();
            // The FORMAT column and the sample columns as the file wrote them, which is what a
            // record's `toString` prints while its genotypes are still lazy.
            let raw = field[8..].join("\t");
            (
                gatk_tools::variants_to_table::Record {
                    contig: field[0].to_string(),
                    start: field[1].parse().unwrap_or(0),
                    id: field[2].to_string(),
                    reference: field[3].to_string(),
                    alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                    qual: field[5].parse().ok(),
                    filters: match field[6] {
                        "." | "PASS" => Vec::new(),
                        names => names.split(';').map(|name| name.to_string()).collect(),
                    },
                    info,
                    genotypes,
                },
                raw,
            )
        })
        .collect()
}

/// `CompareBaseQualities`, which is no GATK tool at all.
///
/// It extends `PicardCommandLineProgram`, so its namespace is Picard's argument set rather than the
/// engine's: no read filter, no interval, no sequence-dictionary validation, and its two SAM files
/// arrive as POSITIONAL arguments rather than under `--input`. Both readers are opened by hand and
/// each is wrapped in a `SecondaryOrSupplementarySkippingIterator`, so the skipping happens per file
/// rather than over the pair.
///
/// The tool RETURNS its verdict: `hasNonDiagonalElements() ? 1 : 0`, which the dispatcher prints and
/// which is not an exit status. `--throw-on-diff` turns the same fact into a refusal instead.
pub fn compare_base_qualities(parser: &Parser) -> Outcome {
    let files = parser.positional_values();
    // The parser has already refused any count but two, so this is a read of what it collected.
    let (first, second) = (files[0].clone(), files[1].clone());

    let read = |path: &str| -> Result<Vec<htsjdk_bam::record::BamRecord>, Thrown> {
        let source =
            gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        source.iter_all().map_err(reads_traversal_error)
    };
    let left = read(&first)?;
    let right = read(&second)?;

    let arguments = gatk_tools::compare_base_qualities::CompareArguments {
        static_quantization_quals: arguments(parser, "static-quantized-quals")
            .iter()
            .filter_map(|value| value.parse().ok())
            .collect(),
        round_down: flag(parser, "round-down-quantized"),
        throw_on_diff: flag(parser, "throw-on-diff"),
    };
    let result =
        gatk_tools::compare_base_qualities::compare_base_qualities(&left, &right, &arguments)
            .map_err(|refusal| match refusal {
                // `--round-down-quantized` alone is the PARSER's refusal and not the tool's, so it carries
                // the argument's name and the bad value rather than a banner of its own.
                gatk_tools::compare_base_qualities::CompareError::RoundDownAlone => {
                    Thrown::command_line(refusal.message())
                }
                gatk_tools::compare_base_qualities::CompareError::QualitiesDiffer => {
                    Thrown::user(refusal.message())
                }
                other => bad_input(other.message()),
            })?;

    // The report is written where `printOutResults` writes it: the file `-O` names, or stdout.
    write_report(&argument(parser, "output"), &result.report)?;
    Ok(Some(result.exit_code.to_string()))
}

/// `VariantContext.toString()`, which a refusal prints whole.
///
/// The lazy branch, `toStringUnparsedGenotypes`, because a record read from a file keeps its
/// genotype text until something decodes it: what the string carries is the FORMAT column and the
/// sample columns as they were written, tabs and all. Measured on row 9 of `VariantsToTable`'s
/// array, which is the only place this port prints one.
fn variant_context_to_string(
    record: &gatk_tools::variants_to_table::Record,
    raw_genotypes: &str,
) -> String {
    let stop = record.start + record.reference.len() as i32 - 1;
    let position = if stop == record.start {
        format!("{}:{}", record.contig, record.start)
    } else {
        format!("{}:{}-{}", record.contig, record.start, stop)
    };
    // `hasLog10PError()`, which is false for a QUAL the file wrote as `.`.
    let qual = match record.qual {
        Some(value) => format!("{value:.2}"),
        None => ".".to_string(),
    };
    // `ParsingUtils.sortList(getAlleles())`: the reference allele carries a `*`, and the list is
    // sorted by `Allele.compareTo`, which puts the reference first and the rest by their bases.
    let mut alleles: Vec<String> = record.alternates.clone();
    alleles.sort();
    let alleles = std::iter::once(format!("{}*", record.reference))
        .chain(alleles)
        .collect::<Vec<String>>()
        .join(", ");
    // `ParsingUtils.sortedString(getAttributes())`: a `TreeMap`'s `toString`, so `{}` when empty and
    // `{k=v, k=v}` sorted by key otherwise.
    let mut attributes: Vec<String> = record
        .info
        .iter()
        .map(|(key, value)| format!("{key}={}", attribute_to_string(value)))
        .collect();
    attributes.sort();
    let attributes = format!("{{{}}}", attributes.join(", "));
    format!(
        "[VC Unknown @ {position} Q{qual} of type={} alleles=[{alleles}] attr={attributes} GT={} filters={}",
        variant_type_name(&record.reference, &record.alternates),
        raw_genotypes,
        record.filters.join(",")
    )
}

/// An attribute's `toString`, which for a list is Java's `[a, b]` and not the comma join a column
/// carries.
fn attribute_to_string(value: &gatk_tools::variants_to_table::Value) -> String {
    match value {
        gatk_tools::variants_to_table::Value::One(one) => one.clone(),
        gatk_tools::variants_to_table::Value::Many(many) => format!("[{}]", many.join(", ")),
    }
}

/// `VariantContext.getType()`, by the same rule [`gatk_tools::remove_nearby_indels`] ports: length
/// decides, and two alternates that disagree are MIXED.
fn variant_type_name(reference: &str, alternates: &[String]) -> &'static str {
    if alternates.is_empty() || (alternates.len() == 1 && alternates[0] == ".") {
        return "NO_VARIATION";
    }
    let mut kind: Option<&'static str> = None;
    for alternate in alternates {
        let this =
            if alternate.starts_with('<') || alternate.contains('[') || alternate.contains(']') {
                "SYMBOLIC"
            } else if alternate.len() == reference.len() {
                if reference.len() == 1 {
                    "SNP"
                } else {
                    "MNP"
                }
            } else {
                "INDEL"
            };
        match kind {
            None => kind = Some(this),
            Some(seen) if seen == this => {}
            Some(_) => return "MIXED",
        }
    }
    kind.unwrap_or("NO_VARIATION")
}

/// The records a variant walker's traversal reaches, which is a QUERY by interval and not a filter.
///
/// Shared by the tools that write a VCF, because getting it wrong is invisible in a file that has
/// one contig: the records outside the window are simply absent from the reference's output.
fn variants_in_traversal<'a>(
    records: &'a [htsjdk_vcf::variant::VariantContext],
    intervals: Option<&[gatk_engine::interval::SimpleInterval]>,
    input: &str,
) -> Result<Vec<&'a htsjdk_vcf::variant::VariantContext>, Thrown> {
    let spans: Vec<Locus> = records
        .iter()
        .map(|record| Locus {
            contig: record.contig.clone(),
            start: record.start as i32,
            stop: record.stop as i32,
        })
        .collect();
    if gatk_engine::variant_source::intervals_for_traversal(intervals).is_some()
        && !has_feature_index(input)
    {
        return Err(Thrown::user(
            gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                path: input.to_string(),
            }
            .message(),
        ));
    }
    Ok(gatk_engine::variant_source::traverse(&spans, intervals)
        .into_iter()
        .map(|span| {
            let index = spans
                .iter()
                .position(|other| std::ptr::eq(other, span))
                .expect("a span of this list");
            &records[index]
        })
        .collect())
}

/// The two lines `getDefaultToolVCFHeaderLines` adds, which `--add-output-vcf-command-line` gates.
fn default_tool_vcf_header_lines(
    parser: &Parser,
    tool: &str,
) -> Vec<htsjdk_vcf::header::HeaderLine> {
    if !flag(parser, "add-output-vcf-command-line") {
        return Vec::new();
    }
    vec![
        htsjdk_vcf::header::HeaderLine::Unstructured {
            key: "source".to_string(),
            value: tool.to_string(),
        },
        htsjdk_vcf::header::HeaderLine::Structured {
            key: "GATKCommandLine".to_string(),
            fields: vec![
                ("ID".to_string(), tool.to_string()),
                (
                    "CommandLine".to_string(),
                    crate::command_line::expanded(tool, parser),
                ),
                ("Version".to_string(), crate::TOOLKIT_VERSION.to_string()),
            ],
        },
    ]
}

/// `RemoveNearbyIndels`, the first tool here that writes a VCF of its own records.
///
/// The buffer holds at most one indel and remembers one it has already thrown away, so three indels
/// in a row lose all three; `onTraversalSuccess` keeps the last one on a reference comparison rather
/// than an equality. Both are the port's business; what the runner adds is the header the output
/// carries and the traversal that decides which records reach the buffer at all.
pub fn remove_nearby_indels(parser: &Parser) -> Outcome {
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "RemoveNearbyIndels")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    // `--min-indel-spacing` is REQUIRED on this tool, which is unusual for a numeric argument with
    // a default: the annotation says `optional = false`, so the parser refuses a command line
    // without it however sensible the default looks.
    let spacing = number_or(parser, "min-indel-spacing", 1);

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let records: Vec<htsjdk_vcf::variant::VariantContext> = kept.into_iter().cloned().collect();

    let mut header = file.header.clone();
    for line in default_tool_vcf_header_lines(parser, "RemoveNearbyIndels") {
        header.lines.push(line);
    }
    let emitted = gatk_tools::remove_nearby_indels::remove_nearby_indels(&records, spacing);
    let written: Vec<htsjdk_vcf::variant::VariantContext> = emitted
        .into_iter()
        .map(|index| records[index].clone())
        .collect();
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> =
        written.into_iter().filter(|record| keep(record)).collect();
    apply_sites_only(parser, &mut header, &mut written);
    let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;

    write_variant_output(parser, &output, &out)?;
    // `onTraversalSuccess` returns the word, which `handleResult` prints.
    Ok(Some("SUCCESS".to_string()))
}

/// The predicate `--variant-output-filtering` builds, which every emitted record passes through.
type VariantOutputFilter<'a> = Box<dyn Fn(&htsjdk_vcf::variant::VariantContext) -> bool + 'a>;

/// `--variant-output-filtering`, which WRAPS the VCF writer rather than the traversal.
///
/// `IntervalFilteringVcfWriter` tests every record the tool emits against the user intervals, so a
/// record the traversal reached can still be kept out of the file: `STARTS_IN` looks at the start
/// position alone, `ENDS_IN` at the end, `CONTAINED` needs one overlapping interval to hold the
/// whole record, and `ANYWHERE` is the default and no filter at all. Measured on rows 10, 13 and 22
/// of `RemoveNearbyIndels`' array, where an indel at 1000 spanning four bases starts inside
/// `chr1:1-1000` and ends outside it.
///
/// A mode other than `ANYWHERE` with no `-L` is refused, after the dictionaries are validated.
fn variant_output_filter<'a>(
    parser: &Parser,
    intervals: Option<&'a [gatk_engine::interval::SimpleInterval]>,
) -> Result<VariantOutputFilter<'a>, Thrown> {
    let mode = scalar(parser, "variant-output-filtering").unwrap_or_else(|| "ANYWHERE".to_string());
    if mode == "ANYWHERE" {
        return Ok(Box::new(|_| true));
    }
    let Some(intervals) = intervals.filter(|list| !list.is_empty()) else {
        return Err(Thrown::command_line(
            "Argument -L or -XL was missing: Intervals are required if --variant-output-filtering \
             was specified or if the tool uses interval filtering.",
        ));
    };
    Ok(Box::new(move |record| {
        let start = record.start as i32;
        let stop = record.stop as i32;
        intervals.iter().any(|interval| {
            interval.contig == record.contig
                && match mode.as_str() {
                    "STARTS_IN" => interval.start <= start && start <= interval.end,
                    "ENDS_IN" => interval.start <= stop && stop <= interval.end,
                    "CONTAINED" => interval.start <= start && stop <= interval.end,
                    // `OVERLAPS`, and anything the parser would have refused before this.
                    _ => interval.start <= stop && start <= interval.end,
                }
        })
    }))
}

/// `--sites-only-vcf-output`, which builds the writer with `DO_NOT_WRITE_GENOTYPES`.
///
/// A header with no samples is what that writes: the `#CHROM` line stops at INFO and no record
/// carries a FORMAT column. The `##FORMAT` declarations stay, because the option drops the
/// genotypes and not the lines that describe them. Measured on rows 0 and 5 of
/// `RemoveNearbyIndels`' array, where the port wrote a genotype column the reference did not.
fn apply_sites_only(
    parser: &Parser,
    header: &mut htsjdk_vcf::header::VcfHeader,
    records: &mut [htsjdk_vcf::variant::VariantContext],
) {
    if !flag(parser, "sites-only-vcf-output") {
        return;
    }
    header.samples.clear();
    for record in records {
        record.genotypes = Vec::new().into();
    }
}

/// The `##contig` lines of a header, as the pairs an index's `DICT:` properties want.
fn sequence_dictionary_of(header: &htsjdk_vcf::header::VcfHeader) -> Vec<(String, i32)> {
    header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Contig { fields, .. } => {
                let id = fields.iter().find(|(key, _)| key == "ID")?.1.clone();
                let length = fields
                    .iter()
                    .find(|(key, _)| key == "length")
                    .and_then(|(_, value)| value.parse().ok())
                    .unwrap_or(0);
                Some((id, length))
            }
            _ => None,
        })
        .collect()
}

/// `UpdateVCFSequenceDictionary`, which replaces a header's dictionary and passes the records
/// through.
///
/// Its `getBestAvailableSequenceDictionary` is OVERRIDDEN, so the dictionary every caller sees is
/// the new one and not the VCF's own: that is what makes the index written beside the output carry
/// the source's contigs. The four kinds of file `--source-dictionary` accepts are read here, because
/// `SAMSequenceDictionaryExtractor` takes a dictionary out of a variant, an alignment, a `.dict` or
/// a reference.
///
/// The records are written AS THEY GO, so a refusal leaves the ones before it on disk. Reproduced,
/// because a refused run's output file is part of what a row compares.
pub fn update_vcf_sequence_dictionary(parser: &Parser) -> Outcome {
    // `getBestAvailableSequenceDictionary` is overridden, and the FIRST thing to call it is
    // `initializeIntervals`, which runs before `validateSequenceDictionaries`. So a command line
    // that names both dictionaries is refused before the two it names are compared to anything, and
    // the refusal is a `CommandLineException` rather than a `UserException`: status one, not two.
    // Measured on rows 0 and 4 of this tool's array.
    if argument(parser, "source-dictionary").is_some()
        && argument(parser, "sequence-dictionary").is_some()
    {
        return Err(Thrown::command_line(
            gatk_tools::update_vcf_sequence_dictionary::UpdateDictionaryError::TwoDictionaries
                .message(),
        ));
    }
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "UpdateVCFSequenceDictionary")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;

    let refuse =
        |error: gatk_tools::update_vcf_sequence_dictionary::UpdateDictionaryError| -> Thrown {
            // Four of the seven are `CommandLineException`s, which exit at ONE; the one-argument
            // `BadArgumentValue` also carries the `Illegal argument value: ` its constructor
            // prefixes. Measured on rows 5 and 7 of this tool's array.
            if error.java_class().starts_with("org.broadinstitute.barclay") {
                let message = error.message();
                return Thrown::command_line(
                    if error
                        .java_class()
                        .ends_with("CommandLineException$BadArgumentValue")
                    {
                        format!("Illegal argument value: {message}")
                    } else {
                        message
                    },
                );
            }
            Thrown {
                failure: Failure::User,
                exception: error.java_class(),
                message: Some(error.message()),
            }
        };

    // `SAMSequenceDictionaryExtractor.extractDictionary`, whose four kinds this reads by name: a
    // `.dict` and a SAM header are the same text, a VCF's is its `##contig` lines, and a FASTA's is
    // the `.dict` beside it.
    let source = match argument(parser, "source-dictionary") {
        None => None,
        Some(path) => Some((path.clone(), dictionary_from_any(&path)?)),
    };
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    let dictionary = gatk_tools::update_vcf_sequence_dictionary::best_available_dictionary(
        source
            .as_ref()
            .map(|(name, sequences)| (name.as_str(), sequences.as_slice())),
        master.as_ref().map(|header| header.sequences.as_slice()),
        reference.as_ref().map(|header| header.sequences.as_slice()),
        !flag(parser, "disable-sequence-dictionary-validation"),
    )
    .map_err(refuse)?;

    // The input's OWN dictionary, read from its header rather than from the engine: the engine
    // would dig one out of an index, and the check is about what the file says.
    let own: Vec<htsjdk_bam::header::SequenceRecord> = sequence_dictionary_of(&file.header)
        .into_iter()
        .map(|(name, length)| htsjdk_bam::header::SequenceRecord::new(&name, length))
        .collect();
    // The traversal's own refusal comes FIRST: `-L` over an input with no random access is refused
    // while the data source is bounded, which is before `onTraversalStart` runs at all. Measured on
    // row 21 of this tool's array, where the port answered the dictionary check.
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    gatk_tools::update_vcf_sequence_dictionary::check_replace(&own, flag(parser, "replace"))
        .map_err(refuse)?;
    let records: Vec<htsjdk_vcf::variant::VariantContext> = kept.into_iter().cloned().collect();

    // `outputHeader.setSequenceDictionary(sourceDictionary)`: the contig lines are REPLACED, and
    // they carry the index they had in the dictionary rather than the one they had in the file.
    let mut header = file.header.clone();
    header
        .lines
        .retain(|line| !matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }));
    for line in default_tool_vcf_header_lines(parser, "UpdateVCFSequenceDictionary") {
        header.lines.push(line);
    }
    for (index, record) in dictionary.iter().enumerate() {
        header.lines.push(htsjdk_vcf::header::HeaderLine::contig(
            &record.name,
            i64::from(record.length),
            index as i32,
        ));
    }

    let (written, refusal) =
        gatk_tools::update_vcf_sequence_dictionary::update_dictionary(&dictionary, &records);
    let emitted: Vec<htsjdk_vcf::variant::VariantContext> = written
        .into_iter()
        .map(|index| records[index].clone())
        .collect();
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let mut emitted: Vec<htsjdk_vcf::variant::VariantContext> =
        emitted.into_iter().filter(|record| keep(record)).collect();
    apply_sites_only(parser, &mut header, &mut emitted);
    let out = htsjdk_vcf::vcf_file::write_vcf(&header, &emitted)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    write_variant_output(parser, &output, &out)?;
    match refusal {
        Some(error) => Err(refuse(error)),
        None => Ok(None),
    }
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`, over the file kinds the corpus can hold.
///
/// A `.dict` and a SAM header are the same lines; a BAM carries them in its header; a VCF declares
/// them as `##contig`; and a FASTA has none of its own, so the dictionary is the `.dict` beside it,
/// which is what `ReferenceSequenceFileFactory` names.
fn dictionary_from_any(path: &str) -> Result<Vec<htsjdk_bam::header::SequenceRecord>, Thrown> {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
        let text = std::fs::read_to_string(path)
            .map_err(|error| Thrown::user(format!("{path}: {error}")))?;
        return Ok(vcf_dictionary(&text).sequences);
    }
    if lower.ends_with(".fasta") || lower.ends_with(".fa") || lower.ends_with(".fna") {
        let beside = dictionary_path(path);
        let text = std::fs::read_to_string(&beside)
            .map_err(|error| Thrown::user(format!("{beside}: {error}")))?;
        return Ok(htsjdk_bam::reader::parse_header_text(&text).sequences);
    }
    if lower.ends_with(".bam") {
        let source =
            gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        return Ok(source.header().sequences.clone());
    }
    let text =
        std::fs::read_to_string(path).map_err(|error| Thrown::user(format!("{path}: {error}")))?;
    Ok(htsjdk_bam::reader::parse_header_text(&text).sequences)
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
    let VariantWalkerStart {
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
    // `UserException$BadInput` puts `Bad input: ` in front of its message, which the port's
    // `message()` leaves to the caller. Measured on six rows of this tool's array.
    .map_err(|refusal| {
        if refusal.java_class().ends_with("UserException$BadInput") {
            bad_input(refusal.message())
        } else {
            Thrown {
                failure: Failure::User,
                exception: refusal.java_class(),
                message: Some(refusal.message()),
            }
        }
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
            tool_command_line: command_line_header_line(parser, "SelectVariants"),
            samples: if sites_only {
                Vec::new()
            } else {
                selection.samples.clone()
            },
        },
    );
    // `VcfUtils.updateHeaderContigLines`, which every tool that writes a VCF beside a REFERENCE
    // calls and which rewrites two kinds of line: every `##contig` is replaced by one built from
    // the dictionary, carrying `assembly=` when a reference path is known, and the `##reference`
    // line is replaced by that path's URI. Without a reference the dictionary is the driving
    // VCF's own, and a file that declares none keeps the lines it had.
    let header = update_header_contig_lines(parser, header)?;

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

    // Whether any argument on this command line READS a genotype, which is what decides whether the
    // output carries the file's own genotype text or a rebuilt one. It is an access and not a
    // change: a run that decodes and rewrites nothing still writes the sorted form, because htsjdk
    // drops the unparsed block on the first accessor. Every clause here is a gate in the reference
    // that reaches `getGenotypes()`.
    let touches_genotypes = !selection.no_samples_specified
        || subset_arguments.remove_unused_alternates
        || subset_arguments.keep_original_chr_counts
        || subset_arguments.keep_original_depth
        || filter_arguments.exclude_non_variants
        || !filter_arguments.select_genotype_expressions.is_empty()
        || filter_arguments.max_filtered_genotypes != i32::MAX
        || filter_arguments.min_filtered_genotypes != 0
        || filter_arguments.max_fraction_filtered_genotypes != 1.0
        || filter_arguments.min_fraction_filtered_genotypes != 0.0
        || filter_arguments.max_nocall_number != i32::MAX
        || filter_arguments.max_nocall_fraction != 1.0
        || output_arguments.set_filtered_genotypes_to_no_call
        || !output_arguments.genotype_annotations_to_drop.is_empty();

    // The genotype text every record was READ with, taken before anything decodes it: the first
    // accessor drops it, and this port has to decode to decide what to keep.
    let unparsed: Vec<Option<String>> = file
        .records
        .iter()
        .map(|record| record.genotypes.unparsed().map(|text| text.to_string()))
        .collect();
    let sample_names: std::sync::Arc<[String]> = file.header.samples.clone().into();

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
        let mut vc = crate::variant_bridge::from_engine(original, &record);
        // And its genotypes come back from the FILE when nothing on this command line reads one.
        // htsjdk writes an untouched context from the text it was read with, so a whole-cohort run
        // keeps the input's `GT:GQ:DP` where a subsetting run recomputes and sorts to `GT:DP:GQ`.
        // The port decodes to decide, so it puts the text back rather than never taking it out.
        if !touches_genotypes {
            if let Some(text) = &unparsed[located.index] {
                vc.genotypes = htsjdk_vcf::genotypes_context::GenotypesContext::lazy(
                    vc.genotypes.to_vec(),
                    text.clone(),
                    sample_names.clone(),
                );
            }
        }
        pending.add(record, vc);
    }
    for (_, vc) in pending.drain() {
        written.push(vc);
    }

    // `--variant-output-filtering` wraps the writer here as it does on the other two tools that
    // write a VCF: the traversal reached the record and the writer decides whether it lands.
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    written.retain(|record| keep(record));

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

/// `VcfUtils.updateHeaderContigLines`: the contig lines and the `##reference` line, rebuilt.
///
/// The dictionary is the REFERENCE's when there is one and the driving variants' otherwise, and a
/// file that declares neither keeps whatever contig lines it had. `assembly=` is the reference
/// file's NAME, and the `##reference` value is its URI unless `--suppress-reference-path` asks for
/// the bare name with its extension cut. Measured on twenty-three rows of `SelectVariants`' array,
/// where the port wrote the input's contig lines and no reference line at all.
fn update_header_contig_lines(
    parser: &Parser,
    header: htsjdk_vcf::header::VcfHeader,
) -> Result<htsjdk_vcf::header::VcfHeader, Thrown> {
    let reference = argument(parser, "reference");
    let dictionary: Vec<(String, i32)> = match &reference {
        Some(path) => {
            let source =
                gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(path))
                    .map_err(|error| Thrown::user(format!("{error:?}")))?;
            gatk_tools::reference_walker::dictionary(&source)
                .sequences
                .iter()
                .map(|sequence| (sequence.name.clone(), sequence.length))
                .collect()
        }
        None => sequence_dictionary_of(&header),
    };
    if dictionary.is_empty() {
        return Ok(header);
    }

    let mut header = header;
    header.lines.retain(|line| {
        !matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }) && line.key() != "reference"
    });
    let assembly = reference.as_ref().map(|path| {
        std::path::Path::new(path)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.clone())
    });
    for (index, (name, length)) in dictionary.iter().enumerate() {
        let mut fields = vec![
            ("ID".to_string(), name.clone()),
            ("length".to_string(), length.to_string()),
        ];
        if let Some(assembly) = &assembly {
            fields.push(("assembly".to_string(), assembly.clone()));
        }
        header.lines.push(htsjdk_vcf::header::HeaderLine::Contig {
            index: index as i32,
            fields,
        });
    }
    if let Some(path) = &reference {
        // `referencePath.toUri()`, which is an ABSOLUTE URI: the runner is handed whatever the
        // command line wrote, so the path is resolved before it is rendered.
        let value = if flag(parser, "suppress-reference-path") {
            let name = std::path::Path::new(path)
                .file_name()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_else(|| path.clone());
            match name.rfind('.') {
                Some(dot) => name[..dot].to_string(),
                None => name,
            }
        } else {
            let absolute = std::path::Path::new(path)
                .canonicalize()
                .unwrap_or_else(|_| std::path::PathBuf::from(path));
            format!("file://{}", absolute.display())
        };
        header
            .lines
            .push(htsjdk_vcf::header::HeaderLine::Unstructured {
                key: "reference".to_string(),
                value,
            });
    }
    Ok(header)
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
    // Two of the three are not `UserException`s: an expression that does not compile is an
    // `IllegalArgumentException` and one that evaluates to the wrong class is the JVM's own
    // `ClassCastException`, so the handler prints the class and the run ends at three.
    if error.java_class().starts_with("java.lang.") {
        return Thrown::non_user(error.java_class(), error.message());
    }
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
    if refused.is_empty() {
        return Ok(());
    }
    Err(Thrown::non_user(
        PORT_LIMITATION,
        format!(
            "SelectVariants in this port does not implement {}: a pedigree's Mendelian \
             violations, the two random fractions, the concordance tracks and the genotype caller \
             each reach behaviour no measured static reproduces, and answering without them would \
             be a different answer rather than a refusal",
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
fn command_line_header_line(parser: &Parser, tool: &str) -> Option<htsjdk_vcf::header::HeaderLine> {
    if !flag(parser, "add-output-vcf-command-line") {
        return None;
    }
    Some(htsjdk_vcf::header::HeaderLine::Structured {
        key: "GATKCommandLine".to_string(),
        fields: vec![
            ("ID".to_string(), tool.to_string()),
            (
                "CommandLine".to_string(),
                crate::command_line::expanded(tool, parser),
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

    // `--create-output-variant-md5`, which digests the file as it was WRITTEN: a block compressed
    // output is digested compressed, because the digest is taken by the stream the writer wraps.
    if flag(parser, "create-output-variant-md5") {
        write_file(
            &format!("{output}.md5"),
            gatk_tools::gather_bam_files::md5_file(&bytes).as_bytes(),
        )?;
    }

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

/// `PathSeqBuildReferenceTaxonomy.doWork`, with every file read in the order the tool opens it.
///
/// The port's arithmetic is [`gatk_tools::pathseq_taxonomy`] and its file is
/// [`gatk_tools::pathseq_kryo::taxonomy_database_file`], both oracle-backed. What is here is the
/// order the refusals come in, which is the order the tool touches its inputs: the catalog check
/// before anything is opened, then the reference's dictionary, then the reference's names (a taxon
/// id that is not a number is refused there), the RefSeq catalog, the GenBank one, `names.dmp`,
/// `nodes.dmp`, the tree, and last the output. So each file is read only when the tool would have
/// opened it, and a refusal from an earlier one is never masked by a missing later one.
///
/// A catalog is gunzipped only when its PATH ends in `.gz`, which is `makeReaderMaybeGzipped`'s
/// rule: the bytes are not sniffed. A catalog that says `.gz` and is not gzip fails to open,
/// and that refusal is the same one a missing file gets.
pub fn path_seq_build_reference_taxonomy(parser: &Parser) -> Outcome {
    use gatk_tools::pathseq_taxonomy::{self as taxonomy, CatalogFormat, TaxonomyError};

    let refused = |error: TaxonomyError| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(error.message()),
    };
    let refseq = argument(parser, "refseq-catalog");
    let genbank = argument(parser, "genbank-catalog");
    if refseq.is_none() && genbank.is_none() {
        return Err(refused(TaxonomyError::NoCatalog));
    }
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let tax_dump = argument(parser, "tax-dump").ok_or_else(|| {
        Thrown::command_line("Argument tax-dump was missing: Argument 'tax-dump' is required")
    })?;
    let min_length: i64 = scalar(parser, "min-non-virus-contig-length")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0);

    let dictionary = reference_dictionary(parser)?.ok_or_else(|| {
        bad_input(
            "Reference sequence dictionary not found. Please build one using \
             CreateSequenceDictionary."
                .to_string(),
        )
    })?;
    let contigs: Vec<(String, i64)> = dictionary
        .sequences
        .iter()
        .map(|record| (record.name.clone(), record.length as i64))
        .collect();

    let mut properties = gatk_engine::java_hash::JavaHashMap::new();
    let by_accession =
        taxonomy::parse_reference_records(&contigs, &mut properties).map_err(refused)?;
    let mut not_found = None;
    if let Some(path) = &refseq {
        let text = catalog_text(path)?;
        not_found = Some(
            taxonomy::parse_catalog(
                &text,
                CatalogFormat::RefSeq,
                &by_accession,
                &mut properties,
                None,
            )
            .map_err(refused)?,
        );
    }
    if let Some(path) = &genbank {
        let text = catalog_text(path)?;
        taxonomy::parse_catalog(
            &text,
            CatalogFormat::GenBank,
            &by_accession,
            &mut properties,
            not_found.as_ref(),
        )
        .map_err(refused)?;
    }
    let names = tar_gz_entry(&tax_dump, "names.dmp")?;
    taxonomy::parse_names(&names, &mut properties).map_err(refused)?;
    let nodes = tar_gz_entry(&tax_dump, "nodes.dmp")?;
    taxonomy::parse_nodes(&nodes, &mut properties).map_err(refused)?;

    let tree = taxonomy::build_taxonomic_tree(&properties).map_err(refused)?;
    taxonomy::remove_unused_tax_ids(&mut properties, &tree);
    let map = taxonomy::build_accession_to_tax_id(&properties, &tree, min_length);
    let file = gatk_tools::pathseq_kryo::taxonomy_database_file(&tree, &map).map_err(|error| {
        Thrown::non_user(
            PORT_LIMITATION,
            format!("a HashMap order this port has not measured: {error:?}"),
        )
    })?;
    std::fs::write(&output, file).map_err(|_| Thrown {
        failure: Failure::User,
        exception:
            "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile",
        message: Some("Could not serialize objects to file".to_string()),
    })?;
    Ok(None)
}

/// `getBufferedReaderGz`: the file, gunzipped when its name ends in `.gz`, as text.
fn catalog_text(path: &str) -> Result<String, Thrown> {
    let cannot_open = || bad_input(format!("Could not open file {path}"));
    let bytes = std::fs::read(path).map_err(|_| cannot_open())?;
    if !path.ends_with(".gz") {
        return Ok(String::from_utf8_lossy(&bytes).into_owned());
    }
    // `GZIPInputStream` reads its header when it is constructed, inside the same `try`, so a file
    // that is not gzip is refused as one that could not be opened.
    if bytes.len() < 2 || bytes[0] != 0x1f || bytes[1] != 0x8b {
        return Err(cannot_open());
    }
    gunzip(&bytes).map_err(|_| Thrown::user("Error reading from catalog file".to_string()))
}

fn gunzip(bytes: &[u8]) -> std::io::Result<String> {
    use std::io::Read;
    let mut text = Vec::new();
    flate2::read::MultiGzDecoder::new(bytes).read_to_end(&mut text)?;
    Ok(String::from_utf8_lossy(&text).into_owned())
}

/// `getBufferedReaderTarGz`: one entry of a gzipped tarball, found by walking its headers.
fn tar_gz_entry(tar_path: &str, name: &str) -> Result<String, Thrown> {
    let cannot_open = || {
        bad_input(format!(
            "Could not open compressed tarball file {name} in {tar_path}"
        ))
    };
    let bytes = std::fs::read(tar_path).map_err(|_| cannot_open())?;
    if bytes.len() < 2 || bytes[0] != 0x1f || bytes[1] != 0x8b {
        return Err(cannot_open());
    }
    let tar = {
        use std::io::Read;
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(bytes.as_slice())
            .read_to_end(&mut out)
            .map_err(|_| cannot_open())?;
        out
    };
    tar_entry(&tar, name)
        .map(|entry| String::from_utf8_lossy(entry).into_owned())
        .ok_or_else(|| bad_input(format!("Could not find file {name} in tarball {tar_path}")))
}

/// The data of the first entry named `name` in an uncompressed tar stream.
///
/// A ustar header is 512 bytes: the name in the first hundred, the size in octal at 124, the type
/// at 156, and a ustar prefix at 345 that the name is joined to. A GNU `L` entry carries the NEXT
/// entry's long name as its data, and a pax `x` entry may carry it as `path=`. An all-zero block
/// ends the archive.
fn tar_entry<'a>(tar: &'a [u8], name: &str) -> Option<&'a [u8]> {
    let field = |block: &[u8], from: usize, to: usize| -> String {
        let raw = &block[from..to];
        let end = raw.iter().position(|byte| *byte == 0).unwrap_or(raw.len());
        String::from_utf8_lossy(&raw[..end]).into_owned()
    };
    let mut offset = 0;
    let mut long_name: Option<String> = None;
    while offset + 512 <= tar.len() {
        let block = &tar[offset..offset + 512];
        if block.iter().all(|byte| *byte == 0) {
            return None;
        }
        let size = if block[124] & 0x80 != 0 {
            block[125..136]
                .iter()
                .fold(0usize, |acc, byte| (acc << 8) | *byte as usize)
        } else {
            usize::from_str_radix(field(block, 124, 136).trim(), 8).unwrap_or(0)
        };
        let data_start = offset + 512;
        let data_end = (data_start + size).min(tar.len());
        let data = &tar[data_start..data_end];
        let next = data_start + size.div_ceil(512) * 512;
        match block[156] {
            b'L' => {
                long_name = Some(field(data, 0, data.len()));
            }
            b'x' => {
                let text = String::from_utf8_lossy(data);
                for record in text.lines() {
                    if let Some((_, rest)) = record.split_once(' ') {
                        if let Some(path) = rest.strip_prefix("path=") {
                            long_name = Some(path.to_string());
                        }
                    }
                }
            }
            b'g' => {}
            _ => {
                let entry_name = long_name.take().unwrap_or_else(|| {
                    let short = field(block, 0, 100);
                    let prefix = if &block[257..262] == b"ustar" {
                        field(block, 345, 500)
                    } else {
                        String::new()
                    };
                    if prefix.is_empty() {
                        short
                    } else {
                        format!("{prefix}/{short}")
                    }
                });
                if entry_name == name {
                    return Some(data);
                }
            }
        }
        offset = next;
    }
    None
}

/// `IOUtils.writeTarGz`: regular files in a gzipped ustar stream, as commons-compress lays them
/// out. A name of a hundred bytes or more goes ahead of its entry in a GNU `././@LongLink` entry,
/// which is `LONGFILE_GNU`; the stream ends with two zero blocks.
fn tar_gz(entries: &[(String, Vec<u8>)]) -> Vec<u8> {
    use std::io::Write;
    let mtime = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_secs());
    let header = |name: &[u8], size: usize, kind: u8| -> [u8; 512] {
        let mut block = [0u8; 512];
        let name_length = name.len().min(100);
        block[..name_length].copy_from_slice(&name[..name_length]);
        block[100..108].copy_from_slice(b"0000644\0");
        block[108..116].copy_from_slice(b"0000000\0");
        block[116..124].copy_from_slice(b"0000000\0");
        block[124..136].copy_from_slice(format!("{size:011o}\0").as_bytes());
        block[136..148].copy_from_slice(format!("{mtime:011o}\0").as_bytes());
        block[156] = kind;
        block[257..263].copy_from_slice(b"ustar\0");
        block[263..265].copy_from_slice(b"00");
        block[148..156].copy_from_slice(b"        ");
        let sum: u32 = block.iter().map(|byte| *byte as u32).sum();
        block[148..156].copy_from_slice(format!("{sum:06o}\0 ").as_bytes());
        block
    };
    let mut tar: Vec<u8> = Vec::new();
    let push = |tar: &mut Vec<u8>, data: &[u8]| {
        tar.extend_from_slice(data);
        tar.resize(tar.len().div_ceil(512) * 512, 0);
    };
    for (name, data) in entries {
        if name.len() >= 100 {
            let mut long = name.as_bytes().to_vec();
            long.push(0);
            tar.extend_from_slice(&header(b"././@LongLink", long.len(), b'L'));
            push(&mut tar, &long);
        }
        tar.extend_from_slice(&header(name.as_bytes(), data.len(), b'0'));
        push(&mut tar, data);
    }
    tar.resize(tar.len() + 1024, 0);
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder
        .write_all(&tar)
        .expect("writing to a vector does not fail");
    encoder.finish().expect("writing to a vector does not fail")
}

/// `PathSeqBuildKmers.doWork`, up to the file it writes.
///
/// The arithmetic is [`gatk_tools::pathseq_kmers`] and the file is
/// [`gatk_tools::pathseq_kryo::kmer_set_file`], both oracle-backed. What the runner adds is what
/// decides the TABLE the file carries, and it is not the set of k-mers:
///
/// * **the contigs are read in the file's own order**, because `getAllReferenceBases` collects them
///   into a `LinkedHashMap` as `nextSequence` hands them over, and the bases are the file's bytes,
///   neither upper-cased nor flattened;
/// * **every k-mer is ADDED, duplicates included**, one contig's array after another, so the
///   insertion order is the reference's own and a repeated k-mer is added again rather than
///   skipped;
/// * **the set is sized by that total**, `numLongs` counting duplicates, which is the capacity the
///   file's first int reports.
///
/// **The reference is tested before the mask and read after it**: the constructor of
/// `ReferenceFileSparkSource` refuses a path that does not exist, then `parseMask` runs, and only
/// then are the bases loaded. A bad mask over a missing reference is therefore the reference's
/// refusal, and a bad mask over a present one never waits for the file to be read.
///
/// The output name gains `.hss` when it does not already end in it, which is `writeKmerSet`'s own
/// rule and not the caller's.
///
/// `--bloom-false-positive-probability` above zero writes a `PSKmerBloomFilter` instead, which is a
/// container this port does not carry: that is a refusal of the port's own, not one of GATK's.
pub fn path_seq_build_kmers(parser: &Parser) -> Outcome {
    use gatk_engine::hopscotch::LargeLongHopscotchSet;
    use gatk_tools::pathseq_kmers::{self, KmerError};

    let refused = |error: KmerError| Thrown::non_user(error.java_class(), error.message());
    let reference = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let kmer_size: usize = scalar(parser, "kmer-size")
        .and_then(|value| value.parse().ok())
        .unwrap_or(31);
    let spacing: usize = scalar(parser, "kmer-spacing")
        .and_then(|value| value.parse().ok())
        .unwrap_or(1);
    let bloom: f64 = scalar(parser, "bloom-false-positive-probability")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.0);
    let mask_argument = scalar(parser, "kmer-mask").unwrap_or_default();

    // `ReferenceFileSparkSource`'s constructor, which runs before the mask is parsed and tests
    // only that the path exists.
    if !std::path::Path::new(&reference).exists() {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$MissingReference",
            message: Some(format!(
                "The specified fasta file ({reference}) does not exist."
            )),
        });
    }
    let positions = pathseq_kmers::parse_mask(&mask_argument, kmer_size).map_err(refused)?;
    let mask = pathseq_kmers::get_mask(&positions, kmer_size);

    let bases = std::fs::read(&reference).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not read {reference}: {error}"))
    })?;
    if gatk_tools::read_walker_refusal::is_block_compressed(&bases) {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "PathSeqBuildKmers reads a block-compressed FASTA through its .gzi, which is not ported",
        ));
    }
    let contigs = fasta_records(&bases);
    let mut per_contig = Vec::new();
    let mut total: i64 = 0;
    for contig in &contigs {
        let kmers =
            pathseq_kmers::masked_kmers(contig, kmer_size, spacing, mask).map_err(refused)?;
        total += kmers.len() as i64;
        per_contig.push(kmers);
    }
    if bloom > 0.0 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "PathSeqBuildKmers writes a PSKmerBloomFilter for --bloom-false-positive-probability \
             above zero, and the Bloom filter is not ported",
        ));
    }
    if total == 0 {
        // `LargeLongHopscotchSet` refuses to be built from nothing, and it is the container that
        // refuses rather than the tool.
        return Err(refused(KmerError::EmptySet));
    }
    let mut set = LargeLongHopscotchSet::new(total);
    for kmers in &per_contig {
        for value in kmers {
            set.add(*value as i64);
        }
    }
    let file = gatk_tools::pathseq_kryo::kmer_set_file(kmer_size as i32, mask.0 as i64, &set);
    let name = if output.to_lowercase().ends_with(".hss") {
        output.clone()
    } else {
        format!("{output}.hss")
    };
    std::fs::write(&name, file).map_err(|error| {
        Thrown::non_user(PORT_FAILURE, format!("could not write {name}: {error}"))
    })?;
    Ok(None)
}

/// The bases of each record of a FASTA, in the file's order and as its own bytes.
///
/// `nextSequence` hands back what is written: a lower-case base stays lower-case and an IUPAC code
/// stays itself, which is what `SVKmerizer` then refuses a window over.
fn fasta_records(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut records: Vec<Vec<u8>> = Vec::new();
    for line in bytes.split(|byte| *byte == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.first() == Some(&b'>') {
            records.push(Vec::new());
        } else if let Some(current) = records.last_mut() {
            current.extend_from_slice(line);
        }
    }
    records
}

/// `CondenseDepthEvidence`: adjacent depth-evidence bins merged.
///
/// The merge is [`gatk_tools::condense_depth_evidence`], where a golden measures it. What the
/// runner adds is the order the refusals come in, which is the engine's and then the tool's:
///
/// * **the input is opened at startup**, by `FeatureManager.getCodecForFile` once the master
///   dictionary is loaded, so an unreadable file is refused before any argument of the tool's own
///   is looked at;
/// * **then `onTraversalStart`**: a minimum above the maximum, then no codec for the output's
///   name, then a codec for another feature type;
/// * **then the sink is opened and its header written**, from the input's sample names, before a
///   single record is merged.
///
/// Three things are refusals of the port's own rather than GATK's. A block-compressed or binary
/// OUTPUT is written through the GKL deflater and a tabix index or through the `.bci` container,
/// none of which this runner carries. A block-compressed or binary INPUT is the same gap on the
/// read side. And a malformed input is refused by Tribble wrapping the codec's exception, whose
/// message names an iterator by its identity hash, so no byte of it can be reproduced.
pub fn condense_depth_evidence(parser: &Parser) -> Outcome {
    use gatk_tools::condense_depth_evidence as condense;
    use gatk_tools::sv_feature_codecs::{self as codecs, Encoding};

    let input = argument(parser, "depth-evidence").ok_or_else(|| {
        Thrown::command_line(
            "Argument depth-evidence was missing: Argument 'depth-evidence' is required",
        )
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let arguments = condense::Arguments {
        max_interval_length: number_or(parser, "max-interval-size", 1000),
        min_interval_length: number_or(parser, "min-interval-size", 0),
    };
    let refused = |error: condense::CondenseError| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException",
        message: Some(error.message()),
    };

    // The engine's own prelude, the one `PrintReadCounts` measured for the same base class: the
    // read filters resolve while the command line is parsed, and `onStartup` loads a master
    // dictionary before it opens anything else.
    let _ = resolve_read_filters(parser, "CondenseDepthEvidence")?;
    let _ = master_dictionary(parser)?;
    // `getCodecForFile` tests that the path is readable before it asks a codec anything.
    let bytes = std::fs::read(&input).map_err(|_| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: java_absolute_path(&input),
            }
            .message(),
        )
    })?;
    match codecs::find(&input) {
        Some(codec) if codec.feature_type != "DepthEvidence" => {
            return Err(Thrown::user(format!(
                "File {input} contains features of the wrong type."
            )));
        }
        Some(codecs::Codec {
            encoding: Encoding::Text {
                block_compressed: false,
            },
            ..
        }) => {}
        Some(_) => {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                "CondenseDepthEvidence reads a block-compressed or .bci depth-evidence file, \
                 which is not ported",
            ));
        }
        // Every other codec the engine knows would be asked here, and is not ported: the SV
        // evidence codecs are the ones this tool can accept.
        None => {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "no SV evidence codec reads {input}, and the engine's others are not ported"
                ),
            ));
        }
    }
    let text = String::from_utf8_lossy(&bytes);
    let (samples, records) = condense::read(&text).map_err(|problem| {
        Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "{input} is a malformed depth-evidence file ({problem}), which Tribble refuses"
            ),
        )
    })?;

    condense::check_lengths(&arguments).map_err(refused)?;
    condense::check_output(&output).map_err(refused)?;
    if codecs::find(&output).map(|codec| codec.encoding)
        != Some(Encoding::Text {
            block_compressed: false,
        })
    {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "CondenseDepthEvidence writes a block-compressed or .bci depth-evidence file, which is \
             not ported",
        ));
    }
    // A mismatch is raised inside `apply`, after the sink has written its header, and it is an
    // `IllegalArgumentException` rather than a user error. What the half-written file then holds
    // is the buffered writer's, which nothing has measured, so only the header is written here.
    let merged = match condense::condense(&records, &arguments) {
        Ok(merged) => merged,
        Err(error) => {
            write_file(&output, condense::write(&samples, &[]).as_bytes())?;
            return Err(Thrown::non_user(error.java_class(), error.message()));
        }
    };
    write_file(&output, condense::write(&samples, &merged).as_bytes())?;
    Ok(None)
}

/// `PrintSVEvidence`, for depth evidence: several files merged into one, rewritten against one
/// sample list.
///
/// The sample list and the sort merger are [`gatk_tools::print_sv_evidence`] and the walk's order
/// is [`gatk_engine::multi_feature_walker`], both oracle-backed. What the runner adds is where each
/// refusal comes from, which is `onStartup` before `onTraversalStart`:
///
/// * **every input is opened at startup**, by `getCodecForFile`, in the order the command line
///   names them, after the read filters and the master dictionary;
/// * **then the dictionary is chosen**, the master one against the reference's by
///   `betterDictionary`. A depth-evidence header carries none, so a run given neither is refused
///   with `No dictionary found`;
/// * **then `onTraversalStart`**: no codec for the output's name, then an input of another type;
/// * **then the walk**, whose refusals are mid-run: an input going backwards, a sample two files
///   both report at one bin, and a contig the dictionary does not name, which the sort merger's
///   `compareLocatables` refuses once there are two records to compare.
///
/// `--sample-names` is a `LinkedHashSet`, so a name given twice is one column, in the order it was
/// first given.
///
/// Refusals of the port's own: any evidence type but depth, a block-compressed or `.bci` file on
/// either side, and a malformed input. What a run that fails mid-walk leaves in its output is the
/// buffered writer's, which nothing has measured: the runner writes the header the sink opened
/// with and nothing after it.
pub fn print_sv_evidence(parser: &Parser) -> Outcome {
    use gatk_engine::multi_feature_walker as walker;
    use gatk_tools::condense_depth_evidence as depth;
    use gatk_tools::print_sv_evidence as print;
    use gatk_tools::sv_feature_codecs::{self as codecs, Encoding};

    let inputs = distinct_feature_inputs(arguments(parser, "evidence-file"));
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let mut requested: Vec<String> = Vec::new();
    for name in arguments(parser, "sample-names") {
        if !requested.contains(&name) {
            requested.push(name);
        }
    }
    let user = |exception: &'static str, message: String| Thrown {
        failure: Failure::User,
        exception,
        message: Some(message),
    };
    let limitation = |message: String| Thrown::non_user(PORT_LIMITATION, message);

    let _ = resolve_read_filters(parser, "PrintSVEvidence")?;
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;

    let mut files = Vec::new();
    for input in &inputs {
        let bytes = std::fs::read(input).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: java_absolute_path(input),
                }
                .message(),
            )
        })?;
        match codecs::find(input) {
            Some(codecs::Codec {
                feature_type: "DepthEvidence",
                encoding:
                    Encoding::Text {
                        block_compressed: false,
                    },
            }) => {}
            _ => {
                return Err(limitation(format!(
                    "PrintSVEvidence reads {input}, and only plain-text depth evidence is ported"
                )))
            }
        }
        let (samples, records) =
            depth::read(&String::from_utf8_lossy(&bytes)).map_err(|problem| {
                limitation(format!(
                    "{input} is a malformed depth-evidence file ({problem}), which Tribble refuses"
                ))
            })?;
        if records
            .iter()
            .any(|record| record.counts.len() != samples.len())
        {
            return Err(limitation(format!(
                "{input} has records whose counts do not match its header, which the walk refuses \
                 with an index out of bounds"
            )));
        }
        files.push(print::EvidenceFile {
            samples,
            records: records
                .into_iter()
                .map(|record| print::DepthEvidence {
                    contig: record.contig,
                    start: record.start,
                    end: record.end,
                    counts: record.counts,
                })
                .collect(),
        });
    }

    let source = |header: Option<SamHeader>, name: &str| {
        header.map(|header| walker::DictSource {
            contigs: header
                .sequences
                .iter()
                .map(|record| record.name.clone())
                .collect(),
            source: name.to_string(),
        })
    };
    let dictionary = walker::choose_dictionary(
        source(master, "sequence-dictionary"),
        source(reference, "reference"),
    )
    .map_err(|error| {
        user(
            "org.broadinstitute.hellbender.exceptions.UserException",
            error.message(),
        )
    })?;

    print::check_types(&output, &inputs).map_err(|error| {
        user(
            "org.broadinstitute.hellbender.exceptions.UserException",
            error.message(),
        )
    })?;
    if codecs::find(&output).map(|codec| codec.encoding)
        != Some(Encoding::Text {
            block_compressed: false,
        })
    {
        return Err(limitation(
            "PrintSVEvidence writes a block-compressed or .bci file, which is not ported"
                .to_string(),
        ));
    }
    let samples = print::sample_names(&requested, &files);

    // The walk, each record carrying where it came from in `text`.
    let located: Vec<Vec<walker::Located>> = files
        .iter()
        .map(|file| {
            file.records
                .iter()
                .enumerate()
                .map(|(index, record)| walker::Located {
                    contig: record.contig.clone(),
                    start: record.start,
                    end: record.end,
                    text: index.to_string(),
                })
                .collect()
        })
        .collect();
    let header_only = || write_file(&output, print::write(&samples, &[]).as_bytes());
    let walked = match walker::merge_with_sources(&located, &dictionary) {
        Ok(walked) => walked,
        Err(error) => {
            header_only()?;
            return Err(user(
                "org.broadinstitute.hellbender.exceptions.UserException",
                error.message(),
            ));
        }
    };
    if walked.len() > 1
        && walked
            .iter()
            .any(|(_, feature)| dictionary.sequence_index(&feature.contig) == -1)
    {
        header_only()?;
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Can't do comparison because Locatables' contigs not found in sequence dictionary",
        ));
    }
    let merged: Vec<(usize, print::DepthEvidence)> = walked
        .iter()
        .map(|(input, feature)| {
            let index: usize = feature.text.parse().expect("a record index");
            (*input, files[*input].records[index].clone())
        })
        .collect();
    match print::run(&files, &merged, &requested) {
        Ok((samples, written)) => {
            write_file(&output, print::write(&samples, &written).as_bytes())?;
            Ok(None)
        }
        Err(error) => {
            header_only()?;
            Err(user(
                "org.broadinstitute.hellbender.exceptions.UserException",
                error.message(),
            ))
        }
    }
}

/// `SiteDepthtoBAF`: per-sample allele depths at a set of sites turned into B-allele fractions.
///
/// The arithmetic, the site iterator and the value's `DecimalFormat` are
/// [`gatk_tools::site_depth_to_baf`], and the walk's order is
/// [`gatk_engine::multi_feature_walker`], all oracle-backed. The runner adds the order the tool's
/// files are opened in:
///
/// * **the depth files at startup**, then the dictionary, as for every `MultiFeatureWalker`. A
///   `.sd.txt` has no header at all, so the dictionary is the master one or the reference's;
/// * **then `onTraversalStart`**: the sites VCF is opened, its dictionary must be the walk's
///   (`assertSameDictionary`), and only then is the output's name asked for a codec, refused when
///   there is none or when it names another feature type;
/// * **the sink writes no header**: a `.baf.txt` is records and nothing else.
///
/// Refusals of the port's own: a sites VCF whose dictionary differs from the walk's, since the
/// reference raises an `AssertionError` whose message prints `SAMSequenceRecord.toString`; a sites
/// VCF with no contig lines, whose missing dictionary the reference dereferences; a depth or BAF
/// file that is block compressed or `.bci`; and a malformed depth file. A run refused mid-walk
/// leaves the file the sink created, empty, since what the buffered writer had flushed is
/// unmeasured.
pub fn site_depth_to_baf(parser: &Parser) -> Outcome {
    use gatk_engine::multi_feature_walker as walker;
    use gatk_tools::site_depth_to_baf as baf;
    use gatk_tools::sv_feature_codecs::{self as codecs, Encoding};
    use std::io::Read;

    let inputs = distinct_feature_inputs(arguments(parser, "site-depth"));
    let sites_path = argument(parser, "baf-sites-vcf").ok_or_else(|| {
        Thrown::command_line(
            "Argument baf-sites-vcf was missing: Argument 'baf-sites-vcf' is required",
        )
    })?;
    let output = argument(parser, "baf-evidence-output").ok_or_else(|| {
        Thrown::command_line(
            "Argument baf-evidence-output was missing: Argument 'baf-evidence-output' is required",
        )
    })?;
    let double = |name: &str, default: f64| {
        scalar(parser, name)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(default)
    };
    let arguments = baf::Arguments {
        max_std_dev: double("max-std", 0.2),
        min_total_depth: number_or(parser, "min-total-depth", 10),
        min_het_probability: double("min-het-probability", 0.5),
    };
    let user = |message: String| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException",
        message: Some(message),
    };
    let limitation = |message: String| Thrown::non_user(PORT_LIMITATION, message);
    let unreadable = |path: &str| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: java_absolute_path(path),
            }
            .message(),
        )
    };

    let _ = resolve_read_filters(parser, "SiteDepthtoBAF")?;
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;

    // `SiteDepthCodec.decode` over every line, a header included: it has none to skip.
    let unsigned = |field: &str| field.parse::<u32>().map(|value| value as i32).ok();
    let mut files: Vec<Vec<baf::SiteDepth>> = Vec::new();
    for input in &inputs {
        let bytes = std::fs::read(input).map_err(|_| unreadable(input))?;
        match codecs::find(input) {
            Some(codecs::Codec {
                feature_type: "SiteDepth",
                encoding:
                    Encoding::Text {
                        block_compressed: false,
                    },
            }) => {}
            _ => {
                return Err(limitation(format!(
                    "SiteDepthtoBAF reads {input}, and only plain-text site depth is ported"
                )))
            }
        }
        let mut records = Vec::new();
        for line in String::from_utf8_lossy(&bytes).lines() {
            let columns: Vec<&str> = line.split('\t').collect();
            let parsed = (columns.len() == 7)
                .then(|| {
                    Some(baf::SiteDepth {
                        contig: columns[0].to_string(),
                        position: unsigned(columns[1])?.wrapping_add(1),
                        sample: columns[2].to_string(),
                        counts: [
                            unsigned(columns[3])?,
                            unsigned(columns[4])?,
                            unsigned(columns[5])?,
                            unsigned(columns[6])?,
                        ],
                    })
                })
                .flatten();
            match parsed {
                Some(record) => records.push(record),
                None => {
                    return Err(limitation(format!(
                        "{input} is a malformed site-depth file, which Tribble refuses"
                    )))
                }
            }
        }
        files.push(records);
    }

    let lengths = |header: &Option<SamHeader>| {
        header.as_ref().map(|header| {
            header
                .sequences
                .iter()
                .map(|record| (record.name.clone(), record.length as i64))
                .collect::<Vec<_>>()
        })
    };
    let master_lengths = lengths(&master);
    let reference_lengths = lengths(&reference);
    let source = |pairs: &Option<Vec<(String, i64)>>, name: &str| {
        pairs.as_ref().map(|pairs| walker::DictSource {
            contigs: pairs.iter().map(|(contig, _)| contig.clone()).collect(),
            source: name.to_string(),
        })
    };
    let dictionary = walker::choose_dictionary(
        source(&master_lengths, "sequence-dictionary"),
        source(&reference_lengths, "reference"),
    )
    .map_err(|error| user(error.message()))?;
    // The chosen one is whichever `betterDictionary` kept, and it is the larger.
    let chosen: Vec<(String, i64)> = [master_lengths, reference_lengths]
        .into_iter()
        .flatten()
        .find(|pairs| {
            pairs.iter().map(|(contig, _)| contig).collect::<Vec<_>>()
                == dictionary.contigs.iter().collect::<Vec<_>>()
        })
        .unwrap_or_default();

    // `onTraversalStart`: the sites VCF, then its dictionary against the walk's.
    let raw = std::fs::read(&sites_path).map_err(|_| unreadable(&sites_path))?;
    let text = if codecs::has_block_compressed_extension(&sites_path) {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(raw.as_slice())
            .read_to_string(&mut text)
            .map_err(|_| limitation(format!("{sites_path} could not be inflated")))?;
        text
    } else {
        String::from_utf8_lossy(&raw).into_owned()
    };
    let contigs: Vec<(String, i64)> = text
        .lines()
        .take_while(|line| line.starts_with("##"))
        .filter_map(|line| line.strip_prefix("##contig=<"))
        .filter_map(|fields| {
            let fields = fields.trim_end_matches('>');
            let mut id = None;
            let mut length = 0i64;
            for field in fields.split(',') {
                match field.split_once('=') {
                    Some(("ID", value)) => id = Some(value.to_string()),
                    Some(("length", value)) => length = value.parse().unwrap_or(0),
                    _ => {}
                }
            }
            id.map(|id| (id, length))
        })
        .collect();
    if contigs.is_empty() {
        return Err(limitation(format!(
            "{sites_path} declares no contigs, and the reference dereferences its null dictionary"
        )));
    }
    if contigs != chosen {
        return Err(limitation(format!(
            "{sites_path}'s dictionary differs from the walk's, which the reference refuses with \
             an AssertionError the port does not reproduce"
        )));
    }
    match codecs::find(&output) {
        None => return Err(user(codecs::no_output_codec(&output))),
        Some(codec) if codec.feature_type != "BafEvidence" => {
            return Err(user(format!(
                "We're intending to write BafEvidence, but the feature type associated with the \
                 output file expects features of type {}",
                codec.feature_type
            )))
        }
        Some(codecs::Codec {
            encoding: Encoding::Text {
                block_compressed: false,
            },
            ..
        }) => {}
        Some(_) => {
            return Err(limitation(
                "SiteDepthtoBAF writes a block-compressed or .bci file, which is not ported"
                    .to_string(),
            ))
        }
    }

    let sites = baf::read_sites(&text);
    let located: Vec<Vec<walker::Located>> = files
        .iter()
        .map(|records| {
            records
                .iter()
                .enumerate()
                .map(|(index, record)| walker::Located {
                    contig: record.contig.clone(),
                    start: record.position,
                    end: record.position,
                    text: index.to_string(),
                })
                .collect()
        })
        .collect();
    let empty = || write_file(&output, b"");
    let walked = match walker::merge_with_sources(&located, &dictionary) {
        Ok(walked) => walked,
        Err(error) => {
            empty()?;
            return Err(user(error.message()));
        }
    };
    let depths: Vec<baf::SiteDepth> = walked
        .iter()
        .map(|(input, feature)| {
            files[*input][feature.text.parse::<usize>().expect("a record index")].clone()
        })
        .collect();
    match baf::run(&depths, &sites, &arguments) {
        Ok(written) => {
            write_file(&output, baf::write(&written).as_bytes())?;
            Ok(None)
        }
        Err(error) => {
            empty()?;
            Err(user(error.message()))
        }
    }
}

/// The refusal each Mutect table reader raises for a file it cannot open: the `IOException` is
/// caught and replaced, and the message names the file as the command line gave it.
fn table_unreadable(path: &str) -> Thrown {
    Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException",
        message: Some(format!(
            "Encountered an IO exception while reading from {path}."
        )),
    }
}

/// A gather's error, whose class is borrowed from the error: a `Thrown` holds a `&'static str`.
fn gather_error(error: gatk_tools::mutect_gathers::GatherError) -> Thrown {
    let class = error.java_class();
    Thrown {
        failure: if class.starts_with("org.broadinstitute.hellbender.exceptions.UserException") {
            Failure::User
        } else {
            Failure::Other
        },
        exception: class,
        message: Some(error.message()),
    }
}

/// `MergeMutectStats`: every shard's statistics summed.
///
/// The sum is [`gatk_tools::mutect_gathers::merge_stats`], where a golden measures it. The runner
/// adds that `--stats` is a `LinkedHashSet<File>`, so a file named twice is read once, and that
/// the files are all read before any statistic is checked: an unreadable shard is refused before
/// an unknown statistic in an earlier one.
pub fn merge_mutect_stats(parser: &Parser) -> Outcome {
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let mut paths: Vec<String> = Vec::new();
    for path in arguments(parser, "stats") {
        if !paths.contains(&path) {
            paths.push(path);
        }
    }
    let mut shards = Vec::new();
    for path in &paths {
        shards.push(std::fs::read_to_string(path).map_err(|_| table_unreadable(path))?);
    }
    let texts: Vec<&str> = shards.iter().map(String::as_str).collect();
    let merged = gatk_tools::mutect_gathers::merge_stats(&texts).map_err(gather_error)?;
    write_file(&output, merged.as_bytes())?;
    Ok(Some("SUCCESS".to_string()))
}

/// `GatherPileupSummaries`: the non-empty shards sorted by their first record, then concatenated.
///
/// The gather is [`gatk_tools::mutect_gathers::gather_pileup_summaries`]. The runner adds what
/// `onStartup` does first: `--sequence-dictionary` is loaded before any shard is read, and one
/// without an `@SQ` line is refused as a malformed file named by its ABSOLUTE path, since
/// `loadFastaDictionary(File)` rewraps it that way. The value `doWork` returns counts the shards
/// kept, which is what the tool prints.
pub fn gather_pileup_summaries(parser: &Parser) -> Outcome {
    use gatk_engine::pileup_summary;

    let dictionary_path = argument(parser, "sequence-dictionary").ok_or_else(|| {
        Thrown::command_line(
            "Argument sequence-dictionary was missing: Argument 'sequence-dictionary' is required",
        )
    })?;
    let output = argument(parser, "O")
        .ok_or_else(|| Thrown::command_line("Argument O was missing: Argument 'O' is required"))?;
    let inputs = arguments(parser, "I");
    let text = std::fs::read_to_string(&dictionary_path).map_err(|_| {
        Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "{dictionary_path} could not be read, and the reference's refusal carries the \
                 message of an IOException the port does not reproduce"
            ),
        )
    })?;
    let header = htsjdk_bam::reader::parse_header_text(&text);
    if header.sequences.is_empty() {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile",
            message: Some(format!(
                "Unknown file is malformed: Could not read sequence dictionary from given fasta \
                 file {}",
                java_absolute_path(&dictionary_path)
            )),
        });
    }
    let dictionary: Vec<String> = header
        .sequences
        .iter()
        .map(|record| record.name.clone())
        .collect();
    let mut shards = Vec::new();
    for path in &inputs {
        shards.push((
            std::fs::read_to_string(path).map_err(|_| table_unreadable(path))?,
            path.clone(),
        ));
    }
    let pairs: Vec<(&str, &str)> = shards
        .iter()
        .map(|(text, path)| (text.as_str(), path.as_str()))
        .collect();
    let mut kept = 0;
    for (text, source) in &pairs {
        let (_, records) = pileup_summary::read_from_file(text, source).map_err(|error| {
            gather_error(gatk_tools::mutect_gathers::GatherError::PileupSummary(
                error,
            ))
        })?;
        if !records.is_empty() {
            kept += 1;
        }
    }
    let gathered = gatk_tools::mutect_gathers::gather_pileup_summaries(&pairs, &dictionary)
        .map_err(gather_error)?;
    write_file(&output, gathered.as_bytes())?;
    Ok(Some(format!("Successfully merged {kept} samples")))
}

/// `GatherNormalArtifactData`: every shard's records in the order given, under one header.
///
/// The concatenation is [`gatk_tools::mutect_gathers::gather_normal_artifact_data`]. The writer is
/// opened before the first shard is read, so a run refused over an unreadable shard still leaves
/// the file with its header and the records of the shards before it.
pub fn gather_normal_artifact_data(parser: &Parser) -> Outcome {
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let inputs = arguments(parser, "input");
    let mut shards = Vec::new();
    for path in &inputs {
        match std::fs::read_to_string(path) {
            Ok(text) => shards.push(text),
            Err(_) => {
                let texts: Vec<&str> = shards.iter().map(String::as_str).collect();
                write_file(
                    &output,
                    gatk_tools::mutect_gathers::gather_normal_artifact_data(&texts).as_bytes(),
                )?;
                return Err(table_unreadable(path));
            }
        }
    }
    let texts: Vec<&str> = shards.iter().map(String::as_str).collect();
    write_file(
        &output,
        gatk_tools::mutect_gathers::gather_normal_artifact_data(&texts).as_bytes(),
    )?;
    Ok(Some("SUCCESS".to_string()))
}

/// `AnnotatedIntervalCollection.create`, up to the records: the file must be readable, then its
/// NAME must be one the codec claims, then the codec reads it.
///
/// The codec claims `.seg`, `.maf` and `.maf.annotated` and nothing else, so a table named `.tsv`
/// is refused as a file that could not be parsed, before a line of it is read. A reader error
/// names the input by its URI, which is how tribble reports its source.
pub(crate) fn annotated_intervals(
    path: &str,
) -> Result<gatk_tools::annotated_interval::AnnotatedIntervalCollection, Thrown> {
    let uri = format!("file://{}", java_absolute_path(path));
    let meta = std::fs::metadata(path);
    let problem = match &meta {
        Err(_) => Some("It doesn't exist."),
        Ok(meta) if !meta.is_file() => Some("It isn't a regular file"),
        Ok(_) => None,
    };
    if let Some(problem) = problem {
        return Err(Thrown {
            failure: Failure::User,
            exception:
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile",
            message: Some(format!("Couldn't read file {uri}. Error was: {problem}")),
        });
    }
    if !(path.ends_with(".seg") || path.ends_with(".maf") || path.ends_with(".maf.annotated")) {
        return Err(bad_input(format!("Could not parse xsv file: {uri}")));
    }
    let text = std::fs::read_to_string(path)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
    gatk_tools::annotated_interval::read(&text).map_err(|error| Thrown {
        failure: Failure::Other,
        exception: error.java_class(),
        message: Some(error.message_with_source(&uri)),
    })
}

/// `getBestAvailableSequenceDictionary` for a tool that requires a reference: the master
/// dictionary when one was given, the reference's otherwise.
fn best_dictionary(parser: &Parser) -> Result<Vec<String>, Thrown> {
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    Ok(master
        .or(reference)
        .map(|header| {
            header
                .sequences
                .iter()
                .map(|record| record.name.clone())
                .collect()
        })
        .unwrap_or_default())
}

/// `MergeAnnotatedRegions`: overlapping regions of a segment file merged.
///
/// The merge and the file format are [`gatk_tools::annotated_interval`], where a golden measures
/// them. The runner adds the engine's startup before `traverse` reads the segments, and the
/// reading itself, whose refusals depend on the file's name before its content.
pub fn merge_annotated_regions(parser: &Parser) -> Outcome {
    use gatk_tools::annotated_interval::{merge_regions, DEFAULT_SEPARATOR};

    let segments = argument(parser, "segments").ok_or_else(|| {
        Thrown::command_line("Argument segments was missing: Argument 'segments' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let _ = resolve_read_filters(parser, "MergeAnnotatedRegions")?;
    let dictionary = best_dictionary(parser)?;
    let mut collection = annotated_intervals(&segments)?;
    gatk_tools::annotated_interval::check_sortable(&collection.records, &dictionary)
        .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?;
    collection.records = merge_regions(&collection.records, &dictionary, DEFAULT_SEPARATOR);
    write_file(&output, collection.write().as_bytes())?;
    Ok(None)
}

/// `getBestAvailableSequenceDictionary` for a tool whose reference is optional: `None` where
/// neither a master dictionary nor a reference was given.
fn optional_best_dictionary(parser: &Parser) -> Result<Option<Vec<String>>, Thrown> {
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    Ok(master.or(reference).map(|header| {
        header
            .sequences
            .iter()
            .map(|record| record.name.clone())
            .collect()
    }))
}

/// `MergeAnnotatedRegionsByAnnotation`: neighbouring regions merged when they agree on the named
/// annotations and lie within the distance.
///
/// The merge is [`gatk_tools::annotated_interval::merge_regions_by_annotation`]. The runner adds
/// the three checks `traverse` makes, in its order: every named annotation is in the file (the
/// missing ones listed in the iteration order of a `HashSet` of the names given), a dictionary is
/// available, which with no reference here comes only from `--sequence-dictionary`, and the
/// distance is not negative. The output is written from the first merged region's annotations, so
/// a file of no regions is refused only then.
pub fn merge_annotated_regions_by_annotation(parser: &Parser) -> Outcome {
    use gatk_engine::java_hash::JavaHashMap;
    use gatk_tools::annotated_interval::{
        merge_regions_by_annotation, write_without_header, DEFAULT_SEPARATOR,
    };

    let segments = argument(parser, "segments").ok_or_else(|| {
        Thrown::command_line("Argument segments was missing: Argument 'segments' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let names = arguments(parser, "annotations-to-match");
    let max_distance: i64 = scalar(parser, "max-merge-distance")
        .and_then(|value| value.parse().ok())
        .unwrap_or(1_000_000);
    let column =
        |name: &str, default: &str| argument(parser, name).unwrap_or_else(|| default.to_string());
    let illegal = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);

    let _ = resolve_read_filters(parser, "MergeAnnotatedRegionsByAnnotation")?;
    let dictionary = optional_best_dictionary(parser)?;
    let collection = annotated_intervals(&segments)?;
    let mut wanted: JavaHashMap<String, ()> =
        JavaHashMap::with_capacity(JavaHashMap::<String, ()>::copy_capacity(names.len()));
    for name in &names {
        wanted.insert(name.clone(), ());
    }
    let missing: Vec<String> = wanted
        .keys()
        .filter(|name| !collection.annotations.contains(name))
        .cloned()
        .collect();
    if !missing.is_empty() {
        return Err(illegal(format!(
            "Input file did not have all of the specified annotations.  Missing annotations were: \
             {}",
            missing.join(", ")
        )));
    }
    let dictionary = dictionary.ok_or_else(|| {
        illegal(
            "Sequence dictionary not available in the input file nor specified in a reference \
             parameter.  Please specify a reference with the -R parameter for this input file."
                .to_string(),
        )
    })?;
    if max_distance < 0 {
        return Err(illegal(
            "Cannot have a negative value for distance.".to_string(),
        ));
    }
    gatk_tools::annotated_interval::check_sortable(&collection.records, &dictionary)
        .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?;
    let merged = merge_regions_by_annotation(
        &collection.records,
        &dictionary,
        &names,
        DEFAULT_SEPARATOR,
        max_distance,
    );
    let text = write_without_header(
        &merged,
        &column("output-contig-column", "CONTIG"),
        &column("output-start-column", "START"),
        &column("output-end-column", "END"),
    )
    .map_err(|error| Thrown::non_user(error.java_class(), error.message()))?;
    write_file(&output, text.as_bytes())?;
    Ok(None)
}

/// `TagGermlineEvents`: tumour segments tagged where a called normal segment matches them.
///
/// The tagging is [`gatk_tools::tag_germline_events::tag_tumour_segments`], where a golden
/// measures it and its refusals. The runner adds the engine's startup, both files read in
/// `traverse` (the tumour first), and the output written from the tumour file's own header with
/// `POSSIBLE_GERMLINE` sorted into its annotations.
pub fn tag_germline_events(parser: &Parser) -> Outcome {
    use gatk_tools::tag_germline_events::tag_tumour_segments;

    let tumour_path = argument(parser, "segments").ok_or_else(|| {
        Thrown::command_line("Argument segments was missing: Argument 'segments' is required")
    })?;
    let normal_path = argument(parser, "called-matched-normal-seg-file").ok_or_else(|| {
        Thrown::command_line(
            "Argument called-matched-normal-seg-file was missing: Argument \
             'called-matched-normal-seg-file' is required",
        )
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let padding = number_or(parser, "endpoint-padding", 1000);
    let call = argument(parser, "input-call-header").unwrap_or_else(|| "CALL".to_string());
    let threshold: f64 = scalar(parser, "reciprocal-threshold")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.75);

    let _ = resolve_read_filters(parser, "TagGermlineEvents")?;
    let dictionary = best_dictionary(parser)?;
    let mut tumour = annotated_intervals(&tumour_path)?;
    let normal = annotated_intervals(&normal_path)?;
    let tagged = tag_tumour_segments(
        &tumour.records,
        &normal.records,
        &call,
        &dictionary,
        "POSSIBLE_GERMLINE",
        padding,
        threshold,
    )
    .map_err(|error| Thrown {
        failure: if error.java_class().contains("UserException") {
            Failure::User
        } else {
            Failure::Other
        },
        exception: error.java_class(),
        message: Some(error.message()),
    })?;
    tumour.records = tagged;
    tumour.annotations.push("POSSIBLE_GERMLINE".to_string());
    tumour.annotations.sort();
    write_file(&output, tumour.write().as_bytes())?;
    Ok(None)
}

/// `GatherTranches`: the tranches of a scattered VQSR run pooled by VQSLOD and cut again at the
/// requested sensitivities.
///
/// The pooling and the cut are [`gatk_tools::gather_tranches`], where a golden measures them. The
/// runner adds `doWork`'s order: htsjdk's `assertFileIsReadable` over EVERY input before any is
/// read, a `SAMException` naming the file by its URI; then each shard read in turn, a malformed
/// one named by the path as the command line gave it, into an output already created; then the
/// file written. `doWork` returns
/// `0`, which the tool prints.
pub fn gather_tranches(parser: &Parser) -> Outcome {
    use gatk_engine::tranches::{Mode, TrancheError};
    use gatk_tools::gather_tranches as tranches;

    let inputs = arguments(parser, "input");
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let mode = match scalar(parser, "mode").as_deref() {
        Some("SNP") => Mode::Snp,
        Some("INDEL") => Mode::Indel,
        Some("BOTH") => Mode::Both,
        _ => {
            return Err(Thrown::command_line(
                "Argument mode was missing: Argument 'mode' is required",
            ))
        }
    };
    let requested = arguments(parser, "truth-sensitivity-tranche");
    let levels: Vec<f64> = if requested.is_empty() {
        vec![100.0, 99.9, 99.0, 90.0]
    } else {
        requested
            .iter()
            .map(|value| value.parse().unwrap_or(f64::NAN))
            .collect()
    };

    for input in &inputs {
        let uri = format!("file://{}", java_absolute_path(input));
        let problem = match std::fs::metadata(input) {
            Err(_) => Some("Cannot read non-existent file: "),
            Ok(meta) if meta.is_dir() => Some("Cannot read file because it is a directory: "),
            Ok(_) => None,
        };
        if let Some(problem) = problem {
            let uri = if problem.contains("directory") {
                format!("{uri}/")
            } else {
                uri
            };
            return Err(Thrown::non_user(
                "htsjdk.samtools.SAMException",
                format!("{problem}{uri}"),
            ));
        }
    }
    // The output stream is opened before the first shard is read, so a refused shard leaves it
    // created and empty.
    write_file(&output, b"")?;
    let mut all = Vec::new();
    for input in &inputs {
        let text = std::fs::read_to_string(input)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{input}: {error}")))?;
        let shard = tranches::read_vqslod_tranches(&text).map_err(|mut error| {
            // The reader was given text, so a malformed-file refusal is named here, by the path
            // as the command line gave it: `MalformedFile(GATKPath, ...)` prints the raw string.
            match &mut error {
                TrancheError::HeaderLength { file, .. } | TrancheError::RowLength { file, .. } => {
                    *file = input.clone()
                }
                _ => {}
            }
            let class = error.class();
            Thrown {
                failure: if class.contains("UserException") {
                    Failure::User
                } else {
                    Failure::Other
                },
                exception: class,
                message: Some(error.message()),
            }
        })?;
        all.extend(shard);
    }
    let gathered = tranches::merge_and_convert(&all, &levels, mode);
    let text = format!(
        "{}{}",
        tranches::print_header(),
        tranches::tranches_string(&gathered)
    );
    write_file(&output, text.as_bytes())?;
    Ok(Some("0".to_string()))
}

/// `AS_FilterStatus` as `getAttributeAsString` reads it: a list the codec split on commas is
/// joined back with them, and an absent attribute is `None`.
fn as_filter_status(record: &htsjdk_vcf::variant::VariantContext) -> Option<String> {
    use htsjdk_vcf::variant::Value;
    record
        .attributes
        .iter()
        .find(|(key, _)| key == "AS_FilterStatus")
        .map(|(_, value)| match value {
            Value::List(items) => items
                .iter()
                .map(|item| item.format().unwrap_or_default())
                .collect::<Vec<_>>()
                .join(","),
            other => other.format().unwrap_or_default(),
        })
}

/// What a Mutect allele filter changes on a record: `VariantContextBuilder.filter` adds the name to
/// the filters the record had, a set, and `attribute` replaces `AS_FilterStatus`.
fn apply_allele_filter(
    record: &mut htsjdk_vcf::variant::VariantContext,
    filtered: bool,
    name: &str,
    status: Option<String>,
) {
    use htsjdk_vcf::variant::Value;
    if filtered {
        let mut filters = record.filters.clone().unwrap_or_default();
        if !filters.iter().any(|filter| filter == name) {
            filters.push(name.to_string());
        }
        record.filters = Some(filters);
    }
    if let Some(status) = status {
        match record
            .attributes
            .iter_mut()
            .find(|(key, _)| key == "AS_FilterStatus")
        {
            Some((_, value)) => *value = Value::Str(status),
            None => record
                .attributes
                .push(("AS_FilterStatus".to_string(), Value::Str(status))),
        }
    }
}

/// `NuMTFilterTool`: mitochondrial calls whose alternate depth a nuclear insertion could explain.
///
/// The cutoff and the per-allele decision are [`gatk_tools::numt_filter`], where a golden measures
/// them. The runner adds the walk and the writer: the header gains the `possible_numt` FILTER line
/// and nothing else, the cutoff is computed once the header is written, a refused record leaves
/// the file closed over the records before it, and every record's
/// genotypes are DECODED, since `getDataByAllele` reads each one's `AD`, so every record is
/// re-encoded rather than copied from the file.
pub fn numt_filter_tool(parser: &Parser) -> Outcome {
    use gatk_tools::numt_filter::{self as numt, Record};

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "NuMTFilterTool")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let double = |name: &str, default: f64| {
        scalar(parser, name)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(default)
    };
    let arguments = numt::Arguments {
        median_autosomal_coverage: double("autosomal-coverage", 0.0),
        max_numt_autosomal_copies: double("max-numt-autosomal-copies", 4.0),
    };
    let illegal = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut header = file.header.clone();
    if !header.has_filter_line(numt::FILTER_NAME) {
        header.lines.push(htsjdk_vcf::header::HeaderLine::Filter {
            id: numt::FILTER_NAME.to_string(),
            description: "Allele depth is below expected coverage of NuMT in autosome".to_string(),
        });
    }
    let cutoff = numt::cutoff_for(&arguments).map_err(|error| illegal(format!("{error:?}")))?;
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let mut written = Vec::new();
    for original in kept {
        let mut record = original.clone();
        let allele_depths: Vec<Option<Vec<i32>>> = record
            .genotypes
            .iter()
            .map(|genotype| genotype.ad.clone())
            .collect();
        let reduced = Record {
            alleles: record
                .alleles
                .iter()
                .map(|allele| allele.display_string())
                .collect(),
            allele_depths,
            filters: record.filters.clone().unwrap_or_default(),
            as_filter_status: as_filter_status(&record),
        };
        let applied = match numt::apply(&reduced, cutoff) {
            Ok(applied) => applied,
            Err(error) => {
                // `closeTool` runs in a `finally`, so the writer is closed over the header and the
                // records it was already given.
                apply_sites_only(parser, &mut header, &mut written);
                let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?;
                write_variant_output(parser, &output, &out)?;
                return Err(illegal(error.message()));
            }
        };
        let filtered = applied.filters.len() > reduced.filters.len();
        let status = (applied.as_filter_status != reduced.as_filter_status)
            .then(|| applied.as_filter_status.clone())
            .flatten();
        apply_allele_filter(&mut record, filtered, numt::FILTER_NAME, status);
        if keep(&record) {
            written.push(record);
        }
    }
    apply_sites_only(parser, &mut header, &mut written);
    let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    write_variant_output(parser, &output, &out)?;
    Ok(None)
}

/// A genotype's `AF`, as `getAttributeAsDoubleArray` reads it: `None` when the genotype has none,
/// and a `.` inside the list as a missing entry.
fn genotype_allele_fractions(genotype: &htsjdk_vcf::variant::Genotype) -> Option<Vec<Option<f64>>> {
    use htsjdk_vcf::variant::Value;
    let value = genotype
        .extended
        .iter()
        .find(|(key, _)| key == "AF")
        .map(|(_, value)| value)?;
    let items: Vec<String> = match value {
        Value::Missing => return None,
        Value::List(items) => items
            .iter()
            .map(|item| item.format().unwrap_or_default())
            .collect(),
        other => other
            .format()
            .unwrap_or_default()
            .split(',')
            .map(str::to_string)
            .collect(),
    };
    if items.len() == 1 && items[0] == "." {
        return None;
    }
    Some(items.iter().map(|item| item.parse::<f64>().ok()).collect())
}

/// `MTLowHeteroplasmyFilterTool`: every low-heteroplasmy allele filtered, but only once more of
/// them than the allowance passed every other filter.
///
/// The two passes are [`gatk_tools::mt_low_heteroplasmy::run`], where a golden measures them,
/// `--low-het-threshold` included: the field is a compile-time constant the command line never
/// reaches. The runner adds the walk. The header gains the `mt_many_low_hets` FILTER line; the
/// second pass reads the file afresh, so a record is re-encoded only when the file failed and its
/// genotypes were read, and is otherwise copied from the file. A first pass refused on a genotype
/// without `AF` leaves the file closed over its header alone.
pub fn mt_low_heteroplasmy_filter_tool(parser: &Parser) -> Outcome {
    use gatk_tools::mt_low_heteroplasmy::{self as low_het, Record};

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "MTLowHeteroplasmyFilterTool")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let arguments = low_het::Arguments {
        max_allowed_low_hets: number_or(parser, "max-allowed-low-hets", 3),
    };

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut header = file.header.clone();
    if !header.has_filter_line(low_het::FILTER_NAME) {
        header.lines.push(htsjdk_vcf::header::HeaderLine::Filter {
            id: low_het::FILTER_NAME.to_string(),
            description: "All low heteroplasmy sites are filtered when at least x low het sites \
                          pass all other filters"
                .to_string(),
        });
    }
    // The records as the module reads them, from clones: the first pass's decoding does not reach
    // the second pass, which reads the file again.
    let reduced: Vec<Record> = kept
        .iter()
        .map(|record| {
            let copy = (*record).clone();
            Record {
                alternates: copy
                    .alleles
                    .iter()
                    .skip(1)
                    .map(|allele| allele.display_string())
                    .collect(),
                allele_fractions: copy
                    .genotypes
                    .iter()
                    .map(genotype_allele_fractions)
                    .collect(),
                filters: copy.filters.clone().unwrap_or_default(),
                as_filter_status: as_filter_status(&copy),
            }
        })
        .collect();
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let results = match low_het::run(&reduced, &arguments) {
        Ok(results) => results,
        Err(error) => {
            let mut none: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
            apply_sites_only(parser, &mut header, &mut none);
            let out = htsjdk_vcf::vcf_file::write_vcf(&header, &none)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            write_variant_output(parser, &output, &out)?;
            let class = match error {
                low_het::LowHetError::NoAlleleFraction => "java.lang.NullPointerException",
                low_het::LowHetError::Filter(_) => "java.lang.IllegalArgumentException",
            };
            return Err(Thrown::non_user(class, error.message()));
        }
    };
    // The file failed exactly when a record changed: failing takes more unfiltered low sites than
    // the allowance, and each of those is then an artifact.
    let failed = results != reduced;
    let mut written = Vec::new();
    for ((original, before), after) in kept.iter().zip(&reduced).zip(&results) {
        let mut record = (*original).clone();
        if failed {
            // `areAllelesArtifacts` reads every record's genotypes once the file has failed.
            record.genotypes.decode();
            let filtered = after.filters != before.filters;
            let status = (after.as_filter_status != before.as_filter_status)
                .then(|| after.as_filter_status.clone())
                .flatten();
            apply_allele_filter(&mut record, filtered, low_het::FILTER_NAME, status);
        }
        if keep(&record) {
            written.push(record);
        }
    }
    apply_sites_only(parser, &mut header, &mut written);
    let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    write_variant_output(parser, &output, &out)?;
    Ok(None)
}

/// `FeatureManager`'s inputs: a `LinkedHashMap` keyed by `FeatureInput`, whose equality is the raw
/// argument, so a file named twice, on the command line or through a `.list`, is one input, in
/// the place it was first named.
fn distinct_feature_inputs(values: Vec<String>) -> Vec<String> {
    let mut distinct: Vec<String> = Vec::new();
    for value in values {
        if !distinct.contains(&value) {
            distinct.push(value);
        }
    }
    distinct
}

/// `ExampleMultiFeatureWalker`: every feature of every input, merged, printed as it is handed over.
///
/// The walk is [`gatk_engine::multi_feature_walker`], oracle-backed through this very tool. The
/// runner adds its startup, which is the SV evidence tools' own (the inputs opened, then the
/// dictionary chosen), and prints each feature's `toString` as the walk reaches it, so a walk
/// refused part way has already printed what came before. Only depth evidence is read; any other
/// feature type is a refusal of the port's own.
pub fn example_multi_feature_walker(parser: &Parser) -> Outcome {
    use gatk_engine::multi_feature_walker as walker;
    use gatk_tools::condense_depth_evidence as depth;
    use gatk_tools::sv_feature_codecs::{self as codecs, Encoding};
    use std::io::Write;

    let inputs = distinct_feature_inputs(arguments(parser, "feature"));
    let _ = resolve_read_filters(parser, "ExampleMultiFeatureWalker")?;
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    let mut located: Vec<Vec<walker::Located>> = Vec::new();
    for input in &inputs {
        let bytes = std::fs::read(input).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: java_absolute_path(input),
                }
                .message(),
            )
        })?;
        if codecs::find(input)
            != Some(codecs::Codec {
                feature_type: "DepthEvidence",
                encoding: Encoding::Text {
                    block_compressed: false,
                },
            })
        {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!("ExampleMultiFeatureWalker reads {input}, and only plain-text depth evidence is ported"),
            ));
        }
        let (_, records) = depth::read(&String::from_utf8_lossy(&bytes)).map_err(|problem| {
            Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "{input} is a malformed depth-evidence file ({problem}), which Tribble refuses"
                ),
            )
        })?;
        located.push(
            records
                .iter()
                .map(|record| walker::Located {
                    contig: record.contig.clone(),
                    start: record.start,
                    end: record.end,
                    text: std::iter::once(format!(
                        "{}\t{}\t{}",
                        record.contig, record.start, record.end
                    ))
                    .chain(record.counts.iter().map(|count| count.to_string()))
                    .collect::<Vec<_>>()
                    .join("\t"),
                })
                .collect(),
        );
    }
    let source = |header: Option<SamHeader>, name: &str| {
        header.map(|header| walker::DictSource {
            contigs: header
                .sequences
                .iter()
                .map(|record| record.name.clone())
                .collect(),
            source: name.to_string(),
        })
    };
    let user = |message: String| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException",
        message: Some(message),
    };
    let dictionary = walker::choose_dictionary(
        source(master, "sequence-dictionary"),
        source(reference, "reference"),
    )
    .map_err(|error| user(error.message()))?;
    // The merge is computed whole, so what a refused walk had printed is the features the heap
    // handed over before the one that went backwards: the prefix of the order up to the refusal.
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    match walker::merge(&located, &dictionary) {
        Ok(features) => {
            for feature in &features {
                let _ = writeln!(out, "{}", feature.text);
            }
            Ok(None)
        }
        Err(error) => {
            let printed = walker::merge_prefix(&located, &dictionary);
            for feature in &printed {
                let _ = writeln!(out, "{}", feature.text);
            }
            Err(user(error.message()))
        }
    }
}

/// `CalculateAverageCombinedAnnotations`: GenomicsDB's summed annotations divided by the number of
/// samples that carry a variant.
///
/// The header and the division are [`gatk_tools::calculate_average_combined_annotations`], where a
/// golden measures them. The runner adds the walk and the writer: a record without `RAW_GT_COUNT`
/// is refused where the walk reaches it, and `closeTool` closes the writer in a `finally`, so the
/// file is left with the header and the records before it. A record the tool does not rebuild
/// keeps its genotypes as the file wrote them.
pub fn calculate_average_combined_annotations(parser: &Parser) -> Outcome {
    use gatk_tools::calculate_average_combined_annotations as average;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "CalculateAverageCombinedAnnotations")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let annotations = arguments(parser, "summed-annotation-to-divide");
    let refused = |error: average::AverageError| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException",
        message: Some(error.message()),
    };

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut header = average::header_with_averages(&file.header, &annotations).map_err(refused)?;
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let mut written = Vec::new();
    let mut failure = None;
    for record in kept {
        match average::apply(record, &annotations) {
            Ok(out) => {
                if keep(&out) {
                    written.push(out);
                }
            }
            Err(error) => {
                failure = Some(refused(error));
                break;
            }
        }
    }
    apply_sites_only(parser, &mut header, &mut written);
    let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    write_variant_output(parser, &output, &out)?;
    match failure {
        Some(thrown) => Err(thrown),
        None => Ok(None),
    }
}

/// `GATKVariantContextUtils.isAlleleInList`: whether an alternate of one record is among the
/// alternates of another at the same start, once the shorter reference is extended to the longer.
///
/// Equal references compare the alternates as they are. Otherwise the longer reference is the
/// common one (`determineReferenceAllele`), and the other record's alleles are EXTENDED by the
/// bases the longer reference has past the shorter's length; two different references of one
/// length are `Err`, which is the `IllegalStateException` the caller turns into its own refusal.
fn is_allele_in_list(
    reference: &str,
    alternate: &str,
    other_reference: &str,
    other_alternates: &[String],
) -> Result<bool, ()> {
    let symbolic = |allele: &str| {
        allele.starts_with('<') || allele.contains('[') || allele.contains(']') || allele == "*"
    };
    let extend = |allele: &str, tail: &str| {
        if symbolic(allele) {
            allele.to_string()
        } else {
            format!("{allele}{tail}")
        }
    };
    if reference == other_reference {
        return Ok(other_alternates.iter().any(|other| other == alternate));
    }
    if reference.len() == other_reference.len() {
        return Err(());
    }
    if reference.len() > other_reference.len() {
        let tail = &reference[other_reference.len()..];
        Ok(other_alternates
            .iter()
            .any(|other| extend(other, tail) == alternate))
    } else {
        let tail = &other_reference[reference.len()..];
        let extended = extend(alternate, tail);
        Ok(other_alternates.contains(&extended))
    }
}

/// `FilterVariantTranches`: a CNN-scored VCF filtered at the scores its own resource sites reach.
///
/// The cutoffs, the band names and the header lines are [`gatk_tools::filter_variant_tranches`],
/// where a golden measures them. The runner adds the two passes over the input and what the
/// resources contribute:
///
/// * **the tranches are validated first**, SNP then indel, then the writer is opened and the
///   header written: every input line in sorted order, the FILTER lines dropped under
///   `--invalidate-previous-filters`, a line per tranche band, the samples as a `TreeSet`, and the
///   tool's default lines. An input without the score's INFO line is refused there;
/// * **the first pass** counts the records carrying the score by type, and takes a record's OWN
///   score once, at the first resource record overlapping it that shares its start and one of its
///   alternates (`isAlleleInList`). A record that is not a SNP counts as an indel there;
/// * **the second pass** sets the band on a filtered record and `PASS` on every record left with
///   no filter, so a record that came in as `.` goes out as `PASS`.
///
/// A resource is queried by interval, so one without an index is refused at its first query. A
/// refusal after the header leaves the file with the header alone, since `closeTool` closes the
/// writer in a `finally`. A score written as a list is a `ClassCastException` in the reference and
/// a refusal of the port's own here.
pub fn filter_variant_tranches(parser: &Parser) -> Outcome {
    use gatk_tools::filter_variant_tranches as tranches;
    use gatk_tools::remove_nearby_indels::{variant_type, VariantType};
    use htsjdk_vcf::variant::Value;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "FilterVariantTranches")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let info_key = argument(parser, "info-key").unwrap_or_else(|| "CNN_2D".to_string());
    let remove_old_filters = flag(parser, "invalidate-previous-filters");
    let doubles = |name: &str, default: f64| -> Vec<f64> {
        let values = arguments(parser, name);
        if values.is_empty() {
            vec![default]
        } else {
            values
                .iter()
                .map(|value| value.parse().unwrap_or(f64::NAN))
                .collect()
        }
    };
    let refused = |error: tranches::FilterVariantTranchesError| {
        let class = error.class();
        Thrown {
            failure: if class.contains("CommandLineException") {
                Failure::CommandLine
            } else {
                Failure::User
            },
            exception: class,
            message: Some(error.message()),
        }
    };
    let bad_input_at = |message: String| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(format!("Bad input: {message}")),
    };

    // The resources are feature inputs, opened at startup.
    let mut resources = Vec::new();
    for path in distinct_feature_inputs(arguments(parser, "resource")) {
        let bytes = std::fs::read(&path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: java_absolute_path(&path),
                }
                .message(),
            )
        })?;
        let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
            let mut inflated = String::new();
            std::io::Read::read_to_string(
                &mut flate2::read::MultiGzDecoder::new(bytes.as_slice()),
                &mut inflated,
            )
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
            inflated
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        let file = htsjdk_vcf::reader::read_vcf(&text)
            .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
        resources.push((path.clone(), file.records));
    }

    let snp_tranches =
        tranches::validate_tranches(&doubles("snp-tranche", 99.95)).map_err(refused)?;
    let indel_tranches =
        tranches::validate_tranches(&doubles("indel-tranche", 99.4)).map_err(refused)?;

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut header = file.header.clone();
    let has_info_key = header.lines.iter().any(|line| {
        matches!(line, htsjdk_vcf::header::HeaderLine::Compound { key, id, .. }
            if key == "INFO" && *id == info_key)
    });
    let header_only = |header: &htsjdk_vcf::header::VcfHeader| -> Result<(), Thrown> {
        let out = htsjdk_vcf::vcf_file::write_vcf(header, &[])
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        write_variant_output(parser, &output, &out)
    };
    if !has_info_key {
        // The writer was opened first, and the header was never written: an empty file.
        write_file(&output, b"")?;
        return Err(refused(
            tranches::FilterVariantTranchesError::InfoKeyNotInHeader(info_key.clone()),
        ));
    }
    if remove_old_filters {
        // A parsed FILTER line is `Structured`, a built one `Filter`: both go.
        header.lines.retain(|line| match line {
            htsjdk_vcf::header::HeaderLine::Filter { .. } => false,
            htsjdk_vcf::header::HeaderLine::Structured { key, .. } => key != "FILTER",
            _ => true,
        });
    }
    for (class, levels) in [
        (tranches::SNP_STRING, &snp_tranches),
        ("INDEL", &indel_tranches),
    ] {
        for (id, description) in tranches::tranche_header_lines(&info_key, class, levels) {
            if !header.has_filter_line(&id) {
                header
                    .lines
                    .push(htsjdk_vcf::header::HeaderLine::Filter { id, description });
            }
        }
    }
    header.samples.sort();
    header.lines.extend(default_tool_vcf_header_lines(
        parser,
        "FilterVariantTranches",
    ));

    let score_of = |record: &htsjdk_vcf::variant::VariantContext| -> Result<Option<f64>, Thrown> {
        match record.attributes.iter().find(|(key, _)| *key == info_key) {
            None => Ok(None),
            Some((_, Value::List(_))) => Err(Thrown::non_user(
                PORT_LIMITATION,
                format!("the {info_key} score is a list, which the reference casts to a String"),
            )),
            Some((_, value)) => {
                let text = value.format().unwrap_or_default();
                text.parse::<f64>().map(Some).map_err(|_| {
                    Thrown::non_user(
                        "java.lang.NumberFormatException",
                        format!("For input string: \"{text}\""),
                    )
                })
            }
        }
    };

    // The first pass.
    let (mut scored_snps, mut scored_indels) = (0usize, 0usize);
    let (mut snp_scores, mut indel_scores): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    for record in &kept {
        let has_score = record.attributes.iter().any(|(key, _)| *key == info_key);
        if !has_score {
            continue;
        }
        let kind = variant_type(record);
        match kind {
            VariantType::Snp => scored_snps += 1,
            VariantType::Indel => scored_indels += 1,
            _ => {}
        }
        let reference = record.reference().display_string();
        'resources: for (path, resource) in &resources {
            // `FeatureDataSource` asks for the index at a source's FIRST query, which is the first
            // scored record that reaches this resource.
            if !has_feature_index(path) {
                header_only(&header)?;
                return Err(Thrown::user(format!(
                    "Input {path} must support random access to enable queries by interval. If \
                     it's a file, please index it using the bundled tool IndexFeatureFile"
                )));
            }
            for other in resource.iter().filter(|other| {
                other.contig == record.contig
                    && other.start <= record.stop
                    && other.stop >= record.start
            }) {
                let other_reference = other.reference().display_string();
                let other_alternates: Vec<String> = other
                    .alternate_alleles()
                    .iter()
                    .map(|allele| allele.display_string())
                    .collect();
                for alternate in record.alternate_alleles() {
                    let matched = record.start == other.start
                        && match is_allele_in_list(
                            &reference,
                            &alternate.display_string(),
                            &other_reference,
                            &other_alternates,
                        ) {
                            Ok(matched) => matched,
                            Err(()) => {
                                header_only(&header)?;
                                return Err(bad_input_at(format!(
                                    "The provided variant file(s) have inconsistent references for the same position(s) at {}:{}, {} in input vs. {} in resource",
                                    other.contig,
                                    other.start,
                                    java_allele(record.reference()),
                                    java_allele(other.reference())
                                )));
                            }
                        };
                    if matched {
                        let score = score_of(record)?.unwrap_or(f64::NAN);
                        if kind == VariantType::Snp {
                            snp_scores.push(score);
                        } else {
                            indel_scores.push(score);
                        }
                        break 'resources;
                    }
                }
            }
        }
    }
    let (snp_cutoffs, indel_cutoffs) = match tranches::cutoffs(
        &snp_scores,
        &indel_scores,
        scored_snps,
        scored_indels,
        &snp_tranches,
        &indel_tranches,
        &info_key,
    ) {
        Ok(cutoffs) => cutoffs,
        Err(error) => {
            header_only(&header)?;
            return Err(refused(error));
        }
    };

    // The second pass.
    let keep = variant_output_filter(parser, intervals.as_deref())?;
    let mut written = Vec::new();
    for original in &kept {
        let mut record = (*original).clone();
        if remove_old_filters {
            record.filters = None;
        }
        if let Some(score) = score_of(&record)? {
            let kind = variant_type(&record);
            let band = if kind == VariantType::Snp
                && !snp_cutoffs.is_empty()
                && tranches::is_tranche_filtered(score, &snp_cutoffs)
            {
                Some(tranches::filter_string_from_score(
                    &info_key,
                    tranches::SNP_STRING,
                    score,
                    &snp_tranches,
                    &snp_cutoffs,
                ))
            } else if kind == VariantType::Indel
                && !indel_cutoffs.is_empty()
                && tranches::is_tranche_filtered(score, &indel_cutoffs)
            {
                Some(tranches::filter_string_from_score(
                    &info_key,
                    "INDEL",
                    score,
                    &indel_tranches,
                    &indel_cutoffs,
                ))
            } else {
                None
            };
            if let Some(band) = band {
                let mut filters = record.filters.clone().unwrap_or_default();
                if !filters.contains(&band) {
                    filters.push(band);
                }
                record.filters = Some(filters);
            }
        }
        if record
            .filters
            .as_ref()
            .is_none_or(|filters| filters.is_empty())
        {
            record.filters = Some(Vec::new());
        }
        if keep(&record) {
            written.push(record);
        }
    }
    apply_sites_only(parser, &mut header, &mut written);
    let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
    write_variant_output(parser, &output, &out)?;
    Ok(None)
}

/// `Allele.toString()`: the bases, with a `*` after a reference allele.
fn java_allele(allele: &htsjdk_vcf::allele::Allele) -> String {
    let text = allele.display_string();
    if allele.is_reference() {
        format!("{text}*")
    } else {
        text
    }
}

/// `DownsampleByDuplicateSet`: whole molecules kept or dropped, one seeded draw per molecule.
///
/// The grouping, the three rejection rules and the draws are
/// [`gatk_tools::downsample_by_duplicate_set`], where a golden measures them. The runner adds the
/// read walker's startup and filters, and the writer. A traversal that hands over no read at all
/// is the reference's own `NullPointerException`; a read without a well-formed `MI` is a refusal of
/// the port's own, since what the reference throws depends on how the tag is malformed.
pub fn downsample_by_duplicate_set(parser: &Parser) -> Outcome {
    use gatk_tools::downsample_by_duplicate_set as downsample;

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "DownsampleByDuplicateSet")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let fraction: f64 = scalar(parser, "fraction-to-keep")
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| {
            Thrown::command_line(
                "Argument fraction-to-keep was missing: Argument 'fraction-to-keep' is required",
            )
        })?;
    let arguments = downsample::Arguments {
        fraction_to_keep: fraction,
        minimum_reads: number_or(parser, "min-reads", 1).max(0) as usize,
        minimum_reads_per_strand: number_or(parser, "min-per-strand-reads", 0).max(0) as usize,
    };

    let filter = read_filter(parser, &filters, &header)?;
    let command_line = crate::command_line::expanded("DownsampleByDuplicateSet", parser);
    let options = gatk_tools::sam_output::Options {
        intervals: intervals.clone(),
        create_output_bam_index: flag(parser, "create-output-bam-index"),
        add_output_sam_program_record: flag(parser, "add-output-sam-program-record"),
        command_line: &command_line,
        version: crate::TOOLKIT_VERSION,
    };
    let (level, deflater) = output_compression(parser);
    match downsample::downsample_with(&source, &options, &filter, level, deflater, &arguments)
        .map_err(reads_traversal_error)?
    {
        None => Err(Thrown::non_user(
            PORT_LIMITATION,
            "DownsampleByDuplicateSet met a read without a well-formed MI tag, which the reference \
             fails on with a JVM exception",
        )),
        // `processLastReadSet` asks the set it never started for its reads: a helpful NPE naming
        // the field, which is deterministic.
        Some(Err(downsample::DownsampleError::NoReads)) => Err(Thrown::non_user(
            "java.lang.NullPointerException",
            "Cannot invoke \"org.broadinstitute.hellbender.tools.walkers.consensus.ReadsWithSameUMI.getReads()\" because \"this.currentReadsWithSameUMI\" is null",
        )),
        Some(Err(error)) => Err(Thrown::user(error.message())),
        Some(Ok((bytes, bai))) => {
            write_bam(parser, &output, &bytes, bai)?;
            // `onTraversalSuccess` returns the word, which `handleResult` prints.
            Ok(Some("SUCCESS".to_string()))
        }
    }
}

/// `ASEReadCounter`: reference and alternate counts at the heterozygous sites of a VCF.
///
/// The counting cascade, the overlapping-fragment handling and the line are
/// [`gatk_tools::ase_read_counter`], where a golden measures them. The runner adds the locus
/// walk and what `apply` asks of each locus, in its order:
///
/// * **the reference base first**, before anything else: with no `--reference` the context holds
///   no bases and `getBase()` indexes an empty array, which is the reference's own
///   `ArrayIndexOutOfBoundsException`;
/// * **then the variants at the locus**, a feature query that needs each input's index, refused
///   at the first locus when one has none, and more than one variant at a locus is refused;
/// * a site that is not biallelic, or where no genotype is heterozygous, is skipped with only a
///   warning, and one whose alternate has no bases is refused.
///
/// The header is written when the traversal starts and every line as its site is reached, so a
/// refusal mid-walk leaves the header and the lines before it, on stdout when no `--output` is
/// given.
pub fn ase_read_counter(parser: &Parser) -> Outcome {
    use gatk_tools::ase_read_counter as ase;
    use htsjdk_vcf::genotype_type::{determine_type, GenotypeType};

    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "ASEReadCounter")?;
    let output = argument(parser, "output");
    let variant_paths = distinct_feature_inputs(arguments(parser, "variant"));
    let format = match scalar(parser, "output-format").as_deref() {
        Some("TABLE") => ase::OutputFormat::Table,
        Some("CSV") => ase::OutputFormat::Csv,
        _ => ase::OutputFormat::RTable,
    };
    let count_type = match scalar(parser, "count-overlap-reads-handling").as_deref() {
        Some("COUNT_READS") => ase::CountType::CountReads,
        Some("COUNT_FRAGMENTS") => ase::CountType::CountFragments,
        _ => ase::CountType::CountFragmentsRequireSameBase,
    };
    let minimum_depth = number_or(parser, "min-depth-of-non-filtered-base", -1);
    let minimum_mapping_quality = number_or(parser, "min-mapping-quality", 0);
    let minimum_base_quality = number_or(parser, "min-base-quality", 0).clamp(0, 255) as u8;

    let mut variants = Vec::new();
    for path in &variant_paths {
        let file =
            htsjdk_vcf::reader::read_vcf(&feature_text(path)?).map_err(|failure| Thrown {
                failure: Failure::User,
                exception: failure.error.class(),
                message: Some(failure.error.message()),
            })?;
        variants.push((path.clone(), file.records));
    }
    let mut reference = match argument(parser, "reference") {
        Some(path) => Some(
            gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?,
        ),
        None => None,
    };

    let filter = read_filter(parser, &filters, &header)?;
    let records = gatk_tools::read_walker::traverse(&source, &intervals, &|_| true)
        .map_err(reads_traversal_error)?;
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
    .map_err(locus_traversal_error)?;

    // `onTraversalStart` opens the stream and writes the header before the first locus.
    let mut text = ase::header(format);
    text.push('\n');
    let finish = |text: &str| -> Result<(), Thrown> {
        match &output {
            Some(path) => write_file(path, text.as_bytes()),
            None => {
                print!("{text}");
                Ok(())
            }
        }
    };
    let refuse = |text: &str, thrown: Thrown| -> Outcome {
        finish(text)?;
        Err(thrown)
    };
    let mut queried = false;
    for one in &applied {
        let contig = &one.context.contig;
        let position = one.context.position;
        let base = match reference.as_mut() {
            Some(source) => match source.query(contig, position, position) {
                Ok(bases) => bases[0],
                Err(error) => return refuse(&text, Thrown::user(format!("{error:?}"))),
            },
            None => {
                return refuse(
                    &text,
                    Thrown::non_user(
                        "java.lang.ArrayIndexOutOfBoundsException",
                        "Index 0 out of bounds for length 0",
                    ),
                )
            }
        };
        if !queried {
            queried = true;
            if let Some((path, _)) = variants.iter().find(|(path, _)| !has_feature_index(path)) {
                return refuse(
                    &text,
                    Thrown::user(format!(
                        "Input {path} must support random access to enable queries by interval. \
                         If it's a file, please index it using the bundled tool IndexFeatureFile"
                    )),
                );
            }
        }
        let at: Vec<&htsjdk_vcf::variant::VariantContext> = variants
            .iter()
            .flat_map(|(_, records)| records.iter())
            .filter(|record| {
                record.contig == *contig
                    && record.start <= i64::from(position)
                    && i64::from(position) <= record.stop
            })
            .collect();
        if at.len() > 1 {
            return refuse(
                &text,
                Thrown::user(format!(
                    "More then one variant context at position: {contig}:{position}"
                )),
            );
        }
        let Some(site) = at.first() else {
            continue;
        };
        if site.alleles.len() != 2 {
            continue;
        }
        let hets = site
            .genotypes
            .iter()
            .filter(|genotype| determine_type(genotype) == GenotypeType::Het)
            .count();
        if hets < 1 {
            continue;
        }
        let alternate = site.alleles[1].display_string();
        if alternate.is_empty() || site.alleles[1].is_symbolic() {
            return refuse(
                &text,
                Thrown::user(
                    "The file of variant sites must contain heterozygous sites and cannot be a \
                     GVCF file containing <NON_REF> alleles."
                        .to_string(),
                ),
            );
        }
        let alternate = alternate.as_bytes()[0];
        let pileup = ase::filter_pileup(&one.context.pileup, count_type);
        let counts = ase::count_site(
            &pileup,
            base,
            alternate,
            minimum_mapping_quality,
            minimum_base_quality,
        );
        if let Some(line) = ase::line(
            contig,
            position,
            &site.id,
            base,
            alternate,
            counts,
            minimum_depth,
            format,
        ) {
            text.push_str(&line);
            text.push('\n');
        }
    }
    finish(&text)?;
    Ok(None)
}

/// `VariantContext.toString()` for a record read from a VCF, which a walker's wrapper prints whole.
///
/// The same string [`variant_context_to_string`] builds for `VariantsToTable`'s own record type, over
/// htsjdk-rs's instead. The source is the caller's: `MultiVariantDataSource` names each record by
/// the path it came from, where a plain variant walker's records say `Unknown`. The span is
/// `getEnd()` and so reads `END`, the
/// alleles are sorted reference first, the attributes are a `TreeMap`'s `toString`, and the
/// genotypes are the unparsed text while they are still lazy and `GenotypesContext.toString()`
/// otherwise, which for a file with no samples is `[]`.
fn java_variant_context_string(
    record: &htsjdk_vcf::variant::VariantContext,
    source: &str,
) -> String {
    use htsjdk_vcf::variant::Value;
    let position = if record.start == record.stop {
        format!("{}:{}", record.contig, record.start)
    } else {
        format!("{}:{}-{}", record.contig, record.start, record.stop)
    };
    let qual = if record.has_log10_p_error() {
        format!("{:.2}", record.phred_scaled_qual())
    } else {
        ".".to_string()
    };
    let reference = record.reference().display_string();
    let alternates: Vec<String> = record
        .alternate_alleles()
        .iter()
        .map(|allele| allele.display_string())
        .collect();
    let mut sorted = alternates.clone();
    sorted.sort();
    let alleles = std::iter::once(format!("{reference}*"))
        .chain(sorted)
        .collect::<Vec<String>>()
        .join(", ");
    fn java(value: &Value) -> String {
        match value {
            Value::Missing => "null".to_string(),
            Value::Bool(flag) => flag.to_string(),
            Value::Str(text) => text.clone(),
            Value::List(items) => format!(
                "[{}]",
                items.iter().map(java).collect::<Vec<String>>().join(", ")
            ),
            other => other.format().unwrap_or_default(),
        }
    }
    let mut attributes: Vec<(String, String)> = record
        .attributes
        .iter()
        .map(|(key, value)| (key.clone(), java(value)))
        .collect();
    attributes.sort();
    let attributes = attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<String>>()
        .join(", ");
    let genotypes = match record.genotypes.unparsed() {
        Some(text) => text.to_string(),
        None if record.genotypes.is_empty() => "[]".to_string(),
        None => "[...]".to_string(),
    };
    format!(
        "[VC {source} @ {position} Q{qual} of type={} alleles=[{alleles}] attr={{{attributes}}} GT={genotypes} filters={}",
        variant_type_name(&reference, &alternates),
        record.filters.clone().unwrap_or_default().join(",")
    )
}

/// `getAttributeAsString(key, default)`: a list is joined at `,`, whatever it held.
fn attribute_as_string(record: &htsjdk_vcf::variant::VariantContext, key: &str) -> Option<String> {
    use htsjdk_vcf::variant::Value;
    record
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| match value {
            Value::Str(text) => text.clone(),
            Value::List(items) => items
                .iter()
                .map(|item| item.format().unwrap_or_default())
                .collect::<Vec<String>>()
                .join(","),
            Value::Bool(flag) => flag.to_string(),
            other => other.format().unwrap_or_default(),
        })
}

/// A recal file's record, as `ApplyVQSR` reads one.
struct RecalVariant<'a>(&'a htsjdk_vcf::variant::VariantContext);

impl gatk_tools::apply_vqsr::RecalRecord for RecalVariant<'_> {
    fn start(&self) -> i32 {
        self.0.start as i32
    }
    fn end(&self) -> i32 {
        self.0.stop as i32
    }
    fn lod_string(&self) -> Option<String> {
        attribute_as_string(self.0, gatk_tools::apply_vqsr::VQS_LOD_KEY)
    }
    fn culprit(&self) -> Option<String> {
        attribute_as_string(self.0, gatk_tools::apply_vqsr::CULPRIT_KEY)
    }
    fn has_positive_label(&self) -> bool {
        self.0
            .attributes
            .iter()
            .any(|(key, _)| key == gatk_tools::apply_vqsr::POSITIVE_LABEL_KEY)
    }
    fn has_negative_label(&self) -> bool {
        self.0
            .attributes
            .iter()
            .any(|(key, _)| key == gatk_tools::apply_vqsr::NEGATIVE_LABEL_KEY)
    }
}

impl gatk_tools::apply_vqsr::AllelicRecalRecord for RecalVariant<'_> {
    fn first_alternate(&self) -> String {
        self.0
            .alternate_alleles()
            .first()
            .map(|allele| allele.display_string())
            .unwrap_or_default()
    }
}

/// `ApplyVQSR`: a VQSR recal file and its tranches applied to the variants they were built from.
///
/// The cut, the filter names, the header lines and both filtering paths are
/// [`gatk_tools::apply_vqsr`], each measured by a golden of its own. The runner is the walker around
/// them, and it is the first `MultiVariantWalker` here:
///
/// * **one `--variant`**: the collection is read, and more than one input is the port's limitation
///   rather than a merge, since a single input's merged header is its own header;
/// * **`onTraversalStart` refuses before the writer exists**: the tranches file (a missing
///   `--tranches-file` is the reference's `NullPointerException`), a previous run's malformed filter
///   name, the two cutoffs given together, and a level no tranche reaches all leave no file;
/// * **the header is a `HashSet` of lines**: the input's, the four VQSR `INFO` lines, `END`, the
///   `PASS` filter line, the three allele-specific lines under `-AS`, the tranche or `LOW_VQSLOD`
///   lines and the tool's own, with an identical line collapsing and nothing else;
/// * **every record queries the recal file first**, so a recal file with no index is refused at the
///   first record, and like every other failure inside `apply` it reaches the user as a
///   `GATKException` naming the locus and the whole record;
/// * **a record of the other mode, or filtered and not ignored, is written untouched**, and
///   `--exclude-filtered` drops only a recalibrated record whose new filter is neither `PASS` nor
///   `.`.
pub fn apply_vqsr(parser: &Parser) -> Outcome {
    use gatk_engine::tranches::{Mode, TruthSensitivityTranche};
    use gatk_tools::apply_vqsr as vqsr;
    use gatk_tools::remove_nearby_indels::{variant_type, VariantType};
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    use htsjdk_vcf::variant::Value;

    let _ = resolve_read_filters(parser, "ApplyVQSR")?;
    let inputs = arguments(parser, "variant");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --variant is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup_over(parser, input)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    // The recal file is a feature input, opened at startup and queried at every record.
    let recal_path = argument(parser, "recal-file").ok_or_else(|| {
        Thrown::command_line("Argument recal-file was missing: Argument 'recal-file' is required")
    })?;
    let recal_bytes = std::fs::read(&recal_path).map_err(|_| {
        Thrown::user(
            index_feature_file::Refusal::CouldNotReadInputFile {
                path: java_absolute_path(&recal_path),
            }
            .message(),
        )
    })?;
    let recal_text = if gatk_tools::read_walker_refusal::is_block_compressed(&recal_bytes) {
        let mut inflated = String::new();
        std::io::Read::read_to_string(
            &mut flate2::read::MultiGzDecoder::new(recal_bytes.as_slice()),
            &mut inflated,
        )
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{recal_path}: {error}")))?;
        inflated
    } else {
        String::from_utf8_lossy(&recal_bytes).into_owned()
    };
    let recal_file = htsjdk_vcf::reader::read_vcf(&recal_text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let recal_indexed = has_feature_index(&recal_path);

    let use_as = flag(parser, "use-allele-specific-annotations");
    let exclude_filtered = flag(parser, "exclude-filtered");
    let ignore_all_filters = flag(parser, "ignore-all-filters");
    let mode = match scalar(parser, "mode").as_deref() {
        Some("INDEL") => Mode::Indel,
        Some("BOTH") => Mode::Both,
        _ => Mode::Snp,
    };
    let number = |name: &str| -> Option<f64> {
        // A boxed `Double` is not a string to `argument`, which reads it as absent.
        scalar(parser, name).map(|value| value.parse().unwrap_or(f64::NAN))
    };
    let level = number("truth-sensitivity-filter-level");
    let cutoff = number("lod-score-cutoff");
    let mut ignored: Vec<String> = arguments(parser, "ignore-filter");
    ignored.sort();
    ignored.dedup();

    let refused = |error: vqsr::ApplyVqsrError| Thrown {
        failure: if error.class().contains("UserException") {
            Failure::User
        } else {
            Failure::Other
        },
        exception: error.class(),
        message: Some(error.message()),
    };

    // `onTraversalStart`: the tranches first, sorted by sensitivity, before anything else is read.
    let tranches: Vec<TruthSensitivityTranche> = match level {
        None => Vec::new(),
        Some(_) => {
            let Some(path) = argument(parser, "tranches-file") else {
                return Err(Thrown::non_user(
                    "java.lang.NullPointerException",
                    "Cannot invoke \"org.broadinstitute.hellbender.engine.GATKPath.toPath()\" because \"f\" is null",
                ));
            };
            let text = std::fs::read_to_string(&path).map_err(|_| {
                Thrown::user(format!(
                    "Couldn't read file {path}. Error was: Can't read tranches file with exception: {path}"
                ))
            })?;
            gatk_engine::tranches::read_tranches(&path, &text).map_err(|error| Thrown {
                failure: if error.class().contains("UserException") {
                    Failure::User
                } else {
                    Failure::Other
                },
                exception: error.class(),
                message: Some(error.message()),
            })?
        }
    };
    let kept = level.map(|level| vqsr::keep(&tranches, level));

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let mut header = file.header.clone();
    let compound = |id: &str, number: Cardinality, line_type: LineType, description: &str| {
        HeaderLine::Compound {
            key: "INFO".to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: description.to_string(),
            extra: Vec::new(),
        }
    };
    let mut added = vec![
        compound(
            "END",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Stop position of the interval",
        ),
        compound(
            vqsr::VQS_LOD_KEY,
            Cardinality::Fixed(1),
            LineType::Float,
            "Log odds of being a true variant versus being false under the trained gaussian mixture model",
        ),
        compound(
            vqsr::CULPRIT_KEY,
            Cardinality::Fixed(1),
            LineType::String,
            "The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out",
        ),
        // A `Flag` line built with a count of one is rewritten to zero by htsjdk's constructor.
        compound(
            vqsr::POSITIVE_LABEL_KEY,
            Cardinality::Fixed(0),
            LineType::Flag,
            "This variant was used to build the positive training set of good variants",
        ),
        compound(
            vqsr::NEGATIVE_LABEL_KEY,
            Cardinality::Fixed(0),
            LineType::Flag,
            "This variant was used to build the negative training set of bad variants",
        ),
        HeaderLine::Filter {
            id: vqsr::PASSES_FILTERS.to_string(),
            description: "Site contains at least one allele that passes filters".to_string(),
        },
    ];
    if use_as {
        added.push(compound(
            vqsr::AS_FILTER_STATUS_KEY,
            Cardinality::A,
            LineType::String,
            "Filter status for each allele, as assessed by ApplyVQSR. Note that the VCF filter field will reflect the most lenient/sensitive status across all alleles.",
        ));
        added.push(compound(
            vqsr::AS_CULPRIT_KEY,
            Cardinality::A,
            LineType::String,
            "For each alt allele, the annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out",
        ));
        added.push(compound(
            vqsr::AS_VQS_LOD_KEY,
            Cardinality::A,
            LineType::String,
            "For each alt allele, the log odds of being a true variant versus being false under the trained gaussian mixture model",
        ));
    }

    // `checkForPreviousApplyRecalRun`, over the INPUT's filter lines.
    let filter_ids: Vec<String> = file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Filter { id, .. } => Some(id.clone()),
            HeaderLine::Structured { key, fields } if key == "FILTER" => fields
                .iter()
                .find(|(name, _)| name == "ID")
                .map(|(_, value)| value.clone()),
            _ => None,
        })
        .collect();
    let runs = vqsr::previous_runs(&filter_ids).map_err(refused)?;

    let cut = match &kept {
        Some(kept) => {
            if cutoff.is_some() {
                return Err(refused(vqsr::ApplyVqsrError::MutuallyExclusiveCutoffs));
            }
            let lines =
                vqsr::tranche_filter_lines(kept, level.unwrap_or_default()).map_err(refused)?;
            for line in lines {
                added.push(HeaderLine::Filter {
                    id: line.id,
                    description: line.description,
                });
            }
            vqsr::Cut::Tranches(kept.clone())
        }
        None => {
            let line = vqsr::low_vqslod_filter_line(cutoff);
            added.push(HeaderLine::Filter {
                id: line.id,
                description: line.description,
            });
            vqsr::Cut::Lod(cutoff.unwrap_or(vqsr::DEFAULT_VQSLOD_CUTOFF))
        }
    };
    // A `HashSet`: an identical line collapses, one differing in anything stays beside it.
    let same = |a: &HeaderLine, b: &HeaderLine| a.render() == b.render();
    for line in added {
        if !header.lines.iter().any(|existing| same(existing, &line)) {
            header.lines.push(line);
        }
    }
    header.samples.sort();
    header.samples.dedup();
    header
        .lines
        .extend(default_tool_vcf_header_lines(parser, "ApplyVQSR"));

    let kept_records = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let is_of_mode = |kind: VariantType| match mode {
        Mode::Snp => matches!(kind, VariantType::Snp | VariantType::Mnp),
        Mode::Indel => matches!(
            kind,
            VariantType::Indel | VariantType::Mixed | VariantType::Symbolic
        ),
        Mode::Both => true,
    };
    let both_modes_were_run =
        (mode == Mode::Snp && runs.indel) || (mode == Mode::Indel && runs.snp);

    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut written = written.to_vec();
        apply_sites_only(parser, &mut header, &mut written);
        let out = htsjdk_vcf::vcf_file::write_vcf(&header, &written)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        write_variant_output(parser, &output, &out)
    };
    for record in kept_records {
        // `MultiVariantWalker.traverse` catches whatever `apply` throws and rethrows it as a
        // `GATKException` naming the locus and the whole record, so the cause is not what prints.
        let wrapped = || {
            Thrown::non_user(
                "org.broadinstitute.hellbender.exceptions.GATKException",
                format!(
                    "Exception thrown at {}:{} {}",
                    record.contig,
                    record.start,
                    java_variant_context_string(record, &input)
                ),
            )
        };
        // `featureContext.getValues(recal, vc.getStart())` runs first, whatever the record is.
        if !recal_indexed {
            finish(&written)?;
            return Err(wrapped());
        }
        let recals: Vec<RecalVariant> = recal_file
            .records
            .iter()
            .filter(|other| {
                other.contig == record.contig
                    && other.start == record.start
                    && other.start <= record.stop
                    && other.stop >= record.start
            })
            .map(RecalVariant)
            .collect();
        let kind = variant_type(record);
        let of_mode = is_of_mode(kind);
        let filters = record.filters.clone().unwrap_or_default();
        let evaluate = use_as || of_mode;
        if !vqsr::recalibrates(evaluate, &filters, ignore_all_filters, &ignored) {
            written.push(record.clone());
            continue;
        }
        let mut out = record.clone();
        let set = |out: &mut htsjdk_vcf::variant::VariantContext, key: &str, value: Value| match out
            .attributes
            .iter_mut()
            .find(|(name, _)| name == key)
        {
            Some(slot) => slot.1 = value,
            None => out.attributes.push((key.to_string(), value)),
        };
        let description = java_variant_context_string(record, &input);
        let filter = if !use_as {
            let (annotation, lod) = match vqsr::site_specific_filtering(
                record.start as i32,
                record.stop as i32,
                &recals,
                &description,
            ) {
                Ok(found) => found,
                Err(_) => {
                    finish(&written)?;
                    return Err(wrapped());
                }
            };
            set(&mut out, vqsr::VQS_LOD_KEY, Value::Str(annotation.vqslod));
            set(&mut out, vqsr::CULPRIT_KEY, Value::Str(annotation.culprit));
            if annotation.positive_label {
                set(&mut out, vqsr::POSITIVE_LABEL_KEY, Value::Bool(true));
            }
            if annotation.negative_label {
                set(&mut out, vqsr::NEGATIVE_LABEL_KEY, Value::Bool(true));
            }
            cut.filter(lod)
        } else {
            let previous = (runs.snp || runs.indel).then(|| {
                vqsr::PreviousAlleleLists::from_attributes(
                    &attribute_as_string(record, vqsr::AS_CULPRIT_KEY).unwrap_or_default(),
                    &attribute_as_string(record, vqsr::AS_VQS_LOD_KEY).unwrap_or_default(),
                    &attribute_as_string(record, vqsr::AS_FILTER_STATUS_KEY).unwrap_or_default(),
                )
            });
            let reference = record.reference().display_string();
            let alternates: Vec<String> = record
                .alternate_alleles()
                .iter()
                .map(|allele| allele.display_string())
                .collect();
            let site = vqsr::AlleleSpecificSite {
                start: record.start as i32,
                end: record.stop as i32,
                reference: &reference,
                alternates: &alternates,
                record: &description,
            };
            let (annotations, best) = match vqsr::allele_specific_filtering_with(
                &site,
                &recals,
                mode,
                &cut,
                previous.as_ref(),
            ) {
                Ok(found) => found,
                Err(_) => {
                    finish(&written)?;
                    return Err(wrapped());
                }
            };
            if annotations.positive_label {
                set(&mut out, vqsr::POSITIVE_LABEL_KEY, Value::Bool(true));
            }
            if annotations.negative_label {
                set(&mut out, vqsr::NEGATIVE_LABEL_KEY, Value::Bool(true));
            }
            for (key, list) in [
                (vqsr::AS_FILTER_STATUS_KEY, &annotations.filter_status),
                (vqsr::AS_VQS_LOD_KEY, &annotations.vqslod),
                (vqsr::AS_CULPRIT_KEY, &annotations.culprit),
            ] {
                if !list.is_empty() {
                    set(&mut out, key, Value::Str(list.join(vqsr::LIST_DELIMITER)));
                }
            }
            let previous_status = attribute_as_string(record, vqsr::AS_FILTER_STATUS_KEY);
            vqsr::site_filter_from_alleles_with(
                kind == VariantType::Mixed,
                of_mode,
                both_modes_were_run,
                previous_status.as_deref(),
                best,
                &cut,
            )
        };
        out.filters = match filter.as_str() {
            vqsr::PASSES_FILTERS => Some(Vec::new()),
            vqsr::UNFILTERED => None,
            other => Some(vec![other.to_string()]),
        };
        if vqsr::writes_out(true, &filter, exclude_filtered) {
            written.push(out);
        }
    }
    finish(&written)?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `SVStratify`: every SV record labelled with the stratum its type, size and track overlap put it
/// in, written to one file or to one file per stratum.
///
/// The engine is [`gatk_tools::sv_stratify`] and the record [`gatk_tools::sv_call_record`]; the
/// runner is the configuration loading and the writers, in `onTraversalStart`'s order:
///
/// * **the master dictionary is required**, and `--sequence-dictionary` is the only thing that
///   supplies it, a reference included;
/// * **the tracks are loaded before the table**, each through `IntervalUtils.loadIntervals`, and a
///   name given twice is refused after its file was read;
/// * **the table's columns are checked when its header is read**, the rows parsed one by one, a
///   stratum validated as it is built and a duplicate name refused as it is added;
/// * **the writers exist before the first record**: one per stratum and a `default` one under
///   `--split-output`, each a BGZF VCF named `<prefix>.<stratum>.vcf.gz` in the output directory,
///   created in the engine's `HashMap` order, so a refusal inside `apply` leaves every file behind
///   with its header;
/// * **every failure inside `apply` is the walker's `GATKException`**, the multiple-match refusal
///   included, since the traversal wraps whatever `apply` throws.
pub fn sv_stratify(parser: &Parser) -> Outcome {
    use gatk_tools::sv_stratify as stratify;

    let _ = resolve_read_filters(parser, "SVStratify")?;
    let inputs = arguments(parser, "variant");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --variant is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup_over(parser, input)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let config = argument(parser, "stratify-config").ok_or_else(|| {
        Thrown::command_line(
            "Argument stratify-config was missing: Argument 'stratify-config' is required",
        )
    })?;
    let prefix = argument(parser, "output-prefix");
    let split_output = flag(parser, "split-output");
    let allow_multiple = flag(parser, "allow-multiple-matches");
    let thresholds = stratification_thresholds(parser);

    let illegal = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);
    let gatk = |message: String| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            message,
        )
    };
    // `onTraversalStart`.
    let Some(dictionary) = master_dictionary(parser)? else {
        return Err(illegal(
            "Reference dictionary is required; please specify with --sequence-dictionary"
                .to_string(),
        ));
    };
    let sequences: Vec<(String, i32)> = dictionary
        .sequences
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();

    let engine = load_stratification_engine(parser, &config, &dictionary)?;

    // `initializeWriters`.
    let prefix = if split_output {
        let Some(prefix) = prefix else {
            return Err(illegal(
                "Argument --output-prefix required if using --split-output".to_string(),
            ));
        };
        if !std::path::Path::new(&output).is_dir() {
            return Err(illegal(
                "Argument --output must be a directory if using split-output".to_string(),
            ));
        }
        Some(prefix)
    } else {
        None
    };
    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let mut header = file.header.clone();
    let strat_line = htsjdk_vcf::header::HeaderLine::Compound {
        key: "INFO".to_string(),
        id: "STRAT".to_string(),
        number: htsjdk_vcf::header::Cardinality::Fixed(1),
        line_type: htsjdk_vcf::header::LineType::String,
        description: "Stratum ID".to_string(),
        extra: Vec::new(),
    };
    if !header
        .lines
        .iter()
        .any(|line| same_compound_id(line, &strat_line))
    {
        header.lines.push(strat_line);
    }
    header.samples.sort();
    header.samples.dedup();
    // One writer per file, in the order they were created: `default` first, then each stratum.
    let mut writers: Vec<(String, String, Vec<htsjdk_vcf::variant::VariantContext>)> = match &prefix
    {
        Some(prefix) => stratify::split_output_files(&engine, prefix)
            .into_iter()
            .zip(
                std::iter::once(stratify::DEFAULT_STRATUM.to_string())
                    .chain(engine.strata.iter().map(|stratum| stratum.name.clone())),
            )
            .map(|(file, stratum)| {
                (
                    stratum,
                    std::path::Path::new(&output)
                        .join(file)
                        .display()
                        .to_string(),
                    Vec::new(),
                )
            })
            .collect(),
        None => vec![(
            stratify::DEFAULT_STRATUM.to_string(),
            output.clone(),
            Vec::new(),
        )],
    };
    let close = |writers: &[(String, String, Vec<htsjdk_vcf::variant::VariantContext>)]| {
        for (_, path, records) in writers {
            let mut header = header.clone();
            let mut records = records.clone();
            apply_sites_only(parser, &mut header, &mut records);
            let out = htsjdk_vcf::vcf_file::write_vcf(&header, &records)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            write_variant_output(parser, path, &out)?;
        }
        Ok::<(), Thrown>(())
    };

    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    for record in kept {
        let wrapped = || {
            gatk(format!(
                "Exception thrown at {}:{} {}",
                record.contig,
                record.start,
                java_variant_context_string(record, &input)
            ))
        };
        let mut bare = record.clone();
        bare.genotypes = Vec::new().into();
        let sv = match gatk_tools::sv_call_record::create(&bare, &sequences) {
            Ok(sv) => sv,
            Err(_) => {
                close(&writers)?;
                return Err(wrapped());
            }
        };
        let written = match stratify::apply(
            &engine,
            &sv.stratify_record(),
            thresholds,
            allow_multiple,
            split_output,
        ) {
            Ok(written) => written,
            Err(_) => {
                close(&writers)?;
                return Err(wrapped());
            }
        };
        for entry in written {
            let mut labelled = record.clone();
            match labelled
                .attributes
                .iter_mut()
                .find(|(key, _)| key == "STRAT")
            {
                Some(slot) => slot.1 = htsjdk_vcf::variant::Value::Str(entry.stratum.clone()),
                None => labelled.attributes.push((
                    "STRAT".to_string(),
                    htsjdk_vcf::variant::Value::Str(entry.stratum.clone()),
                )),
            }
            let target = writers
                .iter_mut()
                .find(|(stratum, _, _)| *stratum == entry.file)
                .expect("a writer per stratum");
            target.2.push(labelled);
        }
    }
    close(&writers)?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// An SV record as `SVConcordance` holds it: the converted call, and what the annotator reads.
struct SvConcordanceItem {
    call: gatk_tools::sv_call_record::SvCallRecord,
    record: gatk_tools::sv_concordance::Record,
    variant: htsjdk_vcf::variant::VariantContext,
}

/// `SVConcordance`: every eval SV annotated with the closest truth SV the linkage allows.
///
/// The closest-record choice and every value the annotator computes are
/// [`gatk_tools::sv_concordance`]; the conversion both ways is [`gatk_tools::sv_call_record`]. The
/// runner is the walk:
///
/// * **the master dictionary is required**, a `UserException` of its own, and the two files'
///   dictionaries are then compared with the contig order checked;
/// * **the two files are walked by [`gatk_engine::concordance_walker`]** with every pair counted
///   concordant and no truth filter, and each step adds its truth record before its eval one;
/// * **a truth record is rebuilt with a dictionary**, so its breakpoints are validated where an
///   eval record's are not, and its genotypes keep only their alleles and `CN`;
/// * **the finder is flushed at each new contig**, every eval record taking the closest truth record
///   of the contig the linkage allows, and the output is sorted by `compareCalls`;
/// * **the record written is `getVariantBuilder`'s**, not the input line: the fields the record owns
///   are written again, a passing record comes out unfiltered, and a null annotation is absent.
///
/// Nothing inside `apply` is wrapped: `AbstractConcordanceWalker.traverse` calls it bare, so a record
/// the conversion refuses reaches the user as the conversion's own exception.
pub fn sv_concordance(parser: &Parser) -> Outcome {
    use gatk_tools::sv_call_record as svr;
    use gatk_tools::sv_concordance as conc;
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    use htsjdk_vcf::variant::Value;

    let _ = resolve_read_filters(parser, "SVConcordance")?;
    let truth_path = argument(parser, "truth").ok_or_else(|| {
        Thrown::command_line("Argument truth was missing: Argument 'truth' is required")
    })?;
    let eval_path = argument(parser, "evaluation").ok_or_else(|| {
        Thrown::command_line("Argument evaluation was missing: Argument 'evaluation' is required")
    })?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    let read_vcf = |path: &str| -> Result<(String, htsjdk_vcf::reader::VcfFile), Thrown> {
        let bytes = std::fs::read(path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: path.to_string(),
                }
                .message(),
            )
        })?;
        let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
            htsjdk_bgzf::read::decompress_all(&bytes)
                .ok()
                .and_then(|bytes| String::from_utf8(bytes).ok())
                .ok_or_else(|| {
                    Thrown::non_user(
                        gatk_tools::read_walker_refusal::SAM_FORMAT,
                        format!("{path} is not a block compressed file"),
                    )
                })?
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })?;
        Ok((text, file))
    };
    let (truth_text, truth_file) = read_vcf(&truth_path)?;
    let (eval_text, eval_file) = read_vcf(&eval_path)?;

    // The engine's own validation, as `Concordance` runs it: the master against the features.
    let truth_dictionary = vcf_dictionary(&truth_text);
    let eval_dictionary = vcf_dictionary(&eval_text);
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
            validate_against_master(master, "features", &truth_dictionary.sequences)?;
        }
        if let Some(reference) = &reference {
            gatk_tools::sequence_dictionary::validate(
                "reference",
                &reference.sequences,
                "features",
                &truth_dictionary.sequences,
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

    // `onTraversalStart`.
    let Some(dictionary) = master else {
        return Err(Thrown::user("Reference sequence dictionary required"));
    };
    gatk_tools::sequence_dictionary::validate(
        "eval",
        &eval_dictionary.sequences,
        "truth",
        &truth_dictionary.sequences,
        false,
        true,
    )
    .map_err(|refusal| Thrown {
        failure: Failure::User,
        exception: refusal.java_class(),
        message: Some(refusal.message()),
    })?;
    let sequences: Vec<(String, i32)> = dictionary
        .sequences
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();
    let parameter = |name: &str, default: f64| -> f64 {
        scalar(parser, name)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let window = |name: &str, default: i32| -> i32 {
        scalar(parser, name)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let linkage = gatk_tools::sv_cluster::Linkage {
        depth: gatk_tools::sv_cluster::ClusteringParameters::depth(
            parameter("depth-interval-overlap", 0.8),
            parameter("depth-size-similarity", 0.0),
            window("depth-breakend-window", 10_000_000),
            parameter("depth-sample-overlap", 0.0),
        ),
        mixed: gatk_tools::sv_cluster::ClusteringParameters::mixed(
            parameter("mixed-interval-overlap", 0.8),
            parameter("mixed-size-similarity", 0.0),
            window("mixed-breakend-window", 1000),
            parameter("mixed-sample-overlap", 0.0),
        ),
        pesr: gatk_tools::sv_cluster::ClusteringParameters::pesr(
            parameter("pesr-interval-overlap", 0.5),
            parameter("pesr-size-similarity", 0.0),
            window("pesr-breakend-window", 500),
            parameter("pesr-sample-overlap", 0.0),
        ),
        cluster_del_with_dup: false,
    };
    let common: Vec<String> = eval_file
        .header
        .samples
        .iter()
        .filter(|sample| truth_file.header.samples.contains(sample))
        .cloned()
        .collect();

    // The header: the eval file's, with the tool's lines added to its set.
    let mut header = eval_file.header.clone();
    let compound = |key: &str, id: &str, number: Cardinality, line_type: LineType, text: &str| {
        HeaderLine::Compound {
            key: key.to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: text.to_string(),
            extra: Vec::new(),
        }
    };
    let one = Cardinality::Fixed(1);
    let added = vec![
        compound(
            "FORMAT",
            "CONC_ST",
            Cardinality::Unbounded,
            LineType::String,
            "The genotype concordance contingency state",
        ),
        compound(
            "FORMAT",
            "TRUTH_CN_EQUAL",
            one,
            LineType::Integer,
            "Truth CNV copy state is equal (1=True, 0=False)",
        ),
        compound(
            "INFO",
            "STATUS",
            one,
            LineType::String,
            "Truth status: TP/FP/FN for true positive/false positive/false negative.",
        ),
        compound(
            "INFO",
            "GENOTYPE_CONCORDANCE",
            one,
            LineType::Float,
            "Genotype concordance",
        ),
        compound(
            "INFO",
            "CNV_CONCORDANCE",
            one,
            LineType::Float,
            "CNV copy number concordance",
        ),
        compound(
            "INFO",
            "NON_REF_GENOTYPE_CONCORDANCE",
            one,
            LineType::Float,
            "Non-ref genotype concordance",
        ),
        compound(
            "INFO",
            "HET_PPV",
            one,
            LineType::Float,
            "Heterozygous genotype positive predictive value",
        ),
        compound(
            "INFO",
            "HET_SENSITIVITY",
            one,
            LineType::Float,
            "Heterozygous genotype sensitivity",
        ),
        compound(
            "INFO",
            "HOMVAR_PPV",
            one,
            LineType::Float,
            "Homozygous genotype positive predictive value",
        ),
        compound(
            "INFO",
            "HOMVAR_SENSITIVITY",
            one,
            LineType::Float,
            "Homozygous genotype sensitivity",
        ),
        compound(
            "INFO",
            "VAR_PPV",
            one,
            LineType::Float,
            "Non-ref genotype positive predictive value",
        ),
        compound(
            "INFO",
            "VAR_SENSITIVITY",
            one,
            LineType::Float,
            "Non-ref genotype sensitivity",
        ),
        compound(
            "INFO",
            "VAR_SPECIFICITY",
            one,
            LineType::Float,
            "Non-ref genotype specificity",
        ),
        compound(
            "INFO",
            "TRUTH_VID",
            one,
            LineType::String,
            "Matching truth set variant id",
        ),
        compound(
            "INFO",
            "TRUTH_AC",
            Cardinality::A,
            LineType::Integer,
            "Truth set allele count",
        ),
        compound(
            "INFO",
            "TRUTH_AN",
            Cardinality::A,
            LineType::Integer,
            "Truth set allele number",
        ),
        compound(
            "INFO",
            "TRUTH_AF",
            Cardinality::A,
            LineType::Float,
            "Truth set allele frequency",
        ),
        compound(
            "INFO",
            "AF",
            Cardinality::A,
            LineType::Float,
            "Allele Frequency, for each ALT allele, in the same order as listed",
        ),
        compound(
            "INFO",
            "AC",
            Cardinality::A,
            LineType::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        ),
        compound(
            "INFO",
            "AN",
            one,
            LineType::Integer,
            "Total number of alleles in called genotypes",
        ),
    ];
    for line in added {
        if !header
            .lines
            .iter()
            .any(|existing| same_compound_id(existing, &line))
        {
            header.lines.push(line);
        }
    }

    // The walk.
    let truth: Vec<ConcordanceLocus> = truth_file
        .records
        .iter()
        .enumerate()
        .map(|(index, record)| ConcordanceLocus {
            index,
            contig: record.contig.clone(),
            start: record.start as i32,
            filtered: record.is_filtered(),
        })
        .collect();
    let eval: Vec<ConcordanceLocus> = eval_file
        .records
        .iter()
        .enumerate()
        .map(|(index, record)| ConcordanceLocus {
            index,
            contig: record.contig.clone(),
            start: record.start as i32,
            filtered: record.is_filtered(),
        })
        .collect();
    let walk_dictionary: Vec<String> = truth_dictionary
        .sequences
        .iter()
        .map(|sequence| sequence.name.clone())
        .collect();
    let steps =
        gatk_engine::concordance_walker::concordance(&truth, &eval, &walk_dictionary, |_, _| true);

    let refused = |error: svr::SvRecordError| Thrown {
        failure: if error.is_user() {
            Failure::User
        } else {
            Failure::Other
        },
        exception: error.class(),
        message: Some(error.message()),
    };
    // The concordance side of a record: the call the linkage reads and the genotypes as indices.
    let item = |variant: &htsjdk_vcf::variant::VariantContext,
                call: gatk_tools::sv_call_record::SvCallRecord,
                keep_counts: bool|
     -> SvConcordanceItem {
        let genotypes = variant
            .genotypes
            .iter()
            .map(|genotype| conc::Genotype {
                sample: genotype.sample_name.clone(),
                alleles: genotype
                    .alleles
                    .iter()
                    .map(|allele| {
                        if allele.is_no_call() {
                            None
                        } else {
                            variant
                                .alleles
                                .iter()
                                .position(|known| known == allele)
                                .map(|index| index as i32)
                        }
                    })
                    .collect(),
                copy_number: genotype
                    .extended
                    .iter()
                    .find(|(key, _)| key == "CN")
                    .and_then(|(_, value)| value.format())
                    .and_then(|text| text.parse().ok()),
            })
            .collect();
        let text = |key: &str| {
            call.attributes
                .iter()
                .find(|(name, _)| name == key)
                .and_then(|(_, value)| match value {
                    Value::Missing => None,
                    other => other.format(),
                })
        };
        let allele_counts = if keep_counts {
            match (text("AC"), text("AF"), text("AN")) {
                (Some(ac), Some(af), Some(an)) => Some(conc::AlleleCounts {
                    count: Vec::new(),
                    frequency: Vec::new(),
                    number: 0,
                    verbatim: Some((ac, af, an)),
                }),
                _ => None,
            }
        } else {
            None
        };
        let record = conc::Record {
            call: gatk_tools::sv_cluster::CallRecord {
                id: call.id.clone(),
                sv_type: call.sv_type,
                contig_a: call.contig_a.clone(),
                position_a: call.position_a,
                contig_b: call.contig_b.clone(),
                position_b: call.position_b,
                strand_a: call.strand_a,
                strand_b: call.strand_b,
                length: call.length,
                algorithms: call.algorithms.clone(),
                carriers: Vec::new(),
            },
            genotypes,
            allele_counts,
        };
        SvConcordanceItem {
            call,
            record,
            variant: variant.clone(),
        }
    };

    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    let mut truth_items: Vec<SvConcordanceItem> = Vec::new();
    let mut eval_items: Vec<SvConcordanceItem> = Vec::new();
    let mut current: Option<String> = None;
    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut records = written.to_vec();
        apply_sites_only(parser, &mut header, &mut records);
        let out = htsjdk_vcf::vcf_file::write_vcf(&header, &records)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        write_variant_output(parser, &output, &out)
    };
    let flush = |truth_items: &mut Vec<SvConcordanceItem>,
                 eval_items: &mut Vec<SvConcordanceItem>,
                 written: &mut Vec<htsjdk_vcf::variant::VariantContext>| {
        let truths: Vec<conc::Record> =
            truth_items.iter().map(|item| item.record.clone()).collect();
        let mut annotated: Vec<(
            gatk_tools::sv_call_record::SvCallRecord,
            htsjdk_vcf::variant::VariantContext,
        )> = Vec::new();
        for item in eval_items.iter() {
            let closest = conc::closest(&linkage, &item.record, &truths);
            let annotation = conc::annotate(&item.record, closest, &common);
            annotated.push(sv_concordance_output(item, &annotation, &common));
        }
        annotated.sort_by(|a, b| svr::compare_calls(&a.0, &b.0, &sequences));
        written.extend(annotated.into_iter().map(|(_, variant)| variant));
        truth_items.clear();
        eval_items.clear();
    };
    for step in steps {
        let mut add = |index: usize, is_truth: bool| -> Result<(), Thrown> {
            let variant = if is_truth {
                &truth_file.records[index]
            } else {
                &eval_file.records[index]
            };
            let call = svr::create(variant, &sequences).map_err(refused)?;
            if current.as_deref() != Some(call.contig_a.as_str()) {
                flush(&mut truth_items, &mut eval_items, &mut written);
                current = Some(call.contig_a.clone());
            }
            if is_truth {
                // `minimizeTruthFootprint`: rebuilt WITH the dictionary, genotypes cut to
                // alleles and `CN`.
                svr::validate_coordinates(&call, &sequences).map_err(refused)?;
                let mut stripped = variant.clone();
                let genotypes: Vec<htsjdk_vcf::variant::Genotype> = variant
                    .genotypes
                    .iter()
                    .map(|genotype| {
                        let mut kept = htsjdk_vcf::variant::Genotype::new(
                            &genotype.sample_name,
                            genotype.alleles.clone(),
                        );
                        kept.extended = genotype
                            .extended
                            .iter()
                            .filter(|(key, _)| key == "CN")
                            .cloned()
                            .collect();
                        kept
                    })
                    .collect();
                stripped.genotypes = genotypes.into();
                truth_items.push(item(&stripped, call, true));
            } else {
                eval_items.push(item(variant, call, true));
            }
            Ok(())
        };
        let outcome = step
            .truth
            .map_or(Ok(()), |index| add(index, true))
            .and_then(|()| step.eval.map_or(Ok(()), |index| add(index, false)));
        if let Err(error) = outcome {
            finish(&written)?;
            return Err(error);
        }
    }
    flush(&mut truth_items, &mut eval_items, &mut written);
    finish(&written)?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// One eval record written back: its genotypes and attributes annotated, then `getVariantBuilder`.
fn sv_concordance_output(
    item: &SvConcordanceItem,
    annotation: &gatk_tools::sv_concordance::Annotation,
    common: &[String],
) -> (
    gatk_tools::sv_call_record::SvCallRecord,
    htsjdk_vcf::variant::VariantContext,
) {
    use htsjdk_vcf::variant::Value;
    let is_cnv = item.call.sv_type == gatk_tools::sv_stratify::SvType::Cnv;
    let genotypes: Vec<htsjdk_vcf::variant::Genotype> = item
        .variant
        .genotypes
        .iter()
        .map(|genotype| {
            let mut genotype = genotype.clone();
            if common.contains(&genotype.sample_name) {
                if is_cnv {
                    let value = annotation
                        .truth_copy_number_equal
                        .iter()
                        .find(|(sample, _)| *sample == genotype.sample_name)
                        .and_then(|(_, equal)| *equal)
                        .map(|equal| Value::Int(if equal { 1 } else { 0 }))
                        .unwrap_or(Value::Missing);
                    genotype
                        .extended
                        .push(("TRUTH_CN_EQUAL".to_string(), value));
                } else if let Some(Some(state)) = annotation
                    .contingency
                    .iter()
                    .find(|(sample, _)| *sample == genotype.sample_name)
                    .map(|(_, state)| state.clone())
                {
                    genotype
                        .extended
                        .push(("CONC_ST".to_string(), Value::Str(state)));
                }
            }
            genotype
        })
        .collect();

    let mut call = item.call.clone();
    let mut put = |key: &str, value: Option<Value>| {
        call.attributes.retain(|(name, _)| name != key);
        if let Some(value) = value {
            call.attributes.push((key.to_string(), value));
        }
    };
    let double = |value: f64| (!value.is_nan()).then_some(Value::Double(value));
    put(
        "TRUTH_VID",
        annotation.truth_variant_id.clone().map(Value::Str),
    );
    put("STATUS", Some(Value::Str(annotation.status.to_string())));
    if is_cnv {
        put(
            "CNV_CONCORDANCE",
            annotation.copy_number_concordance.map(Value::Double),
        );
    } else if let Some(metrics) = &annotation.metrics {
        put("GENOTYPE_CONCORDANCE", double(metrics.genotype_concordance));
        put(
            "NON_REF_GENOTYPE_CONCORDANCE",
            double(metrics.non_ref_genotype_concordance),
        );
        put("HET_PPV", double(metrics.het_ppv));
        put("HET_SENSITIVITY", double(metrics.het_sensitivity));
        put("HOMVAR_PPV", double(metrics.homvar_ppv));
        put("HOMVAR_SENSITIVITY", double(metrics.homvar_sensitivity));
        put("VAR_PPV", double(metrics.var_ppv));
        put("VAR_SENSITIVITY", double(metrics.var_sensitivity));
        put("VAR_SPECIFICITY", double(metrics.var_specificity));
    }
    let counts = |counts: &gatk_tools::sv_concordance::AlleleCounts| -> [Value; 3] {
        match &counts.verbatim {
            Some((ac, af, an)) => [
                Value::Str(ac.clone()),
                Value::Str(af.clone()),
                Value::Str(an.clone()),
            ],
            None => [
                Value::List(counts.count.iter().map(|c| Value::Int(*c as i64)).collect()),
                Value::List(counts.frequency.iter().map(|f| Value::Double(*f)).collect()),
                Value::Int(counts.number as i64),
            ],
        }
    };
    if !is_cnv {
        // The eval record's own counts are written only when it carried none.
        if let Some(own) = &annotation.allele_counts {
            if own.verbatim.is_none() {
                let [ac, af, an] = counts(own);
                put("AC", Some(ac));
                put("AF", Some(af));
                put("AN", Some(an));
            }
        }
        match &annotation.truth_allele_counts {
            None => {
                put("TRUTH_AC", None);
                put("TRUTH_AF", None);
                put("TRUTH_AN", None);
            }
            Some(theirs) => {
                let [ac, af, an] = counts(theirs);
                put("TRUTH_AC", Some(ac));
                put("TRUTH_AF", Some(af));
                put("TRUTH_AN", Some(an));
            }
        }
    }
    let filters = item.variant.filters.clone().unwrap_or_default();
    let variant = gatk_tools::sv_call_record::to_variant(
        &call,
        item.variant.alleles.clone(),
        genotypes,
        &filters,
    );
    (call, variant)
}

/// Two compound header lines under the same key and ID, which is how `VCFHeader.addMetaDataLine`
/// keys them: it keeps the line already there, so an eval file declaring its own `AC` keeps its
/// description where a constructor over a `HashSet` of lines would have kept both. Measured on
/// `SVConcordance`'s array.
fn same_compound_id(
    existing: &htsjdk_vcf::header::HeaderLine,
    line: &htsjdk_vcf::header::HeaderLine,
) -> bool {
    use htsjdk_vcf::header::HeaderLine;
    match (existing, line) {
        (
            HeaderLine::Compound { key, id, .. },
            HeaderLine::Compound {
                key: other_key,
                id: other_id,
                ..
            },
        ) => key == other_key && id == other_id,
        _ => existing.render() == line.render(),
    }
}

/// `ConvertHeaderlessHadoopBamShardToBam.doWork`: a donor's header, the shard's own bytes, and a
/// terminator.
///
/// A `CommandLineProgram` rather than a `GATKTool`, so nothing here is the engine's: the donor is
/// opened by `SamReaderFactory` directly, and its refusals are htsjdk's rather than GATK's. The
/// three steps fail in three ways, and the array measures each:
///
///   - the donor that cannot be opened is a `RuntimeIOException` naming the path AS GIVEN, and the
///     handler prints its class because it is not a `UserException`;
///   - the output and the shard both fail inside one `try`, whose `IOException` becomes
///     `Error writing to <the output's absolute path>` whichever of the two it was;
///   - and a shard that fails leaves the header block behind, because the stream was opened, and
///     written to, before `copyFile` asked for the shard.
///
/// The donor needs no header to speak of: a file that is not a BAM is read as SAM text, and a text
/// with no `@` lines at all is an empty header whose whole text is `@HD VN:1.6`.
pub fn convert_headerless_hadoop_bam_shard_to_bam(parser: &Parser) -> Outcome {
    use gatk_tools::convert_headerless_shard::{header_block, EMPTY_GZIP_BLOCK};
    use std::io::Write;

    let required = |name: &str| {
        argument(parser, name).ok_or_else(|| {
            Thrown::command_line(format!(
                "Argument {name} was missing: Argument '{name}' is required"
            ))
        })
    };
    let shard = required("bam-shard")?;
    let donor = required("bam-with-header")?;
    let output = required("output")?;

    let header = donor_header(&donor)?;

    let writing = || Thrown::user(format!("Error writing to {}", java_absolute_path(&output)));
    let mut out = std::fs::File::create(&output).map_err(|_| writing())?;
    let (level, deflater) = output_compression(parser);
    let block = header_block(&header, level, deflater)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, error.to_string()))?;
    out.write_all(&block).map_err(|_| writing())?;
    // `FileUtils.copyFile(File, OutputStream)` opens the shard only now, so a missing one leaves
    // the header block in the file.
    let bytes = std::fs::read(&shard).map_err(|_| writing())?;
    out.write_all(&bytes).map_err(|_| writing())?;
    out.write_all(&EMPTY_GZIP_BLOCK).map_err(|_| writing())?;
    // `doWork` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `SamReaderFactory.makeDefault().validationStringency(SILENT).open(file).getFileHeader()`.
///
/// Measured on the reference rather than read from htsjdk: a gzip or BGZF stream whose content is
/// not `BAM\1` is read as SAM text, like an uncompressed one, and SAM text's header is its leading
/// `@` lines, so a VCF, an empty file and a bare word all give the empty header.
fn donor_header(donor: &str) -> Result<SamHeader, Thrown> {
    use std::io::Read;

    let unreadable = |reason: &str| {
        Thrown::non_user(
            "htsjdk.samtools.util.RuntimeIOException",
            format!("java.io.FileNotFoundException: {donor} ({reason})"),
        )
    };
    let path = std::path::Path::new(donor);
    if path.is_dir() {
        return Err(unreadable("Is a directory"));
    }
    let bytes = std::fs::read(path).map_err(|_| unreadable("No such file or directory"))?;
    if bytes.starts_with(b"CRAM") {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "A CRAM donor is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let content = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(bytes.as_slice())
            .read_to_end(&mut plain)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{donor}: {error}")))?;
        plain
    } else {
        bytes
    };
    if content.starts_with(&gatk_tools::read_walker_refusal::BAM_MAGIC) {
        return htsjdk_bam::reader::BamReader::new(&content)
            .map(|reader| reader.header.text)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{donor}: {error:?}")));
    }
    let text = String::from_utf8_lossy(&content);
    let leading: String = text
        .lines()
        .take_while(|line| line.starts_with('@'))
        .map(|line| format!("{line}\n"))
        .collect();
    // `SAMTextHeaderCodec.decode` fills a `new SAMFileHeader()`, whose constructor has already set
    // `VN:1.6`: a text with no `@HD` keeps that line, and one with an `@HD` overwrites it in place.
    let parsed = htsjdk_bam::reader::parse_header_text(&leading);
    let mut header = SamHeader::default();
    for (key, value) in parsed.attributes.iter() {
        header.attributes.set(key, value);
    }
    Ok(SamHeader {
        attributes: header.attributes,
        ..parsed
    })
}

/// `SVCluster`: SV records grouped by the canonical linkage and each group written as one record.
///
/// Which records group is [`gatk_tools::sv_cluster`], what a group is written as
/// [`gatk_tools::sv_collapser`], and both conversions [`gatk_tools::sv_call_record`]. The runner is
/// `SVClusterWalker` around them:
///
/// * **the reference is required and supplies the dictionary**, the reference allele of every
///   collapsed record, and the contig lines of the output, which replace the input's;
/// * **the ploidy table fills in every sample a record lacks**, with its ploidy as `ECN`, as `CN` as
///   well for a copy-number record when the header declares `CN`, and as that many reference
///   alleles, or no-calls under `--default-no-call`;
/// * **`--fast-mode` keeps only the carriers' genotypes** of a record that is not a CNV;
/// * **`--flag-field-logic` is declared and never read**: the factory builds the collapser with `OR`;
/// * **the output is sorted by contig and start**, stably, and each record renamed
///   `<prefix><8 hex digits>` in the order it was built when `--variant-prefix` is given.
///
/// `DEFRAGMENT_CNV` is another linkage and another collapser configuration, and is the port's
/// limitation here. A failure inside `apply` is the walker's `GATKException`; a collapse is placed
/// at the end of its contig rather than at the record that completed it, which only a refused
/// collapse can observe.
pub fn sv_cluster(parser: &Parser) -> Outcome {
    let mut walker = SvClusterWalker::start(parser, "SVCluster", Vec::new())?;
    let linkage = sv_cluster_linkage(parser);
    let records = walker.records.clone();
    let mut members: Vec<gatk_tools::sv_collapser::Member> = Vec::new();
    let mut current: Option<String> = None;
    for record in &records {
        let member = match walker.member(record) {
            Ok(member) => member,
            Err(error) => return walker.refuse(error),
        };
        if current.as_deref() != Some(member.call.contig_a.as_str()) {
            if let Err(error) = walker.cluster_and_build(&mut members, &linkage) {
                return walker.refuse(error);
            }
            current = Some(member.call.contig_a.clone());
        }
        members.push(member);
    }
    if let Err(error) = walker.cluster_and_build(&mut members, &linkage) {
        return walker.refuse(error);
    }
    walker.finish()?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `GroupedSVCluster`: each record clustered under the thresholds of the one stratum it falls in.
///
/// The stratification is `SVStratify`'s ([`load_stratification_engine`]) over the REFERENCE's
/// dictionary, and the clustering and writing `SVCluster`'s ([`SvClusterWalker`]); the tool is the
/// wiring, in [`gatk_tools::grouped_sv_cluster`]:
///
/// * **`onTraversalStart` refuses after the writer exists**: no strata, a clustering table whose
///   groups do not match the strata one for one, or a stratum it does not name, each an
///   `IllegalStateException`, leave the header behind;
/// * **a record matching two strata is refused**, inside `apply` and therefore wrapped;
/// * **a record matching none is written straight out**, not clustered, with `MEMBERS` its own ID
///   and `STRAT` `default`;
/// * **a matched record carries `STRAT` into its cluster**, whose representative writes it.
///
/// Each stratum's engine flushes at a new contig of its own, and the engines are flushed at the end
/// in the `HashMap` order the strata already have.
pub fn grouped_sv_cluster(parser: &Parser) -> Outcome {
    use gatk_tools::grouped_sv_cluster as grouped;
    use htsjdk_vcf::variant::Value;

    let strat_line = htsjdk_vcf::header::HeaderLine::Compound {
        key: "INFO".to_string(),
        id: "STRAT".to_string(),
        number: htsjdk_vcf::header::Cardinality::Fixed(1),
        line_type: htsjdk_vcf::header::LineType::String,
        description: "Stratum ID".to_string(),
        extra: Vec::new(),
    };
    let mut walker = SvClusterWalker::start(parser, "GroupedSVCluster", vec![strat_line])?;
    let config = argument(parser, "stratify-config").ok_or_else(|| {
        Thrown::command_line(
            "Argument stratify-config was missing: Argument 'stratify-config' is required",
        )
    })?;
    let clustering = argument(parser, "clustering-config").ok_or_else(|| {
        Thrown::command_line(
            "Argument clustering-config was missing: Argument 'clustering-config' is required",
        )
    })?;
    let thresholds = stratification_thresholds(parser);
    let illegal_state =
        |message: String| Thrown::non_user("java.lang.IllegalStateException", message);

    let engine = match load_stratification_engine(parser, &config, &walker.dictionary) {
        Ok(engine) => engine,
        Err(error) => return walker.refuse(error),
    };
    if engine.strata.is_empty() {
        return walker.refuse(illegal_state(grouped::GroupedError::NoStrata.message()));
    }
    let parameters = match read_clustering_config(&clustering) {
        Ok(parameters) => parameters,
        Err(error) => return walker.refuse(error),
    };
    let engines = grouped::Engines::new(&parameters);
    if let Err(error) = grouped::validate(&engine, &engines) {
        return walker.refuse(illegal_state(error.message()));
    }
    let enable_cnv = flag(parser, "enable-cnv");
    let mut buckets: Vec<(
        String,
        Option<String>,
        Vec<gatk_tools::sv_collapser::Member>,
    )> = engine
        .strata
        .iter()
        .map(|stratum| (stratum.name.clone(), None, Vec::new()))
        .collect();
    let linkage_of = |name: &str| {
        grouped::linkage_for(engines.get(name).expect("a validated group"), enable_cnv)
    };

    let records = walker.records.clone();
    for record in &records {
        let mut member = match walker.member(record) {
            Ok(member) => member,
            Err(error) => return walker.refuse(error),
        };
        let matches = match engine.matches(&member.call.stratify_record(), thresholds) {
            Ok(matches) => matches,
            Err(_) => return walker.refuse(walker.wrapped(record, true)),
        };
        let set = |member: &mut gatk_tools::sv_collapser::Member, key: &str, value: Value| {
            member.call.attributes.retain(|(name, _)| name != key);
            member.call.attributes.push((key.to_string(), value));
        };
        match matches.len() {
            0 => {
                let id = member.call.id.clone();
                set(
                    &mut member,
                    gatk_tools::sv_collapser::CLUSTER_MEMBER_IDS_KEY,
                    Value::List(vec![Value::Str(id)]),
                );
                set(
                    &mut member,
                    "STRAT",
                    Value::List(vec![Value::Str(
                        gatk_tools::sv_stratify::DEFAULT_STRATUM.to_string(),
                    )]),
                );
                let built = walker.build(
                    member.call,
                    member.alleles,
                    member.genotypes,
                    &member.filters,
                );
                if let Err(error) = built {
                    return walker.refuse(error);
                }
            }
            1 => {
                let name = matches[0].name.clone();
                set(
                    &mut member,
                    "STRAT",
                    Value::List(vec![Value::Str(name.clone())]),
                );
                let bucket = buckets
                    .iter_mut()
                    .find(|(stratum, _, _)| *stratum == name)
                    .expect("a bucket per stratum");
                if bucket.1.as_deref() != Some(member.call.contig_a.as_str()) {
                    let flushed = walker.cluster_and_build(&mut bucket.2, &linkage_of(&name));
                    if let Err(error) = flushed {
                        return walker.refuse(error);
                    }
                    bucket.1 = Some(member.call.contig_a.clone());
                }
                bucket.2.push(member);
            }
            _ => return walker.refuse(walker.wrapped(record, true)),
        }
    }
    for (name, _, members) in &mut buckets {
        if let Err(error) = walker.cluster_and_build(members, &linkage_of(name)) {
            return walker.refuse(error);
        }
    }
    walker.finish()?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `StratifiedClusteringTableParser`: the columns checked at the header line, then one row per
/// group with its three fractions and its window.
fn read_clustering_config(
    path: &str,
) -> Result<Vec<gatk_tools::grouped_sv_cluster::StratumParameters>, Thrown> {
    let text = std::fs::read_to_string(path).map_err(|_| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            "IO error while reading config table",
        )
    })?;
    let table_error = |error: gatk_engine::tsv_table::TableError| Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message()),
    };
    let table = gatk_engine::tsv_table::Table::parse(&text, path).map_err(table_error)?;
    let header_line = text
        .lines()
        .position(|line| !line.starts_with('#'))
        .map(|index| index + 1)
        .unwrap_or(0);
    gatk_tools::grouped_sv_cluster::check_columns(&table.columns).map_err(|error| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(format!(
            "Bad input: format error in '{path}' at line {header_line}: {}",
            error.message()
        )),
    })?;
    let mut parameters = Vec::new();
    for (row, line) in table.rows.iter().zip(&table.row_lines) {
        let double = |column: &str| -> Result<f64, Thrown> {
            let value = table.get(row, column).map_err(table_error)?;
            value.parse::<f64>().map_err(|_| Thrown {
                failure: Failure::User,
                exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
                message: Some(format!(
                    "Bad input: format error in '{path}' at line {line}: expected double value for column {column} but found {value}"
                )),
            })
        };
        parameters.push(gatk_tools::grouped_sv_cluster::StratumParameters {
            name: table.get(row, "NAME").map_err(table_error)?.to_string(),
            reciprocal_overlap: double("RECIPROCAL_OVERLAP")?,
            size_similarity: double("SIZE_SIMILARITY")?,
            breakend_window: table
                .get_int(row, "BREAKEND_WINDOW", path, *line)
                .map_err(table_error)?,
            sample_overlap: double("SAMPLE_OVERLAP")?,
        });
    }
    Ok(parameters)
}

/// `SVClusterEngineArgumentsCollection`'s three parameter sets, and `--enable-cnv`.
fn sv_cluster_linkage(parser: &Parser) -> gatk_tools::sv_cluster::Linkage {
    let parameter = |name: &str, default: f64| -> f64 {
        scalar(parser, name)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    let window = |name: &str, default: i32| -> i32 {
        scalar(parser, name)
            .and_then(|value| value.parse().ok())
            .unwrap_or(default)
    };
    gatk_tools::sv_cluster::Linkage {
        depth: gatk_tools::sv_cluster::ClusteringParameters::depth(
            parameter("depth-interval-overlap", 0.8),
            parameter("depth-size-similarity", 0.0),
            window("depth-breakend-window", 10_000_000),
            parameter("depth-sample-overlap", 0.0),
        ),
        mixed: gatk_tools::sv_cluster::ClusteringParameters::mixed(
            parameter("mixed-interval-overlap", 0.8),
            parameter("mixed-size-similarity", 0.0),
            window("mixed-breakend-window", 1000),
            parameter("mixed-sample-overlap", 0.0),
        ),
        pesr: gatk_tools::sv_cluster::ClusteringParameters::pesr(
            parameter("pesr-interval-overlap", 0.5),
            parameter("pesr-size-similarity", 0.0),
            window("pesr-breakend-window", 500),
            parameter("pesr-sample-overlap", 0.0),
        ),
        cluster_del_with_dup: flag(parser, "enable-cnv"),
    }
}

/// `SVClusterWalker`: what `SVCluster` and `GroupedSVCluster` share, from `onTraversalStart` to the
/// sorted write.
struct SvClusterWalker<'p> {
    parser: &'p Parser,
    input: String,
    records: Vec<htsjdk_vcf::variant::VariantContext>,
    output: String,
    algorithm: gatk_tools::sv_cluster::Algorithm,
    breakpoints: gatk_tools::sv_collapser::BreakpointSummary,
    alternates: gatk_tools::sv_collapser::AltAlleleSummary,
    fast_mode: bool,
    omit_members: bool,
    default_no_call: bool,
    prefix: Option<String>,
    numbered: usize,
    reference: gatk_engine::reference::ReferenceFileSource,
    dictionary: SamHeader,
    sequences: Vec<(String, i32)>,
    ploidies: PloidyTable,
    header: htsjdk_vcf::header::VcfHeader,
    has_cn_format: bool,
    built: Vec<htsjdk_vcf::variant::VariantContext>,
}

impl<'p> SvClusterWalker<'p> {
    /// The walker's startup and `onTraversalStart`, with `extra` header lines the subclass adds.
    fn start(
        parser: &'p Parser,
        tool: &str,
        extra: Vec<htsjdk_vcf::header::HeaderLine>,
    ) -> Result<Self, Thrown> {
        use gatk_tools::sv_collapser as collapser;
        use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};

        let _ = resolve_read_filters(parser, tool)?;
        let inputs = arguments(parser, "variant");
        if inputs.len() > 1 {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                "More than one --variant is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
            ));
        }
        let input = inputs.into_iter().next().ok_or_else(|| {
            Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
        })?;
        let VariantWalkerStart {
            input,
            text,
            intervals,
            ..
        } = variant_walker_startup_over(parser, input)?;
        let output = argument(parser, "output").ok_or_else(|| {
            Thrown::command_line("Argument output was missing: Argument 'output' is required")
        })?;
        let ploidy_path = argument(parser, "ploidy-table").ok_or_else(|| {
            Thrown::command_line(
                "Argument ploidy-table was missing: Argument 'ploidy-table' is required",
            )
        })?;
        let algorithm = match scalar(parser, "algorithm").as_deref() {
            Some("MAX_CLIQUE") => gatk_tools::sv_cluster::Algorithm::MaxClique,
            Some("DEFRAGMENT_CNV") => {
                return Err(Thrown::non_user(
                    PORT_LIMITATION,
                    "--algorithm DEFRAGMENT_CNV is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
                ))
            }
            _ => gatk_tools::sv_cluster::Algorithm::SingleLinkage,
        };
        let breakpoints = scalar(parser, "breakpoint-summary-strategy")
            .and_then(|name| collapser::BreakpointSummary::value_of(&name))
            .unwrap_or(collapser::BreakpointSummary::Representative);
        let alternates = match scalar(parser, "alt-allele-summary-strategy").as_deref() {
            Some("MOST_SPECIFIC_SUBTYPE") => collapser::AltAlleleSummary::MostSpecificSubtype,
            _ => collapser::AltAlleleSummary::CommonSubtype,
        };

        // `onTraversalStart`: the reference and its dictionary, then the ploidy table.
        let Some(reference_path) = argument(parser, "reference") else {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                "A cluster walker without --reference is refused by the engine before it starts, which this port does not word yet. This message is the port's own and not GATK's.",
            ));
        };
        let Some(dictionary) = reference_dictionary(parser)? else {
            return Err(Thrown::user("Reference sequence dictionary required"));
        };
        let reference = gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(
            &reference_path,
        ))
        .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let sequences: Vec<(String, i32)> = dictionary
            .sequences
            .iter()
            .map(|sequence| (sequence.name.clone(), sequence.length))
            .collect();
        let ploidies = read_ploidy_table(&ploidy_path)?;

        // `createHeader`: the input's lines, the reference's contigs, and the tool's own lines.
        let file = htsjdk_vcf::reader::read_vcf(&text)
            .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
        let mut header = file.header.clone();
        header.samples.sort();
        header.samples.dedup();
        header
            .lines
            .retain(|line| !matches!(line, HeaderLine::Contig { .. }));
        for (index, (name, length)) in sequences.iter().enumerate() {
            header
                .lines
                .push(HeaderLine::contig(name, i64::from(*length), index as i32));
        }
        let compound =
            |key: &str, id: &str, number: Cardinality, line_type: LineType, text: &str| {
                HeaderLine::Compound {
                    key: key.to_string(),
                    id: id.to_string(),
                    number,
                    line_type,
                    description: text.to_string(),
                    extra: Vec::new(),
                }
            };
        let one = Cardinality::Fixed(1);
        let omit_members = flag(parser, "omit-members");
        let mut added = vec![
            compound(
                "INFO",
                "END",
                one,
                LineType::Integer,
                "Stop position of the interval",
            ),
            compound(
                "INFO",
                "SVLEN",
                Cardinality::Unbounded,
                LineType::Integer,
                "Difference in length between REF and ALT alleles",
            ),
            compound(
                "INFO",
                "SVTYPE",
                one,
                LineType::String,
                "Type of structural variant",
            ),
            compound("INFO", "END2", one, LineType::Integer, "Second position"),
            compound("INFO", "CHR2", one, LineType::String, "Second contig"),
            compound(
                "INFO",
                "STRANDS",
                one,
                LineType::String,
                "First and second strands",
            ),
            compound(
                "INFO",
                "ALGORITHMS",
                Cardinality::Unbounded,
                LineType::String,
                "Source algorithms",
            ),
        ];
        if !omit_members {
            added.push(compound(
                "INFO",
                "MEMBERS",
                Cardinality::Unbounded,
                LineType::String,
                "Cluster variant ids",
            ));
        }
        added.push(compound("FORMAT", "GT", one, LineType::String, "Genotype"));
        added.extend(extra);
        for line in added {
            if !header
                .lines
                .iter()
                .any(|existing| same_compound_id(existing, &line))
            {
                header.lines.push(line);
            }
        }
        let has_cn_format = header.lines.iter().any(|line| {
            matches!(line, HeaderLine::Compound { key, id, .. } if key == "FORMAT" && id == "CN")
        });
        let records = variants_in_traversal(&file.records, intervals.as_deref(), &input)?
            .into_iter()
            .cloned()
            .collect();
        Ok(SvClusterWalker {
            parser,
            input,
            records,
            output,
            algorithm,
            breakpoints,
            alternates,
            fast_mode: flag(parser, "fast-mode"),
            omit_members,
            default_no_call: flag(parser, "default-no-call"),
            prefix: argument(parser, "variant-prefix"),
            numbered: 0,
            reference,
            dictionary,
            sequences,
            ploidies,
            header,
            has_cn_format,
            built: Vec::new(),
        })
    }

    /// The traversal's `GATKException` for a record `apply` failed on. `create` decodes the
    /// genotypes as its very last step, so a record it refused prints them as the file carried
    /// them and a record refused after it prints them decoded.
    fn wrapped(&self, record: &htsjdk_vcf::variant::VariantContext, decoded: bool) -> Thrown {
        let text = if decoded {
            java_variant_context_string_decoded(record, &self.input)
        } else {
            java_variant_context_string(record, &self.input)
        };
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            format!(
                "Exception thrown at {}:{} {}",
                record.contig, record.start, text
            ),
        )
    }

    /// `SVClusterWalker.apply` up to `applyRecord`: the conversion, and `--fast-mode`'s carriers.
    fn member(
        &self,
        record: &htsjdk_vcf::variant::VariantContext,
    ) -> Result<gatk_tools::sv_collapser::Member, Thrown> {
        let call = gatk_tools::sv_call_record::create(record, &self.sequences)
            .map_err(|_| self.wrapped(record, false))?;
        let mut genotypes: Vec<htsjdk_vcf::variant::Genotype> =
            record.genotypes.iter().cloned().collect();
        if self.fast_mode && call.sv_type != gatk_tools::sv_stratify::SvType::Cnv {
            let has_alt = record.alleles[1..]
                .iter()
                .any(|a| !a.is_no_call() && !a.is_reference());
            let mut carriers = Vec::new();
            for genotype in genotypes {
                match gatk_tools::sv_collapser::is_carrier(call.sv_type, has_alt, &genotype) {
                    Ok(true) => carriers.push(genotype),
                    Ok(false) => {}
                    Err(_) => return Err(self.wrapped(record, true)),
                }
            }
            genotypes = carriers;
        }
        Ok(gatk_tools::sv_collapser::Member {
            call,
            alleles: record.alleles.clone(),
            genotypes,
            filters: record.filters.clone().unwrap_or_default(),
        })
    }

    /// One engine's items clustered and collapsed, each group built as it comes out.
    ///
    /// The engine flushes at a new contig; a collapse is placed there rather than at the record
    /// that completed its cluster, which only a refused collapse can observe.
    fn cluster_and_build(
        &mut self,
        members: &mut Vec<gatk_tools::sv_collapser::Member>,
        linkage: &gatk_tools::sv_cluster::Linkage,
    ) -> Result<(), Thrown> {
        use gatk_tools::sv_collapser as collapser;
        let needs_carriers = linkage.depth.sample_overlap > 0.0
            || linkage.mixed.sample_overlap > 0.0
            || linkage.pesr.sample_overlap > 0.0;
        let calls: Vec<gatk_tools::sv_cluster::CallRecord> = members
            .iter()
            .map(|member| {
                let carriers = if needs_carriers {
                    let has_alt = member
                        .alleles
                        .iter()
                        .any(|a| !a.is_no_call() && !a.is_reference());
                    member
                        .genotypes
                        .iter()
                        .filter(|g| {
                            collapser::is_carrier(member.call.sv_type, has_alt, g).unwrap_or(false)
                        })
                        .map(|g| g.sample_name.clone())
                        .collect()
                } else {
                    Vec::new()
                };
                gatk_tools::sv_cluster::CallRecord {
                    id: member.call.id.clone(),
                    sv_type: member.call.sv_type,
                    contig_a: member.call.contig_a.clone(),
                    position_a: member.call.position_a,
                    contig_b: member.call.contig_b.clone(),
                    position_b: member.call.position_b,
                    strand_a: member.call.strand_a,
                    strand_b: member.call.strand_b,
                    length: member.call.length,
                    algorithms: member.call.algorithms.clone(),
                    carriers,
                }
            })
            .collect();
        for cluster in gatk_tools::sv_cluster::cluster_indices(&calls, linkage, self.algorithm) {
            let group: Vec<collapser::Member> = cluster
                .iter()
                .map(|index| members[*index].clone())
                .collect();
            let reference = &mut self.reference;
            let collapsed = collapser::collapse(
                &group,
                self.breakpoints,
                self.alternates,
                collapser::FlagFieldLogic::Or,
                &mut |contig, position| {
                    reference
                        .query(contig, position, position)
                        .ok()
                        .and_then(|bases| bases.first().copied())
                },
            )
            .map_err(|error| Thrown::non_user(error.class(), error.message()))?;
            self.build(
                collapsed.call,
                collapsed.alleles,
                collapsed.genotypes,
                &collapsed.filters,
            )?;
        }
        members.clear();
        Ok(())
    }

    /// `buildVariantContext`: every sample filled in from the ploidy table, the new name, the
    /// dictionary's validation, then the variant.
    fn build(
        &mut self,
        mut call: gatk_tools::sv_call_record::SvCallRecord,
        alleles: Vec<htsjdk_vcf::allele::Allele>,
        mut genotypes: Vec<htsjdk_vcf::variant::Genotype>,
        filters: &[String],
    ) -> Result<(), Thrown> {
        use htsjdk_vcf::variant::Value;
        let is_cnv = matches!(
            call.sv_type,
            gatk_tools::sv_stratify::SvType::Del
                | gatk_tools::sv_stratify::SvType::Dup
                | gatk_tools::sv_stratify::SvType::Cnv
        );
        for sample in &self.header.samples {
            if genotypes.iter().any(|g| g.sample_name == *sample) {
                continue;
            }
            let ploidy = ploidy_of(&self.ploidies, sample, &call.contig_a)?;
            let allele = if self.default_no_call {
                htsjdk_vcf::allele::Allele::no_call()
            } else {
                alleles[0].clone()
            };
            let mut genotype =
                htsjdk_vcf::variant::Genotype::new(sample, vec![allele; ploidy.max(0) as usize]);
            genotype
                .extended
                .push(("ECN".to_string(), Value::Int(i64::from(ploidy))));
            if is_cnv && self.has_cn_format {
                genotype
                    .extended
                    .push(("CN".to_string(), Value::Int(i64::from(ploidy))));
            }
            genotypes.push(genotype);
        }
        if let Some(prefix) = &self.prefix {
            call.id = format!("{prefix}{:08x}", self.numbered);
            self.numbered += 1;
        }
        gatk_tools::sv_call_record::validate_coordinates(&call, &self.sequences)
            .map_err(|error| Thrown::non_user(error.class(), error.message()))?;
        if self.omit_members {
            call.attributes
                .retain(|(key, _)| key != gatk_tools::sv_collapser::CLUSTER_MEMBER_IDS_KEY);
        }
        self.built.push(gatk_tools::sv_call_record::to_variant(
            &call, alleles, genotypes, filters,
        ));
        Ok(())
    }

    /// `closeTool` after a refusal: the writer was opened with its header, and the sorting buffer
    /// is never flushed into it.
    fn refuse(&self, error: Thrown) -> Outcome {
        self.write(&[])?;
        Err(error)
    }

    fn write(&self, records: &[htsjdk_vcf::variant::VariantContext]) -> Result<(), Thrown> {
        let mut header = self.header.clone();
        let mut records = records.to_vec();
        apply_sites_only(self.parser, &mut header, &mut records);
        let out = write_vcf_honouring_lenient(self.parser, &header, &records)?;
        write_variant_output(self.parser, &self.output, &out)
    }

    /// `onTraversalSuccess`: the sorting collection, by contig then start, stably.
    fn finish(&mut self) -> Result<(), Thrown> {
        let sequences = &self.sequences;
        let index_of = |contig: &str| {
            sequences
                .iter()
                .position(|(name, _)| name == contig)
                .unwrap_or(usize::MAX)
        };
        let mut built = std::mem::take(&mut self.built);
        built.sort_by_key(|record| (index_of(&record.contig), record.start));
        self.write(&built)
    }
}

/// Each sample's ploidy per contig, in the table's order.
type PloidyTable = Vec<(String, Vec<(String, i32)>)>;

/// `PloidyTable`: the first column names the sample, every other one a contig.
fn read_ploidy_table(path: &str) -> Result<PloidyTable, Thrown> {
    let text = std::fs::read_to_string(path).map_err(|_| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            "IO error while reading ploidy table",
        )
    })?;
    let table_error = |error: gatk_engine::tsv_table::TableError| Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message()),
    };
    let table = gatk_engine::tsv_table::Table::parse(&text, path).map_err(table_error)?;
    let mut ploidies: PloidyTable = Vec::new();
    for (row, line) in table.rows.iter().zip(&table.row_lines) {
        let mut contigs = Vec::new();
        for column in &table.columns[1..] {
            let value = table
                .get_int(row, column, path, *line)
                .map_err(table_error)?;
            contigs.push((column.clone(), value));
        }
        // `Collectors.toMap` would refuse a repeated sample; the last one is kept here.
        ploidies.retain(|(sample, _)| *sample != row[0]);
        ploidies.push((row[0].clone(), contigs));
    }
    Ok(ploidies)
}

/// `PloidyTable.get`.
fn ploidy_of(ploidies: &PloidyTable, sample: &str, contig: &str) -> Result<i32, Thrown> {
    let illegal = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);
    let Some((_, contigs)) = ploidies.iter().find(|(name, _)| name == sample) else {
        return Err(illegal(format!(
            "Sample {sample} not found in ploidy records"
        )));
    };
    contigs
        .iter()
        .find(|(name, _)| name == contig)
        .map(|(_, ploidy)| *ploidy)
        .ok_or_else(|| {
            illegal(format!(
                "No ploidy entry for sample {sample} at contig {contig}"
            ))
        })
}

/// `VCFWriter` over the tool's own header, with `--lenient` deciding what an undeclared key does.
///
/// `GATKTool.createVCFWriter` adds `ALLOW_MISSING_FIELDS_IN_HEADER` under `--lenient`, and
/// without it the encoder refuses a FILTER, INFO or FORMAT key the header does not declare with an
/// `IllegalStateException` naming the record. Measured on `SVCluster`'s array, whose input filters a
/// record `LowQual` without declaring it: the reference wrote the row with `--lenient` and refused
/// every other one.
fn write_vcf_honouring_lenient(
    parser: &Parser,
    header: &htsjdk_vcf::header::VcfHeader,
    records: &[htsjdk_vcf::variant::VariantContext],
) -> Result<String, Thrown> {
    use htsjdk_vcf::encoder::{EncodeError, MissingFields, VcfEncoder};
    let missing = if flag(parser, "lenient") {
        MissingFields::Allow
    } else {
        MissingFields::Refuse
    };
    let encoder = VcfEncoder::new(header).with_missing_fields(missing);
    let mut out = header.write();
    for record in records {
        encoder
            .encode_into(record, &mut out)
            .map_err(|error| match error {
                EncodeError::MissingFromHeader {
                    key,
                    field,
                    contig,
                    start,
                } => Thrown::non_user(
                    "java.lang.IllegalStateException",
                    format!(
                    "Key {key} found in VariantContext field {field} at {contig}:{start} but this \
                     key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF \
                     headers by default."
                ),
                ),
                other => Thrown::user(format!("{other:?}")),
            })?;
        out.push('\n');
    }
    Ok(out)
}

/// [`java_variant_context_string`] once the genotypes have been decoded, which is
/// `toStringDecodeGenotypes`: `GT=` is `GenotypesContext.toString()`, each genotype in sample-name
/// order as `[sample alleles GQ DP AD PL FT {attributes}]`, the alleles of an unphased genotype
/// sorted reference first, and every field that is absent left out.
fn java_variant_context_string_decoded(
    record: &htsjdk_vcf::variant::VariantContext,
    source: &str,
) -> String {
    use htsjdk_vcf::variant::Value;
    let lazy = java_variant_context_string(record, source);
    let mut genotypes: Vec<&htsjdk_vcf::variant::Genotype> = record.genotypes.iter().collect();
    genotypes.sort_by(|a, b| a.sample_name.cmp(&b.sample_name));
    fn java(value: &Value) -> String {
        match value {
            Value::Missing => "null".to_string(),
            Value::Bool(flag) => flag.to_string(),
            Value::Str(text) => text.clone(),
            Value::List(items) => format!(
                "[{}]",
                items.iter().map(java).collect::<Vec<String>>().join(", ")
            ),
            other => other.format().unwrap_or_default(),
        }
    }
    let rendered: Vec<String> = genotypes
        .iter()
        .map(|genotype| {
            let allele = |a: &htsjdk_vcf::allele::Allele| {
                if a.is_no_call() {
                    ".".to_string()
                } else if a.is_reference() {
                    format!("{}*", a.display_string())
                } else {
                    a.display_string()
                }
            };
            let calls = if genotype.alleles.is_empty() {
                "NA".to_string()
            } else {
                let mut alleles: Vec<&htsjdk_vcf::allele::Allele> =
                    genotype.alleles.iter().collect();
                if !genotype.phased {
                    alleles.sort_by(|a, b| {
                        b.is_reference()
                            .cmp(&a.is_reference())
                            .then(a.display_string().cmp(&b.display_string()))
                    });
                }
                alleles
                    .into_iter()
                    .map(allele)
                    .collect::<Vec<String>>()
                    .join(if genotype.phased { "|" } else { "/" })
            };
            let int = |name: &str, value: Option<i32>| {
                value.map(|v| format!(" {name} {v}")).unwrap_or_default()
            };
            let ints = |name: &str, values: &Option<Vec<i32>>| {
                values
                    .as_ref()
                    .map(|v| {
                        format!(
                            " {name} {}",
                            v.iter()
                                .map(i32::to_string)
                                .collect::<Vec<String>>()
                                .join(",")
                        )
                    })
                    .unwrap_or_default()
            };
            let mut extended: Vec<(String, String)> = genotype
                .extended
                .iter()
                .map(|(key, value)| (key.clone(), java(value)))
                .collect();
            extended.sort();
            let extended = if extended.is_empty() {
                String::new()
            } else {
                format!(
                    " {{{}}}",
                    extended
                        .iter()
                        .map(|(key, value)| format!("{key}={value}"))
                        .collect::<Vec<String>>()
                        .join(", ")
                )
            };
            format!(
                "[{} {}{}{}{}{}{}{}]",
                genotype.sample_name,
                calls,
                int("GQ", genotype.gq),
                int("DP", genotype.dp),
                ints("AD", &genotype.ad),
                ints("PL", &genotype.pl),
                genotype
                    .filters
                    .as_ref()
                    .map(|f| format!(" FT {f}"))
                    .unwrap_or_default(),
                extended
            )
        })
        .collect();
    let decoded = format!("[{}]", rendered.join(","));
    // Everything but the genotypes is the lazy form's.
    let (head, tail) = lazy
        .split_once(" GT=")
        .expect("a VariantContext string has a GT field");
    let filters = tail.rsplit_once(" filters=").map(|(_, f)| f).unwrap_or("");
    format!("{head} GT={decoded} filters={filters}")
}

/// `SVStratificationEngineArgumentsCollection`'s three thresholds, at their defaults when absent.
fn stratification_thresholds(parser: &Parser) -> gatk_tools::sv_stratify::Thresholds {
    gatk_tools::sv_stratify::Thresholds {
        overlap_fraction: scalar(parser, "stratify-overlap-fraction")
            .and_then(|value| value.parse().ok())
            .unwrap_or(0.0),
        num_breakpoint_overlaps: scalar(parser, "stratify-num-breakpoint-overlaps")
            .and_then(|value| value.parse().ok())
            .unwrap_or(1),
        num_breakpoint_overlaps_interchrom: scalar(
            parser,
            "stratify-num-breakpoint-overlaps-interchromosomal",
        )
        .and_then(|value| value.parse().ok())
        .unwrap_or(1),
    }
}

/// `SVStratify.loadStratificationConfig`, which `GroupedSVCluster` calls too: the tracks, each
/// loaded through the interval arguments and refused for a name given twice after its file was read,
/// then the table, checked for its columns at its header line and parsed row by row.
///
/// `dictionary` is the tool's: the master one for `SVStratify`, the reference's for
/// `GroupedSVCluster`.
fn load_stratification_engine(
    parser: &Parser,
    config: &str,
    dictionary: &SamHeader,
) -> Result<gatk_tools::sv_stratify::Engine, Thrown> {
    use gatk_tools::sv_stratify as stratify;
    let illegal = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);
    let gatk = |message: String| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            message,
        )
    };
    let bad_input = |message: String| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(format!("Bad input: {message}")),
    };
    // `loadStratificationConfig`: the tracks, then the table.
    let track_names = arguments(parser, "track-name");
    let track_files = arguments(parser, "track-intervals");
    if track_names.len() != track_files.len() {
        return Err(illegal(
            stratify::StratifyError::TrackCountMismatch.message(),
        ));
    }
    let mut names: Vec<String> = Vec::new();
    let mut loaded: Vec<Vec<stratify::Interval>> = Vec::new();
    for (name, path) in track_names.iter().zip(&track_files) {
        let parameters = gatk_engine::interval_arguments::traversal_parameters(
            std::slice::from_ref(path),
            &[],
            dictionary,
            SetRule::Union,
            MergingRule::All,
            0,
            0,
        )
        .map_err(|error| Thrown {
            failure: Failure::User,
            exception: error.java_class(),
            message: Some(error.message()),
        })?;
        if names.contains(name) {
            return Err(bad_input(
                stratify::StratifyError::DuplicateTrack { name: name.clone() }.message(),
            ));
        }
        names.push(name.clone());
        loaded.push(
            parameters
                .intervals
                .into_iter()
                .map(|interval| stratify::Interval {
                    contig: interval.contig,
                    start: interval.start,
                    end: interval.end,
                })
                .collect(),
        );
    }
    let tracks =
        stratify::Tracks::new(&names, &loaded).map_err(|error| bad_input(error.message()))?;

    let table_text = std::fs::read_to_string(config)
        .map_err(|_| gatk("IO error while reading config table".to_string()))?;
    let table =
        gatk_engine::tsv_table::Table::parse(&table_text, config).map_err(|error| Thrown {
            failure: Failure::User,
            exception: error.java_class(),
            message: Some(error.message()),
        })?;
    // The header line's number, which a format error names: the first line that is no comment.
    let header_line = table_text
        .lines()
        .position(|line| !line.starts_with('#'))
        .map(|index| index + 1)
        .unwrap_or(0);
    if table.columns.is_empty() {
        return Err(bad_input(format!(
            "format error in '{config}' at line {}: premature end of table: header line not found",
            table_text.lines().count()
        )));
    }
    stratify::check_columns(&table.columns).map_err(|error| {
        bad_input(format!(
            "format error in '{config}' at line {header_line}: {}",
            error.message()
        ))
    })?;
    let mut strata: Vec<stratify::Stratum> = Vec::new();
    for row in &table.rows {
        let cell = |column: &str| -> String {
            table
                .get(row, column)
                .map(str::to_string)
                .unwrap_or_default()
        };
        let type_name = cell("SVTYPE");
        let sv_type = stratify::SvType::parse(&type_name).ok_or_else(|| {
            illegal(format!(
                "No enum constant org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.StructuralVariantAnnotationType.{type_name}"
            ))
        })?;
        let name = cell("NAME");
        let bound = |column: &str| -> Result<Option<i32>, Thrown> {
            let value = cell(column);
            if stratify::NULL_TABLE_VALUES.contains(&value.as_str()) {
                return Ok(None);
            }
            value.parse::<i32>().map(Some).map_err(|_| {
                Thrown::non_user(
                    "java.lang.NumberFormatException",
                    format!("For input string: \"{value}\""),
                )
            })
        };
        let min_size = bound("MIN_SIZE")?;
        let max_size = bound("MAX_SIZE")?;
        let track_list = stratify::parse_track_string(&cell("TRACKS"), &names)
            .map_err(|error| gatk(error.message()))?;
        let stratum = stratify::Stratum::new(&name, sv_type, min_size, max_size, track_list)
            .map_err(|error| illegal(error.message()))?;
        if strata.iter().any(|existing| existing.name == stratum.name) {
            return Err(gatk(format!("Encountered duplicate name {}", stratum.name)));
        }
        strata.push(stratum);
    }
    stratify::Engine::new(strata, tracks).map_err(|error| bad_input(error.message()))
}

/// `SVAnnotate`: every SV annotated with what it is predicted to do to the protein-coding genes
/// and non-coding elements it reaches.
///
/// The rules are [`gatk_tools::sv_annotate`]; the runner reads the three inputs and writes:
///
/// * **the GTF and the BED are read in `onTraversalStart`, before the writer**, and only their
///   features on a contig the VCF's own dictionary names are kept; every promoter is built there
///   too, so a window the interval refuses leaves no file;
/// * **the SV type comes from the ALT allele**, not `SVTYPE`: a breakend allele is a BND (a CPX when
///   `CPX_INTERVALS` is present), a symbolic one its first symbol, anything else refused;
/// * **nothing inside `apply` is wrapped**, a `VariantWalker` calling it bare;
/// * **the consequences are added to the record's own attributes**, each a sorted list of names,
///   and `PREDICTED_INTERGENIC` is written whenever a GTF was given, as a flag that a `false` drops.
pub fn sv_annotate(parser: &Parser) -> Outcome {
    use gatk_tools::sv_annotate as annotate;
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    use htsjdk_vcf::variant::Value;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "SVAnnotate")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let promoter_window = scalar(parser, "promoter-window-length")
        .and_then(|value| value.parse().ok())
        .unwrap_or(1000);
    let max_breakend_len = scalar(parser, "max-breakend-as-cnv-length")
        .and_then(|value| value.parse().ok())
        .unwrap_or(-1);
    let illegal = |message: String| Thrown::non_user("java.lang.IllegalArgumentException", message);
    let refused = |error: annotate::AnnotateError| match error {
        annotate::AnnotateError::CpxWithoutIntervals
        | annotate::AnnotateError::CpxWithoutType
        | annotate::AnnotateError::CtxWithoutContig2 => Thrown::user(error.message()),
        annotate::AnnotateError::NumberFormat { .. } => {
            Thrown::non_user("java.lang.NumberFormatException", error.message())
        }
        annotate::AnnotateError::InvalidInterval { .. } => illegal(error.message()),
    };

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let contigs: Vec<String> = sequence_dictionary_of(&file.header)
        .into_iter()
        .map(|(name, _)| name)
        .collect();

    // `onTraversalStart`: the GTF, then the BED, then the writer.
    let read = |path: &str| {
        std::fs::read_to_string(path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: java_absolute_path(path),
                }
                .message(),
            )
        })
    };
    let gtf = argument(parser, "protein-coding-gtf");
    let transcripts = match &gtf {
        Some(path) => annotate::transcripts_from_gtf(&read(path)?, &contigs),
        None => Vec::new(),
    };
    for transcript in &transcripts {
        annotate::promoter_interval(transcript, promoter_window).map_err(refused)?;
    }
    let bed = argument(parser, "non-coding-bed");
    let non_coding = match &bed {
        Some(path) => annotate::non_coding_from_bed(&read(path)?, &contigs).map_err(refused)?,
        None => Vec::new(),
    };

    let mut header = file.header.clone();
    let list = |id: &str, text: &str| HeaderLine::Compound {
        key: "INFO".to_string(),
        id: id.to_string(),
        number: Cardinality::Unbounded,
        line_type: LineType::String,
        description: text.to_string(),
        extra: Vec::new(),
    };
    let mut added = vec![
        list(annotate::LOF, "Gene(s) on which the SV is predicted to have a loss-of-function effect."),
        list(annotate::INT_EXON_DUP, "Gene(s) on which the SV is predicted to result in intragenic exonic duplication without breaking any coding sequences."),
        list(annotate::COPY_GAIN, "Gene(s) on which the SV is predicted to have a copy-gain effect."),
        list(annotate::TSS_DUP, "Gene(s) for which the SV is predicted to duplicate the transcription start site."),
        list(annotate::DUP_PARTIAL, "Gene(s) which are partially overlapped by an SV's duplication, but the transcription start site is not duplicated."),
        list(annotate::INTRONIC, "Gene(s) where the SV was found to lie entirely within an intron."),
        list(annotate::PARTIAL_EXON_DUP, "Gene(s) where the duplication SV has one breakpoint in the coding sequence."),
        list(annotate::INV_SPAN, "Gene(s) which are entirely spanned by an SV's inversion."),
        list(annotate::UTR, "Gene(s) for which the SV is predicted to disrupt a UTR."),
        list(annotate::MSV_EXON_OVERLAP, "Gene(s) on which the multiallelic SV would be predicted to have a LOF, INTRAGENIC_EXON_DUP, COPY_GAIN, DUP_PARTIAL, TSS_DUP, or PARTIAL_EXON_DUP annotation if the SV were biallelic."),
        list(annotate::PROMOTER, "Gene(s) for which the SV is predicted to overlap the promoter region."),
        list(annotate::BREAKEND_EXON, "Gene(s) for which the SV breakend is predicted to fall in an exon."),
        HeaderLine::Compound {
            key: "INFO".to_string(),
            id: annotate::INTERGENIC.to_string(),
            number: Cardinality::Fixed(0),
            line_type: LineType::Flag,
            description: "SV does not overlap any protein-coding genes.".to_string(),
            extra: Vec::new(),
        },
        list(annotate::NONCODING_SPAN, "Class(es) of noncoding elements spanned by SV."),
        list(annotate::NONCODING_BREAKPOINT, "Class(es) of noncoding elements disrupted by SV breakpoint."),
        list(annotate::NEAREST_TSS, "Nearest transcription start site to an intergenic variant."),
        list(annotate::PARTIAL_DISPERSED_DUP, "Gene(s) overlapped partially by the duplicated interval involved in a dispersed duplication event in a complex SV."),
    ];
    added.retain(|line| {
        !header
            .lines
            .iter()
            .any(|existing| same_compound_id(existing, line))
    });
    header.lines.extend(added);
    header
        .lines
        .extend(default_tool_vcf_header_lines(parser, "SVAnnotate"));

    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut records = written.to_vec();
        apply_sites_only(parser, &mut header, &mut records);
        let out = write_vcf_honouring_lenient(parser, &header, &records)?;
        write_variant_output(parser, &output, &out)
    };
    let attribute = |record: &htsjdk_vcf::variant::VariantContext, key: &str| {
        record
            .attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.clone())
    };
    let text_of = |value: &Value| match value {
        Value::Str(text) => text.clone(),
        other => other.format().unwrap_or_default(),
    };

    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for record in kept {
        let outcome = (|| -> Result<htsjdk_vcf::variant::VariantContext, Thrown> {
            // `getSVType`.
            let alternates = record.alternate_alleles();
            if alternates.len() > 1 {
                let names: Vec<String> = alternates.iter().map(|a| a.display_string()).collect();
                return Err(illegal(format!(
                    "Expected single ALT allele, found multiple: [{}]",
                    names.join(", ")
                )));
            }
            let alt = alternates[0].display_string();
            let is_breakpoint = alt.len() > 1 && (alt.contains('[') || alt.contains(']'));
            let sv_type = if is_breakpoint {
                if attribute(record, "CPX_INTERVALS").is_some() {
                    annotate::SvType::Cpx
                } else {
                    annotate::SvType::Bnd
                }
            } else if alternates[0].is_symbolic() {
                let symbol = alt.replace(['<', '>'], "");
                let first = symbol.split(':').next().unwrap_or_default().to_string();
                annotate::sv_type_named(&first).ok_or_else(|| {
                    illegal(format!(
                        "No enum constant org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.StructuralVariantAnnotationType.{first}"
                    ))
                })?
            } else {
                return Err(illegal(format!(
                    "Unexpected ALT allele: {alt}. Expected breakpoint or symbolic ALT allele representing a structural variant record."
                )));
            };
            // `getComplexSubtype`.
            let complex_type = match attribute(record, "CPX_TYPE") {
                None => None,
                Some(value) => {
                    let name = text_of(&value);
                    if !gatk_tools::sv_call_record::COMPLEX_SUBTYPES.contains(&name.as_str()) {
                        let error =
                            gatk_tools::sv_call_record::SvRecordError::InvalidComplexSubtype {
                                subtype: name,
                            };
                        return Err(illegal(error.message()));
                    }
                    annotate::complex_subtype(&name)
                }
            };
            let int = |key: &str, default: i32| -> Result<i32, Thrown> {
                match attribute(record, key) {
                    None | Some(Value::Missing) => Ok(default),
                    Some(value) => {
                        let text = text_of(&value);
                        text.parse().map_err(|_| {
                            Thrown::non_user(
                                "java.lang.NumberFormatException",
                                format!("For input string: \"{text}\""),
                            )
                        })
                    }
                }
            };
            let variant = annotate::Variant {
                id: record.id.clone(),
                contig: record.contig.clone(),
                position: record.start as i32,
                end: record.stop as i32,
                sv_type,
                sv_length: int("SVLEN", 0)?,
                contig2: attribute(record, "CHR2").map(|value| text_of(&value)),
                end2: match attribute(record, "END2") {
                    None => None,
                    Some(_) => Some(int("END2", record.start as i32)?),
                },
                strands: attribute(record, "STRANDS").map(|value| text_of(&value)),
                complex_type,
                complex_intervals: match attribute(record, "CPX_INTERVALS") {
                    None => Vec::new(),
                    Some(Value::List(items)) => items.iter().map(text_of).collect(),
                    Some(value) => vec![text_of(&value)],
                },
            };
            let annotation = annotate::annotate_structural_variant(
                &variant,
                &transcripts,
                &non_coding,
                gtf.is_some(),
                bed.is_some(),
                promoter_window,
                max_breakend_len,
            )
            .map_err(refused)?;
            let mut out = record.clone();
            for (consequence, names) in annotation.consequences {
                out.attributes.retain(|(key, _)| *key != consequence);
                out.attributes.push((
                    consequence,
                    Value::List(names.into_iter().map(Value::Str).collect()),
                ));
            }
            if let Some(intergenic) = annotation.intergenic {
                out.attributes
                    .retain(|(key, _)| key != annotate::INTERGENIC);
                out.attributes
                    .push((annotate::INTERGENIC.to_string(), Value::Bool(intergenic)));
            }
            Ok(out)
        })();
        match outcome {
            Ok(annotated) => written.push(annotated),
            Err(error) => {
                finish(&written)?;
                return Err(error);
            }
        }
    }
    finish(&written)?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `ReferenceBlockConcordance`: two GVCFs' reference blocks as three histograms.
///
/// The accumulation and the metrics file are [`gatk_tools::reference_block_concordance`]; the
/// runner reads the two files, keeps only the records whose first genotype is hom-ref on either
/// side (both walker filters), walks them with [`gatk_engine::concordance_walker`], and writes
/// the three files `onTraversalSuccess` writes, each headed by `getMetricsFile`'s two lines: the
/// command line, and the time the run started, which the covering array does not compare.
///
/// The tool returns `SUCCESS`, which `handleResult` prints.
pub fn reference_block_concordance(parser: &Parser) -> Outcome {
    use gatk_tools::reference_block_concordance as rbc;

    let _ = resolve_read_filters(parser, "ReferenceBlockConcordance")?;
    let required = |name: &str| {
        argument(parser, name).ok_or_else(|| {
            Thrown::command_line(format!(
                "Argument {name} was missing: Argument '{name}' is required"
            ))
        })
    };
    let truth_path = required("truth")?;
    let eval_path = required("evaluation")?;
    let truth_histogram = required("truth-block-histogram")?;
    let eval_histogram = required("eval-block-histogram")?;
    let concordance_histogram = required("confidence-concordance-histogram")?;

    let read_vcf = |path: &str| -> Result<htsjdk_vcf::reader::VcfFile, Thrown> {
        let bytes = std::fs::read(path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: path.to_string(),
                }
                .message(),
            )
        })?;
        let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
            htsjdk_bgzf::read::decompress_all(&bytes)
                .ok()
                .and_then(|bytes| String::from_utf8(bytes).ok())
                .ok_or_else(|| {
                    Thrown::non_user(
                        gatk_tools::read_walker_refusal::SAM_FORMAT,
                        format!("{path} is not a block compressed file"),
                    )
                })?
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })
    };
    let truth_file = read_vcf(&truth_path)?;
    let eval_file = read_vcf(&eval_path)?;
    let dictionary: Vec<String> = sequence_dictionary_of(&truth_file.header)
        .into_iter()
        .map(|(name, _)| name)
        .collect();

    // Both walker filters: the first genotype hom-ref.
    // The concordance walker's own `FeatureDataSource`s name their records `Unknown`, where a
    // `MultiVariantDataSource` names them by path.
    let blocks = |file: &htsjdk_vcf::reader::VcfFile| -> Vec<(rbc::Block, bool)> {
        file.records
            .iter()
            .filter_map(|record| {
                let first = record.genotypes.first()?;
                let hom_ref =
                    !first.alleles.is_empty() && first.alleles.iter().all(|a| a.is_reference());
                hom_ref.then(|| {
                    (
                        rbc::Block {
                            contig: record.contig.clone(),
                            start: record.start as i32,
                            end: record.stop as i32,
                            gq: first.gq.unwrap_or(-1),
                            is_hom_ref: true,
                            genotypes: record.genotypes.len(),
                            rendered: java_variant_context_string_decoded(record, "Unknown"),
                        },
                        record.is_filtered(),
                    )
                })
            })
            .collect()
    };
    let truth = blocks(&truth_file);
    let eval = blocks(&eval_file);
    let loci = |blocks: &[(rbc::Block, bool)]| -> Vec<ConcordanceLocus> {
        blocks
            .iter()
            .enumerate()
            .map(|(index, (block, filtered))| ConcordanceLocus {
                index,
                contig: block.contig.clone(),
                start: block.start,
                filtered: *filtered,
            })
            .collect()
    };
    let steps: Vec<(Option<usize>, Option<usize>)> = gatk_engine::concordance_walker::concordance(
        &loci(&truth),
        &loci(&eval),
        &dictionary,
        |_, _| true,
    )
    .into_iter()
    .map(|step| (step.truth, step.eval))
    .collect();
    let truth_blocks: Vec<rbc::Block> = truth.into_iter().map(|(block, _)| block).collect();
    let eval_blocks: Vec<rbc::Block> = eval.into_iter().map(|(block, _)| block).collect();
    let histograms = rbc::accumulate(&truth_blocks, &eval_blocks, &steps)
        .map_err(|error| Thrown::non_user("java.lang.IllegalStateException", error.message()))?;

    // `getMetricsFile`: the command line, then the time the run started.
    let headers = vec![
        crate::command_line::expanded("ReferenceBlockConcordance", parser),
        format!("Started on: {}", java_display_now()),
    ];
    write_file(
        &truth_histogram,
        rbc::write_histogram(&histograms.truth_blocks, &headers).as_bytes(),
    )?;
    write_file(
        &eval_histogram,
        rbc::write_histogram(&histograms.eval_blocks, &headers).as_bytes(),
    )?;
    write_file(
        &concordance_histogram,
        rbc::write_histogram(&histograms.confidence_concordance, &headers).as_bytes(),
    )?;
    Ok(Some("SUCCESS".to_string()))
}

/// `Utils.getDateTimeForDisplay(ZonedDateTime.now())`, in UTC: the one value a metrics header
/// carries that no comparison reads.
fn java_display_now() -> String {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_secs() as i64)
        .unwrap_or(0);
    let days = seconds.div_euclid(86_400);
    let of_day = seconds.rem_euclid(86_400);
    // Civil date from days since the epoch (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
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
    let (hour, minute, second) = (of_day / 3600, (of_day / 60) % 60, of_day % 60);
    let (twelve, meridiem) = match hour {
        0 => (12, "AM"),
        1..=11 => (hour, "AM"),
        12 => (12, "PM"),
        _ => (hour - 12, "PM"),
    };
    format!(
        "{} {day}, {year} at {twelve}:{minute:02}:{second:02} {meridiem} UTC",
        MONTHS[(month - 1) as usize]
    )
}

/// `CombineSegmentBreakpoints`: two segment files cut at every breakpoint either carries, each
/// piece annotated from both.
///
/// The cutting and the annotation are [`gatk_tools::combine_segment_breakpoints`]; the runner reads
/// the two collections restricted to `--columns-of-interest`, builds the output header the way
/// `SamFileHeaderMerger` does, and writes the combined collection:
///
/// * **the dictionaries are compared with their order checked**: common contigs of different
///   lengths, or a different order, are refused as bad input;
/// * **the merged header is `@HD VN:1.6 GO:none SO:coordinate`** over the merged dictionary, and a
///   merge with no sequences takes the best available one, the reference's, or is refused when
///   there is none;
/// * **a column of interest neither file has is refused**, after both are read;
/// * **an annotation both files carry takes its file's label as a suffix**, and the columns are
///   written sorted.
///
/// Two inputs with DIFFERENT non-empty dictionaries are merged by `SamFileHeaderMerger`'s own
/// rules, which this does not port: the first input's dictionary is kept.
pub fn combine_segment_breakpoints(parser: &Parser) -> Outcome {
    use gatk_tools::sequence_dictionary::{compare, Compatibility};

    let _ = resolve_read_filters(parser, "CombineSegmentBreakpoints")?;
    let segments = arguments(parser, "segments");
    // The declaration's default is `[1, 2]`, and `--labels null` empties it.
    let labels = arguments(parser, "labels");
    let columns: Vec<String> = arguments(parser, "columns-of-interest");
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;

    // `AnnotatedIntervalCollection.create(path, columnsOfInterest)`: only the columns asked for.
    let read = |path: &str| -> Result<gatk_tools::annotated_interval::AnnotatedIntervalCollection, Thrown> {
        let mut collection = annotated_intervals(path)?;
        collection.annotations.retain(|name| columns.contains(name));
        for record in &mut collection.records {
            record.annotations.retain(|name, _| columns.contains(name));
        }
        Ok(collection)
    };
    let first = read(&segments[0])?;
    let second = read(&segments[1])?;

    // `createOutputSamFileHeader`.
    let dictionary_of =
        |collection: &gatk_tools::annotated_interval::AnnotatedIntervalCollection| {
            let text: String = collection
                .header_lines
                .iter()
                .map(|line| format!("{line}\n"))
                .collect();
            htsjdk_bam::reader::parse_header_text(&text).sequences
        };
    let (dictionary1, dictionary2) = (dictionary_of(&first), dictionary_of(&second));
    match compare(&dictionary1, &dictionary2, true) {
        Compatibility::UnequalCommonContigs => {
            return Err(bad_input(
                "Input files had common contigs with different lengths in the sequence dictionaries.  Were these segment files generated with the same reference?".to_string(),
            ))
        }
        Compatibility::OutOfOrder => {
            return Err(bad_input(
                "Input files have different sequence dictionary ordering.  The risk of errors downstream is too high to continue.".to_string(),
            ))
        }
        _ => {}
    }
    let mut merged = if dictionary1.is_empty() {
        dictionary2
    } else {
        dictionary1
    };
    if merged.is_empty() {
        let best = master_dictionary(parser)?.or(reference_dictionary(parser)?);
        match best {
            Some(best) => merged = best.sequences,
            None => {
                return Err(bad_input(
                    "Cannot assemble a reference dictionary.  In order to use this tool, one of the following conditions must be satisfied:  1)  One or both input files have a SAM File header ... 2)  A reference is provided (-R)".to_string(),
                ))
            }
        }
    }

    let mut seen: Vec<String> = first.annotations.clone();
    for name in &second.annotations {
        if !seen.contains(name) {
            seen.push(name.clone());
        }
    }
    let unused: Vec<String> = columns
        .iter()
        .filter(|name| !seen.contains(name))
        .cloned()
        .collect();
    if !unused.is_empty() {
        return Err(bad_input(format!(
            "Some columns of interest specified by the user were not seen in any input files: {}",
            unused.join(", ")
        )));
    }

    // A column both files carry is suffixed with each file's label, read with `get(0)` and
    // `get(1)`: an emptied `--labels` is refused there, and only when some column is shared.
    let columns_of = |collection: &gatk_tools::annotated_interval::AnnotatedIntervalCollection| {
        collection
            .records
            .first()
            .map(|record| record.annotations.keys().cloned().collect::<Vec<String>>())
            .unwrap_or_default()
    };
    let shared = columns_of(&first)
        .iter()
        .any(|name| columns_of(&second).contains(name));
    if shared && labels.len() < 2 {
        return Err(Thrown::non_user(
            "java.lang.IndexOutOfBoundsException",
            format!(
                "Index {} out of bounds for length {}",
                labels.len(),
                labels.len()
            ),
        ));
    }
    let label = |index: usize| labels.get(index).map(String::as_str).unwrap_or("");
    let names: Vec<String> = merged
        .iter()
        .map(|sequence| sequence.name.clone())
        .collect();
    let records = gatk_tools::combine_segment_breakpoints::combine(
        &first.records,
        &second.records,
        &names,
        [label(0), label(1)],
        &columns,
    );
    let mut header = SamHeader::default();
    header.attributes.set("GO", "none");
    header.attributes.set("SO", "coordinate");
    header.sequences = merged;
    let mut annotations: Vec<String> = records
        .first()
        .map(|record| record.annotations.keys().cloned().collect())
        .unwrap_or_default();
    annotations.sort();
    let collection = gatk_tools::annotated_interval::AnnotatedIntervalCollection {
        header_lines: header.encode().lines().map(str::to_string).collect(),
        comments: Vec::new(),
        annotations,
        records,
        contig_column: "CONTIG".to_string(),
        start_column: "START".to_string(),
        end_column: "END".to_string(),
    };
    write_file(&output, collection.write().as_bytes())?;
    Ok(None)
}

/// `MergeMutect2CallsWithMC3`: Mutect2 calls merged into the MC3 call set, one record per step of
/// the concordance walk.
///
/// The states and what each writes are [`gatk_tools::merge_mutect2_mc3`]; the runner builds the
/// records from the two files' own variants rather than from text:
///
/// * **the tumour sample is the eval header's `##tumor_sample`**, and a header without one is the
///   reference's `NullPointerException`;
/// * **the header is the TRUTH file's lines** with the standard `GT` and `AD` lines, the tool's own
///   and `M2_FILTERS`, over the one tumour sample;
/// * **a true positive and a filtered false negative keep the MC3 record** and add `M2` to its
///   `CENTERS`, the latter with the M2 filters; **a false negative** is written unchanged; **a false
///   positive** is rebuilt from the M2 site and alleles alone; **a filtered true negative** is
///   dropped;
/// * **every record's one genotype carries EVERY allele of the site** and the M2 depths, or MC3's
///   `NREF`/`NALT` where M2 has no record.
pub fn merge_mutect2_calls_with_mc3(parser: &Parser) -> Outcome {
    use gatk_engine::concordance_walker::ConcordanceState;
    use gatk_tools::merge_mutect2_mc3 as mc3;
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    use htsjdk_vcf::variant::Value;

    let _ = resolve_read_filters(parser, "MergeMutect2CallsWithMC3")?;
    let required = |name: &str| {
        argument(parser, name).ok_or_else(|| {
            Thrown::command_line(format!(
                "Argument {name} was missing: Argument '{name}' is required"
            ))
        })
    };
    let truth_path = required("truth")?;
    let eval_path = required("evaluation")?;
    let output = required("output")?;
    let read_vcf = |path: &str| -> Result<htsjdk_vcf::reader::VcfFile, Thrown> {
        let bytes = std::fs::read(path).map_err(|_| {
            Thrown::user(
                index_feature_file::Refusal::CouldNotReadInputFile {
                    path: path.to_string(),
                }
                .message(),
            )
        })?;
        let text = if gatk_tools::read_walker_refusal::is_block_compressed(&bytes) {
            htsjdk_bgzf::read::decompress_all(&bytes)
                .ok()
                .and_then(|bytes| String::from_utf8(bytes).ok())
                .ok_or_else(|| {
                    Thrown::non_user(
                        gatk_tools::read_walker_refusal::SAM_FORMAT,
                        format!("{path} is not a block compressed file"),
                    )
                })?
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })
    };
    let truth_file = read_vcf(&truth_path)?;
    let eval_file = read_vcf(&eval_path)?;
    let dictionary: Vec<String> = sequence_dictionary_of(&truth_file.header)
        .into_iter()
        .map(|(name, _)| name)
        .collect();

    // `onTraversalStart`.
    let tumor = eval_file
        .header
        .lines
        .iter()
        .find_map(|line| match line {
            HeaderLine::Unstructured { key, value } if key == "tumor_sample" => Some(value.clone()),
            _ => None,
        })
        .ok_or_else(|| {
            Thrown::non_user(
                "java.lang.NullPointerException",
                "Cannot invoke \"htsjdk.variant.vcf.VCFHeaderLine.getValue()\" because the return value of \"htsjdk.variant.vcf.VCFHeader.getMetaDataLine(String)\" is null",
            )
        })?;
    let mut header = truth_file.header.clone();
    let compound = |key: &str, id: &str, number: Cardinality, line_type: LineType, text: &str| {
        HeaderLine::Compound {
            key: key.to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: text.to_string(),
            extra: Vec::new(),
        }
    };
    let mut added = vec![
        compound(
            "FORMAT",
            "GT",
            Cardinality::Fixed(1),
            LineType::String,
            "Genotype",
        ),
        compound(
            "FORMAT",
            "AD",
            Cardinality::R,
            LineType::Integer,
            "Allelic depths for the ref and alt alleles in the order listed",
        ),
    ];
    added.extend(default_tool_vcf_header_lines(
        parser,
        "MergeMutect2CallsWithMC3",
    ));
    added.push(compound(
        "INFO",
        mc3::M2_FILTERS_KEY,
        Cardinality::Unbounded,
        LineType::String,
        "M2 filters applied to variant.",
    ));
    for line in added {
        // A `HashSet` of lines: an identical one collapses, nothing else does.
        if !header
            .lines
            .iter()
            .any(|existing| existing.render() == line.render())
        {
            header.lines.push(line);
        }
    }
    header.samples = vec![tumor.clone()];

    // The walk, with no filter on either side.
    let loci = |file: &htsjdk_vcf::reader::VcfFile| -> Vec<ConcordanceLocus> {
        file.records
            .iter()
            .enumerate()
            .map(|(index, record)| ConcordanceLocus {
                index,
                contig: record.contig.clone(),
                start: record.start as i32,
                filtered: record.is_filtered(),
            })
            .collect()
    };
    let steps = gatk_engine::concordance_walker::concordance(
        &loci(&truth_file),
        &loci(&eval_file),
        &dictionary,
        |truth, eval| {
            let truth = &truth_file.records[truth.index];
            let eval = &eval_file.records[eval.index];
            match truth.alternate_alleles().first() {
                None => false,
                Some(alternate) => {
                    truth.reference() == eval.reference()
                        && eval.alternate_alleles().contains(alternate)
                }
            }
        },
    );

    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut records = written.to_vec();
        apply_sites_only(parser, &mut header, &mut records);
        let out = write_vcf_honouring_lenient(parser, &header, &records)?;
        write_variant_output(parser, &output, &out)
    };
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for step in steps {
        let truth = step.truth.map(|index| &truth_file.records[index]);
        let eval = step.eval.map(|index| &eval_file.records[index]);
        let int_attribute = |record: &htsjdk_vcf::variant::VariantContext, key: &str| {
            record
                .attributes
                .iter()
                .find(|(name, _)| name == key)
                .and_then(|(_, value)| value.format())
                .and_then(|text| text.parse::<i32>().ok())
                .unwrap_or(0)
        };
        let depths = match eval {
            Some(eval) => {
                let Some(genotype) = eval.genotypes.iter().find(|g| g.sample_name == tumor) else {
                    finish(&written)?;
                    return Err(Thrown::non_user(
                        "java.lang.NullPointerException",
                        "Cannot invoke \"htsjdk.variant.variantcontext.Genotype.getAD()\" because the return value of \"htsjdk.variant.variantcontext.VariantContext.getGenotype(String)\" is null",
                    ));
                };
                genotype.ad.clone()
            }
            None => {
                let truth = truth.expect("a step with no eval has truth");
                Some(vec![
                    int_attribute(truth, mc3::MC3_REF_COUNT_KEY),
                    int_attribute(truth, mc3::MC3_ALT_COUNT_KEY),
                ])
            }
        };
        let site = truth.or(eval).expect("a step has at least one record");
        let mut genotype = htsjdk_vcf::variant::Genotype::new(&tumor, site.alleles.clone());
        genotype.ad = depths;
        let with_center = |mc3_record: &htsjdk_vcf::variant::VariantContext| {
            let mut out = mc3_record.clone();
            let mut centers: Vec<Value> = match out
                .attributes
                .iter()
                .find(|(name, _)| name == mc3::CENTERS_KEY)
            {
                Some((_, Value::List(items))) => items.clone(),
                Some((_, value)) => vec![value.clone()],
                None => Vec::new(),
            };
            centers.push(Value::Str(mc3::M2_CENTER_NAME.to_string()));
            out.attributes.retain(|(name, _)| name != mc3::CENTERS_KEY);
            out.attributes
                .push((mc3::CENTERS_KEY.to_string(), Value::List(centers)));
            out
        };
        let record = match step.state {
            ConcordanceState::TruePositive => Some(with_center(truth.expect("truth"))),
            ConcordanceState::FalsePositive => {
                let m2 = eval.expect("eval");
                let mut out = htsjdk_vcf::variant::VariantContext::new(
                    &m2.contig,
                    m2.start,
                    m2.alleles.clone(),
                );
                out.stop = m2.stop;
                out.attributes.push((
                    mc3::CENTERS_KEY.to_string(),
                    Value::Str(mc3::M2_CENTER_NAME.to_string()),
                ));
                Some(out)
            }
            ConcordanceState::FalseNegative => Some(truth.expect("truth").clone()),
            ConcordanceState::FilteredTrueNegative => None,
            ConcordanceState::FilteredFalseNegative => {
                let mut out = with_center(truth.expect("truth"));
                let filters = eval.expect("eval").filters.clone().unwrap_or_default();
                out.attributes.push((
                    mc3::M2_FILTERS_KEY.to_string(),
                    Value::List(filters.into_iter().map(Value::Str).collect()),
                ));
                Some(out)
            }
        };
        if let Some(mut record) = record {
            record.genotypes = vec![genotype].into();
            written.push(record);
        }
    }
    finish(&written)?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `FilterFuncotations`: a Funcotated VCF read twice and marked with the clinical-significance
/// filters its funcotations match.
///
/// The five filters are [`gatk_tools::filter_funcotations`]; the runner is the two passes and the
/// reading of the `FUNCOTATION` field:
///
/// * **the keys come from the header's `FUNCOTATION` line**, and a header without one is refused
///   before the writer exists;
/// * **each alternate allele's value splits into transcripts at `]#[`** and each transcript into
///   the header's keys at `|`, the value count checked against the key count, and an alternate
///   count that differs from the value count refused;
/// * **the first pass collects the compound het variants**, the second writes every record with
///   `CLINSIG` naming the matched filters, `PASS` when any matched and `NOT_CLINSIG` added to its
///   filters when none did.
pub fn filter_funcotations(parser: &Parser) -> Outcome {
    use gatk_tools::filter_funcotations as ff;
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    use htsjdk_vcf::variant::Value;

    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "FilterFuncotations")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference = match scalar(parser, "ref-version").as_deref() {
        Some("hg38") => ff::Reference::Hg38,
        Some("hg19") => ff::Reference::Hg19,
        _ => ff::Reference::B37,
    };
    let source = match scalar(parser, "allele-frequency-data-source").as_deref() {
        Some("gnomad") => ff::AlleleFrequencySource::Gnomad,
        _ => ff::AlleleFrequencySource::Exac,
    };

    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let description = file.header.lines.iter().find_map(|line| match line {
        HeaderLine::Compound {
            key,
            id,
            description,
            ..
        } if key == "INFO" && id == "FUNCOTATION" => Some(description.clone()),
        _ => None,
    });
    let Some(description) = description else {
        return Err(bad_input(
            "Could not extract Funcotation keys from FUNCOTATION field in input VCF header."
                .to_string(),
        ));
    };
    let keys =
        gatk_tools::funcotator::extract_funcotator_keys_from_header_description(&description);
    let mut header = file.header.clone();
    if !header.has_filter_line(ff::NOT_CLINSIG_FILTER) {
        header.lines.push(HeaderLine::Filter {
            id: ff::NOT_CLINSIG_FILTER.to_string(),
            description: "Filter for clinically insignificant variants.".to_string(),
        });
    }
    let clinsig = HeaderLine::Compound {
        key: "INFO".to_string(),
        id: ff::CLINSIG_INFO_KEY.to_string(),
        number: Cardinality::Fixed(1),
        line_type: LineType::String,
        description:
            "Rule(s) which caused this annotation to be flagged as clinically significant."
                .to_string(),
        extra: Vec::new(),
    };
    if !header
        .lines
        .iter()
        .any(|line| same_compound_id(line, &clinsig))
    {
        header.lines.push(clinsig);
    }
    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let mut header = header.clone();
        let mut records = written.to_vec();
        apply_sites_only(parser, &mut header, &mut records);
        let out = write_vcf_honouring_lenient(parser, &header, &records)?;
        write_variant_output(parser, &output, &out)
    };

    // `createAlleleToFuncotationMapFromFuncotationVcfAttribute`, one list per alternate, then the
    // transcripts of all of them.
    let transcripts_of =
        |record: &htsjdk_vcf::variant::VariantContext| -> Result<Vec<ff::Funcotations>, Thrown> {
            let values: Vec<String> = match record
                .attributes
                .iter()
                .find(|(name, _)| name == "FUNCOTATION")
                .map(|(_, value)| value)
            {
                None => Vec::new(),
                Some(Value::List(items)) => items
                    .iter()
                    .map(|item| item.format().unwrap_or_default())
                    .collect(),
                Some(value) => vec![value.format().unwrap_or_default()],
            };
            if values.len() != record.alternate_alleles().len() {
                return Err(Thrown::non_user(
                "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException",
                "Could not parse FUNCOTATION field properly.",
            ));
            }
            let mut transcripts = Vec::new();
            for value in values {
                for transcript in value.split("]#[").filter(|part| !part.is_empty()) {
                    let mut fields: Vec<&str> = transcript.split('|').collect();
                    if let Some(first) = fields.first_mut() {
                        *first = first.strip_prefix('[').unwrap_or(first);
                    }
                    if let Some(last) = fields.last_mut() {
                        *last = last.strip_suffix(']').unwrap_or(last);
                    }
                    if fields.len() != keys.len() {
                        return Err(Thrown::non_user(
                        "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException",
                        format!(
                            "Cannot parse the funcotation attribute.  Num values: {}   Num keys: {}",
                            fields.len(),
                            keys.len()
                        ),
                    ));
                    }
                    let pairs: Vec<(&str, &str)> = keys
                        .iter()
                        .map(String::as_str)
                        .zip(fields.iter().copied())
                        .collect();
                    transcripts.push(ff::Funcotations::new(&pairs));
                }
            }
            Ok(transcripts)
        };
    let variant_of = |record: &htsjdk_vcf::variant::VariantContext| {
        let het = |g: &htsjdk_vcf::variant::Genotype| {
            g.alleles.len() > 1
                && g.alleles.iter().all(|a| !a.is_no_call())
                && g.alleles.iter().any(|a| *a != g.alleles[0])
        };
        let hom_var = |g: &htsjdk_vcf::variant::Genotype| {
            !g.alleles.is_empty()
                && g.alleles
                    .iter()
                    .all(|a| !a.is_no_call() && !a.is_reference())
                && g.alleles.iter().all(|a| *a == g.alleles[0])
        };
        ff::Variant {
            contig: record.contig.clone(),
            start: record.start as i32,
            end: record.stop as i32,
            reference_allele: record.reference().display_string(),
            alternate_alleles: record
                .alternate_alleles()
                .iter()
                .map(|a| a.display_string())
                .collect(),
            het_count: record.genotypes.iter().filter(|g| het(g)).count() as i32,
            hom_var_count: record.genotypes.iter().filter(|g| hom_var(g)).count() as i32,
        }
    };

    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    // The first pass: every record's transcripts, for the compound het rule.
    let mut first_pass: Vec<(ff::Variant, Vec<ff::Funcotations>)> = Vec::new();
    for record in &kept {
        match transcripts_of(record) {
            Ok(transcripts) => first_pass.push((variant_of(record), transcripts)),
            Err(error) => {
                finish(&[])?;
                return Err(error);
            }
        }
    }
    let compound = ff::compound_het_variants(&first_pass, reference);

    // The second pass.
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for (record, (variant, transcripts)) in kept.iter().zip(&first_pass) {
        let matching =
            match ff::matching_filters(transcripts, variant, reference, source, &compound) {
                Ok(matching) => matching,
                Err(error) => {
                    finish(&written)?;
                    return Err(Thrown::non_user(
                        "java.lang.NumberFormatException",
                        error.message(),
                    ));
                }
            };
        let (value, significant) = ff::clinsig(&matching);
        let mut out = (*record).clone();
        out.attributes
            .retain(|(name, _)| name != ff::CLINSIG_INFO_KEY);
        out.attributes
            .push((ff::CLINSIG_INFO_KEY.to_string(), Value::Str(value)));
        out.filters = if significant {
            Some(Vec::new())
        } else {
            let mut filters = out.filters.clone().unwrap_or_default();
            if !filters.iter().any(|f| f == ff::NOT_CLINSIG_FILTER) {
                filters.push(ff::NOT_CLINSIG_FILTER.to_string());
            }
            Some(filters)
        };
        written.push(out);
    }
    finish(&written)?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}

/// `AnalyzeCovariates`: up to three BQSR reports folded into the csv the plotting script reads.
///
/// The csv is [`gatk_tools::analyze_covariates`]; the runner is `checkArgumentsValues` in its order
/// (each report checked as a file, then that there is one, then each output's location, then that
/// an output was asked for) and the write. The plots are R's and the reference only draws them when
/// `--plots-report-file` is given, which this port does not do: asking for them is the port's
/// limitation. The tool returns `Optional.empty()`, which `handleResult` prints.
pub fn analyze_covariates(parser: &Parser) -> Outcome {
    use gatk_tools::analyze_covariates as ac;

    let bqsr = argument(parser, "bqsr-recal-file");
    let before = argument(parser, "before-report-file");
    let after = argument(parser, "after-report-file");
    let plots = argument(parser, "plots-report-file");
    let csv = argument(parser, "intermediate-csv-file");
    let bad_value = |name: &str, message: String| {
        Thrown::command_line(format!("Argument {name} has a bad value: {message}"))
    };
    for (name, value) in [("BQSR", &bqsr), ("before", &before), ("after", &after)] {
        let Some(path) = value else { continue };
        let meta = std::fs::metadata(path);
        match meta {
            Err(_) => {
                return Err(bad_value(
                    name,
                    format!("input report '{path}' does not exist or is unreachable"),
                ))
            }
            Ok(meta) if !meta.is_file() => {
                return Err(bad_value(
                    name,
                    format!("input report '{path}' is not a regular file"),
                ))
            }
            Ok(_) => {}
        }
    }
    if bqsr.is_none() && before.is_none() && after.is_none() {
        return Err(Thrown::user(ac::AnalyzeCovariatesError::NoReport.message()));
    }
    for (name, value) in [("plots", &plots), ("csv", &csv)] {
        let Some(path) = value else { continue };
        let target = std::path::Path::new(path);
        if target.exists() && !target.is_file() {
            return Err(bad_value(
                name,
                format!("the output file location '{path}' exists as not a file"),
            ));
        }
        let Some(parent) = target
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        else {
            continue;
        };
        if !parent.exists() {
            return Err(bad_value(
                name,
                format!(
                    "the output file parent directory '{}' does not exists or is unreachable",
                    parent.display()
                ),
            ));
        }
        if !parent.is_dir() {
            return Err(bad_value(
                name,
                format!(
                    "the output file parent directory '{}' is not a directory",
                    parent.display()
                ),
            ));
        }
    }
    if plots.is_none() && csv.is_none() {
        return Err(Thrown::user(ac::AnalyzeCovariatesError::NoOutput.message()));
    }

    let mut parsed_reports: Vec<(&str, gatk_engine::recalibration_report::RecalibrationReport)> =
        Vec::new();
    for (role, value) in [("BQSR", &bqsr), ("Before", &before), ("After", &after)] {
        let Some(path) = value else { continue };
        let text = std::fs::read_to_string(path)
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
        let report = gatk_engine::recalibration_report::RecalibrationReport::parse(&text)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        parsed_reports.push((role, report));
    }
    let roles: Vec<ac::RoleReport> = parsed_reports
        .iter()
        .map(|(role, report)| ac::RoleReport { role, report })
        .collect();
    let text = ac::analyze_covariates(&roles, csv.is_some()).map_err(|error| Thrown {
        failure: Failure::User,
        exception: error.java_class(),
        message: Some(error.message()),
    })?;
    if let Some(path) = &csv {
        write_file(path, text.as_bytes())?;
    }
    if plots.is_some() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--plots-report-file draws the plots through R, which this port does not carry. This message is the port's own and not GATK's.",
        ));
    }
    Ok(Some("Optional.empty".to_string()))
}
