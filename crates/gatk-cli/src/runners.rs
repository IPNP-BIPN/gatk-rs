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
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::index_feature_file::{self, Refusal, Source};
use gatk_tools::main_entry::{Thrown, PORT_FAILURE, PORT_LIMITATION};
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
    // Every refusal the port makes here is a UserException in the reference, EXCEPT the one the
    // reference does not make: a tabix index this port cannot write is the port's own gap, and it
    // is reported as one rather than as a status the reference would have exited with.
    let refused = |refusal: Refusal| match refusal {
        Refusal::TabixIsNotWritten { .. } => Thrown::non_user(PORT_LIMITATION, refusal.message()),
        _ => Thrown::user(refusal.message()),
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
    let index = index_feature_file::build(&text, &source, &input).map_err(refused)?;
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
fn read_filter<'a>(
    parser: &'a Parser,
    tool: &str,
    header: &'a SamHeader,
) -> Result<Filter<'a>, Thrown> {
    let mut names: Vec<String> = Vec::new();
    if !flag(parser, "disable-tool-default-read-filters") {
        names.extend(
            gatk_tools::plugin_ownership::default_filters(tool)
                .unwrap_or(&[])
                .iter()
                .map(|name| (*name).to_string()),
        );
    }
    names.extend(arguments(parser, "read-filter"));

    let mut plain: Vec<gatk_readfilter::ReadFilter> = Vec::new();
    let mut wellformed = false;
    let mut parameterized: Vec<gatk_readfilter::Parameterized> = Vec::new();
    for name in &names {
        if name == "WellformedReadFilter" {
            wellformed = true;
        } else if let Some(filter) = gatk_readfilter::by_name(name) {
            plain.push(filter);
        } else if name == "MappingQualityReadFilter" {
            let minimum = scalar(parser, "minimum-mapping-quality")
                .and_then(|text| text.parse::<i32>().ok())
                .unwrap_or(10);
            let maximum =
                scalar(parser, "maximum-mapping-quality").and_then(|text| text.parse::<i32>().ok());
            parameterized.push(gatk_readfilter::Parameterized::MappingQuality {
                min: minimum,
                max: maximum,
            });
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
        if wellformed && !gatk_readfilter::with_header::wellformed(read, header) {
            return false;
        }
        if !plain.iter().all(|filter| filter(read)) {
            return false;
        }
        parameterized
            .iter()
            .all(|filter| filter.decide(read).unwrap_or(false))
    }))
}

/// `CountReads.doWork`, with the input read and the output written.
///
/// Three things the `count-reads-plumbing` golden pins and this reproduces: the tool RETURNS the
/// count, so `handleResult` prints a number; `-O` receives that number and nothing else, with no
/// trailing newline, because the reference writes it with `print`; and `-O` does not suppress the
/// return, so the file is written AND the value comes back.
pub fn count_reads(parser: &Parser) -> Outcome {
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
    // What a walker makes of the file is the READER's answer and not the tool's, and it is three
    // different answers with two different statuses: `read-walker-refusals` measured them.
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
    let intervals_given = !arguments(parser, "intervals").is_empty();
    if let Some(refusal) = gatk_tools::read_walker_refusal::refusal(
        &input,
        path.exists(),
        path.is_dir(),
        decompressed.as_deref(),
        compressed,
        intervals_given,
    ) {
        // The refusal names the exception the reference throws, and the non-user handler PRINTS
        // that class: a walker refusing an interval over a stream with no dictionary answers
        // `java.lang.IllegalArgumentException: ...` and not a banner (#1020).
        return Err(if refusal.is_user() {
            Thrown::user(refusal.message())
        } else {
            Thrown::non_user(refusal.exception(), refusal.message())
        });
    }
    // The index is the one htsjdk's own search finds, not the one a single `with_extension` call
    // guesses. `reads.bam.bai` was the only name asked for here, and htsjdk writes `reads.bai` at
    // least as often and looks for it FIRST, so an interval query over a file indexed the other
    // way found no index and answered zero rather than refusing (#1020).
    let source = match htsjdk_bam::sam_files::find_index(path) {
        Some(index) => ReadsDataSource::open(path, &index),
        None => ReadsDataSource::open_unindexed(path),
    }
    .map_err(|error| Thrown::user(format!("{error:?}")))?;

    let header = source.header().clone();
    let mut intervals = Vec::new();
    for query in arguments(parser, "intervals") {
        let interval = gatk_engine::interval::parse_interval(&query, &header)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        intervals.push(interval);
    }
    let filter = read_filter(parser, "CountReads", &header)?;
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
