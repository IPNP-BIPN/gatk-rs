//! The Picard tools GATK bundles, run the way GATK runs them.
//!
//! `Main` wraps a Picard `CommandLineProgram` in `PicardCommandLineProgramExecutor`, whose
//! `instanceMain` is Picard's own: Picard parses the command line with its Barclay parser, prints a
//! refusal itself and RETURNS 1 rather than throwing, and `doWork` returns a code. Any code but
//! zero becomes a `PicardNonZeroExitException`, which `mainEntry` turns into status 4 after
//! `handleResult` has printed `Tool returned:` and the code; a zero is printed the same way and is
//! a success. What `doWork` throws reaches `mainEntry` like any GATK tool's exception.
//!
//! None of that is the GATK tools' path, which refuses a command line with the tool's usage, a
//! banner and status 2, so these tools are dispatched here before [`crate::run`] looks at them.
//!
//! Ported from `org.broadinstitute.hellbender.cmdline.PicardCommandLineProgramExecutor` and
//! `picard.cmdline.CommandLineProgram.instanceMain`.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{self, Failure, Thrown, PORT_FAILURE};

use crate::runners::{flag, scalar};

mod add_comments_to_bam;
mod build_bam_index;
mod check_terminator_block;
mod filter_vcf;
mod fix_vcf_header;
mod gather_bam_files;
mod gather_vcfs;
mod make_sites_only_vcf;
mod merge_vcfs;
mod rename_sample_in_vcf;
mod set_nm_and_uq_tags;
mod set_nm_md_and_uq_tags;
mod sort_vcf;
mod split_vcfs;
mod update_vcf_sequence_dictionary;
mod vcf_format_converter;
mod vcf_to_interval_list;

/// A Picard `doWork`: its return code, or what it threw. What it prints to standard error before
/// either goes into the string, in order.
pub(crate) type Work = fn(&Parser, &mut String) -> Result<i32, Thrown>;

/// A tool's own `customCommandLineValidation` message, for the tools that refuse more than the
/// parser does.
fn custom_validation(name: &str, parser: &Parser) -> Option<String> {
    match name {
        "FixVcfHeader" => {
            let header =
                crate::runners::argument(parser, "HEADER").filter(|header| header != "null");
            let check = scalar(parser, "CHECK_FIRST_N_RECORDS")
                .and_then(|value| value.parse::<i32>().ok())
                .unwrap_or(-1);
            (header.is_some() && 0 <= check).then(|| {
                "CHECK_FIRST_N_RECORDS should no be specified when HEADER is specified".to_string()
            })
        }
        _ => None,
    }
}

/// The Picard tools this port runs.
pub fn work(name: &str) -> Option<Work> {
    match name {
        "AddCommentsToBam" => Some(add_comments_to_bam::run),
        "BuildBamIndex" => Some(build_bam_index::run),
        "CheckTerminatorBlock" => Some(check_terminator_block::run),
        "FilterVcf" => Some(filter_vcf::run),
        "FixVcfHeader" => Some(fix_vcf_header::run),
        "GatherBamFiles" => Some(gather_bam_files::run),
        "GatherVcfs" => Some(gather_vcfs::run),
        "MakeSitesOnlyVcf" => Some(make_sites_only_vcf::run),
        "MergeVcfs" => Some(merge_vcfs::run),
        "RenameSampleInVcf" => Some(rename_sample_in_vcf::run),
        "SetNmAndUqTags" => Some(set_nm_and_uq_tags::run),
        "SetNmMdAndUqTags" => Some(set_nm_md_and_uq_tags::run),
        "SortVcf" => Some(sort_vcf::run),
        "SplitVcfs" => Some(split_vcfs::run),
        "UpdateVcfSequenceDictionary" => Some(update_vcf_sequence_dictionary::run),
        "VcfFormatConverter" => Some(vcf_format_converter::run),
        "VcfToIntervalList" => Some(vcf_to_interval_list::run),
        _ => None,
    }
}

/// `PicardCommandLineProgramExecutor.instanceMain` and what `mainEntry` makes of it, for a tool
/// [`work`] knows; `None` for any other.
pub fn run(name: &str, args: &[String]) -> Option<crate::Outcome> {
    let work = work(name)?;
    let mut stdout = String::new();
    let mut stderr = String::new();
    // `--help` makes Barclay print the usage and answer false, which is a return of 1 too, and it
    // answers before it checks a required argument.
    let help = args.iter().any(|arg| arg == "-h" || arg == "--help");
    let returned: Result<i32, Thrown> = match parse(name, args) {
        _ if help => {
            stderr.push_str(&crate::tool_usage(name).unwrap_or_default());
            Ok(1)
        }
        // `parseArgs` prints the tool's usage, then the parser's message, and returns false.
        // shortcut: the usage printed here is `usage(false, ...)`, without the common arguments,
        // and is not composed, only the message; compose it if a golden ever compares stderr.
        Err(message) => {
            stderr.push_str(&message);
            stderr.push('\n');
            Ok(1)
        }
        // `customCommandLineValidation`, inside `parseArgs`: the usage, then each message.
        Ok(parser) if custom_validation(name, &parser).is_some() => {
            stderr.push_str(&custom_validation(name, &parser).unwrap_or_default());
            stderr.push('\n');
            Ok(1)
        }
        Ok(parser) => {
            let _ = gatk_tools::tool_stdout::take();
            let result = work(&parser, &mut stderr);
            // `instanceMain`'s `finally`: a non-zero return earns the FAQ link, the run's last line.
            if matches!(result, Ok(code) if code != 0) && !flag(&parser, "QUIET") {
                stderr.push_str(
                    "To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp\n",
                );
            }
            stdout.push_str(&gatk_tools::tool_stdout::take());
            result
        }
    };
    Some(match returned {
        Ok(code) => {
            stdout
                .push_str(&main_entry::tool_returned(Some(&code.to_string())).unwrap_or_default());
            stdout.push('\n');
            let status = if code == 0 {
                0
            } else {
                main_entry::exit_status(Failure::PicardNonZero)
            };
            crate::Outcome {
                stdout,
                stderr,
                status,
            }
        }
        // An `Error` other than running out of memory is no `Exception`: `mainEntry` catches
        // neither, so the JVM's own handler prints it and exits 1.
        Err(thrown)
            if thrown.exception.starts_with("java.lang.")
                && thrown.exception.ends_with("Error")
                && thrown.exception != "java.lang.OutOfMemoryError" =>
        {
            stderr.push_str(&format!(
                "Exception in thread \"main\" {}",
                main_entry::non_user_exception_report(thrown.exception, thrown.message.as_deref())
            ));
            crate::Outcome {
                stdout,
                stderr,
                status: 1,
            }
        }
        Err(thrown) => {
            stderr.push_str(&thrown.report());
            crate::Outcome {
                stdout,
                stderr,
                status: thrown.status(),
            }
        }
    })
}

/// The tool's Barclay parser over the command line, or the message its refusal printed.
fn parse(name: &str, args: &[String]) -> Result<Parser, String> {
    let list = gatk_tools::tool_declarations::declarations(name)
        .expect("the declarations of a Picard tool with a runner");
    let mut parser = crate::parser_for(name, list);
    let borrowed: Vec<&str> = args.iter().map(String::as_str).collect();
    parser
        .parse_arguments(&borrowed)
        .map(|()| parser)
        .map_err(|error| error.message)
}

/// An exception by the class name a library reports it under.
pub(crate) fn thrown(class: &str, message: String) -> Thrown {
    // The class names are a handful per run, and the process ends with the run.
    Thrown::non_user(Box::leak(class.to_string().into_boxed_str()), message)
}

/// `FileExtensions.VCF_LIST`.
const VCF_LIST: &[&str] = &[".vcf", ".vcf.gz", ".bcf", ".vcf.bgz"];

/// `IOUtil.unrollPaths`: a name with a VCF extension is an input, any other is a file listing
/// more of them, and the order given is the order kept.
pub(crate) fn unroll_paths(inputs: &[String]) -> Result<Vec<String>, Thrown> {
    let mut stack: Vec<String> = inputs.to_vec();
    let mut output = Vec::new();
    while let Some(path) = stack.pop() {
        let name = std::path::Path::new(&path)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        if VCF_LIST.iter().any(|extension| name.ends_with(extension)) {
            output.push(path);
            continue;
        }
        let text = std::fs::read_to_string(&path).map_err(|_| {
            thrown(
                "java.lang.IllegalArgumentException",
                format!("had trouble reading from file://{}", absolute(&path)),
            )
        })?;
        stack.extend(
            text.lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(str::to_string),
        );
    }
    output.reverse();
    Ok(output)
}

/// A line of htsjdk's `Log`: level, `yyyy-MM-dd HH:mm:ss` in UTC (the pinned container's zone),
/// the logging class's simple name, and the message.
pub(crate) fn log_line(level: &str, class: &str, message: &str) -> String {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |since| since.as_secs() as i64);
    // Days to a civil date (Howard Hinnant's algorithm).
    let (days, time) = (seconds.div_euclid(86_400), seconds.rem_euclid(86_400));
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{level}\t{year:04}-{month:02}-{day:02} {:02}:{:02}:{:02}\t{class}\t{message}\n",
        time / 3600,
        (time % 3600) / 60,
        time % 60
    )
}

/// `File.getAbsolutePath()`, which Picard's messages print.
pub(crate) fn absolute(path: &str) -> String {
    crate::runners::java_absolute_path(path)
}

/// `IOUtil.assertFileIsReadable`.
pub(crate) fn assert_readable(path: &str) -> Result<(), Thrown> {
    let file = std::path::Path::new(path);
    let refuse = |message: String| Err(Thrown::non_user("htsjdk.samtools.SAMException", message));
    if !file.exists() {
        return refuse(format!(
            "Cannot read non-existent file: file://{}",
            absolute(path)
        ));
    }
    if file.is_dir() {
        return refuse(format!(
            "Cannot read file because it is a directory: file://{}/",
            absolute(path).trim_end_matches('/')
        ));
    }
    if std::fs::File::open(file).is_err() {
        return refuse(format!(
            "File exists but is not readable: {}",
            absolute(path)
        ));
    }
    Ok(())
}

/// `IOUtil.assertFileIsWritable`.
pub(crate) fn assert_writable(path: &str) -> Result<(), Thrown> {
    let file = std::path::Path::new(path);
    let refuse = |message: String| Err(Thrown::non_user("htsjdk.samtools.SAMException", message));
    if file.exists() {
        if file.is_dir() {
            return refuse(format!(
                "Cannot write file because it is a directory: {}",
                absolute(path)
            ));
        }
        if std::fs::OpenOptions::new().append(true).open(file).is_err() {
            return refuse(format!(
                "File exists but is not writable: {}",
                absolute(path)
            ));
        }
        return Ok(());
    }
    let absolute = absolute(path);
    let parent = std::path::Path::new(&absolute)
        .parent()
        .unwrap_or(std::path::Path::new("/"));
    if !parent.exists() {
        return refuse(format!(
            "Cannot write file: {absolute}. Neither file nor parent directory exist."
        ));
    }
    if !parent.is_dir() {
        return refuse(format!(
            "Cannot write file: {absolute}. File does not exist and parent is not a directory."
        ));
    }
    if parent
        .metadata()
        .is_ok_and(|meta| meta.permissions().readonly())
    {
        return refuse(format!(
            "Cannot write file: {absolute}. File does not exist and parent directory is not writable.."
        ));
    }
    Ok(())
}

/// A VCF's whole text, block-compressed or not, as `VCFFileReader` reads it.
pub(crate) fn vcf_text(path: &str) -> Result<String, Thrown> {
    let bytes = std::fs::read(path)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}")))?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let inflated = htsjdk_bgzf::read::decompress_all(&bytes).map_err(|error| {
            Thrown::non_user(
                "htsjdk.tribble.TribbleException",
                format!("{path}: {error:?}"),
            )
        })?;
        return Ok(String::from_utf8_lossy(&inflated).into_owned());
    }
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

/// `VCFHeader.setSequenceDictionary` on a header that has no contig lines: the dictionary's, in
/// front of the column line, as `VCFContigHeaderLine(record, record.getAssembly())` renders them.
pub(crate) fn with_contigs(text: &str, sequences: &[htsjdk_bam::header::SequenceRecord]) -> String {
    let lines: String = sequences
        .iter()
        .map(|record| match record.attributes.get("AS") {
            Some(assembly) => format!(
                "##contig=<ID={},length={},assembly={assembly}>\n",
                record.name, record.length
            ),
            None => format!("##contig=<ID={},length={}>\n", record.name, record.length),
        })
        .collect();
    match text.find("#CHROM") {
        Some(at) => format!("{}{lines}{}", &text[..at], &text[at..]),
        None => text.to_string(),
    }
}

/// A VCF written as `VariantContextWriterBuilder` writes it for a Picard tool: block compressed
/// for a `.gz` name, and indexed on the fly when `index` says so: `CREATE_INDEX` for a tool that
/// hands it to the builder, always for one that keeps htsjdk's default `INDEX_ON_THE_FLY`.
pub(crate) fn write_vcf(
    parser: &Parser,
    output: &str,
    text: &str,
    index: bool,
) -> Result<(), Thrown> {
    write_vcf_with_dictionary(parser, output, text, index, None)
}

/// `write_vcf` for a writer built with `setReferenceDictionary(dictionary)`: the Tribble index
/// takes its `DICT:` properties from that dictionary rather than from the header's contig lines.
pub(crate) fn write_vcf_with_dictionary(
    parser: &Parser,
    output: &str,
    text: &str,
    index: bool,
    dictionary: Option<&[htsjdk_bam::header::SequenceRecord]>,
) -> Result<(), Thrown> {
    let block_compressed = output.ends_with(".gz") || output.ends_with(".bgz");
    let has_dictionary = text
        .lines()
        .take_while(|line| line.starts_with("##"))
        .any(|line| line.starts_with("##contig="));
    // `VariantContextWriterBuilder.build()`: a Tribble index needs a dictionary, a tabix one does
    // not.
    if index && !block_compressed && !has_dictionary {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "A reference dictionary is required for creating Tribble indices on the fly",
        ));
    }
    if index && block_compressed {
        let (level, deflater) = crate::runners::output_compression(parser);
        let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
        std::io::Write::write_all(&mut writer, text.as_bytes())
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
        let compressed = writer
            .into_inner()
            .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
        if let Some(refusal) = tabix_order_refusal(text, &compressed) {
            return Err(refusal);
        }
    }
    let dictionary = dictionary.map(|sequences| {
        sequences
            .iter()
            .map(|sequence| (sequence.name.clone(), sequence.length))
            .collect()
    });
    crate::runners::write_variant_file(parser, output, text, index, false, dictionary)
}

/// `TabixIndexCreator.addFeature`'s order checks, which the on-the-fly index makes as the writer
/// adds each record, at the stream's BGZF virtual offset BEFORE the record. The reference index of
/// a feature is the order its contig first appeared in, not the dictionary's.
fn tabix_order_refusal(text: &str, compressed: &[u8]) -> Option<Thrown> {
    // (compressed address, first uncompressed byte, uncompressed length) of every block.
    let mut blocks = Vec::new();
    let (mut address, mut uncompressed) = (0usize, 0usize);
    while address + 18 <= compressed.len() {
        let size = usize::from(u16::from_le_bytes([
            compressed[address + 16],
            compressed[address + 17],
        ])) + 1;
        let end = address + size;
        if end > compressed.len() || size < 26 {
            break;
        }
        let length = u32::from_le_bytes(compressed[end - 4..end].try_into().ok()?) as usize;
        blocks.push((address, uncompressed, length));
        address = end;
        uncompressed += length;
    }
    let virtual_offset = |offset: usize| -> i64 {
        blocks
            .iter()
            .find(|(_, first, length)| offset >= *first && offset < first + length)
            .map(|(block, first, _)| ((*block as i64) << 16) | (offset - first) as i64)
            .unwrap_or((address as i64) << 16)
    };
    let mut names: Vec<&str> = Vec::new();
    let mut previous: Option<(usize, i64, i64, i64)> = None;
    let mut offset = 0usize;
    for line in text.split_inclusive('\n') {
        let at = offset;
        offset += line.len();
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
        let (contig, start) = (fields[0], fields.get(1)?.parse::<i64>().ok()?);
        // `getEnd()`: the INFO `END` when there is one, else the last base of the REF allele.
        let end = fields
            .get(7)
            .and_then(|info| info.split(';').find_map(|pair| pair.strip_prefix("END=")))
            .and_then(|value| value.parse().ok())
            .unwrap_or(start + fields.get(3).map_or(1, |reference| reference.len() as i64) - 1);
        let index = if names.last() == Some(&contig) {
            names.len() - 1
        } else {
            if names.contains(&contig) {
                // shortcut: htsjdk prints the whole `VariantContext.toString()` here, which this
                // does not reproduce; port it if a row ever reaches this before the check below.
                return Some(Thrown::non_user(
                    "java.lang.IllegalArgumentException",
                    format!("Sequence {contig}:{start} added out sequence of order"),
                ));
            }
            names.push(contig);
            names.len() - 1
        };
        let position = virtual_offset(at);
        if let Some((last_index, last_start, last_end, last_position)) = previous {
            if (last_index, last_start) > (index, start) {
                return Some(Thrown::non_user(
                    "java.lang.IllegalArgumentException",
                    format!(
                        "Features added out of order: previous (TabixFeature{{referenceIndex={last_index}, start={last_start}, end={last_end}, featureStartFilePosition={last_position}, featureEndFilePosition=-1}}) > next (TabixFeature{{referenceIndex={index}, start={start}, end={end}, featureStartFilePosition={position}, featureEndFilePosition=-1}})"
                    ),
                ));
            }
        }
        previous = Some((index, start, end, position));
    }
    None
}

/// `COMPRESSION_LEVEL` when the command line set it, which `output_compression` reads in place of
/// GATK's configured level: the field's default is `Defaults.COMPRESSION_LEVEL`, which is the same
/// configured level, so an unset one is no different from a GATK tool's.
pub(crate) fn compression_level(parser: &Parser) -> Option<u32> {
    parser
        .definitions()
        .iter()
        .find(|definition| definition.long_name() == "COMPRESSION_LEVEL")
        .filter(|definition| definition.has_been_set())
        .and_then(|_| scalar(parser, "COMPRESSION_LEVEL"))
        .and_then(|value| value.parse().ok())
}

#[cfg(test)]
mod tests {
    #[test]
    fn every_picard_tool_parses_its_whole_argument_surface() {
        for tool in [
            "AddCommentsToBam",
            "BuildBamIndex",
            "CheckTerminatorBlock",
            "FilterVcf",
            "FixVcfHeader",
            "GatherBamFiles",
            "GatherVcfs",
            "MakeSitesOnlyVcf",
            "MergeVcfs",
            "RenameSampleInVcf",
            "SetNmAndUqTags",
            "SetNmMdAndUqTags",
            "SortVcf",
            "SplitVcfs",
            "UpdateVcfSequenceDictionary",
            "VcfFormatConverter",
            "VcfToIntervalList",
        ] {
            assert!(super::work(tool).is_some(), "{tool} has no runner");
            assert!(
                crate::parseable(tool),
                "{tool} has an argument the parser cannot convert"
            );
        }
    }

    #[test]
    fn a_refused_command_line_returns_one_and_exits_four() {
        let outcome = super::run("SortVcf", &["--NO_SUCH".to_string(), "1".to_string()]).unwrap();
        assert_eq!(outcome.status, 4);
        assert_eq!(outcome.stdout, "Tool returned:\n1\n");
        assert!(outcome
            .stderr
            .ends_with("NO_SUCH is not a recognized option\n"));
    }
}
