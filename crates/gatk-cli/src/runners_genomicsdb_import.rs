//! `GenomicsDBImport`: the checks around the native importer, and the workspace it writes.
//!
//! Every file of a workspace is the reference's byte for byte; only the names of the fragments and
//! of the `genomicsdb_meta_*.json` files, which hold a random UUID, a thread id and the wall clock,
//! and the `Date` of the command line in `vcfheader.vcf` are this run's own, as they are from one
//! run of the reference to the next. `gatk_tools::genomicsdb` says how each file is made.
//!
//! What the port has not measured is refused with its own message rather than written. A refusal
//! once the workspace exists removes what was written of it, so nothing is left that reads as the
//! reference's output. So does a failure the reference has words for, where the reference leaves
//! the part it had written, which depends on its buffers and threads: the words are the same and
//! the leftovers are not reproduced.

use super::*;
use gatk_tools::genomicsdb as gdb;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use std::path::{Path, PathBuf};

fn limitation(what: &str) -> Thrown {
    Thrown::non_user(
        PORT_LIMITATION,
        format!("{what} This message is the port's own and not GATK's."),
    )
}

/// One input as `SampleNameMap` holds it, with what the reader will open.
struct Sample {
    name: String,
    path: String,
    index: Option<String>,
}

/// `new URI(path)`, whose only failure a local path can reach is an illegal character.
fn check_uri(path: &str) -> Result<(), Thrown> {
    if path.contains("://") {
        return Err(limitation(
            "An input or workspace on a remote file system needs a client this port does not \
             carry.",
        ));
    }
    if path.chars().any(|character| {
        character.is_whitespace()
            || matches!(
                character,
                '"' | '<' | '>' | '\\' | '^' | '`' | '{' | '|' | '}'
            )
    }) {
        return Err(limitation(
            "A path that java.net.URI refuses fails with a message quoting URI's own parser, \
             which this port does not carry.",
        ));
    }
    Ok(())
}

/// The bytes of a feature file, decompressed, and the parsed header of its header lines alone,
/// which is all `getHeaderFromPath` reads.
fn read_header(path: &str) -> Result<(String, htsjdk_vcf::reader::VcfFile), Thrown> {
    if !Path::new(path).exists() {
        let absolute = java_absolute_path(path);
        return Err(Thrown::non_user(
            "htsjdk.tribble.TribbleException$MalformedFeatureFile",
            format!(
                "Unable to parse header with error: {absolute}, for input source: file://{absolute}"
            ),
        ));
    }
    let (_, text) = open_feature_input(path)?;
    let header_end = header_length(&text);
    let file = htsjdk_vcf::reader::read_vcf(&text[..header_end]).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;
    Ok((text, file))
}

/// Where the header ends: the byte after the `#CHROM` line, or the whole text without one.
fn header_length(text: &str) -> usize {
    let mut position = 0;
    for line in text.split_inclusive('\n') {
        position += line.len();
        if line.starts_with("#CHROM") {
            return position;
        }
        if !line.starts_with('#') {
            return position - line.len();
        }
    }
    position
}

/// `SampleNameMap(Path, boolean)`.
fn load_sample_name_map(path: &str, check_compressed: bool) -> Result<Vec<Sample>, Thrown> {
    let text = std::fs::read_to_string(path).map_err(|_| {
        Thrown::user(format!(
            "Couldn't read file file://{}. Error was: exception while reading sample->filename \
             mapping file",
            java_absolute_path(path)
        ))
    })?;
    // `Files.readAllLines`.
    let lines: Vec<&str> = text
        .split('\n')
        .map(|line| line.strip_suffix('\r').unwrap_or(line))
        .collect();
    let lines = match lines.split_last() {
        Some((&"", rest)) => rest,
        _ => &lines[..],
    };
    if lines.is_empty() {
        return Err(Thrown::user(
            "Bad input: At least 1 sample is required but none were found in the sample mapping \
             file",
        ));
    }
    let mut samples: Vec<Sample> = Vec::new();
    for line in lines {
        let split: Vec<&str> = line.split('\t').collect();
        if split.len() != 2 && split.len() != 3 {
            return Err(Thrown::user(format!(
                "Bad input: Sample name map file must have 2 or 3 fields per line in the \
                 format:\nSample\tFile\nor:\nSample\tFile\tIndex\nbut found line: \"{line}\" with \
                 {} fields",
                split.len()
            )));
        }
        let legal = !split[0].trim().is_empty() && split[0].trim() == split[0];
        if !legal || split[1].trim().is_empty() {
            return Err(Thrown::user(format!(
                "Bad input: Sample name map file must have lines in the \
                 format:\nSample\tFile\nor:\nSample\tFile\tIndex\n but found line: '{line}'\nValid \
                 sample names must be non-empty strings that cannot begin or end with whitespace \
                 and valid file names must be non-empty and not all whitespace"
            )));
        }
        let vcf = split[1].trim().to_string();
        let index = match split.get(2) {
            Some(index) if index.trim().is_empty() => {
                return Err(Thrown::user(format!(
                    "Bad input: Found a line in the sample name map file with an empty or \
                     all-whitespace value for the index:\n\"{line}\""
                )))
            }
            Some(index) => Some(index.trim().to_string()),
            None => None,
        };
        check_uri(&vcf)?;
        if let Some(existing) = samples.iter().find(|sample| sample.name == split[0]) {
            return Err(Thrown::user(format!(
                "Bad input: Found two mappings for the same sample: {}\n{vcf}\n{}",
                split[0], existing.path
            )));
        }
        if check_compressed {
            assert_compressed_and_indexed(&vcf, index.as_deref())?;
        }
        samples.push(Sample {
            name: split[0].to_string(),
            path: vcf,
            index,
        });
    }
    samples.sort_by(|left, right| left.name.cmp(&right.name));
    Ok(samples)
}

/// `GATKGenomicsDBUtils.assertVariantFileIsCompressedAndIndexed`.
fn assert_compressed_and_indexed(path: &str, index: Option<&str>) -> Result<(), Thrown> {
    if !path.to_lowercase().ends_with(".vcf.gz") {
        return Err(Thrown::user(format!(
            "Input variant files must be block compressed vcfs when using bypass-feature-reader, \
             but {path} does not end with the standard file extension .vcf.gz"
        )));
    }
    let index_path = match index {
        Some(index) => index.to_string(),
        None => format!("{path}.tbi"),
    };
    if !Path::new(&index_path).is_file() {
        return Err(limitation(
            "A missing index fails in IOUtils.assertFileIsReadable with a message this port has \
             not measured.",
        ));
    }
    Ok(())
}

/// `GATKGenomicsDBUtils.genomicsDBGetAbsolutePath`, which TileDB's `real_dir` answers without
/// touching the file system.
fn workspace_absolute_path(path: &str) -> String {
    let absolute = java_absolute_path(path);
    let mut parts: Vec<&str> = Vec::new();
    for part in absolute.split('/') {
        match part {
            "" | "." => {}
            ".." => {
                parts.pop();
            }
            other => parts.push(other),
        }
    }
    format!("/{}", parts.join("/"))
}

/// A random identifier in the shape of `uuid_generate`, for the names that hold one.
fn uuid_like(counter: &mut u64) -> String {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_nanos() as u64)
        .unwrap_or(0);
    *counter += 1;
    let mix = |seed: u64| {
        let mut value = seed.wrapping_add(0x9e37_79b9_7f4a_7c15);
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    };
    let high = mix(nanos ^ (std::process::id() as u64) << 32 ^ *counter);
    let low = mix(high ^ nanos.rotate_left(17));
    format!(
        "{:08x}-{:04x}-4{:03x}-{:04x}-{:012x}",
        (high >> 32) as u32,
        (high >> 16) as u16,
        (high & 0x0fff) as u16,
        ((low >> 48) as u16 & 0x3fff) | 0x8000,
        low & 0xffff_ffff_ffff
    )
}

fn write_file(path: &Path, bytes: &[u8]) -> Result<(), Thrown> {
    std::fs::write(path, bytes).map_err(|error| {
        limitation(&format!(
            "Writing {} failed ({error}), which the native importer reports in words this port \
             has not measured.",
            path.display()
        ))
    })
}

/// A run's progress on disk, removed unless it ran to the end.
struct Workspace {
    path: PathBuf,
    finished: bool,
}

impl Drop for Workspace {
    fn drop(&mut self) {
        if !self.finished {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }
}

pub fn genomics_db_import(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "GenomicsDBImport")?;

    // `initializeWorkspaceAndToolMode`.
    let update_workspace = argument(parser, "genomicsdb-update-workspace-path")
        .filter(|workspace| !workspace.is_empty());
    let interval_list_output =
        argument(parser, "output-interval-list-to-file").filter(|output| !output.is_empty());
    let incremental = update_workspace.is_some();
    let intervals_from_workspace = incremental && interval_list_output.is_some();
    let variants = arguments(parser, "variant");
    let sample_name_map = argument(parser, "sample-name-map");
    let bypass = flag(parser, "bypass-feature-reader");
    let avoid_nio = flag(parser, "avoid-nio");
    let header_override = argument(parser, "header");
    let overwrite = flag(parser, "overwrite-existing-genomicsdb-workspace");

    // `assertVariantPathsOrSampleNameFileWasSpecified`.
    if variants.is_empty() && sample_name_map.is_none() && !intervals_from_workspace {
        return Err(Thrown::command_line(
            "Argument variant was missing: One of --variant or --sample-name-map must be specified",
        ));
    }
    // `assertOverwriteWorkspaceAndIncrementalImportMutuallyExclusive`.
    if overwrite && incremental {
        return Err(Thrown::command_line(
            "overwrite-existing-genomicsdb-workspace cannot be set to true when \
             genomicsdb-update-workspace-path is set",
        ));
    }
    // `assertAvoidNioConditionsAreValid`.
    if avoid_nio && (!bypass || header_override.is_none()) {
        let mut missing = Vec::new();
        if !bypass {
            missing.push("bypass-feature-reader");
        }
        if header_override.is_none() {
            missing.push("header");
        }
        return Err(Thrown::command_line(format!(
            "Argument {} was missing: If --avoid-nio is set then --bypass-feature-reader and \
             --header must also be specified.",
            missing.join(" and ")
        )));
    }
    // An update reads the existing workspace back through the native library, which needs the
    // reader this port does not have.
    if incremental {
        return Err(limitation(
            "--genomicsdb-update-workspace-path reads an existing workspace's arrays back, which \
             this port does not carry.",
        ));
    }

    if header_override
        .as_deref()
        .is_some_and(|path| !path.starts_with('/'))
    {
        return Err(limitation(
            "A relative --header is recorded absolute in the command line of vcfheader.vcf, which \
             the shared rendering of the command line does not do.",
        ));
    }

    // `initializeHeaderAndSampleMappings`.
    let mut samples: Vec<Sample> = Vec::new();
    let mut merged: Vec<HeaderLine>;
    if !variants.is_empty() {
        let mut headers = Vec::with_capacity(variants.len());
        if let Some(path) = &header_override {
            check_uri(path)?;
            let (_, file) = read_header(path)?;
            merged = file.meta_data_in_input_order();
        } else {
            for path in &variants {
                check_uri(path)?;
                if bypass {
                    assert_compressed_and_indexed(path, None)?;
                }
                let (_, file) = read_header(path)?;
                let count = file.header.samples.len();
                if count != 1 {
                    return Err(Thrown::user(format!(
                        "Input GVCF: {path} was expected to contain a single sample but actually \
                         contained {count} samples."
                    )));
                }
                let name = file.header.samples[0].clone();
                if let Some(existing) = samples.iter().find(|sample| sample.name == name) {
                    return Err(Thrown::user(format!(
                        "Bad input: Duplicate sample: {name}. Sample was found in both {path} and \
                         {}.",
                        existing.path
                    )));
                }
                samples.push(Sample {
                    name,
                    path: path.clone(),
                    index: None,
                });
                headers.push(file);
            }
            let versions: Vec<Option<String>> = headers
                .iter()
                .map(|file| {
                    file.header_version
                        .map(|version| version.version_string().to_string())
                })
                .collect();
            let sources: Vec<htsjdk_vcf::merge::Source> = headers
                .iter()
                .zip(&versions)
                .map(|(file, version)| htsjdk_vcf::merge::Source {
                    header: &file.header,
                    version: version.as_deref(),
                })
                .collect();
            let contigs = |file: &htsjdk_vcf::reader::VcfFile| -> Vec<HeaderLine> {
                file.header
                    .lines
                    .iter()
                    .filter(|line| matches!(line, HeaderLine::Contig { .. }))
                    .cloned()
                    .collect()
            };
            if headers
                .iter()
                .any(|file| contigs(file) != contigs(&headers[0]))
            {
                return Err(limitation(
                    "Inputs whose sequence dictionaries differ leave contigs out of vcfheader.vcf \
                     in a way this port has not measured.",
                ));
            }
            merged = htsjdk_vcf::merge::smart_merge_headers(&sources, true)
                .map_err(|error| Thrown::non_user(error.class(), error.message()))?
                .0;
        }
        // `sampleNameToVcfPath` is a `TreeMap`, so the rows are in the order of the names.
        samples.sort_by(|left, right| left.name.cmp(&right.name));
        if header_override.is_some() {
            return Err(limitation(
                "--header with --variant imports samples this port takes from the headers it does \
                 not then read, which it has not measured.",
            ));
        }
    } else if let Some(map) = &sample_name_map {
        samples = load_sample_name_map(map, bypass && !avoid_nio)?;
        let (_, file) = if let Some(path) = &header_override {
            read_header(path)?
        } else {
            read_header(&samples[0].path)?
        };
        merged = file.meta_data_in_input_order();
    } else {
        unreachable!("the update path is refused above");
    }
    // `mergedHeaderLines.addAll(getDefaultToolVCFHeaderLines())`: a set, so a line already there
    // keeps its place.
    for line in default_tool_vcf_header_lines(parser, "GenomicsDBImport") {
        if !merged.contains(&line) {
            merged.push(line);
        }
    }

    // `new VCFHeader(mergedHeaderLines).getSequenceDictionary()`.
    let has_contig = merged
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }));
    if !has_contig {
        return Err(Thrown::user(
            "The merged vcf header has no sequence dictionary. Please provide a header that \
             contains a sequence dictionary.",
        ));
    }
    if samples.iter().any(|sample| sample.index.is_some()) && bypass {
        return Err(Thrown::user(
            "Indices were specified for some VCFs in the sample name map file, but \
             --bypass-feature-reader was also specified. Specifying explicit indices is not \
             supported when running with --bypass-feature-reader",
        ));
    }
    let dictionary_text: String = merged
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .map(|line| format!("##{}\n", line.render()))
        .collect();
    let dictionary = vcf_dictionary(&dictionary_text);

    // `initializeIntervals`.
    if arguments(parser, "intervals").is_empty()
        && arguments(parser, "exclude-intervals").is_empty()
    {
        return Err(Thrown::user("No intervals specified"));
    }
    let parameters = interval_arguments(parser, &dictionary)?;
    let mut intervals = parameters
        .map(|parameters| parameters.intervals)
        .unwrap_or_default();
    if flag(parser, "merge-input-intervals") {
        // `IntervalUtils.getSpanningIntervals`.
        let mut spanning: Vec<gatk_engine::interval::SimpleInterval> = Vec::new();
        for interval in intervals {
            match spanning
                .iter_mut()
                .find(|existing| existing.contig == interval.contig)
            {
                Some(existing) => {
                    existing.start = existing.start.min(interval.start);
                    existing.end = existing.end.max(interval.end);
                }
                None => spanning.push(interval),
            }
        }
        intervals = spanning;
    }

    // The rest of `GATKTool.onStartup`.
    let master = master_dictionary(parser)?;
    let reference = reference_dictionary(parser)?;
    let inputs = arguments(parser, "input");
    let reads: Vec<Vec<htsjdk_bam::header::SequenceRecord>> = inputs
        .iter()
        .map(|path| reads_dictionary(parser, path, inputs.len()))
        .collect::<Result<_, _>>()?;
    if !flag(parser, "disable-sequence-dictionary-validation") {
        if let Some(master) = &master {
            for one in &reads {
                validate_against_master(master, "reads", one)?;
            }
            if let Some(reference) = &reference {
                validate_against_master(master, "reference", &reference.sequences)?;
            }
        }
        if let Some(reference) = &reference {
            for one in &reads {
                gatk_tools::sequence_dictionary::validate(
                    "reference",
                    &reference.sequences,
                    "reads",
                    one,
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
    }

    // What the native library does for these, this port has not carried.
    if bypass {
        return Err(limitation(
            "--bypass-feature-reader has the native library read the inputs itself, through \
             htslib's parser rather than htsjdk's, which this port does not carry.",
        ));
    }
    if number(parser, "merge-contigs-into-num-partitions") != 0 {
        return Err(limitation(
            "--merge-contigs-into-num-partitions lays several contigs into one array, which this \
             port has not measured.",
        ));
    }
    let batch_size = number(parser, "batch-size");
    let batch_size = if batch_size == 0 {
        samples.len()
    } else {
        batch_size as usize
    };
    if batch_size < samples.len() && flag(parser, "consolidate") {
        return Err(limitation(
            "--consolidate over several batches rewrites the fragments through TileDB's \
             consolidation, which this port does not carry.",
        ));
    }

    // The vid map and the schema come first: a refusal here leaves nothing behind.
    let vid = gdb::vid_map(&merged).map_err(|reason| limitation(&reason))?;
    let schema = gdb::Schema::new(&vid, &merged).map_err(|reason| limitation(&reason))?;

    // `overwriteCreateOrCheckWorkspace`.
    let workspace_argument = argument(parser, "genomicsdb-workspace-path")
        .expect("the parser requires one of the two workspace arguments");
    check_uri(&workspace_argument)?;
    let workspace_path = PathBuf::from(workspace_absolute_path(&workspace_argument));
    create_workspace(&workspace_path, &workspace_argument, overwrite)?;
    let mut workspace = Workspace {
        path: workspace_path.clone(),
        finished: false,
    };
    write_file(&workspace_path.join("__tiledb_workspace.tdb"), &[])?;

    let names: Vec<String> = samples.iter().map(|sample| sample.name.clone()).collect();
    write_file(
        &workspace_path.join("callset.json"),
        gdb::callset_json(&names).as_bytes(),
    )?;
    write_file(
        &workspace_path.join("vidmap.json"),
        vid.to_json().as_bytes(),
    )?;
    let header = VcfHeader {
        lines: merged.clone(),
        samples: Vec::new(),
    };
    write_file(
        &workspace_path.join("vcfheader.vcf"),
        header.write().as_bytes(),
    )?;

    // `executeImport`: batch after batch, each partition in turn, each sample in sorted order.
    let opened = open_samples(&samples)?;
    let mut counter = 0u64;
    let mut last_stamp = 0u128;
    let last_row = samples.len() as i64 - 1;
    let mut lower = 0usize;
    while lower < samples.len() {
        let upper = (lower + batch_size).min(samples.len());
        for sample in &samples[lower..upper] {
            assert_index_exists(sample)?;
        }
        for interval in &intervals {
            let contig = vid
                .contigs
                .iter()
                .find(|contig| contig.name == interval.contig)
                .expect("the interval was validated against this dictionary");
            let begin = contig.offset + i64::from(interval.start) - 1;
            let end = contig.offset + i64::from(interval.end) - 1;
            let array_name = format!("{}${}${}", interval.contig, interval.start, interval.end);
            let array = workspace_path.join(&array_name);
            let meta = array.join("genomicsdb_meta_dir");
            if lower == 0 {
                std::fs::create_dir_all(&meta).map_err(|error| {
                    limitation(&format!("Creating the array failed ({error})."))
                })?;
                write_file(&array.join(".__consolidation_lock"), &[])?;
                write_file(
                    &array.join("__array_schema.tdb"),
                    &schema.serialize(&array.display().to_string()),
                )?;
                write_file(
                    &meta.join("genomicsdb_column_bounds.json"),
                    gdb::column_bounds_json(begin, end).as_bytes(),
                )?;
            }
            let row_bound = ((upper as i64) - 1).min(last_row);
            write_file(
                &meta.join(format!("genomicsdb_meta_{}.json", uuid_like(&mut counter))),
                gdb::row_bounds_json(0, row_bound).as_bytes(),
            )?;

            let mut rows = Vec::with_capacity(upper - lower);
            for (row, sample) in samples.iter().enumerate().take(upper).skip(lower) {
                rows.push(sample_cells(
                    &schema,
                    &opened[row],
                    sample,
                    interval,
                    row as i64,
                    contig.offset,
                )?);
            }
            let cells =
                gdb::load(rows, samples.len(), begin).map_err(|reason| limitation(&reason))?;
            let files = gdb::fragment_files(&schema, &cells);
            if !files.is_empty() {
                let stamp = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map(|elapsed| elapsed.as_millis())
                    .unwrap_or(0)
                    .max(last_stamp + 1);
                last_stamp = stamp;
                let fragment = array.join(format!(
                    "__{}{}_{stamp}",
                    uuid_like(&mut counter),
                    140_000_000_000_000u64 + u64::from(std::process::id())
                ));
                std::fs::create_dir_all(&fragment).map_err(|error| {
                    limitation(&format!("Creating the fragment failed ({error})."))
                })?;
                for (name, bytes) in files {
                    write_file(&fragment.join(name), &bytes)?;
                }
            }
        }
        lower = upper;
    }

    workspace.finished = true;
    Ok(Some("true".to_string()))
}

/// `TileDBUtils::initialize_workspace`, with the reference's words for what it returns.
fn create_workspace(path: &Path, argument: &str, overwrite: bool) -> Result<(), Thrown> {
    let not_created = || Thrown::user(format!("Error creating GenomicsDB workspace: {argument}"));
    if path.is_file() {
        return Err(Thrown::user(format!(
            "Error creating GenomicsDB workspace: {argument} already exists and is not a directory"
        )));
    }
    if is_tiledb_directory(path, "__tiledb_workspace.tdb") {
        if !overwrite {
            return Err(Thrown::user(format!(
                "Error creating GenomicsDB workspace: {argument} already exists"
            )));
        }
        std::fs::remove_dir_all(path).map_err(|_| not_created())?;
    }
    // `StorageManager::workspace_create`.
    if let Some(parent) = path.parent() {
        for marker in [
            "__tiledb_workspace.tdb",
            "__tiledb_group.tdb",
            "__array_schema.tdb",
            "__metadata_schema.tdb",
        ] {
            if is_tiledb_directory(parent, marker) {
                return Err(not_created());
            }
        }
    }
    std::fs::create_dir(path).map_err(|_| not_created())
}

fn is_tiledb_directory(path: &Path, marker: &str) -> bool {
    path.is_dir() && path.join(marker).is_file()
}

/// What `FeatureDataSource` refuses when it opens an input without an index.
fn assert_index_exists(sample: &Sample) -> Result<(), Thrown> {
    let has = |ending: &str| Path::new(&format!("{}{ending}", sample.path)).exists();
    let indexed = if sample.path.ends_with(".gz") {
        has(".tbi")
    } else {
        has(".idx")
    };
    if !indexed
        && sample
            .index
            .as_deref()
            .is_some_and(|index| Path::new(index).exists())
    {
        return Err(limitation(
            "An input read through an index named in the sample map, with none beside it, has \
             not been measured.",
        ));
    }
    if !indexed {
        let source = format!("file://{}", java_absolute_path(&sample.path));
        return Err(Thrown::user(format!(
            "Failed to create reader from {source} because of the following error:\n\tAn index is \
             required, but none found with file ending .idx, for input source: {source}"
        )));
    }
    Ok(())
}

/// The text of each input, already held open for every batch that needs it.
fn open_samples(samples: &[Sample]) -> Result<Vec<htsjdk_vcf::reader::VcfFile>, Thrown> {
    samples
        .iter()
        .map(|sample| {
            let (text, _) = read_header(&sample.path)?;
            if let Some(key) = raw_sample_value_the_reader_loses(&text) {
                return Err(limitation(&format!(
                    "{} has a value of {key} that reaches the workspace as the raw text of its \
                     sample column, which this port's reader does not keep.",
                    sample.path
                )));
            }
            htsjdk_vcf::reader::read_vcf(&text).map_err(|_| {
                limitation(&format!(
                    "{} has a record this port's reader refuses, and the reference only refuses \
                     it, in its own words, if an interval reaches it.",
                    sample.path
                ))
            })
        })
        .collect()
}

/// The FORMAT key of a value the parsed genotype cannot give back. The stream's sample column is the
/// raw text, so a `.` in a String field arrives as a `.`, where the parsed genotype has nothing, and
/// `AD` or `PL` with a `.` inside the list arrives as that list, where the parsed genotype has none.
fn raw_sample_value_the_reader_loses(text: &str) -> Option<String> {
    let mut strings: Vec<&str> = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##FORMAT=<ID=") {
            if line.contains("Type=String") || line.contains("Type=Character") {
                strings.extend(rest.split(',').next());
            }
        } else if !line.starts_with('#') {
            let mut columns = line.split('\t').skip(8);
            let (Some(format), Some(sample)) = (columns.next(), columns.next()) else {
                continue;
            };
            for (key, value) in format.split(':').zip(sample.split(':')) {
                let lost = match key {
                    "GT" => false,
                    "AD" | "PL" => value != "." && value.split(',').any(|item| item == "."),
                    _ => value == "." && (key == "FT" || strings.contains(&key)),
                };
                if lost {
                    return Some(key.to_string());
                }
            }
        }
    }
    None
}

/// `reader.query(contig, start, end)` through the `NoMnpIterator`, then the conversion.
fn sample_cells(
    schema: &gdb::Schema,
    file: &htsjdk_vcf::reader::VcfFile,
    sample: &Sample,
    interval: &gatk_engine::interval::SimpleInterval,
    row: i64,
    contig_offset: i64,
) -> Result<Vec<gdb::Cell>, Thrown> {
    if file.header.samples.len() != 1 {
        return Err(limitation(&format!(
            "{} does not hold exactly one sample, and a sample map's stream of it is one this \
             port has not measured.",
            sample.path
        )));
    }
    let mut cells = Vec::new();
    for record in &file.records {
        if record.contig != interval.contig
            || record.start > i64::from(interval.end)
            || record.stop < i64::from(interval.start)
        {
            continue;
        }
        if gatk_tools::reference_confidence_merger::is_unmixed_mnp_ignoring_non_ref(record) {
            return Err(Thrown::user(format!(
                "Bad input: GenomicsDBImport does not support GVCFs with MNPs. MNP found at \
                 {}:{:2} in VCF {}",
                record.contig,
                record.start,
                java_absolute_path(&sample.path)
            )));
        }
        cells.push(
            schema
                .cell(record, row, contig_offset)
                .map_err(|reason| limitation(&reason))?,
        );
    }
    Ok(cells)
}
