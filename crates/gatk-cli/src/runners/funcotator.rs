//! The three Funcotator tools that take a folder of data sources: `FuncotatorDataSourceDownloader`,
//! which fetches one, and `FuncotateSegments` and `Funcotator`, which annotate from one.
//!
//! A module of its own under the runners, because the three share the folder's plumbing and none
//! of it is anybody else's. It reaches the runners' private helpers (the walker startup, the VCF
//! writer's arguments) through `super`.

use super::*;

/// `Path.toUri().toString()` for a local path: `file://` in front of the ABSOLUTE path, the
/// characters a URI cannot carry escaped, and a trailing slash when the path is an existing
/// directory.
pub(super) fn java_file_uri(path: &str) -> String {
    let absolute = java_absolute_path(path);
    let mut out = String::from("file://");
    for byte in absolute.bytes() {
        let keep = byte.is_ascii_alphanumeric() || b"/-_.!~*'()@&=+$,;:".contains(&byte);
        if keep {
            out.push(byte as char);
        } else {
            out.push_str(&format!("%{byte:02X}"));
        }
    }
    if std::path::Path::new(&absolute).is_dir() && !out.ends_with('/') {
        out.push('/');
    }
    out
}

/// `FuncotatorDataSourceDownloader`: a data source archive copied, checked and unpacked.
///
/// The statics are [`gatk_tools::funcotator_data_source_downloader`]; this is `onStartup` and
/// `doWork` in their order:
///
/// * **the three startup refusals**, data source kind first, reference second, the testing pair
///   third;
/// * **the checksum is read before the copy** when `--validate-integrity` is set, so a missing or
///   empty checksum file refuses with nothing written;
/// * **the copier refuses an existing destination before it opens the source**, unless
///   `--overwrite-output-file`, which is the one refusal a bucket path reaches without the network;
/// * **a checksum mismatch leaves the copy on disk**, because validation follows the copy;
/// * **`--extract-after-download` unpacks into the destination's PARENT**, refusing any entry whose
///   path exists unless overwriting.
///
/// The bucket itself (gs://) is not reachable from this port: a run that gets past the existing
/// destination with a bucket source is refused as the port's limitation.
pub fn funcotator_data_source_downloader(parser: &Parser) -> Outcome {
    use gatk_tools::funcotator_data_source_downloader as dl;

    let arguments = dl::Arguments {
        somatic: flag(parser, "somatic"),
        germline: flag(parser, "germline"),
        hg38: flag(parser, "hg38"),
        hg19: flag(parser, "hg19"),
        testing_data_sources_path: argument(parser, "testing-override-path-for-datasources"),
        testing_sha256_path: argument(parser, "testing-override-path-for-datasources-sha256"),
    };
    dl::startup(&arguments).map_err(|error| Thrown::user(error.message()))?;
    let validate = flag(parser, "validate-integrity");
    let overwrite = flag(parser, "overwrite-output-file");
    let extract = flag(parser, "extract-after-download");

    // `doWork`: somatic before germline, and within each HG38 before HG19; the testing override is
    // the fallthrough.
    let bucket = |kind: dl::DataSourceKind| {
        let reference = if arguments.hg38 { 38 } else { 19 };
        (
            dl::data_sources_path(kind, reference),
            dl::checksum_path(kind, reference),
        )
    };
    let (source, checksum_source) = if arguments.somatic {
        bucket(dl::DataSourceKind::Somatic)
    } else if arguments.germline {
        bucket(dl::DataSourceKind::Germline)
    } else {
        (
            arguments
                .testing_data_sources_path
                .clone()
                .expect("startup checked the testing path"),
            arguments
                .testing_sha256_path
                .clone()
                .expect("startup checked the testing checksum"),
        )
    };
    let in_bucket = |path: &str| path.starts_with("gs://");
    let bucket_refusal = || {
        Thrown::non_user(
            PORT_LIMITATION,
            "The Broad bucket (gs://) is not reachable from this port. This message is the port's own and not GATK's.",
        )
    };

    // `getOutputLocation`: the given file, or the source's own file name in the working directory.
    let output = argument(parser, "output");
    let destination = dl::output_location(
        output.as_deref().map(std::path::Path::new),
        std::path::Path::new(&source),
    )
    .display()
    .to_string();
    let destination_uri = java_file_uri(&destination);
    let source_uri = if in_bucket(&source) {
        source.clone()
    } else {
        java_file_uri(&source)
    };

    // `readSha256SumFromPath`, before the copy.
    let expected = if validate {
        if in_bucket(&checksum_source) {
            return Err(bucket_refusal());
        }
        let checksum_uri = java_file_uri(&checksum_source);
        let contents = std::fs::read(&checksum_source).map_err(|_| {
            Thrown::user(format!(
                "Could not read in sha256sum from file: {checksum_uri}"
            ))
        })?;
        let text = String::from_utf8_lossy(&contents);
        Some(
            dl::expected_sha256(&java_first_line(&text), &checksum_uri)
                .map_err(|error| Thrown::user(error.message()))?,
        )
    } else {
        None
    };

    // `initiateCopy`: the existing destination first, then the copy.
    if std::path::Path::new(&destination).exists() && !overwrite {
        return Err(Thrown {
            failure: Failure::User,
            exception:
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile",
            message: Some(format!(
                "Couldn't write file {destination_uri} because Download aborted!  Output data sources file already exists!"
            )),
        });
    }
    if in_bucket(&source) {
        return Err(bucket_refusal());
    }
    let copy_refusal = || {
        Thrown::user(format!(
            "Could not copy file: {source_uri} -> {destination_uri}"
        ))
    };
    // The source is opened before the destination, so a missing source leaves nothing behind.
    let bytes = std::fs::read(&source).map_err(|_| copy_refusal())?;
    std::fs::write(&destination, &bytes).map_err(|_| copy_refusal())?;

    if let Some(expected) = expected {
        use sha2::Digest;
        let checksum = dl::print_hex_binary(&sha2::Sha256::digest(&bytes));
        dl::validate_integrity(&checksum, &expected)
            .map_err(|error| Thrown::user(error.message()))?;
    }

    if extract {
        let parent = std::path::Path::new(&destination)
            .parent()
            .map(|parent| parent.to_path_buf())
            .unwrap_or_default();
        extract_tar_gz(&destination, &parent, overwrite)?;
    }
    Ok(Some("true".to_string()))
}

/// `Files.lines(path).findFirst()`: the first line, which is ABSENT for an empty file and empty
/// for a file that starts with a line break.
fn java_first_line(text: &str) -> String {
    text.split(['\n', '\r'])
        .next()
        .unwrap_or_default()
        .to_string()
}

/// `IOUtils.extractTarGz(tarGz, destination, overwrite)`: every entry resolved against the
/// destination, an existing path refused unless overwriting, a directory created, a file written.
fn extract_tar_gz(
    archive: &str,
    destination: &std::path::Path,
    overwrite: bool,
) -> Result<(), Thrown> {
    use std::io::Read;

    let archive_uri = java_file_uri(archive);
    let refusal = || Thrown::user(format!("Could not extract data from: {archive_uri}"));
    let compressed = std::fs::read(archive).map_err(|_| refusal())?;
    let mut tar = Vec::new();
    flate2::read::MultiGzDecoder::new(&compressed[..])
        .read_to_end(&mut tar)
        .map_err(|_| refusal())?;
    for entry in tar_entries(&tar) {
        let path = normalized_join(destination, &entry.name);
        let shown = path.display().to_string();
        if path.exists() && !overwrite {
            return Err(Thrown::user(format!(
                "Output destination already exists: {}",
                java_file_uri(&shown)
            )));
        }
        if entry.directory {
            std::fs::create_dir_all(&path).map_err(|_| refusal())?;
        } else {
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent).map_err(|_| refusal())?;
            }
            std::fs::write(&path, &entry.data).map_err(|_| refusal())?;
        }
    }
    Ok(())
}

/// One entry of a tar stream.
struct TarEntry {
    name: String,
    directory: bool,
    data: Vec<u8>,
}

/// The directory and regular entries of an uncompressed tar stream, GNU long names followed.
fn tar_entries(tar: &[u8]) -> Vec<TarEntry> {
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
        let mut name = || {
            long_name.take().unwrap_or_else(|| {
                let prefix = field(block, 345, 500);
                let name = field(block, 0, 100);
                if prefix.is_empty() {
                    name
                } else {
                    format!("{prefix}/{name}")
                }
            })
        };
        match block[156] {
            b'L' => long_name = Some(field(data, 0, data.len())),
            b'5' => out.push(TarEntry {
                name: name(),
                directory: true,
                data: Vec::new(),
            }),
            b'0' | 0 => {
                let name = name();
                // A pre-POSIX archive marks a directory by its trailing slash alone.
                let directory = name.ends_with('/');
                out.push(TarEntry {
                    name,
                    directory,
                    data: data.to_vec(),
                });
            }
            _ => {}
        }
        offset = data_start + size.div_ceil(512) * 512;
    }
    out
}

/// `destination.resolve(name).normalize()`: `.` dropped and `..` taken back, lexically.
fn normalized_join(destination: &std::path::Path, name: &str) -> std::path::PathBuf {
    use std::path::Component;
    let mut out = std::path::PathBuf::new();
    for component in destination.join(name).components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                out.pop();
            }
            other => out.push(other.as_os_str()),
        }
    }
    out
}

// ================================================================================================
// The data sources.
// ================================================================================================

use gatk_tools::funcotate_segments as fs;
use gatk_tools::funcotator_engine as fe;

/// `Properties.load` over a config file: `key = value` lines, `#` and `!` comments, the separator
/// the first `=`, `:` or whitespace, and the value's leading whitespace dropped.
fn java_properties(text: &str) -> std::collections::BTreeMap<String, String> {
    let mut out = std::collections::BTreeMap::new();
    for raw in text.lines() {
        let line = raw.trim_start_matches([' ', '\t', '\u{c}']);
        if line.is_empty() || line.starts_with('#') || line.starts_with('!') {
            continue;
        }
        let end = line
            .find(['=', ':', ' ', '\t', '\u{c}'])
            .unwrap_or(line.len());
        let key = &line[..end];
        let mut rest = line[end..].trim_start_matches([' ', '\t', '\u{c}']);
        if let Some(stripped) = rest.strip_prefix(['=', ':']) {
            rest = stripped.trim_start_matches([' ', '\t', '\u{c}']);
        }
        out.insert(key.to_string(), rest.to_string());
    }
    out
}

/// A data source found under a folder: its config's path and properties.
struct FoundSource {
    config_path: std::path::PathBuf,
    properties: std::collections::BTreeMap<String, String>,
}

/// `getAndValidateDataSourcesFromPaths`: the folders in natural order, each one's manifest
/// checked, then each source directory holding the reference version read for its config.
fn find_data_sources(
    directories: &[String],
    reference_version: &str,
) -> Result<Vec<FoundSource>, Thrown> {
    let mut sorted = directories.to_vec();
    sorted.sort();
    let mut found: Vec<FoundSource> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut has_gencode = false;
    for directory in &sorted {
        let path = std::path::Path::new(directory);
        if !path.is_dir() {
            return Err(Thrown::user(format!(
                "ERROR: Given data source path is not a valid directory: {}",
                java_file_uri(directory)
            )));
        }
        // `logDataSourcesInfo`: the version line, read with the new pattern only.
        let manifest = path.join("MANIFEST.txt");
        if let Ok(text) = std::fs::read_to_string(&manifest) {
            let version = text
                .lines()
                .filter(|line| line.starts_with("Version:"))
                .find_map(fs::parse_manifest_version);
            if let Some(version) = &version {
                match fs::validate_version(version) {
                    fs::VersionVerdict::Refused => {
                        return Err(Thrown::user(fs::version_refusal(Some(version))));
                    }
                    fs::VersionVerdict::DateUnrepresentable => {
                        // `LocalDate.of` threw inside the reader's catch-all, which logs the
                        // exception and keeps the folder.
                        eprintln!(
                            "WARN  DataSourceUtils - Could not read MANIFEST.txt: unable to log data sources version information.\n{}",
                            local_date_exception(version.year, version.month, version.day)
                        );
                    }
                    fs::VersionVerdict::Acceptable => {}
                }
            }
        }
        let mut sources: Vec<std::path::PathBuf> = std::fs::read_dir(path)
            .map_err(|_| {
                Thrown::non_user(
                    "org.broadinstitute.hellbender.exceptions.GATKException",
                    format!("Unable to read contents of: {}", java_file_uri(directory)),
                )
            })?
            .filter_map(|entry| entry.ok().map(|entry| entry.path()))
            .filter(|p| p.is_dir())
            .collect();
        sources.sort();
        for source in sources {
            let under = source.join(reference_version);
            if !under.is_dir() {
                continue;
            }
            let under_text = under.display().to_string();
            let mut configs: Vec<std::path::PathBuf> = std::fs::read_dir(&under)
                .map_err(|_| {
                    Thrown::user(format!(
                        "Unable to read contents of: {}",
                        java_file_uri(&under_text)
                    ))
                })?
                .filter_map(|entry| entry.ok().map(|entry| entry.path()))
                .filter(|p| p.is_file() && p.extension().is_some_and(|e| e == "config"))
                .collect();
            configs.sort();
            if configs.len() > 1 {
                return Err(Thrown::user(format!(
                    "ERROR: Directory contains more than one config file: {}",
                    java_file_uri(&under_text)
                )));
            }
            let Some(config) = configs.pop() else {
                return Err(Thrown::user(format!(
                    "ERROR: Directory does not contain a config file: {}",
                    java_file_uri(&under_text)
                )));
            };
            let config_text = config.display().to_string();
            let properties = java_properties(&std::fs::read_to_string(&config).unwrap_or_default());
            let kind =
                fs::check_config(&java_file_uri(&config_text), &properties).map_err(|error| {
                    Thrown {
                        failure: Failure::User,
                        exception:
                            "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
                        message: Some(error.message()),
                    }
                })?;
            let name = properties["name"].clone();
            if names.contains(&name) {
                return Err(Thrown {
                    failure: Failure::User,
                    exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
                    message: Some(format!(
                        "Bad input: ERROR: contains more than one dataset of name: {name} - one is: {}",
                        java_file_uri(&config_text)
                    )),
                });
            }
            names.push(name);
            has_gencode |= kind == fs::SourceType::Gencode;
            found.push(FoundSource {
                config_path: config,
                properties,
            });
        }
    }
    if found.is_empty() {
        return Err(Thrown::user(format!(
            "ERROR: Could not find any data sources for given reference: {reference_version}"
        )));
    }
    if !has_gencode {
        return Err(Thrown::user("ERROR: a Gencode datasource is required!"));
    }
    Ok(found)
}

/// `LocalDate.of(year, month, day)`'s refusal, as log4j prints the exception's first line.
fn local_date_exception(year: i32, month: i32, day: i32) -> String {
    let _ = year;
    if !(1..=12).contains(&month) {
        return format!(
            "java.time.DateTimeException: Invalid value for MonthOfYear (valid values 1 - 12): {month}"
        );
    }
    if !(1..=31).contains(&day) {
        return format!(
            "java.time.DateTimeException: Invalid value for DayOfMonth (valid values 1 - 28/31): {day}"
        );
    }
    const MONTHS: [&str; 12] = [
        "JANUARY",
        "FEBRUARY",
        "MARCH",
        "APRIL",
        "MAY",
        "JUNE",
        "JULY",
        "AUGUST",
        "SEPTEMBER",
        "OCTOBER",
        "NOVEMBER",
        "DECEMBER",
    ];
    if day == 29 {
        format!("java.time.DateTimeException: Invalid date 'February 29' as '{year}' is not a leap year")
    } else {
        let name = MONTHS[(month - 1) as usize];
        let title = format!("{}{}", &name[..1], name[1..].to_lowercase());
        format!("java.time.DateTimeException: Invalid date '{title} {day}'")
    }
}

/// The GENCODE factory one source makes, with the genes its GTF holds.
struct GencodeSource {
    factory: fe::Factory,
    genes: Vec<fe::Gene>,
}

/// The settings every factory shares, from the command line.
struct FactorySettings {
    mode: fe::TranscriptSelectionMode,
    user_transcripts: Vec<String>,
    five_prime_flank: i32,
    three_prime_flank: i32,
    splice_site_window: i32,
    prefer_mane: bool,
    segment_funcotation: bool,
    min_bases_for_segment: i32,
    overrides: Vec<(String, String)>,
    severities: fe::Severities,
}

/// `createDataSourceFuncotationFactoriesForDataSources`, over the GENCODE sources this port reads.
fn create_factories(
    found: &[FoundSource],
    settings: &FactorySettings,
) -> Result<Vec<GencodeSource>, Thrown> {
    let mut out = Vec::new();
    for source in found {
        let kind = fs::SourceType::parse(&source.properties["type"]).expect("a checked type");
        if kind != fs::SourceType::Gencode {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "The data source {} is of type {}, which this port does not read. This message is the port's own and not GATK's.",
                    source.properties["name"], source.properties["type"]
                ),
            ));
        }
        let directory = source
            .config_path
            .parent()
            .expect("a config sits in a directory");
        let resolve = |key: &str| directory.join(source.properties[key].trim());
        let gtf_path = resolve("src_file");
        let fasta_path = resolve("gencode_fasta_path");
        let gtf_text = std::fs::read_to_string(&gtf_path)
            .map_err(|_| Thrown::user(format!("Couldn't read file {}", gtf_path.display())))?;
        let genes = fe::parse_gencode_gtf(&gtf_text).map_err(|error| Thrown::user(error.0))?;
        let fasta_text = std::fs::read_to_string(&fasta_path)
            .map_err(|_| Thrown::user(format!("Couldn't read file {}", fasta_path.display())))?;
        let name = source.properties["name"].clone();
        let version = source.properties["version"].clone();
        let mut factory = fe::Factory {
            name,
            version,
            ncbi_build: source.properties["ncbi_build_version"].clone(),
            mode: settings.mode,
            user_transcripts: settings.user_transcripts.clone(),
            five_prime_flank: settings.five_prime_flank,
            three_prime_flank: settings.three_prime_flank,
            splice_site_window: settings.splice_site_window,
            prefer_mane: settings.prefer_mane,
            segment_funcotation: settings.segment_funcotation,
            min_bases_for_segment: settings.min_bases_for_segment,
            overrides: Vec::new(),
            severities: settings.severities.clone(),
            transcripts: fe::TranscriptFasta::parse(&fasta_text),
        };
        // `initializeAnnotationOverrides`: the overrides this source supports.
        let supported = factory.supported_fields();
        factory.overrides = settings
            .overrides
            .iter()
            .filter(|(key, _)| supported.contains(key))
            .cloned()
            .collect();
        out.push(GencodeSource { factory, genes });
    }
    Ok(out)
}

/// `splitAnnotationArgsIntoMap`: `KEY:VALUE` pairs into a `LinkedHashMap`.
fn annotation_map(values: &[String]) -> Result<Vec<(String, String)>, Thrown> {
    let mut out: Vec<(String, String)> = Vec::new();
    for value in values {
        let parts: Vec<&str> = value.split(':').collect();
        if parts.len() != 2 {
            return Err(Thrown {
                failure: Failure::User,
                exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
                message: Some(format!(
                    "Bad input: Argument annotation incorrectly formatted: {value}"
                )),
            });
        }
        match out.iter_mut().find(|(k, _)| k == parts[0]) {
            Some(slot) => slot.1 = parts[1].to_string(),
            None => out.push((parts[0].to_string(), parts[1].to_string())),
        }
    }
    Ok(out)
}

/// `processTranscriptList`: one value is tried as a file of IDs first.
fn transcript_list(values: &[String]) -> Vec<String> {
    let mut unique: Vec<String> = Vec::new();
    for value in values {
        if !unique.contains(value) {
            unique.push(value.clone());
        }
    }
    if unique.len() == 1 {
        if let Ok(text) = std::fs::read_to_string(&unique[0]) {
            let mut out: Vec<String> = Vec::new();
            for line in text.lines() {
                if !out.contains(&line.to_string()) {
                    out.push(line.to_string());
                }
            }
            return out;
        }
    }
    unique
}

/// `setVariantClassificationCustomSeverity`.
fn custom_severities(path: &str) -> Result<fe::Severities, Thrown> {
    let bytes = std::fs::read(path).map_err(|_| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile",
        message: Some(format!(
            "Couldn't read file. Error was: Custom severity file does not exist: {path}"
        )),
    })?;
    let text = String::from_utf8_lossy(&bytes);
    let malformed = |message: String| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile",
        message: Some(format!("The input file is malformed: {message}")),
    };
    let mut severities = fe::Severities::default();
    let mut line_number = 1;
    for line in text.split('\n') {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            return Err(malformed(format!(
                "{path}:{line_number} has {} fields!  Each TSV line must have 2 fields!",
                fields.len()
            )));
        }
        let Ok(severity) = fields[1].parse::<i32>() else {
            return Err(malformed(format!(
                "{path}:{line_number}: severity is not an integer  ({})!  Custom severities must be integer values!",
                fields[1]
            )));
        };
        let Some(vc) = fe::VariantClassification::from_name(fields[0]) else {
            return Err(malformed(format!(
                "{path}:{line_number}: invalid/unknown variant classification specified (possible typo): {}",
                fields[0]
            )));
        };
        severities.custom.insert(vc, severity);
        line_number += 1;
    }
    Ok(severities)
}

/// `getUnaccountedForAnnotations`: the annotations no factory supports.
fn unaccounted(sources: &[GencodeSource], map: &[(String, String)]) -> Vec<(String, String)> {
    map.iter()
        .filter(|(key, _)| {
            !sources
                .iter()
                .any(|s| s.factory.supported_fields().contains(key))
        })
        .cloned()
        .collect()
}

/// The genes a query over `[start, end]` returns: every gene overlapping it, in file order.
fn genes_over<'a>(sources: &'a GencodeSource, contig: &str, span: fe::Span) -> Vec<&'a fe::Gene> {
    sources
        .genes
        .iter()
        .filter(|g| g.contig == contig && g.span.overlaps(&span))
        .collect()
}

/// The reference as the factory reads it, loaded whole.
struct LoadedReference {
    contigs: std::collections::BTreeMap<String, Vec<u8>>,
}

impl LoadedReference {
    fn open(path: &str) -> Result<LoadedReference, Thrown> {
        let mut source =
            gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let mut contigs = std::collections::BTreeMap::new();
        for (name, length) in source.sequences().to_vec() {
            let bases = source
                .query(&name, 1, length as i32)
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
            contigs.insert(name, bases);
        }
        Ok(LoadedReference { contigs })
    }
}

impl fe::Reference for LoadedReference {
    fn bases(&self, contig: &str, start: i64, end: i64) -> Vec<u8> {
        let Some(bases) = self.contigs.get(contig) else {
            return Vec::new();
        };
        let from = (start.max(1) - 1) as usize;
        let to = (end.max(0) as usize).min(bases.len());
        bases[from.min(to)..to].to_vec()
    }
    fn length(&self, contig: &str) -> i64 {
        self.contigs.get(contig).map_or(0, |b| b.len() as i64)
    }
}

fn engine_thrown(error: fe::EngineError) -> Thrown {
    match error {
        fe::EngineError::Limitation(what) => Thrown::non_user(
            PORT_LIMITATION,
            format!("Funcotator: {what} This message is the port's own and not GATK's."),
        ),
        fe::EngineError::User { class, message } => {
            if class.starts_with("org.broadinstitute.hellbender.exceptions.UserException") {
                Thrown {
                    failure: Failure::User,
                    exception: class,
                    message: Some(message),
                }
            } else {
                Thrown::non_user(class, message)
            }
        }
    }
}

/// The date `MafOutputRenderer.writeHeader` stamps: `yyyymmdd'T'hhmmss`, where `mm` is the MINUTE
/// and `hh` the twelve-hour clock, in UTC as the pinned container runs.
fn maf_date() -> String {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|since| since.as_secs() as i64)
        .unwrap_or(0);
    let days = seconds.div_euclid(86_400);
    let time = seconds.rem_euclid(86_400);
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + if month <= 2 { 1 } else { 0 };
    let hour = match (time / 3600) % 12 {
        0 => 12,
        h => h,
    };
    let minute = (time % 3600) / 60;
    let second = time % 60;
    format!("{year:04}{minute:02}{day:02}T{hour:02}{minute:02}{second:02}")
}

/// The MAF file's text, rendered row by row and headed when the first row arrives.
struct MafOutput {
    renderer: fe::MafRenderer,
    input_lines: Vec<String>,
    info: String,
    text: String,
}

impl MafOutput {
    fn header(&mut self, columns: &[String]) {
        let date = maf_date();
        self.text.push_str(&self.renderer.header(
            columns,
            &self.input_lines,
            &[],
            crate::TOOLKIT_VERSION,
            &date,
            &self.info,
        ));
        self.renderer.header_written = true;
    }

    /// `write`: the count fields added, then a row per alternate and transcript.
    fn write(
        &mut self,
        alternates: &[fe::Allele],
        counts: Vec<fe::Funcotation>,
        map: &mut fe::FuncotationMap,
    ) {
        for transcript in map.transcripts() {
            map.add(&transcript, &counts);
        }
        for alternate in alternates {
            if alternate.0 == "*" {
                continue;
            }
            for transcript in map.transcripts() {
                let row = self.renderer.row(alternate, map.get(&transcript));
                if !self.renderer.header_written {
                    let columns: Vec<String> = row.iter().map(|(k, _)| k.clone()).collect();
                    self.header(&columns);
                }
                self.text.push_str(
                    &row.iter()
                        .map(|(_, v)| v.as_str())
                        .collect::<Vec<_>>()
                        .join("\t"),
                );
                self.text.push('\n');
            }
        }
    }

    fn close(&mut self) -> String {
        if !self.renderer.header_written {
            let row = self.renderer.row(&fe::Allele("AT".to_string()), &[]);
            let columns: Vec<String> = row.iter().map(|(k, _)| k.clone()).collect();
            self.header(&columns);
        }
        self.text.clone()
    }
}

/// `CustomMafFuncotationCreator.createCustomMafCountFields`: the count columns from the tumor's
/// and the normal's genotypes, or their empty defaults when there is no pair.
fn maf_count_fields(
    alternates: &[fe::Allele],
    record: Option<&htsjdk_vcf::variant::VariantContext>,
    pair: Option<(&str, &str)>,
) -> Result<Vec<fe::Funcotation>, Thrown> {
    const NAMES: [&str; 7] = [
        "Match_Norm_Seq_Allele1",
        "Match_Norm_Seq_Allele2",
        "t_alt_count",
        "t_ref_count",
        "n_alt_count",
        "n_ref_count",
        "tumor_f",
    ];
    let mut out = Vec::new();
    for (index, alternate) in alternates.iter().enumerate() {
        let values: Vec<String> = match (record, pair) {
            (Some(record), Some((tumor, normal))) => {
                let tumor_g = record.genotypes.iter().find(|g| g.sample_name == tumor);
                let tumor_ad = tumor_g.and_then(|g| g.ad.clone());
                if let Some(ad) = &tumor_ad {
                    if ad.len() < 2 {
                        return Err(Thrown::user(format!(
                            "Allelic Depth (AD field) for Variant[{}:{}_{}->{}] Does not contain both a REF and an ALT value (only one value is present)!",
                            record.contig,
                            record.start,
                            record.reference().display_string() + "*",
                            record.alternate_alleles().first().map(|a| a.display_string()).unwrap_or_default()
                        )));
                    }
                }
                let tumor_af = tumor_g.and_then(|g| {
                    g.extended
                        .iter()
                        .find(|(k, _)| k == "AF")
                        .map(|(_, v)| v.format().unwrap_or_default())
                });
                let normal_g = record.genotypes.iter().find(|g| g.sample_name == normal);
                let normal_ad = normal_g.and_then(|g| g.ad.clone());
                let (m1, m2) = match normal_g {
                    None => ("__UNKNOWN__".to_string(), "__UNKNOWN__".to_string()),
                    Some(g) if g.alleles.len() < 2 => {
                        ("__UNKNOWN__".to_string(), "__UNKNOWN__".to_string())
                    }
                    Some(g) => (g.alleles[0].base_string(), g.alleles[1].base_string()),
                };
                vec![
                    m1,
                    m2,
                    tumor_ad
                        .as_ref()
                        .map(|ad| ad[index + 1].to_string())
                        .unwrap_or_default(),
                    tumor_ad
                        .as_ref()
                        .map(|ad| ad[0].to_string())
                        .unwrap_or_default(),
                    normal_ad
                        .as_ref()
                        .map(|ad| ad[index + 1].to_string())
                        .unwrap_or_default(),
                    normal_ad
                        .as_ref()
                        .map(|ad| ad[0].to_string())
                        .unwrap_or_default(),
                    tumor_af
                        .map(|af| af.split(',').nth(index).unwrap_or("").to_string())
                        .unwrap_or_default(),
                ]
            }
            _ => {
                let mut v = vec!["__UNKNOWN__".to_string(), "__UNKNOWN__".to_string()];
                v.extend(std::iter::repeat_n(String::new(), NAMES.len() - 2));
                v
            }
        };
        let fields: Vec<(String, String)> =
            NAMES.iter().map(|n| n.to_string()).zip(values).collect();
        out.push(fe::table(&fields, alternate, "MAF_COUNT_OUTPUT"));
    }
    Ok(out)
}

// ================================================================================================
// The two tools.
// ================================================================================================

/// The annotation defaults and overrides, as `KEY:VALUE` pairs in command-line order.
type AnnotationPairs = Vec<(String, String)>;

/// The settings both tools read off their shared argument collection.
fn factory_settings(
    parser: &Parser,
    segments: bool,
) -> Result<(FactorySettings, AnnotationPairs, AnnotationPairs), Thrown> {
    let mode = match scalar(parser, "transcript-selection-mode").as_deref() {
        Some("BEST_EFFECT") => fe::TranscriptSelectionMode::BestEffect,
        Some("ALL") => fe::TranscriptSelectionMode::All,
        _ => fe::TranscriptSelectionMode::Canonical,
    };
    let defaults_raw = arguments(parser, "annotation-default");
    let overrides_raw = arguments(parser, "annotation-override");
    let (user_transcripts, defaults, overrides) = if segments {
        let defaults = annotation_map(&defaults_raw)?;
        let overrides = annotation_map(&overrides_raw)?;
        (
            transcript_list(&arguments(parser, "transcript-list")),
            defaults,
            overrides,
        )
    } else {
        let user = transcript_list(&arguments(parser, "transcript-list"));
        let defaults = annotation_map(&defaults_raw)?;
        let overrides = annotation_map(&overrides_raw)?;
        (user, defaults, overrides)
    };
    let number = |name: &str, default: i32| {
        scalar(parser, name)
            .and_then(|v| v.parse::<i32>().ok())
            .unwrap_or(default)
    };
    Ok((
        FactorySettings {
            mode,
            user_transcripts,
            five_prime_flank: if segments {
                0
            } else {
                number("five-prime-flank-size", 5000)
            },
            three_prime_flank: if segments {
                0
            } else {
                number("three-prime-flank-size", 0)
            },
            splice_site_window: number("splice-site-window-size", 2),
            prefer_mane: flag(parser, "prefer-mane-transcripts"),
            segment_funcotation: segments,
            min_bases_for_segment: number("min-num-bases-for-segment-funcotation", 150),
            overrides: overrides.clone(),
            severities: fe::Severities::default(),
        },
        defaults,
        overrides,
    ))
}

/// The data sources, the factories and the custom severities, in `onTraversalStart`'s order.
fn start_engine(
    parser: &Parser,
    settings: &mut FactorySettings,
) -> Result<Vec<GencodeSource>, Thrown> {
    let reference_version = scalar(parser, "ref-version").unwrap_or_default();
    let directories = arguments(parser, "data-sources-path");
    let found = find_data_sources(&directories, &reference_version)?;
    let sources = create_factories(&found, settings)?;
    let mut sources = sources;
    if let Some(path) = argument(parser, "custom-variant-classification-order") {
        let severities = custom_severities(&path)?;
        settings.severities = severities.clone();
        for source in &mut sources {
            source.factory.severities = severities.clone();
        }
    }
    Ok(sources)
}

fn output_path(parser: &Parser) -> Result<String, Thrown> {
    argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })
}

/// The `##Funcotator Version` line and the `FUNCOTATION` declaration, after the header's own.
fn funcotator_header(
    lines: &[htsjdk_vcf::header::HeaderLine],
    samples: &[String],
    tool_lines: Vec<htsjdk_vcf::header::HeaderLine>,
    info: &str,
    fields: &[String],
) -> htsjdk_vcf::header::VcfHeader {
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    let mut all: Vec<HeaderLine> = lines.to_vec();
    all.extend(tool_lines);
    all.push(HeaderLine::Unstructured {
        key: "Funcotator Version".to_string(),
        value: format!("{} | {info}", crate::TOOLKIT_VERSION),
    });
    all.push(HeaderLine::Compound {
        key: "INFO".to_string(),
        id: "FUNCOTATION".to_string(),
        number: Cardinality::A,
        line_type: LineType::String,
        description: format!(
            "Functional annotation from the Funcotator tool.  Funcotation fields are: {}",
            fields.join("|")
        ),
        extra: Vec::new(),
    });
    let mut unique: Vec<HeaderLine> = Vec::new();
    for line in all {
        if !unique.contains(&line) {
            unique.push(line);
        }
    }
    htsjdk_vcf::header::VcfHeader {
        lines: unique,
        samples: samples.to_vec(),
    }
}

/// `FuncotateSegments`: a segment file annotated from a folder of data sources.
///
/// The folder is [`find_data_sources`]; the annotation is [`gatk_tools::funcotator_engine`]'s
/// segment path, which names the genes a segment overlaps and the exons its two ends fall in, the
/// ends placed on their transcripts as one-base variants carrying the START's reference base.
/// The runner is the `FeatureWalker` around it:
///
/// * **the segments are read in file order**, and `-L` over them is refused because a segment
///   file has no index;
/// * **each segment becomes a symbolic variant**: `<DEL>`, `<INS>`, `<COPY_NEUTRAL>` or `<*>` from
///   its call, its annotations renamed through `--alias-to-key-mapping`;
/// * **SEG output is two files**, the segments and `<output>.gene_list.txt`, and a segment at or
///   under the minimum length is refused by the gene list;
/// * **VCF and MAF output render the same funcotations** through the variant renderers.
pub fn funcotate_segments(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "FuncotateSegments")?;
    let segments_path = argument(parser, "segments").ok_or_else(|| {
        Thrown::command_line("Argument segments was missing: Argument 'segments' is required")
    })?;
    let output = output_path(parser)?;
    let reference_path = argument(parser, "reference").expect("a required argument");
    let segments_text = std::fs::read_to_string(&segments_path).map_err(|_| {
        Thrown::user(gatk_tools::read_walker_refusal::cannot_read(
            &segments_path,
            false,
        ))
    })?;
    let reference_header = reference_dictionary(parser)?.unwrap_or_default();
    if interval_arguments(parser, &reference_header)?.is_some() {
        return Err(Thrown::user(format!(
            "Input {segments_path} must support random access to enable traversal by intervals. If it's a file, please index it using the bundled tool IndexFeatureFile"
        )));
    }

    // `onTraversalStart`.
    if scalar(parser, "transcript-selection-mode").as_deref() == Some("ALL") {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "Cannot funcotate segments with the ALL transcript selection mode.  Please select another mode.",
        ));
    }
    let (mut settings, defaults, overrides) = factory_settings(parser, true)?;
    let sources = start_engine(parser, &mut settings)?;
    let mapping_raw = {
        let given = arguments(parser, "alias-to-key-mapping");
        let set = parser
            .definitions()
            .iter()
            .any(|d| d.long_name() == "alias-to-key-mapping" && d.has_been_set());
        if !set {
            vec![
                "MEAN_LOG2_COPY_RATIO:Segment_Mean".to_string(),
                "CALL:Segment_Call".to_string(),
                "sample:Sample".to_string(),
                "sample_id:Sample".to_string(),
                "NUM_POINTS_COPY_RATIO:Num_Probes".to_string(),
            ]
        } else {
            given
        }
    };
    let mapping = annotation_map(&mapping_raw)?;
    let default_unaccounted = unaccounted(&sources, &defaults);
    let override_unaccounted = unaccounted(&sources, &overrides);
    let excluded = arguments(parser, "exclude-field");
    let format = scalar(parser, "output-file-format").unwrap_or_default();
    let info = sources
        .iter()
        .map(|s| s.factory.info_string())
        .collect::<Vec<_>>()
        .join(" | ");
    let reference = LoadedReference::open(&reference_path)?;

    let mut tsv = fe::SimpleTsv::new(
        fe::SEG_ALIASES,
        &default_unaccounted,
        &override_unaccounted,
        &excluded,
        true,
    );
    let mut genes = fe::GeneList::new(
        &default_unaccounted,
        &override_unaccounted,
        &excluded,
        settings.min_bases_for_segment,
    );
    let mut maf = MafOutput {
        renderer: fe::MafRenderer::new(
            &sources
                .iter()
                .flat_map(|s| s.factory.supported_fields())
                .collect::<Vec<_>>(),
            &default_unaccounted,
            &override_unaccounted,
            &excluded,
            &scalar(parser, "ref-version").unwrap_or_default(),
        ),
        input_lines: vec!["fileformat=VCFv4.2".to_string()],
        info: info.clone(),
        text: String::new(),
    };
    let mut records: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    let manual: Vec<(String, String)> = default_unaccounted
        .iter()
        .chain(override_unaccounted.iter())
        .fold(Vec::new(), |mut acc, (k, v)| {
            match acc
                .iter_mut()
                .find(|(key, _): &&mut (String, String)| key == k)
            {
                Some(slot) => slot.1 = v.clone(),
                None => acc.push((k.clone(), v.clone())),
            }
            acc
        });
    let included: Vec<String> = sources
        .iter()
        .flat_map(|s| s.factory.supported_fields())
        .chain(manual.iter().map(|(k, _)| k.clone()))
        .filter(|f| !excluded.contains(f))
        .collect();

    let parsed = parse_segments(&segments_text, &segments_path)?;
    let mut failure: Option<Thrown> = None;
    for segment in &parsed {
        let base = fe::Reference::bases(
            &reference,
            &segment.contig,
            segment.start as i64,
            segment.start as i64,
        );
        let reference_allele = fe::Allele(String::from_utf8_lossy(&base).to_string());
        let alternate =
            fe::Allele(fs::allele_of(fs::call_of(&segment.annotations_map())).to_string());
        // The attributes: END, then every annotation renamed through the mapping.
        let mut attributes: Vec<(String, String)> =
            vec![("END".to_string(), segment.end.to_string())];
        for (key, value) in &segment.annotations {
            let renamed = mapping
                .iter()
                .find(|(alias, _)| alias == key)
                .map(|(_, to)| to.clone())
                .unwrap_or_else(|| key.clone());
            attributes.retain(|(k, _)| *k != renamed);
            attributes.push((renamed, value.clone()));
        }
        let variant = fe::Variant {
            contig: segment.contig.clone(),
            start: segment.start,
            end: segment.end,
            reference: reference_allele.clone(),
            alternates: vec![alternate.clone()],
        };
        let source = &sources[0];
        let genes_here = genes_over(
            source,
            &segment.contig,
            fe::Span::new(segment.start, segment.end),
        );
        let funcotations =
            match source
                .factory
                .create_funcotations(&variant, &reference, &genes_here)
            {
                Ok(list) => list,
                Err(error) => {
                    failure = Some(engine_thrown(error));
                    break;
                }
            };
        let mut map = fe::FuncotationMap::no_transcript(&funcotations);
        let attribute = |name: &str| {
            attributes
                .iter()
                .find(|(k, _)| k == name)
                .map(|(_, v)| v.clone())
                .unwrap_or_default()
        };
        for transcript in map.transcripts() {
            let fields: Vec<(String, String)> = [
                "Segment_Mean",
                "Num_Probes",
                "Segment_Call",
                "Sample",
                "build",
            ]
            .iter()
            .map(|name| (name.to_string(), attribute(name)))
            .collect();
            map.add(
                &transcript,
                &[fe::table(&fields, &alternate, "FUNCOTATE_SEGMENTS")],
            );
        }
        let result: Result<(), Thrown> = match format.as_str() {
            "SEG" => (|| {
                tsv.write(&segment.contig, segment.start, segment.end, &mut map)
                    .map_err(engine_thrown)?;
                let vc_string = segment_vc_string(&variant, &attributes);
                genes
                    .write(&variant, &vc_string, &map)
                    .map_err(engine_thrown)
            })(),
            "MAF" => {
                let counts = maf_count_fields(&variant.alternates, None, None)?;
                maf.write(&variant.alternates, counts, &mut map);
                Ok(())
            }
            _ => {
                let value =
                    fe::vcf_funcotation(&variant.alternates, None, &map, &manual, &included);
                let mut record = htsjdk_vcf::variant::VariantContext::new(
                    &segment.contig,
                    i64::from(segment.start),
                    vec![
                        htsjdk_vcf::allele::Allele::from_str(&reference_allele.0, true)
                            .map_err(|e| Thrown::user(format!("{e:?}")))?,
                        htsjdk_vcf::allele::Allele::from_str(&alternate.0, false)
                            .map_err(|e| Thrown::user(format!("{e:?}")))?,
                    ],
                );
                record.stop = i64::from(segment.end);
                let mut pairs = attributes.clone();
                pairs.push(("FUNCOTATION".to_string(), value));
                let order = segment_attribute_order(&segment.annotations, &mapping, &attributes)?;
                record.attributes = order
                    .iter()
                    .map(|key| {
                        let value = pairs
                            .iter()
                            .find(|(k, _)| k == key)
                            .map(|(_, v)| v.clone())
                            .unwrap_or_default();
                        let value = if key == "END" {
                            htsjdk_vcf::variant::Value::Int(value.parse().unwrap_or(0))
                        } else {
                            htsjdk_vcf::variant::Value::Str(value)
                        };
                        (key.clone(), value)
                    })
                    .collect();
                records.push(record);
                Ok(())
            }
        };
        if let Err(error) = result {
            failure = Some(error);
            break;
        }
    }
    if let Some(error) = failure {
        return Err(error);
    }

    // `closeTool`.
    match format.as_str() {
        "SEG" => {
            let segments_text = tsv.close().map_err(engine_thrown)?;
            write_file(&output, segments_text.as_bytes())?;
            let genes_text = genes.close().map_err(engine_thrown)?;
            write_file(
                &format!("{}.gene_list.txt", java_absolute_path(&output)),
                genes_text.as_bytes(),
            )?;
        }
        "MAF" => {
            let text = maf.close();
            write_file(&output, text.as_bytes())?;
        }
        _ => {
            let header = funcotator_header(
                &[],
                &[],
                default_tool_vcf_header_lines(parser, "FuncotateSegments"),
                &info,
                &included,
            );
            let mut header = header;
            apply_sites_only(parser, &mut header, &mut records);
            let rendered = write_vcf_honouring_lenient(parser, &header, &records)?;
            write_variant_output(parser, &output, &rendered)?;
        }
    }
    Ok(Some("true".to_string()))
}

/// The order the VCF encoder walks a segment's attributes in, which is a `HashMap`'s three times
/// over: the converter's builder (END, then the columns), the `Collectors.toMap` that renames them,
/// and the renderer's builder, which copies that map into one sized for it before it adds
/// `FUNCOTATION`. The encoder refuses the first undeclared key it meets in that order.
fn segment_attribute_order(
    annotations: &[(String, String)],
    mapping: &[(String, String)],
    attributes: &[(String, String)],
) -> Result<Vec<String>, Thrown> {
    use gatk_engine::java_hash::JavaHashMap;
    let crowded = |error| Thrown::non_user(PORT_LIMITATION, format!("{error:?}"));
    let mut converted: JavaHashMap<String, ()> = JavaHashMap::new();
    converted.insert("END".to_string(), ());
    for (key, _) in annotations {
        converted.insert(key.clone(), ());
    }
    converted.check().map_err(crowded)?;
    let mut renamed: JavaHashMap<String, ()> = JavaHashMap::new();
    for key in converted.keys() {
        let to = mapping
            .iter()
            .find(|(alias, _)| alias == key)
            .map(|(_, to)| to.clone())
            .unwrap_or_else(|| key.clone());
        renamed.insert(to, ());
    }
    renamed.check().map_err(crowded)?;
    // `new HashMap<>(map)`: the table sized for the map's own size over the load factor.
    let size = attributes.len();
    let mut copy: JavaHashMap<String, ()> =
        JavaHashMap::with_capacity((size as f32 / 0.75f32 + 1.0f32) as usize);
    for key in renamed.keys() {
        copy.insert(key.clone(), ());
    }
    copy.insert("FUNCOTATION".to_string(), ());
    copy.check().map_err(crowded)?;
    Ok(copy.keys().cloned().collect())
}

/// `VariantContext.toString()` for a segment: no source, no quality, sorted attributes.
fn segment_vc_string(variant: &fe::Variant, attributes: &[(String, String)]) -> String {
    let mut sorted = attributes.to_vec();
    sorted.sort();
    let position = if variant.start == variant.end {
        format!("{}:{}", variant.contig, variant.start)
    } else {
        format!("{}:{}-{}", variant.contig, variant.start, variant.end)
    };
    let mut alternates: Vec<String> = variant.alternates.iter().map(|a| a.0.clone()).collect();
    alternates.sort();
    let alleles = std::iter::once(format!("{}*", variant.reference.0))
        .chain(alternates)
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "[VC null @ {position} Q. of type=SYMBOLIC alleles=[{alleles}] attr={{{}}} GT=[] filters=",
        sorted
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// One segment of the file: its locus and its annotations, in column order.
struct Segment {
    contig: String,
    start: i32,
    end: i32,
    annotations: Vec<(String, String)>,
}

impl Segment {
    fn annotations_map(&self) -> std::collections::BTreeMap<String, String> {
        self.annotations.iter().cloned().collect()
    }
}

/// `AnnotatedIntervalCodec` over a segment file: `@` lines skipped, the first other line the
/// columns, the locus columns found by the default config's names.
fn parse_segments(text: &str, path: &str) -> Result<Vec<Segment>, Thrown> {
    const CONTIG: [&str; 12] = [
        "CONTIG",
        "contig",
        "Chromosome",
        "chrom",
        "chromosome",
        "Chrom",
        "seqname",
        "seqnames",
        "CHROM",
        "target_contig",
        "segment_contig",
        "chr",
    ];
    const START: [&str; 14] = [
        "START",
        "start",
        "Start",
        "Start_Position",
        "start_position",
        "chromStart",
        "segment_start",
        "Start_position",
        "target_start",
        "Position",
        "position",
        "pos",
        "POS",
        "segment_start",
    ];
    const END: [&str; 16] = [
        "END",
        "end",
        "End",
        "End_Position",
        "end_position",
        "chromEnd",
        "segment_end",
        "End_position",
        "target_end",
        "stop",
        "Stop",
        "Position",
        "position",
        "pos",
        "POS",
        "segment_end",
    ];
    let mut lines = text.lines().filter(|line| !line.starts_with('@'));
    let Some(header) = lines.next() else {
        return Ok(Vec::new());
    };
    let columns: Vec<&str> = header.split('\t').collect();
    let find = |names: &[&str]| {
        names
            .iter()
            .find_map(|n| columns.iter().position(|c| c == n))
    };
    let (Some(contig), Some(start), Some(end)) = (find(&CONTIG), find(&START), find(&END)) else {
        return Err(Thrown::user(format!(
            "The segment file {path} has no contig, start or end column this port recognises."
        )));
    };
    let mut out = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let number = |i: usize| {
            fields
                .get(i)
                .and_then(|v| v.parse::<i32>().ok())
                .unwrap_or(0)
        };
        out.push(Segment {
            contig: fields.get(contig).unwrap_or(&"").to_string(),
            start: number(start),
            end: number(end),
            annotations: columns
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != contig && *i != start && *i != end)
                .map(|(i, c)| (c.to_string(), fields.get(i).unwrap_or(&"").to_string()))
                .collect(),
        });
    }
    Ok(out)
}

/// `Funcotator`: a VCF annotated from a folder of data sources.
///
/// The folder is [`find_data_sources`]; the annotation is [`gatk_tools::funcotator_engine`]. The
/// runner is `onTraversalStart` and the `VariantWalker` around it:
///
/// * **the reference's dictionary must be a superset of the variants'**, unless validation is off;
/// * **SEG output is refused**, before the data sources are read;
/// * **a VCF that already carries a `Funcotator Version` line is refused** unless reannotating;
/// * **each record's funcotations are its GENCODE ones by transcript, then its own INFO fields**,
///   which the VCF rendering leaves out and the MAF rendering carries as columns;
/// * **filtered records are dropped** under `--remove-filtered-variants`.
pub fn funcotator(parser: &Parser) -> Outcome {
    let output = output_path(parser)?;
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "Funcotator")?;
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let kept = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;

    // `onTraversalStart`.
    if !flag(parser, "disable-sequence-dictionary-validation") {
        let reference = reference_dictionary(parser)?.unwrap_or_default();
        let variants = vcf_dictionary(&text);
        gatk_tools::sequence_dictionary::validate(
            "Reference",
            &reference.sequences,
            "Driving Variants",
            &variants.sequences,
            true,
            false,
        )
        .map_err(|refusal| Thrown {
            failure: Failure::User,
            exception: refusal.java_class(),
            message: Some(refusal.message()),
        })?;
    }
    let format = scalar(parser, "output-file-format").unwrap_or_default();
    if format == "SEG" {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "This tool does not support segment output.  Please see FuncotateSegments.",
        ));
    }
    let (mut settings, defaults, overrides) = factory_settings(parser, false)?;
    let already = file.header.lines.iter().any(|line| {
        matches!(line, htsjdk_vcf::header::HeaderLine::Unstructured { key, .. } if key == "Funcotator Version")
    });
    if !flag(parser, "reannotate-vcf") && already {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            message: Some(format!(
                "Bad input: Given VCF {input} has already been annotated!"
            )),
        });
    }
    let sources = start_engine(parser, &mut settings)?;
    let default_unaccounted = unaccounted(&sources, &defaults);
    let override_unaccounted = unaccounted(&sources, &overrides);
    let excluded = arguments(parser, "exclude-field");
    let info = sources
        .iter()
        .map(|s| s.factory.info_string())
        .collect::<Vec<_>>()
        .join(" | ");
    let manual: Vec<(String, String)> = default_unaccounted
        .iter()
        .chain(override_unaccounted.iter())
        .fold(Vec::new(), |mut acc, (k, v)| {
            match acc
                .iter_mut()
                .find(|(key, _): &&mut (String, String)| key == k)
            {
                Some(slot) => slot.1 = v.clone(),
                None => acc.push((k.clone(), v.clone())),
            }
            acc
        });
    let included: Vec<String> = sources
        .iter()
        .flat_map(|s| s.factory.supported_fields())
        .chain(manual.iter().map(|(k, _)| k.clone()))
        .filter(|f| !excluded.contains(f))
        .collect();
    // `SamplePairExtractor`: the one sample as the tumor with no normal, or TUMOR/NORMAL names.
    let samples = &file.header.samples;
    let pairs: Vec<(String, String)> = if samples.len() == 1 {
        vec![(samples[0].clone(), String::new())]
    } else {
        let tumors: Vec<&String> = samples
            .iter()
            .filter(|s| ["TUMOR", "CASE", "MET"].contains(&s.as_str()))
            .collect();
        let normals: Vec<&String> = samples
            .iter()
            .filter(|s| ["NORMAL", "CONTROL"].contains(&s.as_str()))
            .collect();
        tumors
            .iter()
            .flat_map(|t| normals.iter().map(move |n| ((*t).clone(), (*n).clone())))
            .collect()
    };
    if format == "MAF" && pairs.len() > 1 {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            message: Some(format!(
                "Bad input: Input files with more than one tumor normal pair are currently not supported.  Found: {}",
                pairs.iter().map(|(t, n)| format!("TumorNormalPair{{tumor='{t}', normal='{n}'}}")).collect::<Vec<_>>().join("; ")
            )),
        });
    }
    let input_lines: Vec<String> = text
        .lines()
        .take_while(|line| line.starts_with("##"))
        .map(|line| line[2..].to_string())
        .collect();
    let mut maf = MafOutput {
        renderer: fe::MafRenderer::new(
            &sources
                .iter()
                .flat_map(|s| s.factory.supported_fields())
                .collect::<Vec<_>>(),
            &default_unaccounted,
            &override_unaccounted,
            &excluded,
            &scalar(parser, "ref-version").unwrap_or_default(),
        ),
        input_lines,
        info: info.clone(),
        text: String::new(),
    };
    let reference =
        LoadedReference::open(&argument(parser, "reference").expect("a required argument"))?;
    // The INFO fields the header declares, in its order, which the input funcotations carry.
    let info_fields: Vec<String> = file
        .header
        .lines
        .iter()
        .filter_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Compound { key, id, .. } if key == "INFO" => {
                Some(id.clone())
            }
            _ => None,
        })
        .collect();

    let remove_filtered = flag(parser, "remove-filtered-variants");
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for record in kept {
        if remove_filtered && record.filters.as_ref().is_some_and(|f| !f.is_empty()) {
            continue;
        }
        let variant = fe::Variant {
            contig: record.contig.clone(),
            start: record.start as i32,
            end: record.stop as i32,
            reference: fe::Allele(record.reference().display_string()),
            alternates: record
                .alternate_alleles()
                .iter()
                .map(|a| fe::Allele(a.display_string()))
                .collect(),
        };
        let mut all = Vec::new();
        for source in &sources {
            let span = source
                .factory
                .query_span(fe::Span::new(variant.start, variant.end));
            let genes = genes_over(source, &variant.contig, span);
            all.extend(
                source
                    .factory
                    .create_funcotations(&variant, &reference, &genes)
                    .map_err(engine_thrown)?,
            );
        }
        let mut map = fe::FuncotationMap::from_gencode(&all);
        // `createFuncotations(vc, inputMetadata, INPUT_VCF)`: every attribute must be declared.
        let undeclared: Vec<String> = record
            .attributes
            .iter()
            .map(|(k, _)| k.clone())
            .filter(|k| !info_fields.contains(k))
            .collect();
        if !undeclared.is_empty() {
            return Err(Thrown {
                failure: Failure::User,
                exception: "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile",
                message: Some(format!(
                    "The input file is malformed: Not all attributes in the variant context appear in the metadata: {} .... Please add these attributes to the input metadata (e.g. VCF Header).",
                    undeclared.join(", ")
                )),
            });
        }
        let input_fields: Vec<(String, String)> = info_fields
            .iter()
            .map(|field| {
                let value = record
                    .attributes
                    .iter()
                    .find(|(k, _)| k == field)
                    .map(|(_, v)| match v {
                        htsjdk_vcf::variant::Value::List(items) => items
                            .iter()
                            .map(|item| item.format().unwrap_or_default())
                            .collect::<Vec<_>>()
                            .join(","),
                        other => other.format().unwrap_or_default(),
                    })
                    .unwrap_or_default();
                (field.clone(), value)
            })
            .collect();
        for transcript in map.transcripts() {
            let list: Vec<fe::Funcotation> = variant
                .alternates
                .iter()
                .map(|a| fe::table(&input_fields, a, "INPUT_VCF"))
                .collect();
            map.add(&transcript, &list);
        }
        if format == "MAF" {
            let pair = pairs.first().map(|(t, n)| (t.as_str(), n.as_str()));
            let counts = maf_count_fields(&variant.alternates, Some(record), pair)?;
            maf.write(&variant.alternates, counts, &mut map);
        } else {
            let existing = record
                .attributes
                .iter()
                .find(|(k, _)| k == "FUNCOTATION")
                .and_then(|(_, v)| v.format());
            let value = fe::vcf_funcotation(
                &variant.alternates,
                existing.as_deref(),
                &map,
                &manual,
                &included,
            );
            let mut out = record.clone();
            out.attributes.retain(|(k, _)| k != "FUNCOTATION");
            out.attributes.push((
                "FUNCOTATION".to_string(),
                htsjdk_vcf::variant::Value::Str(value),
            ));
            written.push(out);
        }
    }

    if format == "MAF" {
        let text = maf.close();
        write_file(&output, text.as_bytes())?;
    } else {
        let mut header = funcotator_header(
            &file.header.lines,
            &file.header.samples,
            default_tool_vcf_header_lines(parser, "Funcotator"),
            &info,
            &included,
        );
        let keep = variant_output_filter(parser, intervals.as_deref())?;
        written.retain(|record| keep(record));
        apply_sites_only(parser, &mut header, &mut written);
        let rendered = write_vcf_honouring_lenient(parser, &header, &written)?;
        write_variant_output(parser, &output, &rendered)?;
    }
    Ok(Some("true".to_string()))
}
