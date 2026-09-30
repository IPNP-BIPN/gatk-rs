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
