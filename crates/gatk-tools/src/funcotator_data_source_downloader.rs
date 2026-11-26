//! `FuncotatorDataSourceDownloader`, ported from the tool and the constants it is built on
//! (GATK 4.7.0.0; `FuncotatorDataSourceDownloader.java` and `DataSourceUtils.java` are unchanged
//! since 4.6.2.0).
//!
//! A tar.gz fetched from a bucket, checked against a sha256 beside it, and optionally unpacked.
//! The transport is not ported; every decision around it is: which path is fetched, where it
//! lands, what the checksum file is reduced to before it is compared, and what each refusal says.
//!
//! # With no `-O` the copy lands in the working directory
//!
//! ```java
//! if ( outputFile == null ) {
//!     return IOUtils.getPath(dataSourcesPath.getFileName().toString());
//! }
//! ```
//!
//! That is a RELATIVE path, so the destination has nothing to do with where the source was: it is
//! wherever the process happens to be running.
//!
//! # The computed sum is upper case and the expected one is not
//!
//! `DatatypeConverter.printHexBinary` produces upper-case hex, `cleanExpectedSha256SumString`
//! lower-cases what it read from the file, and `isDestFileValid` lower-cases both before
//! comparing. The failure message prints them as they are, so a mismatch reads as an upper-case
//! sum against a lower-case one.
//!
//! # And the corrupt file is left where it landed
//!
//! `validateIntegrity` runs after `initiateCopy` returns, so a file that fails its checksum has
//! already been written. Nothing removes it.
//!
//! # The cleaner truncates twice, in order
//!
//! ```java
//! if ( cleanString.contains(" ") ) {
//!     cleanString = cleanString.substring(0, cleanString.indexOf(" "));
//! }
//! if ( cleanString.contains("\t")) {
//!     cleanString = cleanString.substring(0, cleanString.indexOf("\t"));
//! }
//! ```
//!
//! The tab test runs on the string the space test already shortened, so `<sum> <name>\t<x>` is cut
//! at the space and the tab is never reached.
//!
//! # The startup checks run in a fixed order
//!
//! The data-source check comes before the reference check, and both come before the
//! all-or-nothing test on the two testing arguments. So passing only the sha256 override is
//! reported as a missing DATA SOURCE rather than as a missing pair, which is not what the argument
//! it named was about.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` is a borrowed, read-only string and `String` an owned one; `format!` is
//!   `String.format`, with `{name}` substituting the variable `name` and `{month:02}` meaning
//!   `%02d`.
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`), which is how a nullable Java
//!   argument such as `outputFile` arrives here.
//! - `Result<T, E>` is "a value or an error"; an `Err` is a `UserException` the reference throws.
//! - `match` must cover every case of what it inspects; the compiler refuses one that does not.

use std::path::{Path, PathBuf};

/// `DataSourceUtils.DATA_SOURCES_BUCKET_PATH`.
pub const DATA_SOURCES_BUCKET_PATH: &str = "gs://gcp-public-data--broad-references/funcotator/";
/// `DataSourceUtils.DATA_SOURCES_NAME_PREFIX`.
pub const DATA_SOURCES_NAME_PREFIX: &str = "funcotator_dataSources";
/// `DataSourceUtils.DS_SOMATIC_NAME_MODIFIER`.
pub const DS_SOMATIC_NAME_MODIFIER: &str = "s";
/// `DataSourceUtils.DS_GERMLINE_NAME_MODIFIER`.
pub const DS_GERMLINE_NAME_MODIFIER: &str = "g";
/// `DataSourceUtils.DS_EXTENSION`.
pub const DS_EXTENSION: &str = ".tar.gz";
/// `DataSourceUtils.DS_CHECKSUM_EXTENSION`.
pub const DS_CHECKSUM_EXTENSION: &str = ".sha256";

/// `MAX_MAJOR_VERSION_NUMBER`, `MAX_MINOR_VERSION_NUMBER` and `MAX_DATE`, as
/// `(major, minor, year, month, day)`: `LocalDate.of(2023, Month.SEPTEMBER, 8)` is the date.
pub const MAX_VERSION: (i32, i32, i32, u32, u32) = (1, 8, 2023, 9, 8);
/// `MIN_MAJOR_VERSION_NUMBER`, `MIN_MINOR_VERSION_NUMBER` and `MIN_DATE`, in the same shape:
/// `LocalDate.of(2019, Month.JANUARY, 24)`.
pub const MIN_VERSION: (i32, i32, i32, u32, u32) = (1, 6, 2019, 1, 24);

/// `getNewDataSourceVersionString`: `v%d.%d.hg%d.%d%02d%02d`, the year, the month and the day
/// coming from `date.getYear()`, `date.getMonthValue()` and `date.getDayOfMonth()`.
pub fn new_data_source_version_string(
    major: i32,
    minor: i32,
    reference: i32,
    year: i32,
    month: u32,
    day: u32,
) -> String {
    format!("v{major}.{minor}.hg{reference}.{year}{month:02}{day:02}")
}

/// `getDataSourceVersionString`: the same without the reference, which is the MINIMUM's shape.
pub fn data_source_version_string(
    major: i32,
    minor: i32,
    year: i32,
    month: u32,
    day: u32,
) -> String {
    format!("v{major}.{minor}.{year}{month:02}{day:02}")
}

/// `getDataSourceMaxVersionString(ref)`.
pub fn max_version_string(reference: i32) -> String {
    let (major, minor, year, month, day) = MAX_VERSION;
    new_data_source_version_string(major, minor, reference, year, month, day)
}

/// `getDataSourceMinVersionString`, which carries no reference at all.
pub fn min_version_string() -> String {
    let (major, minor, year, month, day) = MIN_VERSION;
    data_source_version_string(major, minor, year, month, day)
}

/// Which bundle a run is after: `--somatic` or `--germline`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataSourceKind {
    Somatic,
    Germline,
}

impl DataSourceKind {
    /// The name modifier appended to the base url: `DS_SOMATIC_NAME_MODIFIER` or
    /// `DS_GERMLINE_NAME_MODIFIER`.
    pub fn modifier(self) -> &'static str {
        match self {
            DataSourceKind::Somatic => DS_SOMATIC_NAME_MODIFIER,
            DataSourceKind::Germline => DS_GERMLINE_NAME_MODIFIER,
        }
    }
}

/// The base url a bundle is fetched from, which is where the reference number ends up INSIDE the
/// version string rather than beside it.
///
/// `HG38_BASE_URL` is `DATA_SOURCES_BUCKET_PATH + DATA_SOURCES_NAME_PREFIX + "." +
/// getDataSourceMaxVersionString(38)` (each name qualified by `DataSourceUtils.` in the Java).
/// The Java holds one such constant per reference, each followed by the somatic or germline
/// modifier; this builds the same string for either reference number.
pub fn base_url(kind: DataSourceKind, reference: i32) -> String {
    format!(
        "{DATA_SOURCES_BUCKET_PATH}{DATA_SOURCES_NAME_PREFIX}.{}{}",
        max_version_string(reference),
        kind.modifier()
    )
}

/// The `*_GCLOUD_DATASOURCES_PATH` constants: the base url with `DS_EXTENSION`.
pub fn data_sources_path(kind: DataSourceKind, reference: i32) -> String {
    format!("{}{DS_EXTENSION}", base_url(kind, reference))
}

/// The `*_GCLOUD_DATASOURCES_SHA256_PATH` constants: the base url with `DS_CHECKSUM_EXTENSION`.
pub fn checksum_path(kind: DataSourceKind, reference: i32) -> String {
    format!("{}{DS_CHECKSUM_EXTENSION}", base_url(kind, reference))
}

/// `dataSourceDescription`, which is what the run logs before it starts.
///
/// Any reference other than 38 is described as HG19, which is the Java's `else` branch: the two
/// flags are mutually exclusive and one is required. The testing override's `"TESTING"` is not
/// produced here.
pub fn description(kind: DataSourceKind, reference: i32) -> String {
    let reference = if reference == 38 { "HG38" } else { "HG19" };
    let kind = match kind {
        DataSourceKind::Somatic => "Somatic",
        DataSourceKind::Germline => "Germline",
    };
    format!("{reference}_{kind}")
}

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DownloadError {
    /// Neither `--somatic` nor `--germline` nor the testing override.
    NoDataSource,
    /// Neither `--hg38` nor `--hg19` nor the testing override.
    NoReference,
    /// One testing argument without the other.
    IncompleteTestingArguments,
    /// The checksum file held no first line.
    NoChecksum { uri: String },
    /// The checksums disagreed.
    Corrupt { checksum: String, expected: String },
}

impl DownloadError {
    /// Every refusal here is a plain `UserException`.
    pub fn java_class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    /// The exception's message, as `getMessage()` would return it.
    pub fn message(&self) -> String {
        match self {
            DownloadError::NoDataSource => {
                "Must select either somatic or germline datasources.".to_string()
            }
            DownloadError::NoReference => {
                "Must select either HG19 or HG38 datasources.".to_string()
            }
            DownloadError::IncompleteTestingArguments => {
                "Must specify both a test data sources path and a test data sources sha256sum path."
                    .to_string()
            }
            DownloadError::NoChecksum { uri } => {
                format!("Unable to retrieve expected checksum from: {uri}")
            }
            // Two spaces after the exclamation mark, as the reference writes it.
            DownloadError::Corrupt { checksum, expected } => format!(
                "ERROR: downloaded data sources are corrupt!  Unexpected checksum: \
                 {checksum} != {expected}"
            ),
        }
    }
}

/// The arguments `onStartup` looks at. A testing path that is `None` is the Java's `null`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Arguments {
    pub somatic: bool,
    pub germline: bool,
    pub hg38: bool,
    pub hg19: bool,
    pub testing_data_sources_path: Option<String>,
    pub testing_sha256_path: Option<String>,
}

/// `onStartup`, in its own order: the data source, then the reference, then the pair.
///
/// The order is the behaviour. Passing only `--testing-override-path-for-datasources-sha256`
/// fails the FIRST check, so the message names a data source rather than the argument that was
/// actually incomplete.
pub fn startup(arguments: &Arguments) -> Result<(), DownloadError> {
    let testing = arguments.testing_data_sources_path.is_some();
    if !arguments.somatic && !arguments.germline && !testing {
        return Err(DownloadError::NoDataSource);
    }
    if !arguments.hg38 && !arguments.hg19 && !testing {
        return Err(DownloadError::NoReference);
    }
    // `!=` on the two presence flags is the Java's two-sided "one null and the other not" test.
    if testing != arguments.testing_sha256_path.is_some() {
        return Err(DownloadError::IncompleteTestingArguments);
    }
    Ok(())
}

/// `cleanExpectedSha256SumString`: trim, lower-case, cut at the first space, then cut what is left
/// at the first tab.
///
/// Java's `trim()` removes every character at or below U+0020 and Rust's `trim()` removes Unicode
/// white space; the two agree on spaces, tabs and line ends, which is what a checksum line holds.
pub fn clean_expected_sha256(line: &str) -> String {
    let mut clean = line.trim().to_lowercase();
    // `if let Some(i) = ..` runs its block only when `find` found something: `contains` and
    // `indexOf` in one step. The tab test then runs on the already shortened string.
    if let Some(index) = clean.find(' ') {
        clean = clean[..index].to_string();
    }
    if let Some(index) = clean.find('\t') {
        clean = clean[..index].to_string();
    }
    clean
}

/// `readSha256SumFromPath`: the FIRST line of the file, cleaned, or a refusal naming the path as a
/// URI.
///
/// ```java
/// expectedSha256Sum = Files.lines(sha256Path).findFirst().orElse(null);
/// ```
///
/// The caller passes the file's text; an empty file has no first line and is refused, while a
/// file whose first line is blank yields an empty sum, as in the Java.
pub fn expected_sha256(contents: &str, uri: &str) -> Result<String, DownloadError> {
    match contents.lines().next() {
        Some(line) => Ok(clean_expected_sha256(line)),
        None => Err(DownloadError::NoChecksum {
            uri: uri.to_string(),
        }),
    }
}

/// `getOutputLocation`: the given file, or the SOURCE'S FILE NAME as a relative path.
///
/// The `expect` cannot fire for the paths this tool builds, all of which end in a file name.
pub fn output_location(output_file: Option<&Path>, data_sources_path: &Path) -> PathBuf {
    match output_file {
        Some(path) => path.to_path_buf(),
        None => PathBuf::from(
            data_sources_path
                .file_name()
                .expect("a data sources file name"),
        ),
    }
}

/// `isDestFileValid`, which lower-cases both sides even though only one of them can be upper.
///
/// ```java
/// return checksum.toLowerCase().equals(expectedChecksum.toLowerCase());
/// ```
pub fn is_dest_file_valid(checksum: &str, expected: &str) -> bool {
    checksum.to_lowercase() == expected.to_lowercase()
}

/// `validateIntegrity`, which reports the two sums AS THEY ARE rather than as it compared them:
/// `results.getChecksum() + " != " + results.getExpectedChecksum()`.
pub fn validate_integrity(checksum: &str, expected: &str) -> Result<(), DownloadError> {
    if is_dest_file_valid(checksum, expected) {
        Ok(())
    } else {
        Err(DownloadError::Corrupt {
            checksum: checksum.to_string(),
            expected: expected.to_string(),
        })
    }
}

/// `DatatypeConverter.printHexBinary`, which is upper case with no separator.
///
/// `NioFileCopierWithProgressMeter` applies it to `messageDigest.digest()`. `{byte:02X}` prints
/// each byte as two upper-case hex digits.
pub fn print_hex_binary(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02X}")).collect()
}
