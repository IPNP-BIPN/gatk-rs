//! `FuncotateSegments`: how a copy-number segment file is annotated from a folder of data sources.
//!
//! The folder's SHAPE is as much of the tool as the annotation is, so that is what is ported: the
//! manifest's version gate, the three-level directory walk, the config keys each source type
//! requires, and the two sanity checks at the end. Beside it sits the segment side: which calls
//! become which symbolic allele, and which variant contexts count as segments at all.
//!
//! Reading a GTF, reading a reference and the annotation itself are not ported.
//!
//! Ported from `org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils`,
//! `org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions.DataSourceType`,
//! `AnnotatedIntervalToSegmentVariantContextConverter` and `FuncotatorUtils` (both in
//! `org.broadinstitute.hellbender.tools.funcotator`) in GATK 4.6.2.0. Of these, only
//! `FuncotatorUtils.java` differs in 4.7.0.0, in `getTranscriptIdWithoutVersionNumber`, which this
//! file does not use; `GeneListOutputRenderer.java`, quoted below, is unchanged too.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`). In a function returning an
//!   `Option`, a trailing `?` returns `None` at once when the step yields nothing: a regular
//!   expression that fails to match.
//! - `Result<T, E>` is "a value or an error", here a `UserException` the reference would throw.
//! - `let Some(x) = .. else { .. }` takes the value out or runs the `else` block, which must
//!   leave the function.
//! - `&str` is a borrowed string; `BTreeMap` is Java's `TreeMap`, so its keys iterate sorted.
//! - `match` must cover every case; the compiler refuses one that does not.

use std::collections::BTreeMap;

// ================================================================================================
// The manifest's version.
// ================================================================================================

/// `MIN_MAJOR_VERSION_NUMBER` and its four companions, and the two dates they pair with.
pub const MIN_MAJOR_VERSION: i32 = 1;
/// `MIN_MINOR_VERSION_NUMBER`.
pub const MIN_MINOR_VERSION: i32 = 6;
/// `MIN_DATE`, `LocalDate.of(2019, Month.JANUARY, 24)`, as `(year, month, day)`.
pub const MIN_DATE: (i32, u32, u32) = (2019, 1, 24);
/// `MAX_MAJOR_VERSION_NUMBER`.
pub const MAX_MAJOR_VERSION: i32 = 1;
/// `MAX_MINOR_VERSION_NUMBER`.
pub const MAX_MINOR_VERSION: i32 = 8;
/// `MAX_DATE`, `LocalDate.of(2023, Month.SEPTEMBER, 8)`. A tuple compares field by field, so
/// `(year, month, day)` orders the same way a `LocalDate` does.
pub const MAX_DATE: (i32, u32, u32) = (2023, 9, 8);

/// The two version strings the refusal quotes, which are built from the constants above and not
/// from anything the folder holds.
///
/// `CURRENT_MINIMUM_DATA_SOURCE_VERSION`, `getDataSourceMinVersionString()`.
pub const MINIMUM_VERSION_STRING: &str = "v1.6.20190124";
/// `CURRENT_MAXIMUM_DATA_SOURCE_VERSION`, `getDataSourceMaxVersionString(38)`: hg38 whatever the
/// run's reference.
pub const MAXIMUM_VERSION_STRING: &str = "v1.8.hg38.20230908";

/// `DATA_SOURCES_FTP_PATH`.
const FTP_PATH: &str = "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/";
/// `DATA_SOURCES_BUCKET_PATH`.
const BUCKET_PATH: &str = "gs://gcp-public-data--broad-references/funcotator/";

/// What `NEW_VERSION_PATTERN` pulls out of a `Version:` line, ALREADY MISREAD.
///
/// The pattern has seven groups, because `hg(\d+)` sits between the minor version and the date.
/// The reader takes groups three, four and five as the year, the month and the day, so the `hg`
/// number lands in the year, the year lands in the month, the month lands in the day, and the day
/// becomes the decorator. Every field below is what the reader believes, not what was written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ManifestVersion {
    pub major: i32,
    pub minor: i32,
    /// The `hg` number, believed to be a year.
    pub year: i32,
    /// The written year, believed to be a month.
    pub month: i32,
    /// The written month, believed to be a day.
    pub day: i32,
    /// The written day, believed to be a decorator.
    pub decorator: String,
}

impl ManifestVersion {
    /// `versionMajor + "." + versionMinor + "." + versionYear + versionMonth + versionDay +
    /// versionDecorator`, with no separator after the minor version.
    ///
    /// This is the string the refusal prints back as "Yours", so a manifest reading
    /// `1.2.hg38.20150101` is quoted as `1.2.382015101`.
    pub fn display(&self) -> String {
        format!(
            "{}.{}.{}{}{}{}",
            self.major, self.minor, self.year, self.month, self.day, self.decorator
        )
    }
}

/// `Version:\s+(\d+)\.(\d+)\.hg(\d+)\.(\d\d\d\d)(\d\d)(\d\d)(.*)`, applied to a whole line.
///
/// A line that does not match yields nothing, which is a warning rather than a refusal: the folder
/// is then treated as one with no version at all.
///
/// Written out by hand rather than with a regular expression. Java's `\s` matches six white-space
/// characters; this accepts only the space and the tab, which is what a manifest line carries.
/// A number too large for an `int` makes `Integer.valueOf` throw inside the reader's
/// `catch (final Exception ex)`, and yields nothing here; both leave the folder accepted.
pub fn parse_manifest_version(line: &str) -> Option<ManifestVersion> {
    // `strip_prefix` answers the rest of the string after the prefix, or `None` when the string
    // does not start with it; `?` turns that `None` into "no match".
    let rest = line.strip_prefix("Version:")?;
    let rest = rest.strip_prefix(|c: char| c == ' ' || c == '\t')?;
    let rest = rest.trim_start_matches([' ', '\t']);

    let (major, rest) = leading_digits(rest)?;
    let rest = rest.strip_prefix('.')?;
    let (minor, rest) = leading_digits(rest)?;
    let rest = rest.strip_prefix(".hg")?;
    let (hg, rest) = leading_digits(rest)?;
    let rest = rest.strip_prefix('.')?;
    // `(\d\d\d\d)(\d\d)(\d\d)`: eight digits, cut below as four, two and the rest. Group four
    // (the written year) becomes the MONTH and group five the DAY, as the Java reads them.
    if rest.len() < 8 || !rest.as_bytes()[..8].iter().all(u8::is_ascii_digit) {
        return None;
    }
    Some(ManifestVersion {
        major,
        minor,
        year: hg,
        month: rest[..4].parse().ok()?,
        day: rest[4..6].parse().ok()?,
        decorator: rest[6..].to_string(),
    })
}

/// The greedy `(\d+)` at the front, and what follows it.
///
/// Greedy and never backtracking, which is enough here: in every place it is used the next
/// expected character (a dot) is not a digit, so the regular expression could not have matched by
/// taking fewer digits either.
fn leading_digits(text: &str) -> Option<(i32, &str)> {
    let end = text
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(text.len());
    if end == 0 {
        return None;
    }
    Some((text[..end].parse().ok()?, &text[end..]))
}

/// What `validateVersionInformation` can do with a version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VersionVerdict {
    /// The checks ran and passed.
    Acceptable,
    /// The checks ran and one of them failed.
    Refused,
    /// `LocalDate.of` threw before any date check could run.
    ///
    /// The whole read sits in a `catch (final Exception ex)`, so the throw is SWALLOWED and the
    /// acceptability flag keeps the `true` it was initialised with. A folder that lands here is
    /// therefore accepted, which is why [`acceptable`] treats it as such.
    DateUnrepresentable,
}

/// `validateVersionInformation`: the major and minor checks, then the date built from the misread
/// fields.
///
/// Because the year is the `hg` number and the month is the written year, a well-formed new-style
/// manifest ALWAYS lands on [`VersionVerdict::DateUnrepresentable`]: no real year is a valid month.
/// The date range is therefore unreachable for the manifests it was written for, and only the
/// major and minor numbers can turn a folder away.
pub fn validate_version(version: &ManifestVersion) -> VersionVerdict {
    if version.major < MIN_MAJOR_VERSION || version.major > MAX_MAJOR_VERSION {
        return VersionVerdict::Refused;
    }
    if version.major == MIN_MAJOR_VERSION && version.minor < MIN_MINOR_VERSION {
        return VersionVerdict::Refused;
    }
    if version.major == MAX_MAJOR_VERSION && version.minor > MAX_MINOR_VERSION {
        return VersionVerdict::Refused;
    }
    let Some(date) = as_date(version.year, version.month, version.day) else {
        return VersionVerdict::DateUnrepresentable;
    };
    if date < MIN_DATE || date > MAX_DATE {
        VersionVerdict::Refused
    } else {
        VersionVerdict::Acceptable
    }
}

/// `LocalDate.of(year, month, day)`, which refuses a month outside one to twelve and a day outside
/// the month's own length.
///
/// The leap-year rule is the proleptic Gregorian one `LocalDate` uses. The year's own range check
/// (`LocalDate` refuses years beyond plus or minus 999999999) is not written, because the year here
/// is the `hg` number.
fn as_date(year: i32, month: i32, day: i32) -> Option<(i32, u32, u32)> {
    if !(1..=12).contains(&month) || day < 1 {
        return None;
    }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let length = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        _ if leap => 29,
        _ => 28,
    };
    if day > length {
        return None;
    }
    Some((year, month as u32, day as u32))
}

/// Whether a folder's manifest lets its sources be read at all.
///
/// `None` stands for a manifest that could not be read or did not parse: the version check never
/// happens, and the flag keeps its initial `true`, so the folder is ACCEPTED. The one file that
/// guards the range is the one file that is optional.
pub fn acceptable(version: Option<&ManifestVersion>) -> bool {
    match version {
        None => true,
        Some(version) => validate_version(version) != VersionVerdict::Refused,
    }
}

/// The refusal a folder outside the range produces, message and all.
pub fn version_refusal(version: Option<&ManifestVersion>) -> String {
    let yours = match version {
        Some(version) => version.display(),
        None => "null".to_string(),
    };
    format!(
        "ERROR: Given data source path is too old or too new!  \n\
         \x20      Minimum required version is: {MINIMUM_VERSION_STRING}\n\
         \x20      Maximum allowed version is:  {MAXIMUM_VERSION_STRING}\n\
         \x20      Yours:                       {yours}\n\
         \x20      You must download a compatible version of the data sources from the Broad \
         Institute FTP site: {FTP_PATH}\n\
         \x20      or the Broad Institute Google Bucket: {BUCKET_PATH}\n"
    )
}

// ================================================================================================
// The config file.
// ================================================================================================

/// The `type` key's values, each with the keys it adds to the universal ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceType {
    SimpleXsv,
    LocatableXsv,
    Vcf,
    Gencode,
    Cosmic,
}

impl SourceType {
    /// `getEnum`, which is case-insensitive on the serialised name (`simpleXSV`, `locatableXSV`,
    /// `vcf`, `gencode`, `cosmic`).
    ///
    /// The Java compares with `equalsIgnoreCase`, which also folds the few non-ASCII letters whose
    /// upper case is ASCII; `to_ascii_lowercase` folds only ASCII. A config whose `type` is spelt
    /// with such a letter is accepted there and refused here.
    pub fn parse(text: &str) -> Option<SourceType> {
        let lower = text.to_ascii_lowercase();
        match lower.as_str() {
            "simplexsv" => Some(SourceType::SimpleXsv),
            "locatablexsv" => Some(SourceType::LocatableXsv),
            "vcf" => Some(SourceType::Vcf),
            "gencode" => Some(SourceType::Gencode),
            "cosmic" => Some(SourceType::Cosmic),
            _ => None,
        }
    }

    /// The keys this type adds, in the order they are asserted.
    ///
    /// A VCF and a COSMIC source add none: the source file, the name and the version are all they
    /// need. The lists are each enum constant's `assertConfigFilePropertiesAreValid` in
    /// `FuncotatorArgumentDefinitions.DataSourceType`, in its order.
    pub fn required_keys(self) -> &'static [&'static str] {
        match self {
            SourceType::SimpleXsv => &[
                "xsv_delimiter",
                "xsv_key",
                "xsv_key_column",
                "xsv_permissive_cols",
            ],
            SourceType::LocatableXsv => &[
                "xsv_delimiter",
                "contig_column",
                "start_column",
                "end_column",
            ],
            SourceType::Gencode => &["gencode_fasta_path", "ncbi_build_version"],
            SourceType::Vcf | SourceType::Cosmic => &[],
        }
    }
}

/// The keys EVERY config needs, whatever its type, in the order they are asserted: the six
/// `assertConfigPropertiesContainsKey` calls at the top of `assertConfigFilePropertiesAreValid`.
pub const UNIVERSAL_KEYS: [&str; 6] = [
    "name",
    "version",
    "src_file",
    "origin_location",
    "preprocessing_script",
    "type",
];

/// What `UserException.BadInput` puts in front of its own message.
///
/// It is the EXCEPTION rather than the message that adds it, so the two refusals thrown as plain
/// `UserException` carry no prefix while these do.
pub const BAD_INPUT_PREFIX: &str = "Bad input: ";

/// What reading one config file can go wrong with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConfigError {
    /// A key the config does not carry, named back to the user.
    MissingKey { path: String, key: String },
    /// A `type` value that is not one of the five.
    UnknownType { path: String, value: String },
}

impl ConfigError {
    /// The message as `UserException.BadInput` renders it, prefix included.
    pub fn message(&self) -> String {
        match self {
            ConfigError::MissingKey { path, key } => format!(
                "{BAD_INPUT_PREFIX}Config file for datasource ({path}) does not contain \
                 required key: \"{key}\""
            ),
            ConfigError::UnknownType { path, value } => format!(
                "{BAD_INPUT_PREFIX}ERROR in config file: {path} - Invalid value in \"type\" \
                 field: {value}"
            ),
        }
    }
}

/// The universal keys, then the type, then the keys the type adds.
///
/// The order matters: a config missing both `type` and a locatable XSV's `end_column` is refused
/// for `type`, because the universal keys are checked first.
///
/// Between the universal keys and the type the Java also runs `assertPathFilePropertiesField` on
/// `src_file`, which checks the file exists; that check is not ported. `properties["type"]` below
/// cannot fail: the loop has just checked that the key is present.
pub fn check_config(
    path: &str,
    properties: &BTreeMap<String, String>,
) -> Result<SourceType, ConfigError> {
    for key in UNIVERSAL_KEYS {
        if !properties.contains_key(key) {
            return Err(ConfigError::MissingKey {
                path: path.to_string(),
                key: key.to_string(),
            });
        }
    }
    let value = &properties["type"];
    let Some(source_type) = SourceType::parse(value) else {
        return Err(ConfigError::UnknownType {
            path: path.to_string(),
            value: value.clone(),
        });
    };
    for key in source_type.required_keys() {
        if !properties.contains_key(*key) {
            return Err(ConfigError::MissingKey {
                path: path.to_string(),
                key: (*key).to_string(),
            });
        }
    }
    Ok(source_type)
}

// ================================================================================================
// The folder walk.
// ================================================================================================

/// One config file found under a source's reference-version directory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceConfig {
    /// The path the refusals quote.
    pub path: String,
    /// The reference version whose directory this sits under.
    pub reference: String,
    pub properties: BTreeMap<String, String>,
}

/// One data source folder: a manifest line, if it has one, and the configs beneath it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceFolder {
    pub manifest: Option<ManifestVersion>,
    pub sources: Vec<SourceConfig>,
}

/// What resolving a folder can go wrong with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResolveError {
    /// The version gate, whose message quotes the misread version.
    Version(String),
    /// A config file that could not be read.
    Config(ConfigError),
    /// Two sources under the same reference version with the same `name`.
    DuplicateName { name: String, path: String },
    /// No source had a directory for the requested reference version.
    NoSources { reference: String },
    /// Sources were found, but none of them was a GENCODE one.
    NoGencode,
}

impl ResolveError {
    /// The exception's message. The two `BadInput` refusals carry [`BAD_INPUT_PREFIX`]; the
    /// version refusal and the two sanity checks are plain `UserException`s and do not.
    pub fn message(&self) -> String {
        match self {
            ResolveError::Version(message) => message.clone(),
            ResolveError::Config(error) => error.message(),
            ResolveError::DuplicateName { name, path } => format!(
                "{BAD_INPUT_PREFIX}ERROR: contains more than one dataset of name: {name} - one \
                 is: {path}"
            ),
            ResolveError::NoSources { reference } => {
                format!("ERROR: Could not find any data sources for given reference: {reference}")
            }
            ResolveError::NoGencode => "ERROR: a Gencode datasource is required!".to_string(),
        }
    }
}

/// The sources one folder yields for one reference version.
///
/// The order of the two sanity checks at the end is the order of the messages: an empty result is
/// reported as a missing reference version, and only a NON-EMPTY result without a GENCODE source
/// is reported as a missing GENCODE source.
///
/// `getAndValidateDataSourcesFromPaths`. The Java lists a folder's source directories through
/// `Collectors.toSet()`, so the order in which it meets them is a hash set's; here it is the order
/// of `folder.sources`. That order decides which of two same-named sources is the one the
/// duplicate refusal quotes.
pub fn resolve(folder: &SourceFolder, reference: &str) -> Result<Vec<SourceConfig>, ResolveError> {
    if !acceptable(folder.manifest.as_ref()) {
        return Err(ResolveError::Version(version_refusal(
            folder.manifest.as_ref(),
        )));
    }
    let mut found = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut has_gencode = false;
    for source in &folder.sources {
        // A source with no directory for this reference version is passed over in silence: it is
        // the WHOLE FOLDER coming up empty that is refused, not the individual source.
        if source.reference != reference {
            continue;
        }
        let source_type =
            check_config(&source.path, &source.properties).map_err(ResolveError::Config)?;
        // `check_config` has confirmed `name` is present. `.clone()` copies the string, because the
        // map keeps its own.
        let name = source.properties["name"].clone();
        if names.contains(&name) {
            return Err(ResolveError::DuplicateName {
                name,
                path: source.path.clone(),
            });
        }
        names.push(name);
        has_gencode |= source_type == SourceType::Gencode;
        found.push(source.clone());
    }
    if found.is_empty() {
        return Err(ResolveError::NoSources {
            reference: reference.to_string(),
        });
    }
    if !has_gencode {
        return Err(ResolveError::NoGencode);
    }
    Ok(found)
}

/// The prefix a source puts on its columns: its own `name` key and its own `version` key.
///
/// It is the NAME rather than the directory, so a GENCODE source whose config says `gencode` in
/// lower case produces `gencode_1_genes`, which is not what the renderers look for.
///
/// The two indexings cannot fail on a config that passed [`check_config`], which asserts both
/// keys.
pub fn column_prefix(properties: &BTreeMap<String, String>) -> String {
    format!("{}_{}", properties["name"], properties["version"])
}

// ================================================================================================
// The segments.
// ================================================================================================

/// The three annotation names a call may be written under, tried in this order.
///
/// The FIRST name present decides, even when its value is not a call at all: a file carrying both
/// `CALL` and `Call` is read from `CALL` alone.
///
/// `AnnotatedIntervalToSegmentVariantContextConverter.callAnnotationNames`.
pub const CALL_ANNOTATION_NAMES: [&str; 3] = ["CALL", "Segment_Call", "Call"];

/// `CalledCopyRatioSegment.Call`, by the text a segment file writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Call {
    Deletion,
    Neutral,
    Amplification,
}

impl Call {
    /// `getOutputString()`: the text a segment file writes for the call.
    pub fn output_string(self) -> &'static str {
        match self {
            Call::Deletion => "-",
            Call::Neutral => "0",
            Call::Amplification => "+",
        }
    }

    /// The call whose output string is exactly `text`, or `None`: the `filter` and `findFirst`
    /// inside `retrieveCall`. The comparison is exact, so `+` is an amplification and `+ ` is not.
    pub fn parse(text: &str) -> Option<Call> {
        [Call::Deletion, Call::Neutral, Call::Amplification]
            .into_iter()
            .find(|call| call.output_string() == text)
    }
}

/// The call one segment's annotations carry, or nothing.
///
/// A name that is present with an unrecognised value yields nothing AND STOPS THE SEARCH, so a
/// segment whose `CALL` column reads `NA` is not rescued by a later `Segment_Call` column.
///
/// ```java
/// if (segment.hasAnnotation(callAnnotationName)) {
/// ```
///
/// returns from inside the loop whatever the lookup found, `null` included.
pub fn call_of(annotations: &BTreeMap<String, String>) -> Option<Call> {
    for name in CALL_ANNOTATION_NAMES {
        if let Some(value) = annotations.get(name) {
            return Call::parse(value);
        }
    }
    None
}

/// `<COPY_NEUTRAL>`, and the two the converter shares with structural variants.
pub const COPY_NEUTRAL_ALLELE: &str = "<COPY_NEUTRAL>";
/// `Allele.UNSPECIFIED_ALTERNATE_ALLELE_STRING`, which a segment with no call gets.
pub const UNSPECIFIED_ALLELE: &str = "<*>";

/// The alternate allele a call becomes.
///
/// An amplification becomes `<INS>` rather than `<DUP>`, which is the converter's own choice and
/// not the segment file's.
///
/// ```java
/// case AMPLIFICATION:
///     return Allele.create(SimpleSVType.createBracketedSymbAlleleString("INS"), false);
/// ```
pub fn allele_of(call: Option<Call>) -> &'static str {
    match call {
        None => UNSPECIFIED_ALLELE,
        Some(Call::Deletion) => "<DEL>",
        Some(Call::Amplification) => "<INS>",
        Some(Call::Neutral) => COPY_NEUTRAL_ALLELE,
    }
}

/// `DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT`, which a segment must EXCEED rather than reach.
pub const MIN_BASES_FOR_VALID_SEGMENT: i64 = 150;

/// The alternate alleles a segment may carry: `isSegmentVariantContext`'s `ACCEPTABLE_ALT_ALLELES`.
///
/// The Java builds that list from every `SimpleSVType.SupportedType`, which is
/// `INV, DEL, INS, DUP, DUP_INV`, plus copy-neutral and unspecified. `<DUP_INV>` is missing here,
/// so a variant context whose only alternate is `<DUP_INV>` is a segment to the reference and not
/// to this port. Only membership is asked, so the order does not matter.
pub const ACCEPTABLE_ALTERNATES: [&str; 6] = [
    "<INS>",
    "<DEL>",
    "<DUP>",
    "<INV>",
    COPY_NEUTRAL_ALLELE,
    UNSPECIFIED_ALLELE,
];

/// Whether a variant context could be a copy-number segment.
///
/// Both halves must hold: one acceptable alternate, and a size strictly above the minimum. A
/// segment of exactly the minimum length is not one. The size is `VariantContextUtils.getSize`,
/// `end - start + 1`.
pub fn is_segment(alternates: &[String], start: i64, end: i64, minimum: i64) -> bool {
    let acceptable = alternates
        .iter()
        .any(|allele| ACCEPTABLE_ALTERNATES.contains(&allele.as_str()));
    acceptable && (end - start + 1) > minimum
}

// ================================================================================================
// The gene list beside the output.
// ================================================================================================

/// One row of `<output>.gene_list.txt`: a gene, and the exon within it or nothing.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneRow {
    pub gene: String,
    pub exon: String,
}

/// What one segment contributes to the gene list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SegmentGenes {
    /// Every gene the segment covers, from the `genes` field.
    pub genes: Vec<String>,
    /// The gene and exon the segment STARTS in, if it starts inside one.
    pub start: Option<(String, String)>,
    /// The gene and exon the segment ENDS in, if it ends inside one.
    pub end: Option<(String, String)>,
}

/// The gene list rows, sorted by gene and then by exon.
///
/// A segment covering no gene contributes NOTHING, which is why the file is shorter than the
/// output beside it. Here a gene the segment covers gets a row with an EMPTY exon, and the same
/// gene gets a further row per exon the segment starts or ends in, so one gene can appear twice.
///
/// **This is not the reference's rule.** `GeneListOutputRenderer` keys a `TreeMap` on the start
/// and end `(gene, exon)` pairs, and then adds the covered genes with `""` only when they are
/// neither of those two genes:
///
/// ```java
/// genes.stream().filter(g -> !g.equals(startGene) && !g.equals(endGene))
/// ```
///
/// The golden's `(ALPHA, "")` row is the segment's END pair (its end gene is `ALPHA` and its end
/// exon empty), not a covered-gene row, so the golden cannot tell the two rules apart. A segment
/// covering `A` and `B` that starts in an exon of `A` and ends in an exon of `B` gets two rows
/// there and four here. The map's comparator is also `SimpleNaturalComparator` on both halves, so
/// `G2` sorts before `G10` there; `sort` here compares the strings character by character, which
/// puts `G10` first. The renderer the tool runs, `GeneList` in `funcotator_engine`, follows the
/// reference on both counts; this function is what the suite measures against the golden.
pub fn gene_rows(segments: &[SegmentGenes]) -> Vec<GeneRow> {
    let mut rows: Vec<GeneRow> = Vec::new();
    for segment in segments {
        for gene in &segment.genes {
            push_once(
                &mut rows,
                GeneRow {
                    gene: gene.clone(),
                    exon: String::new(),
                },
            );
        }
        // `.flatten()` over two `Option`s skips the `None`s: a segment that starts or ends outside
        // every gene adds no row for that end.
        for (gene, exon) in [&segment.start, &segment.end].into_iter().flatten() {
            push_once(
                &mut rows,
                GeneRow {
                    gene: gene.clone(),
                    exon: exon.clone(),
                },
            );
        }
    }
    rows.sort();
    rows
}

/// Add `row` unless it is already there: the set half of the reference's `TreeMap` keys.
fn push_once(rows: &mut Vec<GeneRow>, row: GeneRow) {
    if !rows.contains(&row) {
        rows.push(row);
    }
}

/// The text a SEG column with no value carries, which is not the empty string.
///
/// The columns the SEGMENT FILE did not carry are `__UNKNOWN__`; the columns a SOURCE did not fill
/// are empty. The two are different absences and the output keeps them apart.
///
/// `FuncotatorConstants.UNKNOWN_VALUE_STRING`.
pub const UNKNOWN: &str = "__UNKNOWN__";
