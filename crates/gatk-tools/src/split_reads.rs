//! Ported from `org.broadinstitute.hellbender.tools.SplitReads` (GATK 4.6.2.0).
//!
//! The seventh whole tool of the record-transform archetype, and the first whose `-O` is a
//! **directory** and whose run opens more than one writer. Its `apply` is one line; everything
//! worth porting is in which files exist and what header each one got.
//!
//! ```java
//! public void apply(GATKRead read, ReferenceContext r, FeatureContext f) {
//!     outs.computeIfAbsent(getKey(splitters, read), this::createUnknownOutOnDemand).addRead(read);
//! }
//! ```
//!
//! # Every output file of a run has a different header
//!
//! `createSAMWriter` calls `getHeaderForSAMWriter`, which **adds** a `@PG` record to the reads
//! header in place and hands back the same object. A writer serialises the header when it is
//! created, so the nth writer's file carries n records for this one tool, `GATK SplitReads`, then
//! `GATK SplitReads.1`, `.2`. Measured: splitting by library writes four files whose `@PG` lists
//! are one, two, three and four entries long, and the fourth is the on-demand `unknown` writer,
//! created during the traversal rather than before it.
//!
//! Nothing in this archetype's first six tools could show that, because they open one writer.
//! [`crate::sam_output::header_for_sam_writer`] already reproduced the in-place mutation; here it
//! is observable.
//!
//! # The files come from the header, not from the reads
//!
//! `createWriters` takes the cross product of each splitter's `getSplitsBy(header)`, and
//! `ReadGroupSplitter.getSplitsBy` maps `header.getReadGroups()` through its selector. A read group
//! that no read belongs to still gets a file, and that file is a valid empty BAM with an index.
//!
//! The cross product **does not deduplicate**: two read groups with the same sample give the same
//! key twice, so `prepareSAMFileWriter` runs twice on one path and only the second writer is ever
//! closed. The file that survives is the second writer's, and it carries the first writer's `@PG`
//! record too, because that record was added to the shared header before the second writer read
//! it.
//!
//! # A null value is spelled two ways, and on three splitters that aborts the run
//!
//! `addKey` concatenates the value as an `Object`, so a read group with no `LB` contributes the
//! four characters `null` to the key built from the header. `getKey` substitutes
//! `UNKNOWN_OUT_PREFIX` for the same null, so the key built from a read says `unknown`. The read's
//! key is therefore not in the map, and `createUnknownOutOnDemand` accepts exactly one string,
//! `.unknown`:
//!
//! * with one splitter the read's key **is** `.unknown`, so a writer is made on demand and the run
//!   finishes with one empty `.null` file and one `.unknown` file holding the reads;
//! * with three splitters it is `.s2.rg3.unknown`, which is refused:
//!   `ShouldNeverReachHereException`, and the whole run dies.
//!
//! # A read with no read group at all is a null pointer
//!
//! `ReadGroupSplitter.getSplitBy` runs its selector on `getSAMReadGroupRecord(record, header)`,
//! which is null when the read carries no `RG`. What saves the tool is a filter rather than the
//! tool: `WellformedReadFilter`'s `HAS_READ_GROUP` drops the read first. Measured with that filter
//! disabled, the run ends in `java.lang.NullPointerException`.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use gatk_engine::reads::{ReadsDataSource, ReadsError};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK SplitReads";

/// `UNKNOWN_OUT_PREFIX`: what a **read** contributes when its splitter has no value for it.
pub const UNKNOWN_OUT_PREFIX: &str = "unknown";

/// What a **header** contributes for the same missing value, because `addKey` concatenates the
/// value as an `Object` and Java prints a null reference as these four characters.
pub const NULL_FROM_THE_HEADER: &str = "null";

/// The three `ReaderSplitter`s this tool offers, in the order `onTraversalStart` adds them: the
/// order of the key's components is the order of the arguments in the source, not the order they
/// were given on the command line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Splitter {
    /// `--split-sample`, `SampleNameSplitter`.
    Sample,
    /// `--split-read-group`, `ReadGroupIdSplitter`.
    ReadGroupId,
    /// `--split-library-name`, `LibraryNameSplitter`.
    LibraryName,
}

impl Splitter {
    /// The selector a `ReadGroupSplitter` applies to a read group.
    fn value(self, group: &htsjdk_bam::header::ReadGroup) -> Option<String> {
        match self {
            Splitter::Sample => group.attributes.get("SM").map(str::to_string),
            Splitter::ReadGroupId => Some(group.id.clone()),
            Splitter::LibraryName => group.attributes.get("LB").map(str::to_string),
        }
    }

    /// `getSplitsBy(header)`: the selector over **every** read group, duplicates and nulls kept.
    pub fn splits_by(self, header: &SamHeader) -> Vec<Option<String>> {
        header
            .read_groups
            .iter()
            .map(|group| self.value(group))
            .collect()
    }

    /// `getSplitBy(record, header)`: the selector over the read's own read group.
    ///
    /// [`NoReadGroup`] is the null pointer the reference dereferences: the selector is applied to
    /// a null read group rather than guarded against one.
    pub fn split_by(
        self,
        read: &BamRecord,
        header: &SamHeader,
    ) -> Result<Option<String>, NoReadGroup> {
        let id = match read.tags.get(Tag::new(b"RG")) {
            Some(TagValue::Str(id)) => id.clone(),
            _ => return Err(NoReadGroup),
        };
        let group = header
            .read_groups
            .iter()
            .find(|group| group.id == id)
            .ok_or(NoReadGroup)?;
        Ok(self.value(group))
    }
}

/// The read had no read group to ask, and the reference asked it anyway.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoReadGroup;

impl NoReadGroup {
    /// The reference throws out of a method reference applied to null, and the message is empty.
    pub fn message(&self) -> String {
        "null".to_string()
    }

    pub fn class(&self) -> &'static str {
        "java.lang.NullPointerException"
    }
}

/// Why a run stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitError {
    /// A read that cannot answer any splitter.
    NoReadGroup(NoReadGroup),
    /// `createUnknownOutOnDemand` was handed a key it does not recognise, which is any key with a
    /// missing value in it other than the whole key being `.unknown`.
    UnrecognizedKey(String),
}

impl SplitError {
    pub fn message(&self) -> String {
        match self {
            SplitError::NoReadGroup(error) => error.message(),
            SplitError::UnrecognizedKey(key) => {
                format!("Unrecognized attribute value found: {key}")
            }
        }
    }

    pub fn class(&self) -> &'static str {
        match self {
            SplitError::NoReadGroup(error) => error.class(),
            SplitError::UnrecognizedKey(_) => {
                "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException"
            }
        }
    }
}

/// `addKey`: the cross product of the splitters' header values, **in order and with duplicates**.
///
/// One `.value` per splitter, and a missing value spelled [`NULL_FROM_THE_HEADER`]. With no
/// splitters at all this is a single empty key, which is what makes a run with no `--split-*`
/// write one file named after the input.
pub fn keys_from_header(header: &SamHeader, splitters: &[Splitter]) -> Vec<String> {
    let mut keys = vec![String::new()];
    for splitter in splitters {
        let values = splitter.splits_by(header);
        let mut next = Vec::with_capacity(keys.len() * values.len());
        for key in &keys {
            for value in &values {
                next.push(format!(
                    "{key}.{}",
                    value.as_deref().unwrap_or(NULL_FROM_THE_HEADER)
                ));
            }
        }
        keys = next;
    }
    keys
}

/// `getKey`: the same shape for one read, with a missing value spelled [`UNKNOWN_OUT_PREFIX`].
pub fn key_for_read(
    read: &BamRecord,
    header: &SamHeader,
    splitters: &[Splitter],
) -> Result<String, NoReadGroup> {
    let mut key = String::new();
    for splitter in splitters {
        let value = splitter.split_by(read, header)?;
        key.push('.');
        key.push_str(value.as_deref().unwrap_or(UNKNOWN_OUT_PREFIX));
    }
    Ok(key)
}

/// One file the run left behind.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputFile {
    /// `<input base name><key><input extension>`.
    pub name: String,
    pub bam: Vec<u8>,
    pub index: Option<Vec<u8>>,
}

/// `prepareSAMFileWriter`: the file a key is written to.
pub fn file_name(base_name: &str, key: &str, extension: &str) -> String {
    format!("{base_name}{key}{extension}")
}

/// What a run produces: every file in its output directory, or the refusal.
pub type RunResult = Result<Result<Vec<OutputFile>, SplitError>, ReadsError>;

/// `SplitReads`: the reads the traversal reaches, written into one file per key.
///
/// The order of the returned files is the order their writers were created, which is the order
/// their headers grew in: the nth carries n `@PG` records for this tool. A caller that wants them
/// by name sorts them itself, as the reference's own directory listing does.
pub fn split_reads(
    source: &ReadsDataSource,
    options: &Options,
    splitters: &[Splitter],
    base_name: &str,
    extension: &str,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> RunResult {
    let input_header = source.header().clone();
    let records = crate::read_walker::traverse(source, &options.intervals, filter)?;

    // `outs`, a LinkedHashMap: a key put twice keeps its first position and takes the second
    // value, and the writer the first put made is never closed. What survives on disk is the
    // second writer's file, and it carries the @PG record the first writer's creation added.
    let mut keys: Vec<String> = Vec::new();
    let mut headers: Vec<SamHeader> = Vec::new();
    // The reads header the reference mutates in place, once per writer created.
    let mut growing = input_header.clone();
    let mut open_writer = |key: &str, keys: &mut Vec<String>, headers: &mut Vec<SamHeader>| {
        growing = header_for_sam_writer(&growing, TOOL_NAME, options);
        match keys.iter().position(|existing| existing == key) {
            Some(at) => headers[at] = growing.clone(),
            None => {
                keys.push(key.to_string());
                headers.push(growing.clone());
            }
        }
    };

    for key in keys_from_header(&input_header, splitters) {
        open_writer(&key, &mut keys, &mut headers);
    }

    let mut written: Vec<Vec<BamRecord>> = vec![Vec::new(); keys.len()];
    for record in &records {
        let key = match key_for_read(record, &input_header, splitters) {
            Ok(key) => key,
            Err(error) => return Ok(Err(SplitError::NoReadGroup(error))),
        };
        let at = match keys.iter().position(|existing| *existing == key) {
            Some(at) => at,
            None => {
                // `createUnknownOutOnDemand` accepts exactly one string and refuses every other.
                if key != format!(".{UNKNOWN_OUT_PREFIX}") {
                    return Ok(Err(SplitError::UnrecognizedKey(key)));
                }
                open_writer(&key, &mut keys, &mut headers);
                written.push(Vec::new());
                keys.len() - 1
            }
        };
        written[at].push(record.clone());
    }

    let mut files = Vec::with_capacity(keys.len());
    for ((key, header), records) in keys.iter().zip(&headers).zip(&written) {
        let (bam, index) = write_records(header, records, options.create_output_bam_index)?;
        files.push(OutputFile {
            name: file_name(base_name, key, extension),
            bam,
            index,
        });
    }
    Ok(Ok(files))
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::ReadGroup;

    fn header() -> SamHeader {
        let mut header = SamHeader::default();
        for (id, sample, library) in [
            ("rg1", "s1", Some("lib1")),
            ("rg2", "s1", Some("lib2")),
            ("rg3", "s2", None),
        ] {
            let mut group = ReadGroup::new(id);
            group.attributes.set("SM", sample);
            if let Some(library) = library {
                group.attributes.set("LB", library);
            }
            header.read_groups.push(group);
        }
        header
    }

    fn read(group: Option<&str>) -> BamRecord {
        let mut read = BamRecord::default();
        if let Some(group) = group {
            read.tags
                .insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
        }
        read
    }

    #[test]
    fn the_header_decides_the_keys_and_keeps_the_duplicates() {
        let keys = keys_from_header(&header(), &[Splitter::Sample]);
        assert_eq!(keys, [".s1", ".s1", ".s2"], "two read groups, one sample");
    }

    #[test]
    fn a_missing_library_is_null_from_the_header_and_unknown_from_a_read() {
        let header = header();
        assert_eq!(
            keys_from_header(&header, &[Splitter::LibraryName]),
            [".lib1", ".lib2", ".null"]
        );
        assert_eq!(
            key_for_read(&read(Some("rg3")), &header, &[Splitter::LibraryName]).unwrap(),
            ".unknown"
        );
    }

    #[test]
    fn no_splitters_at_all_is_one_empty_key() {
        assert_eq!(keys_from_header(&header(), &[]), [""]);
        assert_eq!(
            key_for_read(&read(Some("rg1")), &header(), &[]).unwrap(),
            ""
        );
        assert_eq!(file_name("plain", "", ".bam"), "plain.bam");
    }

    #[test]
    fn a_read_with_no_read_group_is_the_reference_s_null_pointer() {
        let error = key_for_read(&read(None), &header(), &[Splitter::Sample]).unwrap_err();
        assert_eq!(error.class(), "java.lang.NullPointerException");
        assert_eq!(error.message(), "null");
    }

    /// The key that kills the run: three splitters, and the third has no value.
    #[test]
    fn three_splitters_build_a_key_the_on_demand_writer_refuses() {
        let header = header();
        let splitters = [
            Splitter::Sample,
            Splitter::ReadGroupId,
            Splitter::LibraryName,
        ];
        let key = key_for_read(&read(Some("rg3")), &header, &splitters).unwrap();
        assert_eq!(key, ".s2.rg3.unknown");
        assert!(
            !keys_from_header(&header, &splitters).contains(&key),
            "the header spells the same read group .s2.rg3.null"
        );
        let error = SplitError::UnrecognizedKey(key);
        assert_eq!(
            error.message(),
            "Unrecognized attribute value found: .s2.rg3.unknown"
        );
    }

    #[test]
    fn the_cross_product_is_over_read_groups_not_over_tuples() {
        // Three read groups, three splitters, 27 combinations and 18 distinct keys: the product is
        // taken per splitter, so combinations no read group actually has are keys too.
        let keys = keys_from_header(
            &header(),
            &[
                Splitter::Sample,
                Splitter::ReadGroupId,
                Splitter::LibraryName,
            ],
        );
        assert_eq!(keys.len(), 27);
        let mut distinct = keys.clone();
        distinct.sort();
        distinct.dedup();
        assert_eq!(distinct.len(), 18);
        assert!(distinct.contains(&".s1.rg3.lib1".to_string()));
    }
}
