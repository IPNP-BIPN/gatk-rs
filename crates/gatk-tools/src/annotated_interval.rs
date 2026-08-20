//! The annotated-interval collection, ported from
//! `org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval` and the
//! `annotated_region_default.config` beside it (GATK 4.6.2.0).
//!
//! A tab-separated file of regions with any number of annotation columns, which four copy-number
//! tools read and write.
//!
//! # The locatable columns are found by name from a fixed list
//!
//! The config names twelve spellings of the contig column, sixteen of the start and sixteen of the
//! end, and `Position` (with `position`, `pos` and `POS`) is in BOTH the start list and the end
//! list. A file whose only coordinate column is `Position` therefore parses as one-base regions,
//! its start and its end being the same column.
//!
//! # The output is always this format, whatever the input was
//!
//! A SAM header, then three `@CO` lines naming the locatable columns the writer renamed, then a
//! column line of `CONTIG`, `START`, `END` and the annotations SORTED ALPHABETICALLY. An input's
//! own comments become `@CO` lines ahead of those three; an input's `@SQ` lines are kept.
//!
//! # A file of no rows throws
//!
//! The collection reads its annotation names off the first record, so a file with a column line and
//! nothing under it dies with an `IndexOutOfBoundsException` rather than writing an empty result.
//! That is what the golden's `no-rows` row is, and this port returns it as a refusal.

use std::collections::BTreeMap;

/// `annotated_region_default.config`'s `contig_column`.
pub const CONTIG_COLUMNS: [&str; 12] = [
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

/// `start_column`. `segment_start` appears twice in the config and is listed once here.
pub const START_COLUMNS: [&str; 15] = [
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
    "segment_start",
];

/// `end_column`, which shares `Position`, `position`, `pos` and `POS` with the start list.
pub const END_COLUMNS: [&str; 16] = [
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

/// `MergeAnnotatedRegions.DEFAULT_SEPARATOR`.
pub const DEFAULT_SEPARATOR: &str = "__";

/// One region and its annotations, which are held sorted by name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnnotatedInterval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    /// A `SortedMap`, which is what makes the written columns alphabetical.
    pub annotations: BTreeMap<String, String>,
}

/// A whole file: the SAM header lines it came with, the annotation names, and the records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnnotatedIntervalCollection {
    /// The `@` lines of the input, `@HD` included, in the order they were read.
    pub header_lines: Vec<String>,
    /// The `#` comments of the input, without their marker.
    pub comments: Vec<String>,
    /// Sorted alphabetically, as the collection's constructor sorts them.
    pub annotations: Vec<String>,
    pub records: Vec<AnnotatedInterval>,
    /// The names the locatable columns had in the input, which the output records as comments.
    pub contig_column: String,
    pub start_column: String,
    pub end_column: String,
}

/// What reading a collection refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CollectionError {
    /// No column of the input names a contig, a start or an end.
    NoLocatableColumns,
    /// A column line with no rows under it: the annotations are read off record zero.
    NoRecords,
}

impl CollectionError {
    pub fn java_class(&self) -> &'static str {
        match self {
            CollectionError::NoLocatableColumns => {
                "htsjdk.tribble.TribbleException$MalformedFeatureFile"
            }
            CollectionError::NoRecords => "java.lang.IndexOutOfBoundsException",
        }
    }

    pub fn message(&self) -> String {
        self.message_with_source("")
    }

    /// The message, with the input's URI as tribble appends it.
    ///
    /// The reader is given text rather than a path here, so the source is the caller's to supply;
    /// an empty one leaves the message without tribble's trailing clause.
    pub fn message_with_source(&self, source: &str) -> String {
        match self {
            CollectionError::NoLocatableColumns => {
                let mut text = format!(
                    "Unable to parse header with error: Bad input: Input did not contain any \
                     headers from the list: {}",
                    CONTIG_COLUMNS.join(",")
                );
                if !source.is_empty() {
                    text.push_str(&format!(", for input source: {source}"));
                }
                text
            }
            CollectionError::NoRecords => "Index 0 out of bounds for length 0".to_string(),
        }
    }
}

/// The first column of `columns` that appears in `names`, which is how the codec picks each one.
fn pick<'a>(columns: &[&str], names: &'a [String]) -> Option<&'a String> {
    names.iter().find(|name| columns.contains(&name.as_str()))
}

/// Read a file. `text` is the whole file, lines and all.
pub fn read(text: &str) -> Result<AnnotatedIntervalCollection, CollectionError> {
    let mut header_lines = Vec::new();
    let mut comments = Vec::new();
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with('@') {
            header_lines.push(line.to_string());
        } else if let Some(comment) = line.strip_prefix('#') {
            comments.push(comment.to_string());
        } else {
            rows.push(line);
        }
    }
    let Some(column_line) = rows.first() else {
        return Err(CollectionError::NoLocatableColumns);
    };
    let names: Vec<String> = column_line.split('\t').map(str::to_string).collect();
    let (Some(contig_column), Some(start_column), Some(end_column)) = (
        pick(&CONTIG_COLUMNS, &names),
        pick(&START_COLUMNS, &names),
        pick(&END_COLUMNS, &names),
    ) else {
        return Err(CollectionError::NoLocatableColumns);
    };
    let (contig_column, start_column, end_column) = (
        contig_column.clone(),
        start_column.clone(),
        end_column.clone(),
    );
    let mut records = Vec::new();
    for row in &rows[1..] {
        let fields: Vec<&str> = row.split('\t').collect();
        let field = |wanted: &str| {
            names
                .iter()
                .position(|name| name == wanted)
                .and_then(|index| fields.get(index))
                .copied()
                .unwrap_or("")
        };
        let mut annotations = BTreeMap::new();
        for (index, name) in names.iter().enumerate() {
            if *name == contig_column || *name == start_column || *name == end_column {
                continue;
            }
            annotations.insert(
                name.clone(),
                fields.get(index).copied().unwrap_or("").to_string(),
            );
        }
        records.push(AnnotatedInterval {
            contig: field(&contig_column).to_string(),
            start: field(&start_column).parse().unwrap_or(0),
            end: field(&end_column).parse().unwrap_or(0),
            annotations,
        });
    }
    // The collection takes its annotation names off the first record.
    let Some(first) = records.first() else {
        return Err(CollectionError::NoRecords);
    };
    let annotations: Vec<String> = first.annotations.keys().cloned().collect();
    Ok(AnnotatedIntervalCollection {
        header_lines,
        comments,
        annotations,
        records,
        contig_column,
        start_column,
        end_column,
    })
}

impl AnnotatedIntervalCollection {
    /// The file this collection writes: the SAM header, the three column comments, the column
    /// line, then the records.
    ///
    /// The `@HD` line is written whether or not the input had one, and the input's own comments
    /// come before the three the writer adds.
    pub fn write(&self) -> String {
        let mut text = String::new();
        let mut wrote_hd = false;
        for line in &self.header_lines {
            if line.starts_with("@HD") {
                wrote_hd = true;
            }
            text.push_str(line);
            text.push('\n');
        }
        if !wrote_hd {
            text.insert_str(0, "@HD\tVN:1.6\n");
        }
        for comment in &self.comments {
            text.push_str(&format!("@CO\t{comment}\n"));
        }
        text.push_str(&format!("@CO\t_ContigHeader={}\n", "CONTIG"));
        text.push_str(&format!("@CO\t_StartHeader={}\n", "START"));
        text.push_str(&format!("@CO\t_EndHeader={}\n", "END"));
        let mut columns = vec!["CONTIG".to_string(), "START".to_string(), "END".to_string()];
        columns.extend(self.annotations.iter().cloned());
        text.push_str(&columns.join("\t"));
        text.push('\n');
        for record in &self.records {
            let mut fields = vec![
                record.contig.clone(),
                record.start.to_string(),
                record.end.to_string(),
            ];
            for annotation in &self.annotations {
                fields.push(
                    record
                        .annotations
                        .get(annotation)
                        .cloned()
                        .unwrap_or_default(),
                );
            }
            text.push_str(&fields.join("\t"));
            text.push('\n');
        }
        text
    }
}

/// `IntervalUtils.sortLocatablesBySequenceDictionary`, which puts an unknown contig last rather
/// than refusing it.
pub fn sort_by_dictionary(records: &mut [AnnotatedInterval], dictionary: &[String]) {
    let index_of = |contig: &str| {
        dictionary
            .iter()
            .position(|name| name == contig)
            .map(|index| index as i64)
            .unwrap_or(i64::MAX)
    };
    records.sort_by(|left, right| {
        (index_of(&left.contig), left.start, left.end).cmp(&(
            index_of(&right.contig),
            right.start,
            right.end,
        ))
    });
}

/// `IntervalUtils.overlaps`, which is a real overlap: abutting regions do not.
fn overlaps(left: &AnnotatedInterval, right: &AnnotatedInterval) -> bool {
    left.contig == right.contig && left.start <= right.end && right.start <= left.end
}

/// `renderConflict`: split both values on the separator, deduplicate, sort, rejoin.
///
/// The split is `StringUtils.splitByWholeSeparator`, which answers an EMPTY ARRAY for an empty
/// string rather than an array holding one empty string. So an empty value contributes nothing:
/// merging `` with `7` gives `7`, not `__7`. Rust's own `split` would give the second.
pub fn render_conflict(first: &str, second: &str, separator: &str) -> String {
    let split = |text: &str| -> Vec<String> {
        if text.is_empty() {
            Vec::new()
        } else {
            text.split(separator).map(str::to_string).collect()
        }
    };
    let mut values: Vec<String> = Vec::new();
    for value in split(first).into_iter().chain(split(second)) {
        if !values.contains(&value) {
            values.push(value);
        }
    }
    values.sort();
    values.join(separator)
}

/// `merge`, on two regions already known to overlap.
fn merge(
    first: &AnnotatedInterval,
    second: &AnnotatedInterval,
    separator: &str,
) -> AnnotatedInterval {
    let mut annotations: BTreeMap<String, String> = BTreeMap::new();
    for key in first.annotations.keys().chain(second.annotations.keys()) {
        if annotations.contains_key(key) {
            continue;
        }
        let value = match (first.annotations.get(key), second.annotations.get(key)) {
            (Some(left), Some(right)) => render_conflict(left, right, separator),
            (Some(left), None) => left.clone(),
            (None, Some(right)) => right.clone(),
            (None, None) => continue,
        };
        annotations.insert(key.clone(), value);
    }
    AnnotatedInterval {
        contig: first.contig.clone(),
        start: first.start.min(second.start),
        end: first.end.max(second.end),
        annotations,
    }
}

/// `AnnotatedIntervalUtils.mergeRegions`: one pass with a peek, over the sorted list.
///
/// The merged region is what the next comparison uses, so a chain of overlaps collapses into one
/// region however long it is.
pub fn merge_regions(
    records: &[AnnotatedInterval],
    dictionary: &[String],
    separator: &str,
) -> Vec<AnnotatedInterval> {
    let mut sorted = records.to_vec();
    sort_by_dictionary(&mut sorted, dictionary);
    let mut merged = Vec::new();
    let mut index = 0;
    while index < sorted.len() {
        let mut current = sorted[index].clone();
        index += 1;
        while index < sorted.len() && overlaps(&current, &sorted[index]) {
            current = merge(&current, &sorted[index], separator);
            index += 1;
        }
        merged.push(current);
    }
    merged
}

/// `getDistance`: zero for OVERLAPPING regions, the gap between the nearer endpoints otherwise.
///
/// The method's own comment says overlapping *or abutting* regions answer zero. They do not:
/// `IntervalUtils.overlaps` does not count an abuttal, so `1-100` and `101-200` are one apart and a
/// maximum merge distance of zero merges nothing at all.
///
/// Regions on different contigs are `Long.MAX_VALUE` apart.
pub fn distance(left: &AnnotatedInterval, right: &AnnotatedInterval) -> i64 {
    if left.contig != right.contig {
        return i64::MAX;
    }
    if left.start <= right.end && right.start <= left.end {
        return 0;
    }
    if left.end < right.start {
        i64::from(right.start - left.end)
    } else {
        i64::from(left.start - right.end)
    }
}

/// `mergeRegionsByAnnotation`: neighbours merge when the named annotations agree and the distance
/// is within the maximum.
///
/// The comparison is against the region built so far, so a run of short hops travels however far
/// it goes in total. The annotations NOT named are reconciled the same way an overlap merge
/// reconciles them.
pub fn merge_regions_by_annotation(
    records: &[AnnotatedInterval],
    dictionary: &[String],
    annotation_names: &[String],
    separator: &str,
    max_distance: i64,
) -> Vec<AnnotatedInterval> {
    let mut sorted = records.to_vec();
    sort_by_dictionary(&mut sorted, dictionary);
    let mut merged = Vec::new();
    let mut index = 0;
    while index < sorted.len() {
        let mut current = sorted[index].clone();
        index += 1;
        while index < sorted.len() {
            let next = &sorted[index];
            let agrees = annotation_names
                .iter()
                .all(|name| current.annotations.get(name) == next.annotations.get(name));
            if current.contig != next.contig || distance(&current, next) > max_distance || !agrees {
                break;
            }
            current = merge(&current, next, separator);
            index += 1;
        }
        merged.push(current);
    }
    merged
}

/// The file `MergeAnnotatedRegionsByAnnotation` writes, which goes through the writer rather than
/// the collection: three column names of the caller's choosing, the annotations of the FIRST
/// region, and no SAM header at all.
pub fn write_without_header(
    records: &[AnnotatedInterval],
    contig_column: &str,
    start_column: &str,
    end_column: &str,
) -> Result<String, CollectionError> {
    let first = records.first().ok_or(CollectionError::NoRecords)?;
    let annotations: Vec<String> = first.annotations.keys().cloned().collect();
    let mut columns = vec![
        contig_column.to_string(),
        start_column.to_string(),
        end_column.to_string(),
    ];
    columns.extend(annotations.iter().cloned());
    let mut text = columns.join("\t");
    text.push('\n');
    for record in records {
        let mut fields = vec![
            record.contig.clone(),
            record.start.to_string(),
            record.end.to_string(),
        ];
        for annotation in &annotations {
            fields.push(
                record
                    .annotations
                    .get(annotation)
                    .cloned()
                    .unwrap_or_default(),
            );
        }
        text.push_str(&fields.join("\t"));
        text.push('\n');
    }
    Ok(text)
}
