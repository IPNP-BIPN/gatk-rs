//! `SVStratify` and `SVStatificationEngine`, ported from GATK 4.6.2.0.
//!
//! Structural variants sorted into groups by type, size and reference-track overlap. Reading the
//! VCF is not ported; deciding which strata a record belongs to, and what each refusal says, is.
//!
//! # The strata live in a HashMap, so the configuration's order is not the output's
//!
//! ```java
//! strata = new HashMap<>();
//! ```
//!
//! `getMatches` walks `strata.values()`, so the order a record's matches come back in, the order
//! the refusal message lists them, and the order the split-output files are created are all
//! `java.util.HashMap` iteration order over the stratum NAMES. A configuration declaring
//! `DEL_small_RM` before `DEL_small_both` is reported the other way round.
//!
//! # Without split output every stratum writes to the same file
//!
//! ```java
//! final VariantContextWriter writer = splitOutput ? writers.get(stratum.getName()) : writers.get(DEFAULT_STRATUM);
//! ```
//!
//! So a record matching two strata under `--allow-multiple-matches` is written TWICE into one
//! file, same position and same ID, with two different `STRAT` values.
//!
//! # An insertion ignores both thresholds
//!
//! ```java
//! if (record.getType() == INS) { return matchesTrackBreakpointOverlap(record, 1); }
//! ```
//!
//! Neither `--stratify-overlap-fraction` nor `--stratify-num-breakpoint-overlaps` reaches an
//! insertion: it always asks for exactly one end in a track. BND and CTX take the interchromosomal
//! threshold instead, and only the remaining types run the fraction.
//!
//! # A breakpoint counts once per end, not once per track
//!
//! `countAnyTrackOverlap` returns 1 as soon as ANY named track overlaps, so a stratum naming two
//! tracks cannot reach two from a single end.
//!
//! # Only `-1` is a null bound
//!
//! The null set is `{"-1", "", "NULL", "NA"}`, so `-2` parses as a number and is then refused by
//! the constructor as a negative bound. A null minimum becomes `Integer.MIN_VALUE` and a null
//! maximum `Integer.MAX_VALUE`, which is why a record with no length matches only a stratum with
//! neither bound.
//!
//! # And the column-count message prints the same number twice
//!
//! ```java
//! throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
//! ```
//!
//! Both halves read the same value, so an extra column is reported as `Expected 6 columns but
//! found 6`.

use gatk_engine::java_hash::{hash_map_order, string_hash_code};

/// `SVStratify.DEFAULT_STRATUM`.
pub const DEFAULT_STRATUM: &str = "default";
/// `SVStatificationEngine.NULL_TABLE_VALUES`.
pub const NULL_TABLE_VALUES: &[&str] = &["-1", "", "NULL", "NA"];
/// The five columns the configuration table must have, and only those.
pub const COLUMN_NAMES: &[&str] = &["NAME", "SVTYPE", "MIN_SIZE", "MAX_SIZE", "TRACKS"];

/// `GATKSVVCFConstants.StructuralVariantAnnotationType`, over the types this decides on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvType {
    Del,
    Dup,
    Ins,
    Inv,
    /// A copy number variant, which stands in for both a deletion and a duplication when two
    /// simple CNVs are compared. SVStratify never names it in a configuration; the clustering
    /// linkage does.
    Cnv,
    Cpx,
    Bnd,
    Ctx,
}

impl SvType {
    pub fn parse(text: &str) -> Option<SvType> {
        Some(match text {
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            "INS" => SvType::Ins,
            "INV" => SvType::Inv,
            "CNV" => SvType::Cnv,
            "CPX" => SvType::Cpx,
            "BND" => SvType::Bnd,
            "CTX" => SvType::Ctx,
            _ => return None,
        })
    }

    fn is_interchromosomal(self) -> bool {
        matches!(self, SvType::Bnd | SvType::Ctx)
    }
}

/// A closed interval, as the tracks and the records both use.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Interval {
    pub fn length(&self) -> i64 {
        i64::from(self.end - self.start + 1)
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.end && other.start <= self.end
    }

    /// `SimpleInterval.intersect(...).size()`, which is zero when they do not meet.
    pub fn intersection_length(&self, other: &Interval) -> i64 {
        if !self.overlaps(other) {
            return 0;
        }
        i64::from(self.end.min(other.end) - self.start.max(other.start) + 1)
    }
}

/// `IntervalUtils.sortAndMergeIntervals(..., IntervalMergingRule.ALL)`, which merges abutting
/// intervals as well as overlapping ones.
pub fn sort_and_merge(intervals: &[Interval]) -> Vec<Interval> {
    let mut sorted = intervals.to_vec();
    sorted.sort();
    let mut merged: Vec<Interval> = Vec::new();
    for interval in sorted {
        match merged.last_mut() {
            Some(current)
                if current.contig == interval.contig && interval.start <= current.end + 1 =>
            {
                current.end = current.end.max(interval.end);
            }
            _ => merged.push(interval),
        }
    }
    merged
}

/// One record, reduced to what a stratum reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct CallRecord {
    pub id: String,
    pub sv_type: SvType,
    pub contig_a: String,
    pub position_a: i32,
    pub contig_b: String,
    pub position_b: i32,
    /// `getLength()`, which is absent for the types that have no length of their own.
    pub length: Option<i32>,
}

/// What the configuration and the thresholds refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StratifyError {
    MissingColumn {
        column: String,
    },
    /// The message both halves of which read the SAME count.
    ColumnCount {
        count: usize,
    },
    UnknownTrack {
        name: String,
    },
    ReservedName,
    DuplicateTrack {
        name: String,
    },
    TrackCountMismatch,
    MinGreaterThanMax,
    NegativeMax,
    ReservedMax,
    NegativeMin,
    SizedInterchromosomal {
        name: String,
    },
    ThresholdsBothZero,
    MultipleMatches {
        id: String,
        names: Vec<String>,
    },
}

impl StratifyError {
    pub fn message(&self) -> String {
        match self {
            StratifyError::MissingColumn { column } => format!("Missing column {column}"),
            StratifyError::ColumnCount { count } => {
                format!("Expected {count} columns but found {count}")
            }
            StratifyError::UnknownTrack { name } => {
                format!("Could not find track with name {name}")
            }
            StratifyError::ReservedName => format!(
                "Stratification configuration contains entry with reserved ID \"{DEFAULT_STRATUM}\""
            ),
            StratifyError::DuplicateTrack { name } => {
                format!("Duplicate track name was specified: {name}")
            }
            StratifyError::TrackCountMismatch => {
                "Arguments --track-name and --track-intervals must be specified the same number of \
                 times."
                    .to_string()
            }
            StratifyError::MinGreaterThanMax => {
                "Min size must be strictly less than max size".to_string()
            }
            StratifyError::NegativeMax => "Max size cannot be less than 0".to_string(),
            StratifyError::ReservedMax => {
                format!("Max size {} is reserved", i32::MAX)
            }
            StratifyError::NegativeMin => "Min size cannot be less than 0".to_string(),
            StratifyError::SizedInterchromosomal { name } => {
                format!("BND/CTX categories cannot have min or max size ({name})")
            }
            StratifyError::ThresholdsBothZero => {
                "Overlap fraction and overlapping breakpoints thresholds cannot both be 0"
                    .to_string()
            }
            StratifyError::MultipleMatches { id, names } => format!(
                "Record {id} matched multiple groups: {}. Bypass this error using the \
                 --allow-multiple-matches argument",
                names.join(", ")
            ),
        }
    }
}

/// `parseIntegerMaybeNull`.
pub fn parse_integer_maybe_null(value: &str) -> Option<i32> {
    if NULL_TABLE_VALUES.contains(&value) {
        None
    } else {
        Some(value.parse().expect("an integer bound"))
    }
}

/// `parseTrackString`, which refuses a name no track was registered under.
pub fn parse_track_string(value: &str, tracks: &[String]) -> Result<Vec<String>, StratifyError> {
    if NULL_TABLE_VALUES.contains(&value) {
        return Ok(Vec::new());
    }
    let mut names = Vec::new();
    for name in value.split(',') {
        if !tracks.iter().any(|track| track == name) {
            return Err(StratifyError::UnknownTrack {
                name: name.to_string(),
            });
        }
        if !names.contains(&name.to_string()) {
            names.push(name.to_string());
        }
    }
    Ok(names)
}

/// `tableParser`'s two checks, in its own order: every expected column, then the count.
pub fn check_columns(columns: &[String]) -> Result<(), StratifyError> {
    for column in COLUMN_NAMES {
        if !columns.iter().any(|name| name == column) {
            return Err(StratifyError::MissingColumn {
                column: (*column).to_string(),
            });
        }
    }
    if columns.len() != COLUMN_NAMES.len() {
        // Both halves of the message read the SAME number, so this never says what was expected.
        return Err(StratifyError::ColumnCount {
            count: columns.len(),
        });
    }
    Ok(())
}

/// `SVStatificationEngine.Stratum`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Stratum {
    pub name: String,
    pub sv_type: SvType,
    /// Inclusive, `Integer.MIN_VALUE` when the configuration left it null.
    pub min_size: i32,
    /// Exclusive, `Integer.MAX_VALUE` when the configuration left it null.
    pub max_size: i32,
    /// Sorted, as the constructor sorts them.
    pub track_names: Vec<String>,
}

impl Stratum {
    /// The constructor's validations, in its own order.
    pub fn new(
        name: &str,
        sv_type: SvType,
        min_size: Option<i32>,
        max_size: Option<i32>,
        track_names: Vec<String>,
    ) -> Result<Stratum, StratifyError> {
        if let (Some(max), Some(min)) = (max_size, min_size) {
            if max <= min {
                return Err(StratifyError::MinGreaterThanMax);
            }
        }
        if let Some(max) = max_size {
            if max < 0 {
                return Err(StratifyError::NegativeMax);
            }
            if max == i32::MAX {
                return Err(StratifyError::ReservedMax);
            }
        }
        if let Some(min) = min_size {
            if min < 0 {
                return Err(StratifyError::NegativeMin);
            }
        }
        if sv_type.is_interchromosomal() && (min_size.is_some() || max_size.is_some()) {
            return Err(StratifyError::SizedInterchromosomal {
                name: name.to_string(),
            });
        }
        let mut sorted = track_names;
        sorted.sort();
        Ok(Stratum {
            name: name.to_string(),
            sv_type,
            min_size: min_size.unwrap_or(i32::MIN),
            max_size: max_size.unwrap_or(i32::MAX),
            track_names: sorted,
        })
    }

    pub fn matches_type(&self, record: &CallRecord) -> bool {
        record.sv_type == self.sv_type
    }

    /// A record with no length matches ONLY a stratum with neither bound.
    pub fn matches_size(&self, record: &CallRecord) -> bool {
        match record.length {
            None => self.min_size == i32::MIN && self.max_size == i32::MAX,
            Some(length) => length >= self.min_size && length < self.max_size,
        }
    }
}

/// The tracks a run was given, in the order they were named.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Tracks {
    entries: Vec<(String, Vec<Interval>)>,
}

impl Tracks {
    /// `loadStratificationConfig`'s two argument checks.
    pub fn new(names: &[String], intervals: &[Vec<Interval>]) -> Result<Tracks, StratifyError> {
        if names.len() != intervals.len() {
            return Err(StratifyError::TrackCountMismatch);
        }
        let mut entries: Vec<(String, Vec<Interval>)> = Vec::new();
        for (name, list) in names.iter().zip(intervals) {
            if entries.iter().any(|(existing, _)| existing == name) {
                return Err(StratifyError::DuplicateTrack { name: name.clone() });
            }
            entries.push((name.clone(), list.clone()));
        }
        Ok(Tracks { entries })
    }

    pub fn names(&self) -> Vec<String> {
        self.entries.iter().map(|(name, _)| name.clone()).collect()
    }

    fn intervals(&self, name: &str) -> &[Interval] {
        self.entries
            .iter()
            .find(|(existing, _)| existing == name)
            .map(|(_, list)| list.as_slice())
            .expect("a registered track")
    }

    fn overlaps_any(&self, name: &str, interval: &Interval) -> bool {
        self.intervals(name)
            .iter()
            .any(|entry| entry.overlaps(interval))
    }
}

/// The three thresholds, validated the way `matchesTracks` validates them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Thresholds {
    pub overlap_fraction: f64,
    pub num_breakpoint_overlaps: i32,
    pub num_breakpoint_overlaps_interchrom: i32,
}

impl Thresholds {
    pub fn check(&self) -> Result<(), StratifyError> {
        // The only one of the four validations a measured run can reach: the others are guarded by
        // the argument parser's own ranges.
        if self.overlap_fraction == 0.0 && self.num_breakpoint_overlaps == 0 {
            return Err(StratifyError::ThresholdsBothZero);
        }
        Ok(())
    }
}

/// `matchesTracks`, which routes three ways by type before it looks at a threshold.
pub fn matches_tracks(
    stratum: &Stratum,
    record: &CallRecord,
    tracks: &Tracks,
    thresholds: Thresholds,
) -> Result<bool, StratifyError> {
    thresholds.check()?;
    Ok(match record.sv_type {
        // An insertion asks for one end in a track and nothing else.
        SvType::Ins => matches_breakpoint_overlap(stratum, record, tracks, 1),
        SvType::Bnd | SvType::Ctx => matches_breakpoint_overlap(
            stratum,
            record,
            tracks,
            thresholds.num_breakpoint_overlaps_interchrom,
        ),
        _ => {
            matches_overlap_fraction(stratum, record, tracks, thresholds.overlap_fraction)
                && matches_breakpoint_overlap(
                    stratum,
                    record,
                    tracks,
                    thresholds.num_breakpoint_overlaps,
                )
        }
    })
}

/// `matchesTrackOverlapFraction`: the merged union of every named track, compared with `>=`.
pub fn matches_overlap_fraction(
    stratum: &Stratum,
    record: &CallRecord,
    tracks: &Tracks,
    overlap_fraction: f64,
) -> bool {
    if overlap_fraction <= 0.0 || stratum.track_names.is_empty() {
        return true;
    }
    let interval = Interval {
        contig: record.contig_a.clone(),
        start: record.position_a,
        end: record.position_b,
    };
    let mut overlaps: Vec<Interval> = Vec::new();
    for name in &stratum.track_names {
        overlaps.extend(
            tracks
                .intervals(name)
                .iter()
                .filter(|entry| entry.overlaps(&interval))
                .cloned(),
        );
    }
    let merged = sort_and_merge(&overlaps);
    let overlap_length: i64 = merged
        .iter()
        .map(|entry| interval.intersection_length(entry))
        .sum();
    overlap_length as f64 / interval.length() as f64 >= overlap_fraction
}

/// `matchesTrackBreakpointOverlap`, where an end counts once however many tracks cover it.
pub fn matches_breakpoint_overlap(
    stratum: &Stratum,
    record: &CallRecord,
    tracks: &Tracks,
    required: i32,
) -> bool {
    if required <= 0 || stratum.track_names.is_empty() {
        return true;
    }
    let end_a = Interval {
        contig: record.contig_a.clone(),
        start: record.position_a,
        end: record.position_a,
    };
    let end_b = Interval {
        contig: record.contig_b.clone(),
        start: record.position_b,
        end: record.position_b,
    };
    count_any_track_overlap(stratum, tracks, &end_a)
        + count_any_track_overlap(stratum, tracks, &end_b)
        >= required
}

fn count_any_track_overlap(stratum: &Stratum, tracks: &Tracks, interval: &Interval) -> i32 {
    for name in &stratum.track_names {
        if tracks.overlaps_any(name, interval) {
            return 1;
        }
    }
    0
}

/// The strata a run holds, in the order a `java.util.HashMap` gives them back.
#[derive(Debug, Clone, PartialEq)]
pub struct Engine {
    /// Already in HashMap iteration order.
    pub strata: Vec<Stratum>,
    pub tracks: Tracks,
}

impl Engine {
    /// The strata as `strata.values()` walks them, which is a HashMap over the NAMES.
    pub fn new(strata: Vec<Stratum>, tracks: Tracks) -> Result<Engine, StratifyError> {
        if strata.iter().any(|stratum| stratum.name == DEFAULT_STRATUM) {
            return Err(StratifyError::ReservedName);
        }
        let entries: Vec<(Stratum, i32)> = strata
            .into_iter()
            .map(|stratum| {
                let hash = string_hash_code(&stratum.name);
                (stratum, hash)
            })
            .collect();
        let ordered = hash_map_order(&entries).expect("a hashable stratum set");
        Ok(Engine {
            strata: ordered,
            tracks,
        })
    }

    /// `getMatches`, in that same order.
    pub fn matches(
        &self,
        record: &CallRecord,
        thresholds: Thresholds,
    ) -> Result<Vec<&Stratum>, StratifyError> {
        let mut result = Vec::new();
        for stratum in &self.strata {
            if stratum.matches_type(record)
                && stratum.matches_size(record)
                && matches_tracks(stratum, record, &self.tracks, thresholds)?
            {
                result.push(stratum);
            }
        }
        Ok(result)
    }
}

/// Where one record's copies are written, and under what stratum name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Written {
    /// The file, which is `DEFAULT_STRATUM` for every record unless the output is split.
    pub file: String,
    pub stratum: String,
}

/// `apply`: the record goes to the default group when nothing matched, and once per stratum when
/// something did.
pub fn apply(
    engine: &Engine,
    record: &CallRecord,
    thresholds: Thresholds,
    allow_multiple_matches: bool,
    split_output: bool,
) -> Result<Vec<Written>, StratifyError> {
    let matches = engine.matches(record, thresholds)?;
    if matches.is_empty() {
        return Ok(vec![Written {
            file: DEFAULT_STRATUM.to_string(),
            stratum: DEFAULT_STRATUM.to_string(),
        }]);
    }
    if !allow_multiple_matches && matches.len() > 1 {
        return Err(StratifyError::MultipleMatches {
            id: record.id.clone(),
            names: matches.iter().map(|s| s.name.clone()).collect(),
        });
    }
    Ok(matches
        .into_iter()
        .map(|stratum| Written {
            file: if split_output {
                stratum.name.clone()
            } else {
                DEFAULT_STRATUM.to_string()
            },
            stratum: stratum.name.clone(),
        })
        .collect())
}

/// `generateGroupOutputPath`, and the order `initializeWriters` creates the files in: the default
/// group first, then the strata in the engine's order.
pub fn split_output_files(engine: &Engine, prefix: &str) -> Vec<String> {
    let mut names = vec![format!("{prefix}.{DEFAULT_STRATUM}.vcf.gz")];
    for stratum in &engine.strata {
        names.push(format!("{prefix}.{}.vcf.gz", stratum.name));
    }
    names
}
