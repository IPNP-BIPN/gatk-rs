//! `MarkDuplicatesSpark`: `MarkDuplicatesSparkUtils` and the records it groups, which are GATK's own
//! and not Picard's, over Picard's `OpticalDuplicateFinder`.
//!
//! Two things make this differ from `picard MarkDuplicates`. A read's verdict is its NAME's: the
//! grouping answers with the names that are not duplicates, and every read of such a name, secondary
//! and supplementary ones included, follows it. And the winner of a group is decided by the score
//! and then by the location parsed out of the name, so the order the reads came in only ever shows
//! in which optical duplicate of a set is the one kept.
//!
//! Ported from `org.broadinstitute.hellbender.tools.spark.transforms.markduplicates`,
//! `...utils.read.markduplicates` and its `sparkrecords` (GATK 4.6.2.0), and
//! `picard.sam.markduplicates.util.OpticalDuplicateFinder` and `picard.sam.util.ReadNameParser`.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use gatk_engine::read;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricBean, MetricsFile, Value};

/// `LibraryIdGenerator.UNKNOWN_LIBRARY`.
pub const UNKNOWN_LIBRARY: &str = "Unknown Library";
/// `MarkDuplicates.DUPLICATE_TYPE_TAG`, and its two values.
const DUPLICATE_TYPE: [u8; 2] = *b"DT";
const SEQUENCING: &str = "SQ";
const LIBRARY: &str = "LB";
/// The markers `mark` reads back off a name: never a duplicate, and an optical duplicate.
const NO_OPTICAL_MARKER: i32 = -1;
const OPTICAL_DUPLICATE_MARKER: i32 = -2;
/// `FlowBasedReadUtils.FLOW_BASED_INSIGNIFICANT_END`.
const INSIGNIFICANT_END: i32 = 0;
/// `DEFAULT_MAX_DUPLICATE_SET_SIZE`: a longer list is not looked at.
const MAX_DUPLICATE_SET_SIZE: usize = 300_000;

/// `MarkDuplicatesScoringStrategy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scoring {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
    FlowSumOfBaseQualities,
}

/// `MarkDuplicates.DuplicateTaggingPolicy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tagging {
    DontTag,
    OpticalOnly,
    All,
}

/// `READ_NAME_REGEX`.
#[derive(Debug, Clone)]
pub enum ReadNameRegex {
    /// The optimized capture of the last three colon-separated fields.
    Default,
    /// `null`: no locations at all.
    Off,
    /// The pattern, anchored at both ends because `Matcher.matches` is.
    Pattern(regex::Regex),
}

impl ReadNameRegex {
    /// The argument as given: the default's own text, nothing where it was set to `null`, or a
    /// pattern; one that does not compile is `Pattern.compile`'s failure, whose text it returns.
    pub fn from_arg(value: Option<&str>) -> Result<Self, String> {
        match value {
            Some("<optimized capture of last three ':' separated fields as numeric values>") => {
                Ok(ReadNameRegex::Default)
            }
            None => Ok(ReadNameRegex::Off),
            Some(pattern) => regex::Regex::new(&format!("^(?:{pattern})$"))
                .map(ReadNameRegex::Pattern)
                .map_err(|error| error.to_string()),
        }
    }
}

/// What decides the marking.
#[derive(Debug, Clone)]
pub struct MarkArguments {
    pub scoring: Scoring,
    pub dont_mark_unmapped_mates: bool,
    pub tagging: Tagging,
    pub read_name_regex: ReadNameRegex,
    pub pixel_distance: i32,
    pub flow: FlowArguments,
}

/// The flow half of `MarkDuplicatesSparkArgumentCollection`.
#[derive(Debug, Clone, Default)]
pub struct FlowArguments {
    pub quality_sum_strategy: bool,
    pub end_location_significant: bool,
    pub end_uncertainty: i32,
    pub use_clipped_locations: bool,
    pub skip_start_homopolymers: i32,
    pub q_is_known_end: bool,
}

impl FlowArguments {
    /// `isFlowEnabled`, which `FLOW_Q_IS_KNOWN_END` and the uncertainty alone do not turn on.
    fn enabled(&self) -> bool {
        self.quality_sum_strategy
            || self.end_location_significant
            || self.use_clipped_locations
            || self.skip_start_homopolymers != 0
    }
}

/// A failure inside a Spark task: the job grouping the names, or, `late`, the one zipping their
/// answers back onto the reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaskFailure {
    pub message: String,
    pub late: bool,
}

fn early(message: String) -> TaskFailure {
    TaskFailure {
        message,
        late: false,
    }
}

/// `MarkDuplicatesScoringStrategy.score`: capped so that a pair's two cannot overflow, and a read
/// that failed the vendor check loses half the range.
fn score(record: &BamRecord, scoring: Scoring) -> i16 {
    let raw: i32 = match scoring {
        Scoring::SumOfBaseQualities => record
            .base_qualities
            .iter()
            .map(|quality| i32::from(*quality as i8))
            .filter(|quality| *quality >= 15)
            .sum(),
        Scoring::TotalMappedReferenceLength => {
            if read::is_unmapped(record) {
                0
            } else {
                record.cigar.reference_length() as i32
            }
        }
        Scoring::FlowSumOfBaseQualities => {
            // `FlowBasedReadUtils.flowSumOfBaseQualities`: a homopolymer scores its first quality.
            let mut sum = 0;
            let (mut last, mut effective) = (0u8, 0i32);
            for (index, base) in record.read_bases.iter().enumerate() {
                if *base != last {
                    effective = i32::from(record.base_qualities[index] as i8);
                }
                if effective >= 15 {
                    sum += effective;
                }
                last = *base;
            }
            sum
        }
    };
    let failed = if read::fails_vendor_quality_check(record) {
        i32::from(i16::MIN / 2)
    } else {
        0
    };
    (raw.min(i32::from(i16::MAX / 2)) + failed) as i16
}

/// `ReadUtils.readHasMappedMate`.
fn has_mapped_mate(record: &BamRecord) -> bool {
    read::is_paired(record) && !read::mate_is_unmapped(record)
}

/// `ReadUtils.getStrandedUnclippedStart`.
fn stranded_unclipped_start(record: &BamRecord) -> i32 {
    if read::is_reverse_strand(record) {
        gatk_engine::read_utils::unclipped_end(record)
    } else {
        gatk_engine::read_utils::unclipped_start(record)
    }
}

/// `FlowBasedRead.CLIPPING_TAG_NAME`.
fn clipping_tag(record: &BamRecord) -> Option<&str> {
    match record.tags.get(Tag::new(b"tm")) {
        Some(TagValue::Str(text)) => Some(text),
        _ => None,
    }
}

/// `FlowBasedReadUtils.readEndMarkedUnclipped`.
fn end_marked_unclipped(record: &BamRecord, q_is_known_end: bool) -> bool {
    clipping_tag(record).is_some_and(|tm| tm.contains('A') || (q_is_known_end && tm.contains('Q')))
}

/// `FlowBasedReadUtils.readEndMarkedUncertain`.
fn end_marked_uncertain(record: &BamRecord) -> bool {
    clipping_tag(record).is_some_and(|tm| tm.contains('Q') || tm.contains('Z'))
}

/// `FlowBasedReadUtils.getReadFlowOrder`: the read's group's, else the first group's that has one.
fn flow_order<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a str> {
    read_group_of(record)
        .and_then(|id| header.read_groups.iter().find(|rg| rg.id == id))
        .and_then(|group| group.attributes.get("FO"))
        .or_else(|| {
            header
                .read_groups
                .iter()
                .find_map(|group| group.attributes.get("FO"))
        })
}

/// The length `getMarkDupReadStart` and `getMarkDupReadEnd` skip off the read's start or end: whole
/// homopolymers, with the flows that read nothing counted against the number asked for.
fn skipped_homopolymers(
    record: &BamRecord,
    order: Option<&str>,
    skip: i32,
    from_end: bool,
) -> Result<i32, TaskFailure> {
    let bases = &record.read_bases;
    let at = |index: usize| {
        bases[if from_end {
            bases.len() - 1 - index
        } else {
            index
        }]
    };
    if bases.is_empty() {
        return Err(early(
            "java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0"
                .to_string(),
        ));
    }
    let order = order.map(str::as_bytes);
    // The reference never comes back from a base its flow order does not hold.
    if order.is_some_and(|flows| !bases.iter().all(|base| flows.contains(base))) {
        return Err(early(format!(
            "Read {} holds a base outside its flow order, where the reference never returns. This \
             message is the port's own and not GATK's.",
            record.read_name
        )));
    }
    let mut base = at(0);
    let mut flow = 0usize;
    let mut left = skip;
    if let Some(flows) = order {
        while flows[flow] != base {
            flow = (flow + 1) % flows.len();
            left -= 1;
        }
    }
    let mut size = 1usize;
    while size < bases.len() {
        if at(size) != base {
            left -= 1;
            if left <= 0 {
                break;
            }
            base = at(size);
            if let Some(flows) = order {
                flow = (flow + 1) % flows.len();
                while flows[flow] != base {
                    left -= 1;
                    flow = (flow + 1) % flows.len();
                }
                if left <= 0 {
                    break;
                }
            }
        }
        size += 1;
    }
    Ok(size as i32)
}

/// `FlowBasedReadUtils.getMarkDupReadStart`.
fn flow_read_start(
    record: &BamRecord,
    end_semantics: bool,
    header: &SamHeader,
    flow: &FlowArguments,
) -> Result<i32, TaskFailure> {
    let unclipped = gatk_engine::read_utils::unclipped_start(record);
    if !end_semantics && flow.skip_start_homopolymers != 0 {
        let start = unclipped
            + skipped_homopolymers(
                record,
                flow_order(record, header),
                flow.skip_start_homopolymers,
                false,
            )?;
        return Ok(if flow.use_clipped_locations {
            start.max(gatk_engine::read_utils::start(record))
        } else {
            start
        });
    }
    Ok(if end_marked_unclipped(record, flow.q_is_known_end) {
        unclipped
    } else if end_semantics && end_marked_uncertain(record) {
        INSIGNIFICANT_END
    } else if flow.use_clipped_locations {
        gatk_engine::read_utils::start(record)
    } else {
        unclipped
    })
}

/// `FlowBasedReadUtils.getMarkDupReadEnd`.
fn flow_read_end(
    record: &BamRecord,
    end_semantics: bool,
    header: &SamHeader,
    flow: &FlowArguments,
) -> Result<i32, TaskFailure> {
    let unclipped = gatk_engine::read_utils::unclipped_end(record);
    if !end_semantics && flow.skip_start_homopolymers != 0 {
        let end = unclipped
            - skipped_homopolymers(
                record,
                flow_order(record, header),
                flow.skip_start_homopolymers,
                true,
            )?;
        return Ok(if flow.use_clipped_locations {
            end.min(gatk_engine::read_utils::end(record))
        } else {
            end
        });
    }
    Ok(if end_marked_unclipped(record, flow.q_is_known_end) {
        unclipped
    } else if end_semantics && end_marked_uncertain(record) {
        INSIGNIFICANT_END
    } else if flow.use_clipped_locations {
        gatk_engine::read_utils::end(record)
    } else {
        unclipped
    })
}

/// `FlowModeFragment.getFlowSumOfBaseQualities`: a homopolymer scores its first quality.
fn flow_fragment_sum(record: &BamRecord, start: i32, end: i32) -> i32 {
    let (start, end) = if read::is_reverse_strand(record) {
        (end, start)
    } else {
        (start, end)
    };
    let first = (start - gatk_engine::read_utils::unclipped_start(record)).max(0) as usize;
    let trimmed = (gatk_engine::read_utils::unclipped_end(record) - end).max(0) as usize;
    let mut sum = 0;
    let (mut last, mut effective) = (0u8, 0i32);
    for index in first..record.read_bases.len().saturating_sub(trimmed) {
        let base = record.read_bases[index];
        if base != last {
            effective = i32::from(record.base_qualities[index] as i8);
        }
        if effective >= 15 {
            sum += effective;
        }
        last = base;
    }
    sum
}

/// `ReadsKey.longKeyForFragment`, sign extension included: its mask keeps an int negative when the
/// reference index is, so such a read's key loses its start.
pub fn fragment_key(start: i32, reverse: bool, reference_index: i32, library: u8) -> i64 {
    (i64::from(start) << 32)
        | i64::from((reference_index << 16) & 0xFFFF_0000_u32 as i32)
        | i64::from((i32::from(library as i8) << 8) & 0xFF00)
        | i64::from(reverse)
}

/// `ReadsKey.longKeyForPair`.
fn pair_key(start: i32, reverse: bool, reference_index: i32) -> i64 {
    (i64::from(start) << 32)
        | i64::from((reference_index << 16) & 0xFFFF_0000_u32 as i32)
        | i64::from(reverse)
}

/// `ReadsKey`, without the passthrough's hash: a passthrough is answered wherever it lands.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Key {
    Fragment(i64),
    Pair(i64, i64),
}

/// `TransientFieldPhysicalLocation`, the part `ReadNameParser` fills. `-1` throughout until a name
/// parses, and the pixels are compressed to signed shorts for space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Location {
    tile: i16,
    x: i32,
    y: i32,
}

impl Location {
    const NONE: Location = Location {
        tile: -1,
        x: -1,
        y: -1,
    };

    /// `PhysicalLocation.hasLocation`: `NO_VALUE` is a tile too, whatever the name said.
    fn known(&self) -> bool {
        self.tile != -1
    }

    fn of(tile: i32, x: i32, y: i32) -> Location {
        Location {
            tile: tile as i16,
            x: i32::from(x as i16),
            y: i32::from(y as i16),
        }
    }
}

/// `ReadNameParser.rapidParseInt`: digits until the first thing that is not one, and no digits at
/// all is the `NumberFormatException` the caller turns into no location.
fn rapid_parse_int(text: &str) -> Option<i32> {
    let bytes = text.as_bytes();
    let negative = bytes.first() == Some(&b'-');
    let mut value: i32 = 0;
    let mut digits = false;
    for byte in &bytes[usize::from(negative)..] {
        if !byte.is_ascii_digit() {
            break;
        }
        value = value.wrapping_mul(10).wrapping_add(i32::from(byte - b'0'));
        digits = true;
    }
    digits.then_some(if negative { -value } else { value })
}

/// `Integer.parseInt`.
fn java_parse_int(text: &str) -> Option<i32> {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}

/// `ReadNameParser.readLocationInformation`.
fn location(name: &str, regex: &ReadNameRegex) -> Location {
    match regex {
        ReadNameRegex::Off => Location::NONE,
        ReadNameRegex::Default => {
            let fields: Vec<&str> = name.split(':').collect();
            if fields.len() != 5 && fields.len() != 7 {
                return Location::NONE;
            }
            let parsed: Vec<Option<i32>> = fields[fields.len() - 3..]
                .iter()
                .map(|field| rapid_parse_int(field))
                .collect();
            match parsed[..] {
                [Some(tile), Some(x), Some(y)] => Location::of(tile, x, y),
                _ => Location::NONE,
            }
        }
        ReadNameRegex::Pattern(pattern) => {
            let Some(groups) = pattern.captures(name) else {
                return Location::NONE;
            };
            let field = |index: usize| groups.get(index).and_then(|m| java_parse_int(m.as_str()));
            match (field(1), field(2), field(3)) {
                (Some(tile), Some(x), Some(y)) => Location::of(tile, x, y),
                _ => Location::NONE,
            }
        }
    }
}

/// What a `Fragment` or a `Pair` holds once the name has been parsed.
#[derive(Debug, Clone)]
struct Scored {
    partition: usize,
    name: String,
    score: i16,
    first_reverse: bool,
    second_reverse: bool,
    was_flipped: bool,
    read_group: i32,
    location: Location,
    end: i32,
}

impl Scored {
    /// `Pair.getOrientationForOpticalDuplicates`, as far as `handlePairs` looks: an FR and an RF
    /// are the two it counts apart, and any other is neither.
    fn optical_orientation(&self) -> Option<bool> {
        match (self.first_reverse, self.second_reverse) {
            (true, false) => Some(self.was_flipped),
            (false, true) => Some(!self.was_flipped),
            _ => None,
        }
    }
}

/// `PAIRED_ENDS_SCORE_COMPARATOR`, the location reversed; `max` keeps the later of two equal.
fn better(first: &Scored, second: &Scored) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let located = |a: &Scored, b: &Scored| -> Ordering {
        if a.first_reverse != b.first_reverse {
            return if a.first_reverse {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        a.location
            .tile
            .cmp(&b.location.tile)
            .then(a.location.x.cmp(&b.location.x))
            .then(a.location.y.cmp(&b.location.y))
            .then_with(|| a.name.cmp(&b.name))
    };
    first.score.cmp(&second.score).then(located(second, first))
}

/// `stream.max`: `reduce` over `BinaryOperator.maxBy`, which keeps the left one on a tie.
fn best(group: &[Scored]) -> usize {
    let mut kept = 0;
    for index in 1..group.len() {
        if better(&group[kept], &group[index]).is_lt() {
            kept = index;
        }
    }
    kept
}

/// `OpticalDuplicateFinder.closeEnough`.
fn close_enough(left: &Scored, right: &Scored, distance: i32) -> bool {
    left.location.known()
        && right.location.known()
        && left.read_group == right.read_group
        && left.location.tile == right.location.tile
        && (left.location.x - right.location.x).abs() <= distance
        && (left.location.y - right.location.y).abs() <= distance
}

/// `OpticalDuplicateFinder.findOpticalDuplicates`, the keeper by its index if it is in `list` at all.
///
/// The fast path's threshold is read off the keeper passed in, which here is always one, and not
/// off the one that survives `keeperOrNull`.
fn optical_duplicates(
    list: &[&Scored],
    keeper: Option<usize>,
    regex: &ReadNameRegex,
    distance: i32,
) -> Vec<bool> {
    let mut flags = vec![false; list.len()];
    if matches!(regex, ReadNameRegex::Off) || list.len() < 2 || list.len() > MAX_DUPLICATE_SET_SIZE
    {
        return flags;
    }
    let keeper = keeper.filter(|index| list[*index].location.known());
    if list.len() >= 4 {
        return optical_duplicates_with_graph(list, keeper, distance);
    }
    if let Some(keeper) = keeper {
        for (index, other) in list.iter().enumerate() {
            flags[index] = index != keeper && close_enough(list[keeper], other, distance);
        }
    }
    for left in 0..list.len() {
        if Some(left) == keeper {
            continue;
        }
        for right in (left + 1)..list.len() {
            if Some(right) == keeper || (flags[left] && flags[right]) {
                continue;
            }
            if close_enough(list[left], list[right], distance) {
                let index = if flags[right] { left } else { right };
                flags[index] = true;
            }
        }
    }
    flags
}

/// `getOpticalDuplicatesFlagWithGraph`: union-find over the ends within the distance of each other
/// inside a bucket of `(read group << 16) + tile`, then one read kept per cluster, the keeper in its
/// own and the smallest x, then y, in any other.
fn optical_duplicates_with_graph(
    list: &[&Scored],
    keeper: Option<usize>,
    distance: i32,
) -> Vec<bool> {
    fn find(parent: &mut [usize], index: usize) -> usize {
        let mut root = index;
        while parent[root] != root {
            root = parent[root];
        }
        let mut walk = index;
        while parent[walk] != root {
            let next = parent[walk];
            parent[walk] = root;
            walk = next;
        }
        root
    }

    let mut flags = vec![false; list.len()];
    let mut parent: Vec<usize> = (0..list.len()).collect();
    let mut buckets: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (index, end) in list.iter().enumerate() {
        if end.location.known() {
            let key = (end.read_group << 16).wrapping_add(i32::from(end.location.tile));
            buckets.entry(key).or_default().push(index);
        }
    }
    for members in buckets.values() {
        for (position, left) in members.iter().enumerate() {
            for right in &members[position + 1..] {
                if (list[*left].location.x - list[*right].location.x).abs() <= distance
                    && (list[*left].location.y - list[*right].location.y).abs() <= distance
                {
                    let (a, b) = (find(&mut parent, *left), find(&mut parent, *right));
                    if a != b {
                        parent[a] = b;
                    }
                }
            }
        }
    }
    let cluster: Vec<usize> = (0..list.len())
        .map(|index| find(&mut parent, index))
        .collect();

    let mut representative: Vec<Option<usize>> = vec![None; list.len()];
    if let Some(keeper) = keeper {
        representative[cluster[keeper]] = Some(keeper);
    }
    let keeper_cluster = keeper.map(|index| cluster[index]);
    for index in 0..list.len() {
        let group = cluster[index];
        match representative[group] {
            Some(held) if Some(index) != keeper => {
                let (here, there) = (list[index].location, list[held].location);
                let smaller = here.x < there.x || (here.x == there.x && here.y < there.y);
                if keeper_cluster != Some(group) && smaller {
                    flags[held] = true;
                    representative[group] = Some(index);
                } else {
                    flags[index] = true;
                }
            }
            Some(_) => {}
            None => representative[group] = Some(index),
        }
    }
    flags
}

/// `MarkDuplicatesSparkUtils.getLibraryForRead`, and the two refusals behind it.
fn library_of<'a>(record: &BamRecord, header: &'a SamHeader) -> Result<&'a str, TaskFailure> {
    let group = read_group_of(record);
    match group.and_then(|id| header.read_groups.iter().find(|rg| rg.id == id)) {
        Some(known) => Ok(known.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY)),
        None => Err(early(missing_read_group(record))),
    }
}

fn read_group_of(record: &BamRecord) -> Option<&str> {
    match record.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => Some(id),
        _ => None,
    }
}

/// `UserException.ReadMissingReadGroup` and `HeaderMissingReadGroup`.
fn missing_read_group(record: &BamRecord) -> String {
    let name = &record.read_name;
    match read_group_of(record) {
        None => format!(
            "SAM/BAM/CRAM file (unknown) is malformed: Read {name} is missing the read group (RG) \
             tag, which is required by the GATK.  Please use \
             http://gatkforums.broadinstitute.org/discussion/59/\
             companion-utilities-replacereadgroups to fix this problem"
        ),
        Some(id) => format!(
            "SAM/BAM/CRAM file (unknown) is malformed: Read {name} contains an (RG) tag with the \
             group {id} which is not found in the file header."
        ),
    }
}

/// `constructLibraryIndex`: the libraries in the order the read groups first name them.
pub fn library_index(header: &SamHeader) -> Vec<String> {
    let mut libraries: Vec<String> = Vec::new();
    for group in &header.read_groups {
        let library = group.attributes.get("LB").unwrap_or(UNKNOWN_LIBRARY);
        if !libraries.iter().any(|known| known == library) {
            libraries.push(library.to_string());
        }
    }
    libraries
}

/// One key's records: `splitByType`, an empty fragment counting only for being there.
#[derive(Default)]
struct KeyGroup {
    fragments: Vec<Scored>,
    pairs: Vec<Scored>,
    has_empty_fragment: bool,
}

/// `transformToDuplicateNames`: for each partition, the names that are not duplicates, with what
/// `mark` reads off each. A partition's reads must be grouped by name.
fn names_of_non_duplicates(
    partitions: &[Vec<BamRecord>],
    header: &SamHeader,
    arguments: &MarkArguments,
    optical: bool,
) -> Result<Vec<HashMap<String, i32>>, TaskFailure> {
    let libraries = library_index(header);
    let library_byte = |record: &BamRecord| -> Result<u8, TaskFailure> {
        let library = library_of(record, header)?;
        Ok(libraries.iter().position(|l| l == library).unwrap_or(0) as u8)
    };
    let reference = |record: &BamRecord| record.reference_index;
    let single_key = |record: &BamRecord, library: u8| {
        fragment_key(
            stranded_unclipped_start(record),
            read::is_reverse_strand(record),
            reference(record),
            library,
        )
    };

    let mut groups: HashMap<Key, KeyGroup> = HashMap::new();
    let mut names: Vec<(usize, String, i32)> = Vec::new();
    for (partition, reads) in partitions.iter().enumerate() {
        let mapped: Vec<&BamRecord> = reads.iter().filter(|r| !read::is_unmapped(r)).collect();
        // `spanByKey`: neighbours by name, so a name that comes back later is another group.
        for span in mapped.chunk_by(|a, b| a.read_name == b.read_name) {
            let mut had_non_primary = None;
            let mut primaries: Vec<&BamRecord> = Vec::new();
            for record in span {
                if read::is_secondary_alignment(record) || read::is_supplementary_alignment(record)
                {
                    had_non_primary = Some(*record);
                    continue;
                }
                let library = library_byte(record)?;
                if has_mapped_mate(record) {
                    let key = Key::Fragment(single_key(record, library));
                    groups.entry(key).or_default().has_empty_fragment = true;
                } else {
                    let flow = &arguments.flow;
                    let (start, end, fragment_score) = if flow.enabled() {
                        // `FlowModeFragment`, which keys on a start of its own.
                        let reverse = read::is_reverse_strand(record);
                        let start = if reverse {
                            flow_read_end(record, false, header, flow)?
                        } else {
                            flow_read_start(record, false, header, flow)?
                        };
                        let end = if !flow.end_location_significant {
                            INSIGNIFICANT_END
                        } else if reverse {
                            flow_read_start(record, true, header, flow)?
                        } else {
                            flow_read_end(record, true, header, flow)?
                        };
                        let fragment_score = if end == INSIGNIFICANT_END {
                            -1
                        } else if flow.quality_sum_strategy
                            && crate::flow_based_read::has_flow_tags(record)
                        {
                            flow_fragment_sum(record, start, end).min(i32::from(i16::MAX / 2))
                                as i16
                        } else {
                            score(record, arguments.scoring)
                        };
                        (start, end, fragment_score)
                    } else {
                        (
                            stranded_unclipped_start(record),
                            INSIGNIFICANT_END,
                            score(record, arguments.scoring),
                        )
                    };
                    let key = Key::Fragment(fragment_key(
                        start,
                        read::is_reverse_strand(record),
                        reference(record),
                        library,
                    ));
                    groups.entry(key).or_default().fragments.push(Scored {
                        partition,
                        name: record.read_name.clone(),
                        score: fragment_score,
                        first_reverse: read::is_reverse_strand(record),
                        second_reverse: false,
                        was_flipped: false,
                        read_group: -1,
                        location: Location::NONE,
                        end,
                    });
                }
                primaries.push(record);
            }
            if primaries.is_empty() {
                if let Some(read) = had_non_primary {
                    names.push((partition, read.read_name.clone(), NO_OPTICAL_MARKER));
                }
                continue;
            }
            if primaries.len() > 2 {
                return Err(early(format!(
                    "MarkDuplicatesSpark only supports singleton fragments and pairs. We found the \
                     following group with >2 primary reads: ( {} number of reads). \n{}.",
                    primaries.len(),
                    primaries
                        .iter()
                        .map(|record| format!(
                                "indexpair[{partition},{}]",
                                gatk_engine::assembly_region::render_read(record, header)
                            ))
                        .collect::<Vec<_>>()
                        .join("\n")
                )));
            }
            let mates: Vec<&BamRecord> = primaries
                .iter()
                .copied()
                .filter(|record| has_mapped_mate(record))
                .collect();
            match mates[..] {
                [read1, read2] => {
                    let (key, scored) =
                        new_pair(read1, read2, partition, header, arguments, &library_byte)?;
                    groups.entry(key).or_default().pairs.push(scored);
                }
                [only] => names.push((partition, only.read_name.clone(), NO_OPTICAL_MARKER)),
                _ => {}
            }
        }
    }

    for group in groups.into_values() {
        if !group.fragments.is_empty() && !group.has_empty_fragment {
            if arguments.flow.end_location_significant {
                handle_fragments_with_end(group.fragments, arguments, &mut names);
            } else {
                names.push(handle_fragments(group.fragments, arguments));
            }
        }
        if !group.pairs.is_empty() {
            handle_pairs(group.pairs, arguments, optical, &mut names);
        }
    }

    let mut answers: Vec<HashMap<String, i32>> = vec![HashMap::new(); partitions.len()];
    for (partition, name, value) in names {
        if answers[partition].insert(name.clone(), value).is_some() {
            return Err(TaskFailure {
                message: format!(
                    "Detected multiple mark duplicate records objects corresponding to read with \
                     name '{name}', this could be the result of the file sort order being \
                     incorrect or that a previous tool has let readnames span multiple partitions"
                ),
                late: true,
            });
        }
    }
    Ok(answers)
}

/// `handleFragments`: the best of them, once the names are parsed.
fn handle_fragments(fragments: Vec<Scored>, arguments: &MarkArguments) -> (usize, String, i32) {
    let located: Vec<Scored> = fragments
        .into_iter()
        .map(|fragment| Scored {
            location: location(&fragment.name, &arguments.read_name_regex),
            ..fragment
        })
        .collect();
    let kept = &located[best(&located)];
    (kept.partition, kept.name.clone(), NO_OPTICAL_MARKER)
}

/// `handleFragmentsWithEndPosition`: by end, then runs whose ends lie within the uncertainty of
/// the run so far, a best from each. An insignificant end goes with whichever run is open.
fn handle_fragments_with_end(
    mut fragments: Vec<Scored>,
    arguments: &MarkArguments,
    names: &mut Vec<(usize, String, i32)>,
) {
    if let [only] = &fragments[..] {
        names.push((only.partition, only.name.clone(), 0));
        return;
    }
    fragments.sort_by_key(|fragment| fragment.end);
    let uncertainty = arguments.flow.end_uncertainty;
    let mut run: Vec<Scored> = Vec::new();
    let (mut low, mut high) = (0, 0);
    for fragment in fragments {
        let end = fragment.end;
        if run.is_empty() {
            if end != INSIGNIFICANT_END {
                low = end - uncertainty;
                high = end + uncertainty;
            }
        } else if end == INSIGNIFICANT_END {
            // Simply accumulated.
        } else if low == 0 {
            low = end - uncertainty;
            high = end + uncertainty;
        } else if end >= low && end <= high {
            low = low.min(end - uncertainty);
            high = high.max(end + uncertainty);
        } else {
            names.push(handle_fragments(std::mem::take(&mut run), arguments));
            low = end - uncertainty;
            high = end + uncertainty;
        }
        run.push(fragment);
    }
    if !run.is_empty() {
        names.push(handle_fragments(run, arguments));
    }
}

/// `Pair`'s constructor: the read further left first, with FR forced when both start alike.
fn new_pair(
    read1: &BamRecord,
    read2: &BamRecord,
    partition: usize,
    header: &SamHeader,
    arguments: &MarkArguments,
    library_byte: &dyn Fn(&BamRecord) -> Result<u8, TaskFailure>,
) -> Result<(Key, Scored), TaskFailure> {
    let (start1, start2) = (
        stranded_unclipped_start(read1),
        stranded_unclipped_start(read2),
    );
    let (index1, index2) = (read1.reference_index, read2.reference_index);
    let in_order = if index1 != index2 {
        index1 < index2
    } else {
        start1 <= start2
    };
    let (mut first, mut second) = if in_order {
        (read1, read2)
    } else {
        (read2, read1)
    };
    if index1 == index2
        && start1 == start2
        && read::is_reverse_strand(first)
        && !read::is_reverse_strand(second)
    {
        std::mem::swap(&mut first, &mut second);
    }
    let key = Key::Pair(
        fragment_key(
            stranded_unclipped_start(first),
            read::is_reverse_strand(first),
            first.reference_index,
            library_byte(first)?,
        ),
        pair_key(
            stranded_unclipped_start(second),
            read::is_reverse_strand(second),
            second.reference_index,
        ),
    );
    // The read group is the FIRST read's in the template's order, not the leftmost's.
    let read_group = read_group_of(read1)
        .and_then(|id| header.read_groups.iter().position(|rg| rg.id == id))
        .ok_or_else(|| early(missing_read_group(read1)))?;
    let scored = Scored {
        partition,
        name: read1.read_name.clone(),
        score: score(read1, arguments.scoring).wrapping_add(score(read2, arguments.scoring)),
        first_reverse: read::is_reverse_strand(first),
        second_reverse: read::is_reverse_strand(second),
        was_flipped: read::is_first_of_pair(second),
        read_group: i32::from(read_group as i16),
        location: Location::NONE,
        end: INSIGNIFICANT_END,
    };
    Ok((key, scored))
}

/// `handlePairs`, and `countOpticalDuplicates` under it. One pair is settled before any name is
/// parsed; FR and RF are counted apart only when both are there.
fn handle_pairs(
    mut pairs: Vec<Scored>,
    arguments: &MarkArguments,
    optical: bool,
    names: &mut Vec<(usize, String, i32)>,
) {
    if pairs.len() == 1 {
        names.push((pairs[0].partition, pairs[0].name.clone(), 0));
        return;
    }
    for pair in &mut pairs {
        pair.location = location(&pair.name, &arguments.read_name_regex);
    }
    let kept = best(&pairs);
    let mut count = |members: &[usize]| {
        let list: Vec<&Scored> = members.iter().map(|index| &pairs[*index]).collect();
        let keeper = members.iter().position(|index| *index == kept);
        let flags = optical_duplicates(
            &list,
            keeper,
            &arguments.read_name_regex,
            arguments.pixel_distance,
        );
        let mut total = 0;
        for (position, flag) in flags.iter().enumerate() {
            if *flag {
                total += 1;
                if optical {
                    let flagged = list[position];
                    names.push((
                        flagged.partition,
                        flagged.name.clone(),
                        OPTICAL_DUPLICATE_MARKER,
                    ));
                }
            }
        }
        total
    };
    let orientation = |wanted: bool| -> Vec<usize> {
        (0..pairs.len())
            .filter(|index| pairs[*index].optical_orientation() == Some(wanted))
            .collect()
    };
    let (fr, rf) = (orientation(false), orientation(true));
    let total = if !fr.is_empty() && !rf.is_empty() {
        count(&fr) + count(&rf)
    } else {
        count(&(0..pairs.len()).collect::<Vec<_>>())
    };
    names.push((pairs[kept].partition, pairs[kept].name.clone(), total));
}

/// The marked reads with, for each, the optical count `mark` leaves in its transient attribute.
pub struct Marked {
    pub reads: Vec<BamRecord>,
    optical_counts: Vec<i32>,
}

/// `MarkDuplicatesSpark.mark`, over partitions already ordered and grouped by name.
pub fn mark(
    partitions: Vec<Vec<BamRecord>>,
    header: &SamHeader,
    arguments: &MarkArguments,
) -> Result<Marked, TaskFailure> {
    let optical = arguments.tagging != Tagging::DontTag;
    let mark_unmapped_mates = !arguments.dont_mark_unmapped_mates;
    let mut answers = names_of_non_duplicates(&partitions, header, arguments, optical)?;
    let tag = Tag::new(&DUPLICATE_TYPE);
    let mut reads = Vec::new();
    let mut optical_counts = Vec::new();
    for (partition, records) in partitions.into_iter().enumerate() {
        for mut record in records {
            let answers = &mut answers[partition];
            record.flags &= !DUPLICATE;
            record.tags.remove(tag);
            let mut count = NO_OPTICAL_MARKER;
            match answers.get(&record.read_name).copied() {
                Some(OPTICAL_DUPLICATE_MARKER) => {
                    record.flags |= DUPLICATE;
                    record
                        .tags
                        .insert(tag, TagValue::Str(SEQUENCING.to_string()));
                }
                Some(_) => {
                    if mark_unmapped_mates || !read::is_unmapped(&record) {
                        if let Some(value) =
                            answers.insert(record.read_name.clone(), NO_OPTICAL_MARKER)
                        {
                            count = value;
                        }
                    }
                }
                None => {
                    let unmapped_pair = read::is_unmapped(&record)
                        && (!read::is_paired(&record) || read::mate_is_unmapped(&record));
                    if !unmapped_pair && (mark_unmapped_mates || !read::is_unmapped(&record)) {
                        record.flags |= DUPLICATE;
                        if arguments.tagging == Tagging::All {
                            record.tags.insert(tag, TagValue::Str(LIBRARY.to_string()));
                        }
                    }
                }
            }
            reads.push(record);
            optical_counts.push(count);
        }
    }
    Ok(Marked {
        reads,
        optical_counts,
    })
}

const DUPLICATE: u16 = 0x400;

/// `GATKDuplicationMetrics`.
#[derive(Debug, Clone, Default)]
struct DuplicationMetrics {
    library: String,
    unpaired_reads_examined: i64,
    read_pairs_examined: i64,
    secondary_or_supplementary: i64,
    unmapped_reads: i64,
    unpaired_read_duplicates: i64,
    read_pair_duplicates: i64,
    read_pair_optical_duplicates: i64,
    percent_duplication: f64,
    estimated_library_size: Option<i64>,
}

impl MetricBean for DuplicationMetrics {
    fn class_name(&self) -> &str {
        "org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        &[
            "LIBRARY",
            "UNPAIRED_READS_EXAMINED",
            "READ_PAIRS_EXAMINED",
            "SECONDARY_OR_SUPPLEMENTARY_RDS",
            "UNMAPPED_READS",
            "UNPAIRED_READ_DUPLICATES",
            "READ_PAIR_DUPLICATES",
            "READ_PAIR_OPTICAL_DUPLICATES",
            "PERCENT_DUPLICATION",
            "ESTIMATED_LIBRARY_SIZE",
        ]
    }

    fn values(&self) -> Vec<Value> {
        vec![
            Value::Str(self.library.clone()),
            Value::Long(self.unpaired_reads_examined),
            Value::Long(self.read_pairs_examined),
            Value::Long(self.secondary_or_supplementary),
            Value::Long(self.unmapped_reads),
            Value::Long(self.unpaired_read_duplicates),
            Value::Long(self.read_pair_duplicates),
            Value::Long(self.read_pair_optical_duplicates),
            Value::Double(self.percent_duplication),
            self.estimated_library_size.map_or(Value::Null, Value::Long),
        ]
    }
}

impl DuplicationMetrics {
    /// `DuplicationMetrics.calculateDerivedFields`, which also takes back the zero `generateMetrics`
    /// put where the size has no estimate.
    fn derive(&mut self) {
        self.estimated_library_size = estimate_library_size(
            self.read_pairs_examined - self.read_pair_optical_duplicates,
            self.read_pairs_examined - self.read_pair_duplicates,
        );
        let examined = self.unpaired_reads_examined + self.read_pairs_examined;
        self.percent_duplication = if examined != 0 {
            (self.unpaired_read_duplicates + self.read_pair_duplicates * 2) as f64
                / (self.unpaired_reads_examined + self.read_pairs_examined * 2) as f64
        } else {
            0.0
        };
    }

    /// `calculateRoiHistogram`, a hundred bins from one times the coverage up.
    fn roi_histogram(&self) -> Option<Histogram> {
        let size = self.estimated_library_size?;
        let unique = self.read_pairs_examined - self.read_pair_duplicates;
        let roi = |x: f64| {
            size as f64
                * (1.0
                    - jmath::strict_math::exp(-(x * self.read_pairs_examined as f64) / size as f64))
                / unique as f64
        };
        Some(Histogram {
            bin_label: "BIN".to_string(),
            value_label: "CoverageMult".to_string(),
            key_class: "java.lang.Double".to_string(),
            bins: (1..=100)
                .map(|step| {
                    let x = f64::from(step);
                    (format!("{x:.1}"), roi(x))
                })
                .collect(),
        })
    }
}

/// `DuplicationMetrics.f`.
fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + jmath::strict_math::exp(-n / x)
}

/// `DuplicationMetrics.estimateLibrarySize`: forty bisections, and nothing where there is nothing
/// to estimate from.
fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    if read_pairs <= 0 || read_pairs - unique_read_pairs <= 0 {
        return None;
    }
    let (unique, pairs) = (unique_read_pairs as f64, read_pairs as f64);
    let (mut low, mut high) = (1.0f64, 100.0f64);
    while f(high * unique, unique, pairs) > 0.0 {
        high *= 10.0;
    }
    for _ in 0..40 {
        let middle = (low + high) / 2.0;
        let u = f(middle * unique, unique, pairs);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            low = middle;
        } else if u < 0.0 {
            high = middle;
        }
    }
    Some((unique * (low + high) / 2.0) as i64)
}

/// `generateMetrics` and `saveMetricsRDD`: a row for each library the header names or a read fell
/// in, by name, and a histogram when exactly one library held reads.
pub fn metrics_file(marked: &Marked, header: &SamHeader, headers: &[String]) -> String {
    let mut rows: BTreeMap<String, DuplicationMetrics> = BTreeMap::new();
    for (record, optical) in marked.reads.iter().zip(&marked.optical_counts) {
        let library = read_group_of(record)
            .and_then(|id| header.read_groups.iter().find(|rg| rg.id == id))
            .and_then(|group| group.attributes.get("LB"))
            .unwrap_or(UNKNOWN_LIBRARY);
        let row = rows
            .entry(library.to_string())
            .or_insert_with(|| DuplicationMetrics {
                library: library.to_string(),
                ..Default::default()
            });
        let unmapped = read::is_unmapped(record);
        let other =
            read::is_secondary_alignment(record) || read::is_supplementary_alignment(record);
        let mate = has_mapped_mate(record);
        if unmapped {
            row.unmapped_reads += 1;
        } else if other {
            row.secondary_or_supplementary += 1;
        } else if !mate {
            row.unpaired_reads_examined += 1;
        } else {
            row.read_pairs_examined += 1;
        }
        if read::is_duplicate(record) && !other && !unmapped {
            if mate {
                row.read_pair_duplicates += 1;
            } else {
                row.unpaired_read_duplicates += 1;
            }
        }
        if *optical > -1 {
            row.read_pair_optical_duplicates += i64::from(*optical);
        }
    }
    let held: BTreeSet<String> = rows.keys().cloned().collect();
    for row in rows.values_mut() {
        row.read_pairs_examined /= 2;
        row.read_pair_duplicates /= 2;
    }
    for library in library_index(header) {
        rows.entry(library.clone())
            .or_insert_with(|| DuplicationMetrics {
                library,
                ..Default::default()
            });
    }
    let mut file = MetricsFile::default();
    for line in headers {
        file.add_header(line);
    }
    let mut histogram = None;
    for row in rows.values_mut() {
        row.derive();
        if held.len() == 1 && held.contains(&row.library) {
            histogram = row.roi_histogram();
        }
    }
    for row in rows.into_values() {
        file.add_metric(&row);
    }
    file.histograms.extend(histogram);
    file.write()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_name_parses_only_with_five_or_seven_fields() {
        let default = ReadNameRegex::Default;
        assert_eq!(location("a:b:1:2:3", &default), Location::of(1, 2, 3));
        assert_eq!(location("a:b:c:d:1:2:3", &default), Location::of(1, 2, 3));
        assert_eq!(location("b:1:2:3", &default), Location::NONE);
        assert_eq!(location("plain", &default), Location::NONE);
    }

    #[test]
    fn the_size_has_no_estimate_without_duplicates() {
        assert_eq!(estimate_library_size(10, 10), None);
        assert_eq!(estimate_library_size(0, 0), None);
        assert!(estimate_library_size(10, 8).is_some());
    }
}
