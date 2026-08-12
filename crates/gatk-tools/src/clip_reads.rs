//! Ported from `org.broadinstitute.hellbender.tools.ClipReads` (GATK 4.6.2.0).
//!
//! The ninth whole tool of the record-transform archetype, and the first that writes a second
//! output that is not a BAM. Everything it does to a read goes through
//! [`gatk_engine::clipping::ReadClipper`], which is ported and measured on its own; what is here is
//! the three clippers that build the ops, the representation that decides how they are applied, and
//! the statistics file, which is Java text formatting rather than htsjdk bytes.
//!
//! # The quality clipper reads the read in machine-cycle order
//!
//! ```java
//! for (int i = readLen - 1; i >= 0; i--) {
//!     int baseIndex = read.isReverseStrand() ? readLen - i - 1 : i;
//!     byte qual = quals[baseIndex];
//!     clipSum += (qTrimmingThreshold - qual);
//!     if (clipSum >= 0 && (clipSum >= lastMax)) { lastMax = clipSum; clipPoint = baseIndex; }
//! }
//! ```
//!
//! The loop counts down, but the index it reads is flipped for a reverse-strand read, so the walk
//! goes **up** the array on one strand and **down** it on the other. The op that comes out is
//! flipped to match: `0..clipPoint` on a reverse read and `clipPoint..readLen-1` on a forward one.
//!
//! Measured, at `-QT 10`, on two reads whose qualities are mirrored in array order: the forward one
//! comes back `ACGTANNNNN` and the reverse one `NNNNNACGTA`. A port that walked the array backwards
//! on both would find no clip point at all on the reverse read, and would emit a file that looks
//! healthier and is wrong.
//!
//! # The representation decides whether the writer sorts
//!
//! `presorted` is true only for `WRITE_NS`, `WRITE_NS_Q0S` and `WRITE_Q0S`. The three
//! representations that can move a read get a sorting writer, so the output order is not the
//! traversal order and the index is the index of the sorted order. Observable twice on the
//! measured fixture: under `SOFTCLIP_BASES` a front clip moves a read from 5 to 10, and under
//! `HARDCLIP_BASES` a reverted soft clip moves another from 6 back to 3, past it.
//!
//! This is the second tool of this archetype whose writer sorts, after
//! [`crate::print_distant_mates`], and it is the first where whether it sorts depends on an
//! argument.
//!
//! # `HARDCLIP_BASES` and `REVERT_SOFTCLIPPED_BASES` revert before they clip
//!
//! In `apply`, before the `ReadClipper` is constructed, so every clipper below sees a read whose
//! cigar and start have already changed. Measured: `3S7M` at 6 leaves as `10M` at 3.
//!
//! # The sequence clipper is a regex, and this port is not
//!
//! The reference compiles each `-X` argument and each `-XF` record with `Pattern.compile(seq,
//! Pattern.CASE_INSENSITIVE)` and matches it against `read.getBasesString()`, looping until
//! `find()` fails. **This port searches for the sequence literally**, ASCII-case-insensitively,
//! restarting after each match, which is what `find()` does for a pattern with no metacharacter in
//! it. A sequence containing a regex metacharacter would diverge, and none of the reference's own
//! callers or tests passes one: the argument is documented as "an exact match algorithm". The bound
//! is written here rather than assumed away, and it is the only place in this file where the port
//! is narrower than the reference.
//!
//! The reverse-strand pattern is built with `BaseUtils.simpleReverseComplement`, which **uppercases
//! whatever it complements** (`a` becomes `T`, not `t`), unlike htsjdk's `SequenceUtil.complement`,
//! which preserves case. That difference never reaches the output here, because the statistics are
//! keyed by the forward string, and it is reproduced anyway.
//!
//! # The adapter clipper does not flip for strand
//!
//! `XF` and `XT` are used as they are, on either strand. `XT` is the first base to clip and `XF` is
//! the first base **not** to clip, both 1-based. Both zero clips the whole read, through a
//! `ClippingOp(0, length)` whose stop is one past the end and which survives only because
//! `overwriteFromStartToStop` takes a `Math.min` against the array length.
//!
//! It is also the one place the tool adds something rather than removing bases: `tf` and `tm` are
//! written onto the read, appending `A` to whatever was there. And it counts `xf` rather than
//! `xf - 1` five-prime bases, which is off by one against the op it just built; that is the
//! reference's arithmetic and it is what the statistics file says.
//!
//! # `--read` drops what it does not name
//!
//! The whole of `apply`, including the write, is inside `read.getName().equals(onlyDoRead)`, so a
//! read that is not named is not passed through: it is gone. A `--read` that names nothing produces
//! an empty BAM and a statistics file whose percentages are `NaN`, which is what
//! `String.format("%.2f", 0.0 / 0)` prints in Java and what
//! [`gatk_engine::java_format::format_decimals`] prints here.

use htsjdk_bam::coordinate;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use gatk_engine::clipping::{ClipError, ClippingOp, ClippingRepresentation, ReadClipper};
use gatk_engine::java_format::format_decimals;
use gatk_engine::read;
use gatk_engine::reads::{ReadsDataSource, ReadsError};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK ClipReads";

/// `ClipReads.FIVE_PRIME_TRIMMING_TAG`.
pub const FIVE_PRIME_TRIMMING_TAG: &[u8; 2] = b"tf";
/// `ClipReads.THREE_PRIME_TRIMMING_TAG`.
pub const THREE_PRIME_TRIMMING_TAG: &[u8; 2] = b"tm";
/// `ClipReads.FIVE_PRIME_ADAPTER_LOCATION_TAG`. The first base **not** to clip, 1-based.
pub const FIVE_PRIME_ADAPTER_LOCATION_TAG: &[u8; 2] = b"XF";
/// `ClipReads.THREE_PRIME_ADAPTER_LOCATION_TAG`. The first base **to** clip, 1-based.
pub const THREE_PRIME_ADAPTER_LOCATION_TAG: &[u8; 2] = b"XT";

/// What the tool refuses on, rather than returning a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ClipReadsError {
    /// `RuntimeException("Badly formatted cyclesToClip argument: " + cyclesToClipArg)`.
    BadlyFormattedCycles(String),
    /// A clip the `ReadClipper` refuses; see [`ClipError`].
    Clip(ClipError),
}

impl From<ClipError> for ClipReadsError {
    fn from(error: ClipError) -> Self {
        ClipReadsError::Clip(error)
    }
}

/// `ClipReads.SeqToClip`: one sequence, and the reverse complement used on a reverse-strand read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqToClip {
    /// `CMDLINE-<n>` for a `-X` argument, or the FASTA record's name for a `-XF` one. It reaches
    /// nothing but the reference's log line, and is kept because the reference keeps it.
    pub name: String,
    /// The argument as typed. This is the statistics file's key, so its case is observable.
    pub seq: String,
    /// `BaseUtils.simpleReverseComplement(seq)`, which uppercases.
    pub rev_seq: String,
}

impl SeqToClip {
    pub fn new(name: &str, bases: &[u8]) -> SeqToClip {
        SeqToClip {
            name: name.to_string(),
            seq: String::from_utf8_lossy(bases).into_owned(),
            rev_seq: String::from_utf8_lossy(&simple_reverse_complement(bases)).into_owned(),
        }
    }
}

/// `BaseUtils.simpleComplement` and `BaseUtils.simpleReverseComplement`, which now live in
/// [`gatk_engine::base_utils`] because `ContextCovariate` reverse-complements a read the same way.
/// Re-exported so this tool reads as it did.
pub use gatk_engine::base_utils::{simple_complement, simple_reverse_complement};

/// The arguments that are this tool's own, beside the [`Options`] every tool in the archetype has.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClipArguments {
    /// `-QT`. Below zero is off, and the default is `-1` rather than `0`: at a threshold of zero
    /// the running sum can still reach zero on a zero-quality base, so `0` is a clipper and `-1`
    /// is not.
    pub q_trimming_threshold: i32,
    /// `-CT`, unparsed, because the parse can fail and the failure is observable.
    pub cycles_to_clip: Option<String>,
    /// `-X`, in the order given.
    pub clip_sequences: Vec<String>,
    /// `-XF`, already read: `(name, bases)` per FASTA record, in file order.
    pub clip_sequence_file: Vec<(String, Vec<u8>)>,
    /// `-CR`, whose default is `WRITE_NS`.
    pub clipping_representation: ClippingRepresentation,
    /// `--read`. Everything else is dropped, not passed through.
    pub only_do_read: Option<String>,
    /// `-CA`.
    pub clip_adapter: bool,
    /// `--min-read-length-to-output`, whose default is 0.
    pub min_read_length: i32,
}

impl Default for ClipArguments {
    fn default() -> Self {
        ClipArguments {
            q_trimming_threshold: -1,
            cycles_to_clip: None,
            clip_sequences: Vec::new(),
            clip_sequence_file: Vec::new(),
            clipping_representation: ClippingRepresentation::WriteNs,
            only_do_read: None,
            clip_adapter: false,
            min_read_length: 0,
        }
    }
}

impl ClipArguments {
    /// `onTraversalStart`: the `-X` arguments first, named `CMDLINE-1` upwards, then every record
    /// of the `-XF` file. The order is the statistics file's insertion order, which the `TreeMap`
    /// then throws away, and is the order the ops are built in, which it does not.
    pub fn sequences_to_clip(&self) -> Vec<SeqToClip> {
        let mut sequences = Vec::new();
        for (index, sequence) in self.clip_sequences.iter().enumerate() {
            sequences.push(SeqToClip::new(
                &format!("CMDLINE-{}", index + 1),
                sequence.as_bytes(),
            ));
        }
        for (name, bases) in &self.clip_sequence_file {
            sequences.push(SeqToClip::new(name, bases));
        }
        sequences
    }

    /// `onTraversalStart`: `start1-end1,start2-end2`, 1-based and inclusive, stored 0-based.
    ///
    /// The reference wraps the whole parse in one `try`, so a negative start, a stop before the
    /// start, a missing dash and a non-number all raise the same message.
    pub fn cycles(&self) -> Result<Option<Vec<(i32, i32)>>, ClipReadsError> {
        let Some(argument) = &self.cycles_to_clip else {
            return Ok(None);
        };
        let bad = || ClipReadsError::BadlyFormattedCycles(argument.clone());
        let mut ranges = Vec::new();
        for range in argument.split(',') {
            let mut parts = range.splitn(2, '-');
            let start: i32 = parts
                .next()
                .ok_or_else(bad)?
                .parse::<i32>()
                .map_err(|_| bad())?
                - 1;
            let stop: i32 = parts
                .next()
                .ok_or_else(bad)?
                .parse::<i32>()
                .map_err(|_| bad())?
                - 1;
            if start < 0 || stop < start {
                return Err(bad());
            }
            ranges.push((start, stop));
        }
        Ok(Some(ranges))
    }
}

/// `ReferenceSequenceFileFactory.getReferenceSequenceFile(path)` then `nextSequence()` until null,
/// for the plain FASTA `-XF` takes.
///
/// The name is truncated at the first whitespace, which is that factory's default, and the bases
/// are the sequence lines concatenated with their whitespace removed and their case left alone: the
/// pattern is compiled case-insensitively, so the case only reaches the reverse-complement string.
pub fn parse_clip_sequence_file(text: &str) -> Vec<(String, Vec<u8>)> {
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            records.push((name, Vec::new()));
        } else if let Some((_, bases)) = records.last_mut() {
            bases.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    records
}

/// `ClipReads.ClippingData`: the counters, and the per-sequence counts in a `TreeMap`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ClippingData {
    pub n_total_reads: i64,
    pub n_total_bases: i64,
    pub n_clipped_reads: i64,
    pub n_clipped_bases: i64,
    pub n_q_clipped_bases: i64,
    pub n_range_clipped_bases: i64,
    pub n_seq_clipped_bases: i64,
    pub n_adapter_clipped_bases: i64,
    /// The `TreeMap<String, Long>`, kept sorted by key: `String.compareTo` on ASCII is byte order,
    /// which puts an upper-case argument before a lower-case one.
    pub seq_clip_counts: Vec<(String, i64)>,
}

impl ClippingData {
    /// `new ClippingData(clipSeqs)`: every sequence starts at zero and is therefore printed even
    /// when it matched nothing.
    pub fn new(sequences: &[SeqToClip]) -> ClippingData {
        let mut data = ClippingData::default();
        for sequence in sequences {
            data.put(&sequence.seq, 0);
        }
        data
    }

    fn put(&mut self, key: &str, value: i64) {
        match self
            .seq_clip_counts
            .binary_search_by(|(existing, _)| existing.as_bytes().cmp(key.as_bytes()))
        {
            Ok(index) => self.seq_clip_counts[index].1 = value,
            Err(index) => self.seq_clip_counts.insert(index, (key.to_string(), value)),
        }
    }

    fn get(&self, key: &str) -> Option<i64> {
        self.seq_clip_counts
            .binary_search_by(|(existing, _)| existing.as_bytes().cmp(key.as_bytes()))
            .ok()
            .map(|index| self.seq_clip_counts[index].1)
    }

    pub fn inc_n_q_clipped_bases(&mut self, n: i32) {
        self.n_q_clipped_bases += n as i64;
        self.n_clipped_bases += n as i64;
    }

    pub fn inc_n_range_clipped_bases(&mut self, n: i32) {
        self.n_range_clipped_bases += n as i64;
        self.n_clipped_bases += n as i64;
    }

    pub fn inc_n_adapter_clipped_bases(&mut self, n: i32) {
        self.n_adapter_clipped_bases += n as i64;
        self.n_clipped_bases += n as i64;
    }

    pub fn inc_seq_clipped_bases(&mut self, seq: &str, n: i32) {
        self.n_seq_clipped_bases += n as i64;
        self.n_clipped_bases += n as i64;
        // `seqClipCounts.get(seq) + n`, which throws if the key is absent. It never is: the map is
        // seeded from the same list the caller iterates.
        let previous = self.get(seq).unwrap_or(0);
        self.put(seq, previous + n as i64);
    }

    /// `addData`: every counter summed, and every key of the other map added to this one.
    pub fn add_data(&mut self, other: &ClippingData) {
        self.n_total_reads += other.n_total_reads;
        self.n_total_bases += other.n_total_bases;
        self.n_clipped_reads += other.n_clipped_reads;
        self.n_clipped_bases += other.n_clipped_bases;
        self.n_q_clipped_bases += other.n_q_clipped_bases;
        self.n_range_clipped_bases += other.n_range_clipped_bases;
        self.n_seq_clipped_bases += other.n_seq_clipped_bases;
        self.n_adapter_clipped_bases += other.n_adapter_clipped_bases;
        for (key, count) in &other.seq_clip_counts {
            let total = count + self.get(key).unwrap_or(0);
            self.put(key, total);
        }
    }

    /// `toString`: the statistics file, byte for byte.
    ///
    /// Two eighty-dash rules, nine counters with the label padded out to a fixed width by hand in
    /// the format string, two `%.2f` percentages that are `NaN` when nothing was examined, one
    /// `%8d` row per clip sequence in `TreeMap` order, and an adapter row present only under `-CA`.
    pub fn to_text(&self, clip_adapter: bool) -> String {
        let rule = "-".repeat(80);
        let percent = |numerator: i64, denominator: i64| {
            format_decimals(100.0 * numerator as f64 / denominator as f64, 2)
        };
        let mut text = String::new();
        text.push_str(&rule);
        text.push('\n');
        text.push_str(&format!(
            "Number of examined reads              {}\n",
            self.n_total_reads
        ));
        text.push_str(&format!(
            "Number of clipped reads               {}\n",
            self.n_clipped_reads
        ));
        text.push_str(&format!(
            "Percent of clipped reads              {}\n",
            percent(self.n_clipped_reads, self.n_total_reads)
        ));
        text.push_str(&format!(
            "Number of examined bases              {}\n",
            self.n_total_bases
        ));
        text.push_str(&format!(
            "Number of clipped bases               {}\n",
            self.n_clipped_bases
        ));
        text.push_str(&format!(
            "Percent of clipped bases              {}\n",
            percent(self.n_clipped_bases, self.n_total_bases)
        ));
        text.push_str(&format!(
            "Number of quality-score clipped bases {}\n",
            self.n_q_clipped_bases
        ));
        text.push_str(&format!(
            "Number of range clipped bases         {}\n",
            self.n_range_clipped_bases
        ));
        text.push_str(&format!(
            "Number of sequence clipped bases      {}\n",
            self.n_seq_clipped_bases
        ));
        for (key, count) in &self.seq_clip_counts {
            text.push_str(&format!("  {count:8} clip sites matching {key}\n"));
        }
        if clip_adapter {
            text.push_str(&format!(
                "Number of adapter clipped bases       {}\n",
                self.n_adapter_clipped_bases
            ));
        }
        text.push_str(&rule);
        text.push('\n');
        text
    }
}

/// `strandAwarePositions`: a forward-strand span, turned into this read's span.
fn strand_aware_positions(read: &BamRecord, start: i32, stop: i32) -> (i32, i32) {
    if read::is_reverse_strand(read) {
        let length = read.read_bases.len() as i32;
        (length - stop - 1, length - start - 1)
    } else {
        (start, stop)
    }
}

/// `clipBadQualityScores`: BWA's trimming scan, walked in machine-cycle order.
fn clip_bad_quality_scores(
    read: &BamRecord,
    threshold: i32,
    data: &mut ClippingData,
) -> Option<ClippingOp> {
    let read_len = read.read_bases.len() as i32;
    let reverse = read::is_reverse_strand(read);
    let (mut clip_sum, mut last_max, mut clip_point) = (0i32, -1i32, -1i32);
    for i in (0..read_len).rev() {
        let base_index = if reverse { read_len - i - 1 } else { i };
        let qual = read.base_qualities[base_index as usize] as i32;
        clip_sum += threshold - qual;
        // `>= lastMax`, not `>`: a plateau moves the clip point on, so the longest tail at the
        // maximum wins rather than the first one to reach it.
        if clip_sum >= 0 && clip_sum >= last_max {
            last_max = clip_sum;
            clip_point = base_index;
        }
    }
    if clip_point == -1 {
        return None;
    }
    let (start, stop) = if reverse {
        (0, clip_point)
    } else {
        (clip_point, read_len - 1)
    };
    let op = ClippingOp { start, stop };
    data.inc_n_q_clipped_bases(op.stop - op.start + 1);
    Some(op)
}

/// `clipCycles`: 1-based inclusive ranges, clamped at the end and dropped past it.
fn clip_cycles(
    read: &BamRecord,
    cycles: &[(i32, i32)],
    data: &mut ClippingData,
) -> Vec<ClippingOp> {
    let read_len = read.read_bases.len() as i32;
    let mut ops = Vec::new();
    for (cycle_start, cycle_stop) in cycles {
        if *cycle_start >= read_len {
            continue;
        }
        // "we do tolerate [for convenience) clipping when the stop is beyond the end of the read"
        let cycle_stop = if *cycle_stop >= read_len {
            read_len - 1
        } else {
            *cycle_stop
        };
        let (start, stop) = strand_aware_positions(read, *cycle_start, cycle_stop);
        let op = ClippingOp { start, stop };
        data.inc_n_range_clipped_bases(op.stop - op.start + 1);
        ops.push(op);
    }
    ops
}

/// `clipSequences`: every match of every sequence, in the order the sequences were given.
///
/// See the module docs for what this port does instead of `java.util.regex`.
fn clip_sequences(
    read: &BamRecord,
    sequences: &[SeqToClip],
    data: &mut ClippingData,
) -> Vec<ClippingOp> {
    let mut ops = Vec::new();
    if sequences.is_empty() {
        return ops;
    }
    let reverse = read::is_reverse_strand(read);
    for sequence in sequences {
        let pattern = if reverse {
            sequence.rev_seq.as_bytes()
        } else {
            sequence.seq.as_bytes()
        };
        for (start, stop) in find_all(&read.read_bases, pattern) {
            let op = ClippingOp { start, stop };
            data.inc_seq_clipped_bases(&sequence.seq, op.stop - op.start + 1);
            ops.push(op);
        }
    }
    ops
}

/// `Matcher.find()` in a loop: every non-overlapping match, ASCII-case-insensitive, as inclusive
/// `(start, stop)` offsets.
fn find_all(haystack: &[u8], needle: &[u8]) -> Vec<(i32, i32)> {
    let mut found = Vec::new();
    if needle.is_empty() || needle.len() > haystack.len() {
        return found;
    }
    let mut at = 0;
    while at + needle.len() <= haystack.len() {
        let matched = haystack[at..at + needle.len()]
            .iter()
            .zip(needle)
            .all(|(a, b)| a.eq_ignore_ascii_case(b));
        if matched {
            found.push((at as i32, (at + needle.len() - 1) as i32));
            at += needle.len();
        } else {
            at += 1;
        }
    }
    found
}

/// `read.getAttributeAsInteger(tag)` for the two adapter tags.
fn attribute_as_integer(read: &BamRecord, tag: &[u8; 2]) -> Option<i32> {
    match read.tags.get(Tag::new(tag)) {
        Some(TagValue::Int(value)) => Some(*value as i32),
        _ => None,
    }
}

/// `addAdapterTag`: `A` if the tag is absent, else `A` appended unless one is already there.
fn add_adapter_tag(read: &mut BamRecord, tag: &[u8; 2]) {
    let current = match read.tags.get(Tag::new(tag)) {
        Some(TagValue::Str(value)) => Some(value.clone()),
        _ => None,
    };
    let next = match current {
        None => "A".to_string(),
        Some(value) if !value.contains('A') => format!("{value}A"),
        Some(value) => value,
    };
    read.tags.insert(Tag::new(tag), TagValue::Str(next));
}

/// `clipAdapter`: the ops, and the tags written onto the read.
fn clip_adapter(read: &mut BamRecord, data: &mut ClippingData) -> Vec<ClippingOp> {
    let mut ops = Vec::new();
    let length = read.read_bases.len() as i32;
    let xf = attribute_as_integer(read, FIVE_PRIME_ADAPTER_LOCATION_TAG);
    let xt = attribute_as_integer(read, THREE_PRIME_ADAPTER_LOCATION_TAG);

    if xf == Some(0) && xt == Some(0) {
        // Stop is `read.getLength()`, one past the last offset. `clipRead` shortens it back and
        // `overwriteFromStartToStop` takes a `Math.min`, so the reference never indexes past the
        // array; the op is left as the reference builds it.
        ops.push(ClippingOp {
            start: 0,
            stop: length,
        });
        data.inc_n_adapter_clipped_bases(length);
        return ops;
    }
    if let Some(xt) = xt.filter(|xt| *xt <= length) {
        ops.push(ClippingOp {
            start: xt - 1,
            stop: length,
        });
        add_adapter_tag(read, THREE_PRIME_TRIMMING_TAG);
        data.inc_n_adapter_clipped_bases(length - xt + 1);
    }
    if let Some(xf) = xf.filter(|xf| *xf > 1) {
        // Stop is included, so `xf - 2` is the last base before the first one to keep.
        ops.push(ClippingOp {
            start: 0,
            stop: xf - 2,
        });
        add_adapter_tag(read, FIVE_PRIME_TRIMMING_TAG);
        // `xf`, not `xf - 1`: the reference counts one more base than the op it just built covers.
        data.inc_n_adapter_clipped_bases(xf);
    }
    ops
}

/// What one read's pass through `apply` produced.
pub struct Clipped {
    /// The read as it will be written, after the representation was applied.
    pub read: BamRecord,
    /// This read's own counters, which `accumulate` folds into the run's.
    pub data: ClippingData,
    /// `ReadClipper.wasClipped()`: whether any op was built at all.
    pub was_clipped: bool,
    /// `clipper.getRead().getLength()`, which is the reverted length under the two representations
    /// that revert, and is what the examined-bases counter is made of.
    pub examined_bases: i32,
}

/// `apply` for one read, minus the write and the accumulation.
///
/// The reference builds the `ReadClipper` first and mutates its read through `getRead()` when the
/// adapter clipper writes a tag. Here the ops are built against the read, the tags are written to
/// the copy, and the clipper is constructed from that copy last. The two are the same because no op
/// depends on a tag and no representation reads one: the tags decide the ops, never the reverse.
pub fn clip_one(
    read: &BamRecord,
    header: Option<&SamHeader>,
    sequences: &[SeqToClip],
    cycles: Option<&Vec<(i32, i32)>>,
    arguments: &ClipArguments,
) -> Result<Clipped, ClipReadsError> {
    let mut working = read.clone();
    if matches!(
        arguments.clipping_representation,
        ClippingRepresentation::HardclipBases | ClippingRepresentation::RevertSoftclippedBases
    ) {
        working = gatk_engine::clipping::revert_soft_clipped_bases(&working, header)?;
    }

    let mut data = ClippingData::new(sequences);
    let mut ops = Vec::new();
    // The four run in this order, and `clipRead` applies the ops in the order they were added
    // against a read each previous op has already shortened.
    if let Some(op) = clip_bad_quality_scores(&working, arguments.q_trimming_threshold, &mut data) {
        ops.push(op);
    }
    if let Some(cycles) = cycles {
        ops.extend(clip_cycles(&working, cycles, &mut data));
    }
    ops.extend(clip_sequences(&working, sequences, &mut data));
    if arguments.clip_adapter {
        ops.extend(clip_adapter(&mut working, &mut data));
    }

    let examined_bases = working.read_bases.len() as i32;
    let was_clipped = !ops.is_empty();
    let mut clipper = ReadClipper::new(&working, header);
    for op in ops {
        clipper.add_op(op);
    }
    let clipped = clipper.clip_read(arguments.clipping_representation)?;
    Ok(Clipped {
        read: clipped,
        data,
        was_clipped,
        examined_bases,
    })
}

/// `createSAMWriter(OUTPUT, presorted)`: only the three representations that cannot move a read
/// are written presorted.
pub fn presorted(representation: ClippingRepresentation) -> bool {
    matches!(
        representation,
        ClippingRepresentation::WriteNs
            | ClippingRepresentation::WriteNsQ0s
            | ClippingRepresentation::WriteQ0s
    )
}

/// What a run produces: the output BAM, its index, and the statistics file's text.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>, String), ClipReadsError>, ReadsError>;

/// `ClipReads`: every read the traversal reaches, clipped, and the statistics beside them.
pub fn clip_reads(
    source: &ReadsDataSource,
    options: &Options,
    arguments: &ClipArguments,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> RunResult {
    let cycles = match arguments.cycles() {
        Ok(cycles) => cycles,
        Err(error) => return Ok(Err(error)),
    };
    let sequences = arguments.sequences_to_clip();
    let source_header = source.header();

    let traversed = crate::read_walker::traverse(source, &options.intervals, filter)?;

    let mut accumulator = ClippingData::new(&sequences);
    let mut records = Vec::with_capacity(traversed.len());
    for record in &traversed {
        if let Some(name) = &arguments.only_do_read {
            if &record.read_name != name {
                continue;
            }
        }
        let clipped = match clip_one(
            record,
            Some(source_header),
            &sequences,
            cycles.as_ref(),
            arguments,
        ) {
            Ok(clipped) => clipped,
            Err(error) => return Ok(Err(error)),
        };
        if clipped.read.read_bases.len() as i32 >= arguments.min_read_length {
            records.push(clipped.read);
        }
        accumulator.n_total_reads += 1;
        accumulator.n_total_bases += clipped.examined_bases as i64;
        if clipped.was_clipped {
            accumulator.n_clipped_reads += 1;
            accumulator.add_data(&clipped.data);
        }
    }

    if !presorted(arguments.clipping_representation) {
        // htsjdk's sorting collection, which is stable: two records that compare equal keep the
        // order they were written in.
        records.sort_by(coordinate::compare);
    }

    let header = header_for_sam_writer(source_header, TOOL_NAME, options);
    let (bam, bai) = write_records(&header, &records, options.create_output_bam_index)?;
    Ok(Ok((bam, bai, accumulator.to_text(arguments.clip_adapter))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::text_parse::parse_cigar;

    fn read(
        name: &str,
        start: i32,
        cigar: &str,
        bases: &str,
        quals: &str,
        flags: u16,
    ) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            reference_index: 0,
            alignment_start: start,
            mapping_quality: 60,
            cigar: parse_cigar(cigar).unwrap(),
            read_bases: bases.as_bytes().to_vec(),
            base_qualities: quals.bytes().map(|b| b - 33).collect(),
            ..BamRecord::default()
        }
    }

    fn clip(record: &BamRecord, arguments: &ClipArguments) -> Clipped {
        let sequences = arguments.sequences_to_clip();
        let cycles = arguments.cycles().unwrap();
        clip_one(record, None, &sequences, cycles.as_ref(), arguments).unwrap()
    }

    #[test]
    fn the_quality_scan_clips_the_tail_of_a_forward_read() {
        let record = read("r0", 1, "10M", "ACGTAGGTAC", "IIIII#####", 0);
        let clipped = clip(
            &record,
            &ClipArguments {
                q_trimming_threshold: 10,
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.read_bases, b"ACGTANNNNN");
        assert_eq!(clipped.data.n_q_clipped_bases, 5);
    }

    #[test]
    fn and_the_front_of_a_reverse_one_with_the_same_qualities_mirrored() {
        let record = read("r1", 5, "10M", "GGGGGACGTA", "#####IIIII", 16);
        let clipped = clip(
            &record,
            &ClipArguments {
                q_trimming_threshold: 10,
                ..ClipArguments::default()
            },
        );
        assert_eq!(
            clipped.read.read_bases, b"NNNNNACGTA",
            "the scan walks the array forwards on a reverse-strand read"
        );
    }

    #[test]
    fn a_cycle_range_past_the_end_of_the_read_never_starts() {
        let record = read("r6", 35, "5M", "ACGTA", "IIIII", 0);
        let clipped = clip(
            &record,
            &ClipArguments {
                cycles_to_clip: Some("1-3,8-12".to_string()),
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.read_bases, b"NNNTA");
        assert_eq!(clipped.data.n_range_clipped_bases, 3);
    }

    #[test]
    fn a_sequence_clips_a_read_as_often_as_it_matches() {
        let record = read("r5", 25, "10M", "GGGGGGGGGG", "IIIIIIIIII", 0);
        let clipped = clip(
            &record,
            &ClipArguments {
                clip_sequences: vec!["GGGGG".to_string()],
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.read_bases, b"NNNNNNNNNN");
        assert_eq!(clipped.data.get("GGGGG"), Some(10));
    }

    #[test]
    fn a_reverse_read_is_matched_against_the_reverse_complement() {
        // `ACGT` is its own reverse complement, so it matches twice on a reverse read; `GGGGG`
        // becomes `CCCCC` and does not match at all.
        let record = read("r3", 15, "10M", "ACGTACGTAC", "IIIIIIIIII", 16);
        let clipped = clip(
            &record,
            &ClipArguments {
                clip_sequences: vec!["GGGGG".to_string(), "acgt".to_string()],
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.read_bases, b"NNNNNNNNAC");
        assert_eq!(clipped.data.get("GGGGG"), Some(0));
        assert_eq!(clipped.data.get("acgt"), Some(8), "case-insensitively");
    }

    #[test]
    fn the_per_sequence_counts_are_in_ascii_order_of_the_argument() {
        let arguments = ClipArguments {
            clip_sequences: vec!["acgt".to_string(), "GGGGG".to_string()],
            ..ClipArguments::default()
        };
        let data = ClippingData::new(&arguments.sequences_to_clip());
        let keys: Vec<&str> = data
            .seq_clip_counts
            .iter()
            .map(|(k, _)| k.as_str())
            .collect();
        assert_eq!(keys, vec!["GGGGG", "acgt"], "upper case first");
    }

    #[test]
    fn the_adapter_tags_are_not_flipped_for_strand() {
        let mut record = read("r3", 15, "10M", "ACGTACGTAC", "IIIIIIIIII", 16);
        record
            .tags
            .insert(Tag::new(FIVE_PRIME_ADAPTER_LOCATION_TAG), TagValue::Int(3));
        record
            .tags
            .insert(Tag::new(THREE_PRIME_ADAPTER_LOCATION_TAG), TagValue::Int(8));
        let clipped = clip(
            &record,
            &ClipArguments {
                clip_adapter: true,
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.read_bases, b"NNGTACGNNN");
        assert_eq!(
            clipped.read.tags.get(Tag::new(FIVE_PRIME_TRIMMING_TAG)),
            Some(&TagValue::Str("A".to_string()))
        );
        assert_eq!(
            clipped.read.tags.get(Tag::new(THREE_PRIME_TRIMMING_TAG)),
            Some(&TagValue::Str("A".to_string()))
        );
        // Three bases from the 3' op, and `xf` rather than `xf - 1` from the 5' one.
        assert_eq!(clipped.data.n_adapter_clipped_bases, 3 + 3);
    }

    #[test]
    fn both_adapter_tags_zero_clips_the_whole_read() {
        let mut record = read("r4", 21, "10M", "TTTTTTTTTT", "IIIIIIIIII", 0);
        record
            .tags
            .insert(Tag::new(FIVE_PRIME_ADAPTER_LOCATION_TAG), TagValue::Int(0));
        record
            .tags
            .insert(Tag::new(THREE_PRIME_ADAPTER_LOCATION_TAG), TagValue::Int(0));
        let clipped = clip(
            &record,
            &ClipArguments {
                clip_adapter: true,
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.read_bases, b"NNNNNNNNNN");
        assert_eq!(clipped.data.n_adapter_clipped_bases, 10);
    }

    #[test]
    fn reverting_a_soft_clip_happens_before_anything_is_clipped() {
        let record = read("r2", 6, "3S7M", "TTTGGTACCA", "IIIIIIIIII", 0);
        let clipped = clip(
            &record,
            &ClipArguments {
                q_trimming_threshold: 10,
                clipping_representation: ClippingRepresentation::HardclipBases,
                ..ClipArguments::default()
            },
        );
        assert_eq!(clipped.read.cigar.to_text(), "10M");
        assert_eq!(clipped.read.alignment_start, 3);
    }

    #[test]
    fn only_three_representations_are_written_presorted() {
        assert!(presorted(ClippingRepresentation::WriteNs));
        assert!(presorted(ClippingRepresentation::WriteQ0s));
        assert!(presorted(ClippingRepresentation::WriteNsQ0s));
        assert!(!presorted(ClippingRepresentation::SoftclipBases));
        assert!(!presorted(ClippingRepresentation::HardclipBases));
        assert!(!presorted(ClippingRepresentation::RevertSoftclippedBases));
    }

    #[test]
    fn a_badly_formatted_cycle_argument_is_refused() {
        for argument in ["1", "0-3", "5-2", "a-b", "1-"] {
            let arguments = ClipArguments {
                cycles_to_clip: Some(argument.to_string()),
                ..ClipArguments::default()
            };
            assert_eq!(
                arguments.cycles(),
                Err(ClipReadsError::BadlyFormattedCycles(argument.to_string())),
                "{argument}"
            );
        }
    }

    #[test]
    fn nothing_examined_prints_nan_rather_than_zero() {
        let text = ClippingData::default().to_text(false);
        assert!(
            text.contains("Percent of clipped reads              NaN\n"),
            "{text}"
        );
        assert!(
            text.contains("Percent of clipped bases              NaN\n"),
            "{text}"
        );
        assert!(
            !text.contains("adapter"),
            "the adapter row is present only under -CA"
        );
    }

    #[test]
    fn the_reverse_complement_a_pattern_is_built_from_uppercases() {
        assert_eq!(simple_reverse_complement(b"acgt"), b"ACGT");
        assert_eq!(simple_reverse_complement(b"GGTACC"), b"GGTACC");
        assert_eq!(
            simple_reverse_complement(b"NxN"),
            b"NxN",
            "and passes the rest through"
        );
    }
}
