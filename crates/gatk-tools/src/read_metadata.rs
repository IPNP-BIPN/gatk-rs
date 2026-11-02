//! `ReadMetadata`, `LibraryStatistics` and `IntHistogram.CDF` (GATK 4.6.2.0): what
//! `CalcMetadataSpark` computes, the first step of the SV pipeline.
//!
//! # One partition
//!
//! The file names partitions, and the reference has as many as Disq cut the BAM into. This port runs
//! what a single partition observes, which is X.1's claim everywhere else, so [`ReadMetadata`] holds
//! one [`PartitionBounds`] and the caller refuses a BAM that Disq would have cut.
//!
//! # Every read is reached through the unpaired trap
//!
//! `isTemplateLenTestable` asks about the mate of any primary read that passes the evidence bar, and
//! `mateIsUnmapped` throws on an unpaired read rather than answering. A single-ended BAM with a good
//! alignment therefore stops the tool, unless the bar is what stops the read first: mapping quality
//! under 20 or under 45 aligned bases.
//!
//! # The binary form is what a buffer happened to flush
//!
//! ```java
//! try ( final OutputStream os = new FileOutputStream(outputFile) ) {
//!     new ReadMetadata.Serializer().write(new Kryo(), new Output(os), readMetadata);
//! ```
//!
//! Nothing flushes that `Output`, so `--write-as-binary` leaves the file holding only the bytes that
//! overflowed Kryo's 4096-byte buffer, and the stream's tail is lost. [`FlushedOutput`] follows the
//! buffer, since which bytes survive depends on the size of each write.
//!
//! # Fractions of nothing
//!
//! A library with no testable template divides zero by zero, so its whole CDF is NaN. A binary search
//! over NaNs goes by `floatToIntBits` and lands at index zero, so the median is 0 and both deviations
//! are 0, and the histogram line prints a run of zeroes because `Math.round(NaN)` is 0.

use gatk_engine::java_hash::{HashOrderError, JavaHashCode, JavaHashMap};
use gatk_engine::kryo::Output;
use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// `maxTrackedFragmentLength`, which `CalcMetadataSpark` fixes.
pub const MAX_TRACKED_FRAGMENT_LENGTH: usize = 10000;

/// `FindBreakpointEvidenceSparkArgumentCollection`'s defaults: `CalcMetadataSpark` builds one and
/// sets nothing.
const MIN_EVIDENCE_MAPQ: u8 = 20;
const MIN_EVIDENCE_MATCH_LENGTH: i64 = 45;
const ALLOWED_SHORT_FRAGMENT_OVERHANG: i32 = 10;

/// `PartitionBounds.UNMAPPED`.
const UNMAPPED: i32 = i32::MAX;

const MIN_COVERAGE: f32 = 10.0;

/// What stops the constructor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MetadataFailure {
    /// `mateIsUnmapped` on an unpaired read, inside the task.
    Unpaired,
    /// `nReadBases / nReads` with no mapped primary read.
    NoReads,
    Order(HashOrderError),
}

impl From<HashOrderError> for MetadataFailure {
    fn from(error: HashOrderError) -> MetadataFailure {
        MetadataFailure::Order(error)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PartitionBounds {
    pub first_contig: i32,
    pub first_start: i32,
    pub last_contig: i32,
    pub last_end: i32,
    pub span: i32,
}

/// `IntHistogram.CDF`.
#[derive(Debug, Clone)]
pub struct Cdf {
    pub fractions: Vec<f32>,
    pub counts: i64,
}

impl Cdf {
    fn new(histogram: &[i64], total: i64) -> Cdf {
        let mut sum = 0i64;
        let fractions = histogram
            .iter()
            .map(|count| {
                sum += count;
                sum as f32 / total as f32
            })
            .collect();
        Cdf {
            fractions,
            counts: total,
        }
    }

    /// `internalPopStat`.
    fn pop_stat(&self, fraction: f32) -> i32 {
        match binary_search_floats(&self.fractions, fraction) {
            index if index < 0 => -index - 1,
            index => index,
        }
    }
}

/// `Arrays.binarySearch(float[], float)`, whose ties go by `floatToIntBits`.
fn binary_search_floats(values: &[f32], key: f32) -> i32 {
    let mut low = 0i32;
    let mut high = values.len() as i32 - 1;
    while low <= high {
        let middle = (low + high) >> 1;
        let value = values[middle as usize];
        if value < key {
            low = middle + 1;
        } else if value > key {
            high = middle - 1;
        } else {
            let bits = float_to_int_bits(value);
            let key_bits = float_to_int_bits(key);
            if bits == key_bits {
                return middle;
            } else if bits < key_bits {
                low = middle + 1;
            } else {
                high = middle - 1;
            }
        }
    }
    -(low + 1)
}

/// `Float.floatToIntBits`: one NaN.
fn float_to_int_bits(value: f32) -> i32 {
    if value.is_nan() {
        0x7fc0_0000
    } else {
        value.to_bits() as i32
    }
}

/// `LibraryStatistics`.
#[derive(Debug, Clone)]
pub struct LibraryStatistics {
    pub cdf: Cdf,
    pub median: i32,
    pub negative_mad: f32,
    pub positive_mad: f32,
    pub coverage: f32,
    pub mean_base_quality: f32,
    pub reads: i64,
    pub read_start_frequency: f32,
}

/// A library's name, which may be null twice over: no read group, or a group with no `LB`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LibraryName(pub Option<String>);

impl JavaHashCode for LibraryName {
    fn java_hash_code(&self) -> i32 {
        self.0.as_ref().map_or(0, JavaHashCode::java_hash_code)
    }
}

struct RawStatistics {
    fragment_sizes: Vec<i64>,
    observations: i64,
    reads: i64,
    bases: i64,
    total_base_quality: i64,
}

#[derive(Debug, Clone)]
pub struct ReadMetadata {
    pub reads: i64,
    pub avg_read_len: i32,
    pub ref_bases: i64,
    pub max_reads_in_partition: i64,
    pub coverage: f32,
    pub mean_base_quality: f32,
    pub partition: PartitionBounds,
    pub contigs: Vec<String>,
    pub contig_ids: JavaHashMap<String, i32>,
    pub read_group_to_library: JavaHashMap<String, Option<String>>,
    pub libraries: JavaHashMap<LibraryName, LibraryStatistics>,
}

/// `SVUtils.hashMapCapacity`.
fn hash_map_capacity(elements: usize) -> usize {
    (elements as u64 * 4 / 3) as usize + 1
}

/// `CigarUtils.countAlignedBases`: `CigarOperator.isAlignment`.
fn aligned_bases(record: &BamRecord) -> i64 {
    record
        .cigar
        .elements
        .iter()
        .filter(|element| matches!(element.op, Op::M | Op::Eq | Op::X))
        .map(|element| i64::from(element.length))
        .sum()
}

const PAIRED: u16 = 0x1;
const UNMAPPED_FLAG: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// `SVReadFilter.isMappedPrimary`.
fn is_mapped_primary(record: &BamRecord) -> bool {
    record.flags & (DUPLICATE | QC_FAIL | UNMAPPED_FLAG | SECONDARY | SUPPLEMENTARY) == 0
        && record.reference_index >= 0
        && record.alignment_start != 0
}

/// `SVReadFilter.isTemplateLenTestable`, given a mapped primary read.
fn is_template_length_testable(record: &BamRecord) -> Result<bool, MetadataFailure> {
    if record.mapping_quality < MIN_EVIDENCE_MAPQ
        || aligned_bases(record) < MIN_EVIDENCE_MATCH_LENGTH
    {
        return Ok(false);
    }
    if record.flags & PAIRED == 0 {
        return Err(MetadataFailure::Unpaired);
    }
    let mate_unmapped = record.flags & MATE_UNMAPPED != 0
        || record.mate_reference_index < 0
        || record.mate_alignment_start == 0;
    Ok(!mate_unmapped
        && record.flags & REVERSE == 0
        && record.flags & MATE_REVERSE != 0
        && record.reference_index == record.mate_reference_index
        && record.alignment_start - ALLOWED_SHORT_FRAGMENT_OVERHANG <= record.mate_alignment_start)
}

impl ReadMetadata {
    /// The constructor, over one partition's unfiltered reads.
    pub fn new(header: &SamHeader, records: &[BamRecord]) -> Result<ReadMetadata, MetadataFailure> {
        let mut contig_ids = JavaHashMap::with_capacity(hash_map_capacity(header.sequences.len()));
        for (id, sequence) in header.sequences.iter().enumerate() {
            contig_ids.insert(sequence.name.clone(), id as i32);
        }
        let mut contigs = vec![String::new(); contig_ids.len()];
        for (name, id) in contig_ids.iter() {
            contigs[*id as usize] = name.clone();
        }
        let mut read_group_to_library =
            JavaHashMap::with_capacity(hash_map_capacity(header.read_groups.len()));
        for group in &header.read_groups {
            read_group_to_library.insert(
                group.id.clone(),
                group.attributes.get("LB").map(str::to_string),
            );
        }

        // `PartitionStatistics`: its map is a plain `new HashMap<>()`.
        let mut raw: JavaHashMap<LibraryName, RawStatistics> = JavaHashMap::new();
        let mapped: Vec<&BamRecord> = records
            .iter()
            .filter(|record| is_mapped_primary(record))
            .collect();
        let mut partition = PartitionBounds {
            first_contig: UNMAPPED,
            first_start: -1,
            last_contig: UNMAPPED,
            last_end: -1,
            span: 0,
        };
        if let (Some(first), Some(last)) = (mapped.first(), mapped.last()) {
            let mut span = -first.alignment_start;
            let mut contig = first.reference_index;
            let mut previous_end = 0;
            for (position, record) in mapped.iter().enumerate() {
                if position > 0 && record.reference_index != contig {
                    span += previous_end;
                    contig = record.reference_index;
                    span -= record.alignment_start;
                }
                let group = match record.tags.get(Tag::new(b"RG")) {
                    Some(TagValue::Str(id)) => Some(id.to_string()),
                    _ => None,
                };
                let library = LibraryName(
                    group
                        .and_then(|id| read_group_to_library.get(&id).cloned())
                        .flatten(),
                );
                let testable = is_template_length_testable(record)?;
                if !raw.contains_key(&library) {
                    raw.insert(
                        library.clone(),
                        RawStatistics {
                            fragment_sizes: vec![0; MAX_TRACKED_FRAGMENT_LENGTH + 2],
                            observations: 0,
                            reads: 0,
                            bases: 0,
                            total_base_quality: 0,
                        },
                    );
                }
                let stats = raw.get_mut(&library).expect("just inserted");
                stats.reads += 1;
                stats.bases += aligned_bases(record);
                stats.total_base_quality += record
                    .base_qualities
                    .iter()
                    .map(|quality| i64::from(*quality as i8))
                    .sum::<i64>();
                if testable {
                    let size = record.inferred_insert_size.unsigned_abs() as usize;
                    let last_bin = stats.fragment_sizes.len() - 1;
                    stats.fragment_sizes[size.min(last_bin)] += 1;
                    stats.observations += 1;
                }
                previous_end = record.alignment_end() + 1;
            }
            let last_end = last.alignment_end() + 1;
            partition = PartitionBounds {
                first_contig: first.reference_index,
                first_start: first.alignment_start,
                last_contig: last.reference_index,
                last_end,
                span: span + last_end,
            };
        }
        raw.check()?;

        let reads: i64 = raw.iter().map(|(_, stats)| stats.reads).sum();
        let bases: i64 = raw.iter().map(|(_, stats)| stats.bases).sum();
        if reads == 0 {
            return Err(MetadataFailure::NoReads);
        }
        let ref_bases: i64 = header.sequences.iter().map(|s| i64::from(s.length)).sum();
        let coverage = bases as f32 / ref_bases as f32;
        let total_quality: i64 = raw.iter().map(|(_, stats)| stats.total_base_quality).sum();

        let mut libraries = JavaHashMap::with_capacity(hash_map_capacity(raw.len()));
        for (name, stats) in raw.iter() {
            let cdf = Cdf::new(&stats.fragment_sizes, stats.observations);
            let median = cdf.pop_stat(0.5);
            let negative_mad = (median - cdf.pop_stat(0.25)) as f32;
            let positive_mad = (cdf.pop_stat(0.75) - median) as f32;
            libraries.insert(
                name.clone(),
                LibraryStatistics {
                    cdf,
                    median,
                    negative_mad,
                    positive_mad,
                    coverage: stats.bases as f32 / ref_bases as f32,
                    mean_base_quality: stats.total_base_quality as f32 / stats.bases as f32,
                    reads: stats.reads,
                    read_start_frequency: stats.reads as f32 / ref_bases as f32,
                },
            );
        }
        contig_ids.check()?;
        read_group_to_library.check()?;
        libraries.check()?;

        Ok(ReadMetadata {
            reads,
            avg_read_len: (bases / reads) as i32,
            ref_bases,
            max_reads_in_partition: reads,
            coverage: if coverage >= MIN_COVERAGE {
                coverage
            } else {
                MIN_COVERAGE
            },
            mean_base_quality: total_quality as f32 / bases as f32,
            partition,
            contigs,
            contig_ids,
            read_group_to_library,
            libraries,
        })
    }

    /// `ReadMetadata.writeMetadata`.
    pub fn text(&self) -> String {
        use jmath::to_string::float_to_string as float;
        let mut text = format!(
            "#reads:\t{}\n#partitions:\t1\nmax reads/partition:\t{}\ncoverage:\t{}\nmeanQ:\t{}\n",
            self.reads,
            self.max_reads_in_partition,
            float(self.coverage),
            float(self.mean_base_quality)
        );
        text.push_str("\nLibrary Statistics\n");
        for (name, stats) in self.libraries.iter() {
            text.push_str(&format!(
                "{}:\t{}-{}+{}\t{}\t{}\t{}\t{}\n",
                name.0.as_deref().unwrap_or("NoGroup"),
                stats.median,
                float(stats.negative_mad),
                float(stats.positive_mad),
                float(stats.coverage),
                float(stats.mean_base_quality),
                stats.reads,
                float(stats.read_start_frequency)
            ));
            text.push_str("template size cumulative counts:");
            for fraction in &stats.cdf.fractions {
                // A float times a long is a float, and `Math.round(float)` is an int.
                let cumulative = fraction * stats.cdf.counts as f32;
                text.push_str(&format!(
                    "\t{}",
                    (f64::from(cumulative) + 0.5).floor() as i32
                ));
            }
            text.push('\n');
        }
        text.push_str("\nPartition Boundaries\n");
        let bounds = &self.partition;
        text.push_str(&format!(
            "0\t{}\t{}\t{}\t{}\n",
            bounds.first_contig, bounds.first_start, bounds.last_contig, bounds.last_end
        ));
        text.push_str("contigs map:\n");
        for (id, name) in self.contigs.iter().enumerate() {
            text.push_str(&format!("{id}:{name}\n"));
        }
        text
    }

    /// `Serializer.write` through an `Output` nothing flushes.
    pub fn flushed_binary(&self) -> Vec<u8> {
        let mut output = FlushedOutput::default();
        output.write_int(0); // The cross-contig ignore set, which is empty here.
        output.write_int(self.read_group_to_library.len() as i32);
        for (group, library) in self.read_group_to_library.iter() {
            output.write_string(Some(group));
            output.write_string(library.as_deref());
        }
        output.write_int(self.contig_ids.len() as i32);
        for (name, id) in self.contig_ids.iter() {
            output.write_string(Some(name));
            output.write_int(*id);
        }
        output.write_long(self.reads);
        output.write_int(self.avg_read_len);
        output.write_long(self.ref_bases);
        output.write_long(self.max_reads_in_partition);
        output.write_float(self.coverage);
        output.write_float(self.mean_base_quality);
        output.write_int(1);
        let bounds = &self.partition;
        for value in [
            bounds.first_contig,
            bounds.first_start,
            bounds.last_contig,
            bounds.last_end,
            bounds.span,
        ] {
            output.write_int(value);
        }
        output.write_int(self.libraries.len() as i32);
        for (name, stats) in self.libraries.iter() {
            output.write_string(name.0.as_deref());
            output.write_int(stats.cdf.fractions.len() as i32);
            for fraction in &stats.cdf.fractions {
                output.write_float(*fraction);
            }
            output.write_long(stats.cdf.counts);
            output.write_float(stats.coverage);
            output.write_float(stats.mean_base_quality);
            output.write_long(stats.reads);
            output.write_float(stats.read_start_frequency);
        }
        output.flushed()
    }
}

/// Kryo's `Output` over a stream, as far as what reached the stream.
///
/// `require(n)` flushes when fewer than `n` bytes are free, so the survivors end wherever the last
/// write that did not fit began. An ASCII string fills the buffer to its last byte first
/// (`writeAscii_slow`), and the length form asks for a byte at a time.
#[derive(Default)]
struct FlushedOutput {
    bytes: Vec<u8>,
    flushed: usize,
}

const KRYO_BUFFER: usize = 4096;

impl FlushedOutput {
    fn require(&mut self, needed: usize) {
        if KRYO_BUFFER - (self.bytes.len() - self.flushed) < needed {
            self.flushed = self.bytes.len();
        }
    }

    fn write_int(&mut self, value: i32) {
        self.require(4);
        self.bytes.extend_from_slice(&value.to_be_bytes());
    }

    fn write_long(&mut self, value: i64) {
        self.require(8);
        self.bytes.extend_from_slice(&value.to_be_bytes());
    }

    /// `writeFloat` is `writeInt(floatToIntBits)`.
    fn write_float(&mut self, value: f32) {
        self.write_int(float_to_int_bits(value));
    }

    fn write_string(&mut self, value: Option<&str>) {
        let mut encoded = Output::new();
        encoded.write_string(value);
        let encoded = encoded.bytes();
        let characters = value.map_or(0, |text| text.chars().count());
        if characters > 1 && value.is_some_and(str::is_ascii) {
            let mut remaining = encoded;
            while !remaining.is_empty() {
                let free = KRYO_BUFFER - (self.bytes.len() - self.flushed);
                if free == 0 {
                    self.flushed = self.bytes.len();
                    continue;
                }
                let step = free.min(remaining.len());
                self.bytes.extend_from_slice(&remaining[..step]);
                remaining = &remaining[step..];
            }
        } else {
            for byte in encoded {
                self.require(1);
                self.bytes.push(*byte);
            }
        }
    }

    /// The stream's bytes: what is still in the buffer never gets there.
    fn flushed(mut self) -> Vec<u8> {
        self.bytes.truncate(self.flushed);
        self.bytes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_search_over_nothing_lands_at_zero() {
        let cdf = Cdf::new(&[0, 0, 0, 0], 0);
        assert!(cdf.fractions.iter().all(|fraction| fraction.is_nan()));
        assert_eq!(cdf.pop_stat(0.5), 0);
    }

    #[test]
    fn a_median_is_the_first_bin_at_or_over_the_fraction() {
        let cdf = Cdf::new(&[0, 1, 1, 2], 4);
        assert_eq!(cdf.pop_stat(0.25), 1);
        assert_eq!(cdf.pop_stat(0.5), 2);
        assert_eq!(cdf.pop_stat(0.75), 3);
    }

    #[test]
    fn only_what_overflowed_the_buffer_reaches_the_file() {
        let mut output = FlushedOutput::default();
        for value in 0..1023 {
            output.write_int(value);
        }
        // 4092 bytes held: the next long does not fit, so they all go out and it stays behind.
        output.write_long(7);
        assert_eq!(output.flushed().len(), 4092);
    }
}
