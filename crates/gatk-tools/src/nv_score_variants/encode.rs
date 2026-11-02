//! `scorevariants/encoders.py` and `readers.py`: a variant and its reads as the model's tensors.
//!
//! Each `Err` is a Python exception the reference would have ended on. The tool's Java half turns
//! any of them into `Tool returned: 1` and a success.

use super::vcf::{Definition, InfoEntry, InfoValue, Kind, Number, Record};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use std::collections::{BTreeMap, HashMap};

/// The tensors' edge: reads by positions, and the window of reference.
pub const WINDOW: usize = 128;
const READ_LIMIT: usize = 128;
/// `ReadTensorEncoder._create_channel_map`.
pub const READ_CHANNELS: usize = 15;
const READ_INDEL: usize = 4;
const REFERENCE_BASE: usize = 5;
const REFERENCE_INDEL: usize = 9;
const FLAG_REVERSE: usize = 10;
const FLAG_MATE_REVERSE: usize = 11;
const FLAG_FIRST: usize = 12;
const FLAG_SECOND: usize = 13;
const MAPPING_QUALITY: usize = 14;
const MAPPING_QUALITY_MAX: f64 = 60.0;

const SKIP: u8 = b'~';
const INDEL: u8 = b'*';

/// `SeqIO.to_dict(SeqIO.parse(reference_file, format="fasta"))`, upper-cased where it is used.
pub fn read_fasta(path: &str) -> Result<HashMap<String, Vec<u8>>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("FileNotFoundError: {path}: {e}"))?;
    let mut contigs = HashMap::new();
    let mut current: Option<(String, Vec<u8>)> = None;
    for line in bytes.split(|b| *b == b'\n') {
        if line.first() == Some(&b'>') {
            if let Some((name, sequence)) = current.take() {
                if contigs.insert(name.clone(), sequence).is_some() {
                    return Err(format!("ValueError: Duplicate key '{name}'"));
                }
            }
            let title = String::from_utf8_lossy(&line[1..]);
            let name = title.split_whitespace().next().unwrap_or("").to_string();
            current = Some((name, Vec::new()));
        } else if let Some((_, sequence)) = &mut current {
            sequence.extend(line.iter().filter(|b| **b != b' ' && **b != b'\r'));
        }
    }
    if let Some((name, sequence)) = current.take() {
        if contigs.insert(name.clone(), sequence).is_some() {
            return Err(format!("ValueError: Duplicate key '{name}'"));
        }
    }
    Ok(contigs)
}

/// `get_variant_window`.
fn variant_window(position: i64, reference_size: i64) -> (i64, i64) {
    let half = WINDOW as i64 / 2;
    let start = position - half;
    let mut end = position + half;
    end += (-start).max(0).min(reference_size - 1);
    (start, end)
}

/// `_type_of_biallelic_variant` and `determine_polymorphic_type`, through `VariantTypeEncoder`: 0 a
/// SNP, 1 an indel, and -1 for a multinucleotide change, a mixture and no alternate at all.
pub fn variant_type(record: &Record) -> i64 {
    #[derive(PartialEq, Clone, Copy)]
    enum Biallelic {
        Snp,
        Mnp,
        Indel,
    }
    let mut kind = None;
    for allele in &record.alternates {
        if *allele == record.reference {
            continue;
        }
        let one = if allele.len() == record.reference.len() {
            if allele.len() == 1 {
                Biallelic::Snp
            } else {
                Biallelic::Mnp
            }
        } else {
            Biallelic::Indel
        };
        match kind {
            None => kind = Some(one),
            Some(existing) if existing != one => return -1,
            Some(_) => {}
        }
    }
    match kind {
        Some(Biallelic::Snp) => 0,
        Some(Biallelic::Indel) => 1,
        _ => -1,
    }
}

/// Python's `float(text)`: blanks around it, and underscores between digits.
fn python_float(text: &str) -> Option<f64> {
    let trimmed = text.trim();
    let bytes = trimmed.as_bytes();
    let misplaced = bytes.iter().enumerate().any(|(at, byte)| {
        *byte == b'_'
            && !(at > 0
                && bytes[at - 1].is_ascii_digit()
                && bytes.get(at + 1).is_some_and(u8::is_ascii_digit))
    });
    if misplaced {
        return None;
    }
    trimmed.replace('_', "").parse().ok()
}

/// `encoding[i] = variant.info.get(key, 0.0)`: `bcf_info_get_value`, then what NumPy makes of it.
///
/// `tools/nvscorevariants/probe_edges.py` measures each arm. Only `Number=1`, or no definition at
/// all, is a scalar to pysam; anything else is a tuple, even of one, which NumPy will not take. What
/// is missing is `None`, which it takes as a NaN, and text goes through `float`.
fn annotation(entry: &InfoEntry, definition: Option<&Definition>) -> Result<f64, String> {
    let sequence_error = || "ValueError: setting an array element with a sequence.".to_string();
    if definition.is_some_and(|d| d.kind == Kind::Flag) {
        return Ok(1.0);
    }
    let scalar = definition.is_none_or(|d| d.number == Number::Fixed(1));
    if !scalar {
        return Err(sequence_error());
    }
    match &entry.value {
        InfoValue::Nothing => Ok(f64::NAN),
        InfoValue::Integers(values) => match values[..] {
            [value] => Ok(value.map_or(f64::NAN, f64::from)),
            _ => Err(sequence_error()),
        },
        InfoValue::Floats(values) => match values[..] {
            [value] => Ok(value.map_or(f64::NAN, f64::from)),
            _ => Err(sequence_error()),
        },
        InfoValue::Text(text) => python_float(text)
            .ok_or_else(|| format!("ValueError: could not convert string to float: '{text}'")),
    }
}

/// `AnnotationEncoder` with `BEST_PRACTICES`.
pub fn annotations(
    record: &Record,
    definitions: &HashMap<String, Definition>,
) -> Result<[f32; 7], String> {
    const BEST_PRACTICES: [&str; 7] = ["MQ", "DP", "SOR", "FS", "QD", "MQRankSum", "ReadPosRankSum"];
    let mut out = [0f32; 7];
    for (slot, key) in out.iter_mut().zip(BEST_PRACTICES) {
        if let Some(entry) = record.info.iter().find(|e| e.key == key) {
            *slot = annotation(entry, definitions.get(key))? as f32;
        }
    }
    Ok(out)
}

/// `ReferenceEncoder` with `_LeftStrategy` and `DNAVectorEncoder("N")`: `(window, 4)`.
pub fn reference_tensor(
    contigs: &HashMap<String, Vec<u8>>,
    record: &Record,
) -> Result<Vec<f32>, String> {
    let (start, end) = variant_window(record.pos - 1, record.reference.len() as i64);
    let mut tensor = vec![0f32; WINDOW * 4];
    let contig = contigs
        .get(&record.contig)
        .ok_or_else(|| format!("KeyError: '{}'", record.contig))?;
    let interval_start = start.max(0);
    let interval_end = end.min(contig.len() as i64) + i64::from(start < 0);
    for i in interval_start..interval_end {
        let Some(base) = contig.get(i as usize) else {
            return Err("IndexError: index out of range".to_string());
        };
        let position = (i - interval_start) as usize;
        if position >= WINDOW {
            return Err(format!(
                "IndexError: index {position} is out of bounds for axis 0 with size {WINDOW}"
            ));
        }
        let at = position * 4;
        match base.to_ascii_uppercase() {
            b'A' => tensor[at] = 1.0,
            b'C' => tensor[at + 1] = 1.0,
            b'G' => tensor[at + 2] = 1.0,
            b'T' => tensor[at + 3] = 1.0,
            // `AMBIGUITY_CODES`: all a quarter, through the encoder's getter.
            b'K' | b'M' | b'R' | b'Y' | b'S' | b'W' | b'B' | b'V' | b'H' | b'D' | b'X' | b'N' => {
                tensor[at..at + 4].fill(0.25);
            }
            _ => {}
        }
    }
    Ok(tensor)
}

/// `random_generator.Random`, which is `java.util.Random` by way of NumPy's wrapping integers.
struct JavaRandom {
    seed: i64,
}

const MULTIPLIER: i64 = 0x5DEECE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;

impl JavaRandom {
    fn new(seed: i64) -> JavaRandom {
        JavaRandom {
            seed: (seed ^ MULTIPLIER) & MASK,
        }
    }

    fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND)) & MASK;
        (self.seed >> (48 - bits)) as i32
    }

    fn next_int(&mut self, n: i64) -> i32 {
        if n & -n == n {
            return ((n.wrapping_mul(i64::from(self.next(31)))) >> 31) as i32;
        }
        loop {
            let bits = i64::from(self.next(31));
            let value = bits % n;
            if (bits - value + (n - 1)) as i32 >= 0 {
                return value as i32;
            }
        }
    }
}

/// `pysam`'s `query_alignment_start`.
fn query_alignment_start(read: &BamRecord) -> i64 {
    let mut offset = 0i64;
    for element in &read.cigar.elements {
        match element.op {
            Op::H => {}
            Op::S => offset += i64::from(element.length),
            _ => break,
        }
    }
    offset
}

/// `get_start`.
fn start_of(read: &BamRecord) -> i64 {
    if read.read_unmapped() {
        return -1;
    }
    i64::from(read.alignment_start) - 1 - query_alignment_start(read)
}

/// `read.cigartuples`, as `encoders.py`'s codes.
fn cigar_tuples(read: &BamRecord) -> Vec<(u8, i64)> {
    read.cigar
        .elements
        .iter()
        .map(|e| (e.op as u8, i64::from(e.length)))
        .collect()
}

const CIGAR_M: u8 = 0;
const CIGAR_I: u8 = 1;
const CIGAR_D: u8 = 2;
const CIGAR_S: u8 = 4;
const CODE_TO_CIGAR: &[u8] = b"MIDNSHP=XB";

fn counts(code: u8) -> bool {
    matches!(code, CIGAR_M | CIGAR_I | CIGAR_D | CIGAR_S)
}

/// `utilities.variant_is_snp`.
fn variant_is_snp(record: &Record) -> bool {
    record.reference.len() <= 1
        && !record.alternates.is_empty()
        && record
            .alternates
            .iter()
            .all(|alt| matches!(alt.as_str(), "A" | "C" | "G" | "T" | "N" | "*"))
}

/// `get_base_to_sort_by`.
fn base_to_sort_by(read: &BamRecord, record: &Record, snp: bool) -> u8 {
    if read.read_bases.is_empty() {
        return b'Z';
    }
    let last = read.read_bases.len() as i64 - 1;
    let start = start_of(read);
    if snp {
        return read.read_bases[(record.pos - start - 1).min(last).max(0) as usize];
    }
    let variant_index = record.pos - start;
    let mut index = 0;
    for (code, length) in cigar_tuples(read) {
        index += length;
        if index > variant_index {
            return if code == CIGAR_M {
                read.read_bases[variant_index.min(last).max(0) as usize]
            } else {
                CODE_TO_CIGAR[code as usize]
            };
        }
    }
    b'Y'
}

/// `base_quality_to_p_hot_array`'s two numbers for every quality a byte can hold, as floats.
///
/// `10.0 ** exponent` is glibc's `pow`, and it goes through `1.0 - x` and a rounding to a float
/// afterwards. `StrictMath`'s `pow` comes to the same floats for all 256, which is what the test
/// below pins from `tools/nvscorevariants/probe_edges.py`, and no host library is involved.
fn phred_table() -> Vec<(f32, f32)> {
    (0..=255)
        .map(|quality: i32| {
            let exponent = f64::from(-quality) / 10.0;
            let p = 1.0 - jmath::strict_math::pow(10.0, exponent);
            let not_p = (1.0 - p) / 3.0;
            (p as f32, not_p as f32)
        })
        .collect()
}

/// `AMBIGUITY_CODES`, for a read's base.
fn read_ambiguity(base: u8) -> Option<[f64; 4]> {
    Some(match base {
        b'K' => [0.0, 0.0, 0.5, 0.5],
        b'M' => [0.5, 0.5, 0.0, 0.0],
        b'R' => [0.5, 0.0, 0.0, 0.5],
        b'Y' => [0.0, 0.5, 0.5, 0.0],
        b'S' => [0.0, 0.5, 0.0, 0.5],
        b'W' => [0.5, 0.0, 0.5, 0.0],
        b'B' => [0.0, 0.333, 0.333, 0.334],
        b'V' => [0.333, 0.333, 0.0, 0.334],
        b'H' => [0.333, 0.333, 0.334, 0.0],
        b'D' => [0.333, 0.0, 0.333, 0.334],
        b'X' | b'N' => [0.25, 0.25, 0.25, 0.25],
        _ => return None,
    })
}

/// `_add_insertions` over any sequence. A read's quality and a read's bases go through it in step.
fn add_insertions<T: Clone>(
    sequence: Vec<T>,
    insertions: &BTreeMap<i64, i64>,
    fill: &T,
    beginning: i64,
    padded: Option<&T>,
) -> Vec<T> {
    let mut sequence = sequence;
    for (i, count) in insertions.iter().rev() {
        let insertion = vec![fill.clone(); (*count).max(0) as usize];
        if *i < beginning {
            let mut joined = insertion;
            joined.extend(sequence);
            sequence = joined;
        } else if padded.is_none() || *i as usize <= sequence.len() {
            let cut = (*i as usize).min(sequence.len());
            let tail = sequence.split_off(cut);
            sequence.extend(insertion);
            sequence.extend(tail);
        } else if (*i as usize) < WINDOW {
            if let Some(nothing) = padded {
                sequence.extend(std::iter::repeat_n(nothing.clone(), *i as usize - sequence.len()));
            }
            sequence.extend(insertion);
        }
    }
    sequence
}

/// `ReadTensorEncoder`, with its reservoir, whose state runs across the whole file.
pub struct ReadTensorEncoder {
    random: JavaRandom,
    phred: Vec<(f32, f32)>,
}

impl Default for ReadTensorEncoder {
    fn default() -> Self {
        Self::new()
    }
}

/// What `filter_read` keeps, or the exception one of its `len()`s or the unbound name ends on.
fn keeps(read: &BamRecord) -> Result<bool, String> {
    if read.base_qualities.is_empty() {
        return Err("TypeError: object of type 'NoneType' has no len()".to_string());
    }
    if read.read_bases.is_empty() {
        return Err("TypeError: object of type 'NoneType' has no len()".to_string());
    }
    if !read.read_unmapped() && read.cigar.elements.is_empty() {
        return Err("TypeError: object of type 'NoneType' has no len()".to_string());
    }
    let Some(group) = read_group(read) else {
        return Err(
            "UnboundLocalError: local variable 'read_group' referenced before assignment"
                .to_string(),
        );
    };
    Ok(!group.to_ascii_lowercase().contains("artificial"))
}

fn read_group(read: &BamRecord) -> Option<String> {
    match read.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(group)) => Some(group.clone()),
        // `get_tag` answers with whatever it is, and `.lower()` is what then fails.
        Some(TagValue::Char(group)) => Some((*group as char).to_string()),
        _ => None,
    }
}

impl ReadTensorEncoder {
    pub fn new() -> ReadTensorEncoder {
        ReadTensorEncoder {
            random: JavaRandom::new(47382911),
            phred: phred_table(),
        }
    }

    /// `__call__`: `(reads, window, channels)`. `fetched` is `fetch(contig, start, stop)`.
    pub fn encode(
        &mut self,
        contigs: &HashMap<String, Vec<u8>>,
        record: &Record,
        fetched: Vec<BamRecord>,
    ) -> Result<Vec<f32>, String> {
        let (start, end) = variant_window(record.pos - 1, record.reference.len() as i64);
        // `get_reads`.
        let mut reads: Vec<BamRecord> = Vec::new();
        let mut alignment_number = 0i64;
        for read in fetched {
            if !keeps(&read)? {
                continue;
            }
            alignment_number += 1;
            if alignment_number as usize <= READ_LIMIT {
                reads.push(read);
            } else {
                let slot = self.random.next_int(alignment_number);
                if slot < READ_LIMIT as i32 {
                    // A negative slot would count from the end, as a Python index does.
                    let at = if slot < 0 { reads.len() as i32 + slot } else { slot };
                    reads[at as usize] = read;
                }
            }
        }
        let insertions = get_insertions(&reads, start);
        let snp = variant_is_snp(record);
        reads.sort_by_key(start_of);
        reads.sort_by_key(|read| base_to_sort_by(read, record, snp));

        let mut tensor = vec![0f32; READ_LIMIT * WINDOW * READ_CHANNELS];
        // `get_reference`.
        let contig = contigs
            .get(&record.contig)
            .ok_or_else(|| format!("KeyError: '{}'", record.contig))?;
        let from = start.max(0) as usize;
        let to = ((end + i64::from(start < 0)).max(0) as usize).min(contig.len());
        let slice: Vec<u8> = if from < to {
            contig[from..to].iter().map(u8::to_ascii_uppercase).collect()
        } else {
            Vec::new()
        };
        let reference = add_insertions(slice, &insertions, &INDEL, 0, Some(&0u8));
        // `insert_reference`.
        for (i, base) in reference.iter().take(WINDOW).enumerate() {
            let channel = match base {
                b'A' => Some(REFERENCE_BASE),
                b'C' => Some(REFERENCE_BASE + 1),
                b'G' => Some(REFERENCE_BASE + 2),
                b'T' => Some(REFERENCE_BASE + 3),
                b'*' => Some(REFERENCE_INDEL),
                _ => None,
            };
            for row in 0..READ_LIMIT {
                let cell = (row * WINDOW + i) * READ_CHANNELS;
                if let Some(channel) = channel {
                    tensor[cell + channel] = 1.0;
                } else if read_ambiguity(*base).is_some() {
                    for at in 0..4 {
                        tensor[cell + REFERENCE_BASE + at] = 0.25;
                    }
                }
            }
        }
        // `insert_reads_and_flags`.
        for (row, read) in reads.iter().enumerate() {
            let (sequence, quality) = sequence_and_quality(read, &insertions, start);
            self.fill_row(row, read, &sequence, &quality, &mut tensor);
        }
        Ok(tensor)
    }

    /// `_fill_rows`.
    fn fill_row(
        &self,
        row: usize,
        read: &BamRecord,
        sequence: &[u8],
        quality: &[u8],
        tensor: &mut [f32],
    ) {
        let mut flag_start: i64 = -1;
        let mut flag_end: i64 = 0;
        for (position, (base, qual)) in sequence.iter().zip(quality).enumerate() {
            if position == WINDOW {
                break;
            }
            if *base == SKIP {
                continue;
            }
            if flag_start == -1 {
                flag_start = position as i64;
            } else {
                flag_end = position as i64;
            }
            let cell = (row * WINDOW + position) * READ_CHANNELS;
            match base {
                b'*' => tensor[cell + READ_INDEL] = 1.0,
                b'A' | b'C' | b'G' | b'T' => {
                    let (p, not_p) = self.phred[*qual as usize];
                    let which = match base {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    };
                    for at in 0..4 {
                        tensor[cell + at] = if at == which { p } else { not_p };
                    }
                }
                _ => {
                    if let Some(values) = read_ambiguity(*base) {
                        for (at, value) in values.iter().enumerate() {
                            tensor[cell + at] = *value as f32;
                        }
                    }
                }
            }
        }
        // Python's slice, so a read that never got a second base marks nothing.
        let flag = |bit: u16| if read.flags & bit != 0 { 1.0 } else { 0.0 };
        let mapping_quality = (f64::from(read.mapping_quality) / MAPPING_QUALITY_MAX) as f32;
        let first = flag_start.max(0) as usize;
        for position in first..flag_end.max(0) as usize {
            let cell = (row * WINDOW + position) * READ_CHANNELS;
            tensor[cell + FLAG_REVERSE] = flag(0x10);
            tensor[cell + FLAG_MATE_REVERSE] = flag(0x20);
            tensor[cell + FLAG_FIRST] = flag(0x40);
            tensor[cell + FLAG_SECOND] = flag(0x80);
            tensor[cell + MAPPING_QUALITY] = mapping_quality;
        }
    }
}

/// `get_insertions`.
fn get_insertions(reads: &[BamRecord], interval_start: i64) -> BTreeMap<i64, i64> {
    let mut insertions = BTreeMap::new();
    for read in reads {
        let index_difference = interval_start - start_of(read);
        if index_difference.abs() >= WINDOW as i64 {
            continue;
        }
        let mut index = 0;
        for (code, length) in cigar_tuples(read) {
            if code == CIGAR_I {
                let at = index - index_difference;
                let slot = insertions.entry(at).or_insert(length);
                if *slot < length {
                    *slot = length;
                }
            }
            if counts(code) {
                index += length;
            }
        }
    }
    insertions
}

/// `_get_sequence_and_quality`.
fn sequence_and_quality(
    read: &BamRecord,
    insertions: &BTreeMap<i64, i64>,
    interval_start: i64,
) -> (Vec<u8>, Vec<u8>) {
    let index_difference = interval_start - start_of(read);
    // `_get_sequence_specific_indel_dict`.
    let mut own: BTreeMap<i64, i64> = BTreeMap::new();
    let mut index = 0;
    for (code, length) in cigar_tuples(read) {
        let at = index - index_difference;
        if code == CIGAR_I && insertions.contains_key(&at) {
            own.insert(at, insertions[&at] - length);
        } else if code == CIGAR_D {
            own.insert(at, length);
        }
        if counts(code) {
            index += length;
        }
    }
    for (key, value) in insertions {
        own.entry(*key).or_insert(*value);
    }

    let keep = read.read_bases.len().min(WINDOW);
    let mut sequence = read.read_bases[..keep].to_vec();
    let keep = read.base_qualities.len().min(WINDOW);
    let mut quality = read.base_qualities[..keep].to_vec();
    if index_difference > 0 {
        let cut = index_difference as usize;
        sequence = sequence.get(cut..).unwrap_or(&[]).to_vec();
        quality = quality.get(cut..).unwrap_or(&[]).to_vec();
    } else if index_difference < 0 {
        let pad = (-index_difference) as usize;
        let mut padded = vec![SKIP; pad];
        padded.extend(sequence);
        sequence = padded;
        let mut padded = vec![0u8; pad];
        padded.extend(quality);
        quality = padded;
    }
    let sequence = add_insertions(sequence, &own, &INDEL, 1, None);
    let quality = add_insertions(quality, &own, &0u8, 1, None);
    (sequence, quality)
}

/// `AlignmentFile.fetch`'s test on top of what the index returned: `bam_endpos` against `[beg, end)`.
pub fn overlaps(read: &BamRecord, begin: i64, end: i64) -> bool {
    let start = i64::from(read.alignment_start) - 1;
    let mut length = 0i64;
    if !read.read_unmapped() {
        for element in &read.cigar.elements {
            if matches!(element.op, Op::M | Op::D | Op::N | Op::Eq | Op::X) {
                length += i64::from(element.length);
            }
        }
    }
    if length == 0 {
        length = 1;
    }
    start < end && start + length > begin
}

#[cfg(test)]
mod tests {
    use super::*;
    use sha2::{Digest, Sha256};

    #[test]
    fn the_window_grows_only_at_the_start_of_a_contig() {
        assert_eq!(variant_window(1000, 1), (936, 1064));
        assert_eq!(variant_window(1000, 5), (936, 1064));
        assert_eq!(variant_window(10, 1), (-54, 74));
        assert_eq!(variant_window(10, 4), (-54, 77));
    }

    /// `tools/nvscorevariants/probe_edges.py`, on the reference's own class under NumPy 1.26.
    #[test]
    fn the_reservoir_answers_as_numpy_does() {
        let mut random = JavaRandom::new(47382911);
        let answers: Vec<String> = (129..20129).map(|n| random.next_int(n).to_string()).collect();
        assert_eq!(answers[..8].join(","), "114,25,128,65,80,43,12,82");
        assert_eq!(answers[256 - 129], "54");
        let digest: String = Sha256::digest(answers.join(",").as_bytes())
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        assert_eq!(
            digest,
            "56a7bb1903da461b8362379eba8cc20ab52d47e8d5cc1a320a883232cb3a2022"
        );
    }

    /// The same probe: glibc's `pow` as Python called it, through to the floats the tensor holds.
    #[test]
    fn the_quality_table_is_the_floats_python_gets() {
        let mut bytes = Vec::new();
        for (p, not_p) in phred_table() {
            bytes.extend_from_slice(&p.to_le_bytes());
            bytes.extend_from_slice(&not_p.to_le_bytes());
        }
        let digest: String = Sha256::digest(&bytes).iter().map(|b| format!("{b:02x}")).collect();
        assert_eq!(
            digest,
            "eb86654aff4bb3bcfe0ce8957efb16ef21d264805585ba47771b8fb6608e2df9"
        );
    }

    #[test]
    fn python_float_takes_what_python_takes() {
        assert_eq!(python_float(" 2.5"), Some(2.5));
        assert_eq!(python_float("1_0"), Some(10.0));
        assert_eq!(python_float("7"), Some(7.0));
        assert!(python_float("nan").is_some_and(f64::is_nan));
        assert_eq!(python_float("x"), None);
        assert_eq!(python_float("."), None);
        assert_eq!(python_float("1,2"), None);
        assert_eq!(python_float("_1"), None);
    }
}
