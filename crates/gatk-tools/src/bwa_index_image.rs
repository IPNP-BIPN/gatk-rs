//! The BWA-MEM index image `BwaMemIndexImageCreator` writes, built from the FASTA the way BWA
//! builds it.
//!
//! GATK hands the FASTA to BWA through JNI: `bwa_idx_build` packs the sequence (`bns_fasta2bntseq`),
//! builds the Burrows-Wheeler transform of the forward strand followed by its reverse complement
//! (`bwt_pac2bwt`), interleaves the occurrence counts (`bwt_bwtupdate_core`), samples the suffix
//! array every 32 ranks (`bwt_cal_sa`), and `bwa_idx2mem` then lays the three structures out in one
//! block of memory which is written as it stands. Every byte of that block is a function of the
//! FASTA except the POINTERS the C structures carry, which are addresses of the process that wrote
//! them: see `docs/pointers-that-reach-the-output.md`.
//!
//! This builds the same block. Where BWA writes an address it writes zero, which is what BWA's own
//! loader (`bwa_mem2idx`) overwrites anyway: the image is usable, and comparing it with the
//! reference's means masking the fields [`pointer_fields`] names on both sides.
//!
//! Ported from `bwa/bntseq.c`, `bwa/bwtindex.c`, `bwa/bwt.c`, `bwa/bwa.c` and `bwa/kseq.h` of the
//! BWA that `gatk-bwamem-jni` 1.0.4 ships in GATK 4.6.2.0, each of which carries an MIT header.
//! The suffix array is NOT ported from `is.c`: it is unique, and [`suffix_array`] is an
//! independent implementation of the published induced-sorting algorithm.

/// `sizeof(bwt_t)`: primary, `L2[5]`, `seq_len`, `bwt_size`, the `bwt` pointer, `cnt_table[256]`,
/// `sa_intv` with its padding, `n_sa` and the `sa` pointer.
pub const BWT_STRUCT_BYTES: usize = 1120;
/// `sizeof(bntseq_t)`.
pub const BNS_STRUCT_BYTES: usize = 48;
/// `sizeof(bntann1_t)`.
pub const ANN_STRUCT_BYTES: usize = 40;
/// `sizeof(bntamb1_t)`.
pub const AMB_STRUCT_BYTES: usize = 16;
/// `bwa_idx_build`'s suffix-array sampling interval.
pub const SA_INTERVAL: usize = 32;
/// `OCC_INTERVAL`: one block of four counts per 128 bases of the transform.
const OCC_INTERVAL: usize = 128;
/// The seed `bns_fasta2bntseq` hands `srand48` before replacing ambiguous bases.
pub const BNS_SEED: u32 = 11;

/// One FASTA record as `kseq_read` returns it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub name: Vec<u8>,
    pub comment: Vec<u8>,
    pub sequence: Vec<u8>,
}

/// `kseq_read` over a whole file, FASTA branch.
///
/// A record starts at the first `>` or `@` anywhere, its name runs to the first white-space
/// character, and the comment is what follows it on the line. Sequence lines are read whole until a
/// line starting with `>`, `+` or `@`; empty lines are skipped, anything else, spaces included, is a
/// sequence character. A `+` line ends the record the way a FASTQ quality line would, and what
/// follows it up to the next header is skipped here, which is where kseq's quality reading leaves
/// the stream for a FASTA record.
pub fn read_records(text: &[u8]) -> Vec<Record> {
    let mut records = Vec::new();
    let mut at = 0;
    // Jump to the first header character.
    while at < text.len() && text[at] != b'>' && text[at] != b'@' {
        at += 1;
    }
    while at < text.len() {
        // `at` is on a header character.
        at += 1;
        let start = at;
        while at < text.len() && !is_c_space(text[at]) {
            at += 1;
        }
        let name = text[start..at].to_vec();
        let mut comment = Vec::new();
        if at < text.len() && text[at] != b'\n' {
            at += 1;
            let start = at;
            while at < text.len() && text[at] != b'\n' {
                at += 1;
            }
            comment = strip_carriage_return(&text[start..at]).to_vec();
        }
        if at < text.len() {
            at += 1;
        }
        let mut sequence = Vec::new();
        let mut ended_on = None;
        while at < text.len() {
            let c = text[at];
            if c == b'>' || c == b'+' || c == b'@' {
                ended_on = Some(c);
                break;
            }
            at += 1;
            if c == b'\n' {
                continue;
            }
            let start = at - 1;
            while at < text.len() && text[at] != b'\n' {
                at += 1;
            }
            sequence.extend_from_slice(strip_carriage_return(&text[start..at]));
            if at < text.len() {
                at += 1;
            }
        }
        if ended_on == Some(b'+') {
            // kseq reads a quality string: the rest of the `+` line is skipped, then whole lines
            // are taken until there are as many characters as bases. A different count is an
            // error, and `bns_fasta2bntseq` stops reading at the first one, record dropped.
            while at < text.len() && text[at] != b'\n' {
                at += 1;
            }
            if at < text.len() {
                at += 1;
            }
            let mut quality = 0usize;
            while at < text.len() && quality < sequence.len() {
                let start = at;
                while at < text.len() && text[at] != b'\n' {
                    at += 1;
                }
                quality += strip_carriage_return(&text[start..at]).len();
                if at < text.len() {
                    at += 1;
                }
            }
            if quality != sequence.len() {
                break;
            }
            records.push(Record {
                name,
                comment,
                sequence,
            });
            while at < text.len() && text[at] != b'>' && text[at] != b'@' {
                at += 1;
            }
            continue;
        }
        records.push(Record {
            name,
            comment,
            sequence,
        });
    }
    records
}

/// `ks_getuntil2` with `KS_SEP_LINE` drops a trailing carriage return, so a CRLF file reads as the
/// same bases.
fn strip_carriage_return(line: &[u8]) -> &[u8] {
    if line.len() > 1 && line[line.len() - 1] == b'\r' {
        &line[..line.len() - 1]
    } else {
        line
    }
}

/// C's `isspace` in the "C" locale.
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `nst_nt4_table`: the four bases in either case, and four for everything else.
fn nt4(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// POSIX `drand48`'s generator, which is what `lrand48` draws from.
struct Rand48 {
    state: u64,
}

impl Rand48 {
    /// `srand48(seed)`.
    fn new(seed: u32) -> Self {
        Self {
            state: ((seed as u64) << 16) | 0x330e,
        }
    }

    /// `lrand48()`: the high 31 bits of the next 48-bit state.
    fn next(&mut self) -> u64 {
        self.state =
            (0x5_DEEC_E66D_u64.wrapping_mul(self.state).wrapping_add(0xB)) & ((1u64 << 48) - 1);
        self.state >> 17
    }
}

/// One contig of `bntseq_t`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub offset: i64,
    pub length: i32,
    pub ambiguous_runs: i32,
    pub name: Vec<u8>,
    pub annotation: Vec<u8>,
}

/// One run of ambiguous bases, `bntamb1_t`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hole {
    pub offset: i64,
    pub length: i32,
    pub base: u8,
}

/// What `bns_fasta2bntseq` makes of the records: the bases as 2-bit codes, ambiguous ones replaced
/// by `lrand48() & 3`, and the contigs and holes the `.ann` and `.amb` files hold.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Packed {
    pub bases: Vec<u8>,
    pub annotations: Vec<Annotation>,
    pub holes: Vec<Hole>,
}

/// What BWA refuses once the FASTA has been read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    /// `bns_restore_core` cannot read back an `.amb` line whose base is white space, and calls
    /// `err_fatal`, which exits the process with status one.
    AmbiguousBaseIsWhiteSpace,
}

/// `bns_fasta2bntseq`, forward strand only, followed by the `.ann`/`.amb` round trip
/// `bwa_idx_load_from_disk` makes.
pub fn pack(records: &[Record]) -> Result<Packed, BuildError> {
    let mut random = Rand48::new(BNS_SEED);
    let mut bases = Vec::new();
    let mut annotations = Vec::new();
    let mut holes: Vec<Hole> = Vec::new();
    for record in records {
        let offset = bases.len() as i64;
        let mut runs = 0;
        let mut last = 0u8;
        for (i, &c) in record.sequence.iter().enumerate() {
            let mut code = nt4(c);
            if code >= 4 {
                if last == c {
                    if let Some(hole) = holes.last_mut() {
                        hole.length += 1;
                    }
                } else {
                    holes.push(Hole {
                        offset: offset + i as i64,
                        length: 1,
                        base: c,
                    });
                    runs += 1;
                }
            }
            last = c;
            if code >= 4 {
                code = (random.next() & 3) as u8;
            }
            bases.push(code);
        }
        let (name, annotation) = restored_name_and_annotation(&record.name, &record.comment);
        annotations.push(Annotation {
            offset,
            length: record.sequence.len() as i32,
            ambiguous_runs: runs,
            name,
            annotation,
        });
    }
    if holes.iter().any(|hole| is_c_space(hole.base)) {
        return Err(BuildError::AmbiguousBaseIsWhiteSpace);
    }
    Ok(Packed {
        bases,
        annotations,
        holes,
    })
}

/// The name and comment as `bns_restore_core` reads them back from the `.ann` line `bns_dump`
/// wrote, which is `<gi> <name>` followed by ` <comment>` when there is one, `(null)` standing for
/// none. The name is read with `%s`, which skips white space first and stops at the next, so a
/// header whose name is empty (`>  chr1 c`) comes back as `chr1` with the comment `c`; the comment
/// is the rest of the line without its first character, and ` (null)` is read as empty.
fn restored_name_and_annotation(name: &[u8], comment: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut line = name.to_vec();
    let written: &[u8] = if comment.is_empty() {
        b"(null)"
    } else {
        comment
    };
    line.push(b' ');
    line.extend_from_slice(written);
    let mut at = 0;
    while at < line.len() && is_c_space(line[at]) {
        at += 1;
    }
    let start = at;
    while at < line.len() && !is_c_space(line[at]) {
        at += 1;
    }
    let restored = line[start..at].to_vec();
    let rest = &line[at..];
    let annotation = if rest.len() > 1 && rest != b" (null)" {
        rest[1..].to_vec()
    } else {
        Vec::new()
    };
    (restored, annotation)
}

/// The suffix array of `text`, whose last symbol is a unique smallest sentinel, by induced sorting
/// (Nong, Zhang and Chan, 2009). The array is unique, so any correct construction gives BWA's.
fn suffix_array(text: &[u32], alphabet: usize) -> Vec<usize> {
    let n = text.len();
    if n == 1 {
        return vec![0];
    }
    // S-type is true, L-type false.
    let mut small = vec![false; n];
    small[n - 1] = true;
    for i in (0..n - 1).rev() {
        small[i] = text[i] < text[i + 1] || (text[i] == text[i + 1] && small[i + 1]);
    }
    let is_lms = |i: usize| i > 0 && small[i] && !small[i - 1];
    let mut counts = vec![0usize; alphabet];
    for &c in text {
        counts[c as usize] += 1;
    }
    let heads = |counts: &[usize]| {
        let mut out = vec![0usize; counts.len()];
        let mut sum = 0;
        for (c, &count) in counts.iter().enumerate() {
            out[c] = sum;
            sum += count;
        }
        out
    };
    let tails = |counts: &[usize]| {
        let mut out = vec![0usize; counts.len()];
        let mut sum = 0;
        for (c, &count) in counts.iter().enumerate() {
            sum += count;
            out[c] = sum;
        }
        out
    };
    const EMPTY: usize = usize::MAX;
    let induce = |order: &[usize]| {
        let mut sa = vec![EMPTY; n];
        let mut end = tails(&counts);
        for &p in order.iter().rev() {
            let c = text[p] as usize;
            end[c] -= 1;
            sa[end[c]] = p;
        }
        let mut head = heads(&counts);
        for i in 0..n {
            let j = sa[i];
            if j != EMPTY && j > 0 && !small[j - 1] {
                let c = text[j - 1] as usize;
                sa[head[c]] = j - 1;
                head[c] += 1;
            }
        }
        let mut end = tails(&counts);
        for i in (0..n).rev() {
            let j = sa[i];
            if j != EMPTY && j > 0 && small[j - 1] {
                let c = text[j - 1] as usize;
                end[c] -= 1;
                sa[end[c]] = j - 1;
            }
        }
        sa
    };
    let lms: Vec<usize> = (0..n).filter(|&i| is_lms(i)).collect();
    let sa = induce(&lms);
    let sorted: Vec<usize> = sa.into_iter().filter(|&i| is_lms(i)).collect();
    // Name the LMS substrings in sorted order.
    let mut names = vec![EMPTY; n];
    let mut name = 0usize;
    let mut previous: Option<usize> = None;
    for &p in &sorted {
        if let Some(q) = previous {
            let mut d = 0;
            let same = loop {
                if text[p + d] != text[q + d] || small[p + d] != small[q + d] {
                    break false;
                }
                if d > 0 && (is_lms(p + d) || is_lms(q + d)) {
                    break is_lms(p + d) && is_lms(q + d);
                }
                d += 1;
            };
            if !same {
                name += 1;
            }
        }
        names[p] = name;
        previous = Some(p);
    }
    let order = if name + 1 == sorted.len() {
        sorted
    } else {
        let reduced: Vec<u32> = lms.iter().map(|&p| names[p] as u32).collect();
        let reduced_sa = suffix_array(&reduced, name + 1);
        reduced_sa.into_iter().map(|r| lms[r]).collect()
    };
    induce(&order)
}

/// `bwt_t` as `bwa_idx_load_from_disk` holds it after `bwt_restore_bwt` and `bwt_restore_sa`.
struct Transform {
    primary: u64,
    l2: [u64; 5],
    seq_len: u64,
    /// The interleaved words: occurrence counts every 128 bases, then the bases 16 to a word.
    words: Vec<u32>,
    sa: Vec<u64>,
}

/// `bwt_pac2bwt` with `is_bwt`, `bwt_bwtupdate_core` and `bwt_cal_sa(bwt, 32)`, over the forward
/// strand followed by its reverse complement.
fn transform(forward: &[u8]) -> Transform {
    let mut text: Vec<u8> = forward.to_vec();
    text.extend(forward.iter().rev().map(|&b| 3 - b));
    let n = text.len();
    let mut l2 = [0u64; 5];
    for &c in &text {
        l2[1 + c as usize] += 1;
    }
    for i in 2..=4 {
        l2[i] += l2[i - 1];
    }
    let shifted: Vec<u32> = text
        .iter()
        .map(|&c| c as u32 + 1)
        .chain(std::iter::once(0))
        .collect();
    let sa = suffix_array(&shifted, 5);
    let mut primary = 0usize;
    let mut bwt = Vec::with_capacity(n);
    for (rank, &position) in sa.iter().enumerate() {
        if position == 0 {
            primary = rank;
        } else {
            bwt.push(text[position - 1]);
        }
    }
    let mut packed = vec![0u32; n.div_ceil(16)];
    for (i, &c) in bwt.iter().enumerate() {
        packed[i >> 4] |= (c as u32) << ((15 - (i & 15)) << 1);
    }
    let occurrences = n.div_ceil(OCC_INTERVAL) + 1;
    let mut words = Vec::with_capacity(packed.len() + occurrences * 8);
    let mut counts = [0u64; 4];
    let push_counts = |words: &mut Vec<u32>, counts: &[u64; 4]| {
        for &count in counts {
            words.push(count as u32);
            words.push((count >> 32) as u32);
        }
    };
    for (i, &c) in bwt.iter().enumerate() {
        if i % OCC_INTERVAL == 0 {
            push_counts(&mut words, &counts);
        }
        if i % 16 == 0 {
            words.push(packed[i / 16]);
        }
        counts[c as usize] += 1;
    }
    push_counts(&mut words, &counts);
    let samples = (n + SA_INTERVAL) / SA_INTERVAL;
    let mut sampled = Vec::with_capacity(samples);
    sampled.push(u64::MAX);
    for k in 1..samples {
        sampled.push(sa[k * SA_INTERVAL] as u64);
    }
    Transform {
        primary: primary as u64,
        l2,
        seq_len: n as u64,
        words,
        sa: sampled,
    }
}

/// `bwt_gen_cnt_table`: for each byte of four 2-bit bases, how many of each base it holds.
fn count_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    for (i, entry) in table.iter_mut().enumerate() {
        let mut x = 0u32;
        for j in 0..4usize {
            let hits = ((i & 3) == j) as u32
                + ((i >> 2 & 3) == j) as u32
                + ((i >> 4 & 3) == j) as u32
                + ((i >> 6) == j) as u32;
            x |= hits << (j << 3);
        }
        *entry = x;
    }
    table
}

/// The whole image for a FASTA's bytes: `bwa_idx_build`, `bwa_idx_load_from_disk` and
/// `bwa_idx2mem`, with zero where BWA writes an address.
pub fn build_image(fasta: &[u8]) -> Result<Vec<u8>, BuildError> {
    let packed = pack(&read_records(fasta))?;
    let bwt = transform(&packed.bases);
    let mut out = Vec::new();
    let u64le = |out: &mut Vec<u8>, v: u64| out.extend_from_slice(&v.to_le_bytes());
    let u32le = |out: &mut Vec<u8>, v: u32| out.extend_from_slice(&v.to_le_bytes());
    // bwt_t
    u64le(&mut out, bwt.primary);
    for v in bwt.l2 {
        u64le(&mut out, v);
    }
    u64le(&mut out, bwt.seq_len);
    u64le(&mut out, bwt.words.len() as u64);
    u64le(&mut out, 0); // the `bwt` pointer, which `bwa_idx2mem` clears before the copy
    for v in count_table() {
        u32le(&mut out, v);
    }
    u32le(&mut out, SA_INTERVAL as u32);
    u32le(&mut out, 0); // padding
    u64le(&mut out, bwt.sa.len() as u64);
    u64le(&mut out, 0); // the `sa` pointer
    debug_assert_eq!(out.len(), BWT_STRUCT_BYTES);
    for &w in &bwt.words {
        u32le(&mut out, w);
    }
    for &s in &bwt.sa {
        u64le(&mut out, s);
    }
    // bntseq_t
    let l_pac = packed.bases.len();
    u64le(&mut out, l_pac as u64);
    u32le(&mut out, packed.annotations.len() as u32);
    u32le(&mut out, BNS_SEED);
    u64le(&mut out, 0); // anns
    u32le(&mut out, packed.holes.len() as u32);
    u32le(&mut out, 0); // padding
    u64le(&mut out, 0); // ambs
    u64le(&mut out, 0); // fp_pac
    for hole in &packed.holes {
        u64le(&mut out, hole.offset as u64);
        u32le(&mut out, hole.length as u32);
        out.push(hole.base);
        out.extend_from_slice(&[0, 0, 0]);
    }
    for ann in &packed.annotations {
        u64le(&mut out, ann.offset as u64);
        u32le(&mut out, ann.length as u32);
        u32le(&mut out, ann.ambiguous_runs as u32);
        u32le(&mut out, 0); // gi
        u32le(&mut out, 0); // is_alt
        u64le(&mut out, 0); // name
        u64le(&mut out, 0); // anno
    }
    for ann in &packed.annotations {
        out.extend_from_slice(&ann.name);
        out.push(0);
        out.extend_from_slice(&ann.annotation);
        out.push(0);
    }
    // The forward-only pac, `l_pac / 4 + 1` bytes.
    let mut pac = vec![0u8; l_pac / 4 + 1];
    for (i, &b) in packed.bases.iter().enumerate() {
        pac[i >> 2] |= b << ((!i & 3) << 1);
    }
    out.extend_from_slice(&pac);
    Ok(out)
}

/// The byte ranges of an image that hold an address or a structure's padding, which are the bytes
/// no two runs of the reference agree on. `None` when the bytes are not an image this layout reads.
pub fn pointer_fields(image: &[u8]) -> Option<Vec<std::ops::Range<usize>>> {
    let read64 = |at: usize| -> Option<u64> {
        image
            .get(at..at + 8)
            .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
    };
    let read32 = |at: usize| -> Option<u32> {
        image
            .get(at..at + 4)
            .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
    };
    let bwt_size = read64(56)? as usize;
    let n_sa = read64(1104)? as usize;
    let mut fields = vec![64..72, 1100..1104, 1112..1120];
    let bns = BWT_STRUCT_BYTES
        .checked_add(bwt_size.checked_mul(4)?)?
        .checked_add(n_sa.checked_mul(8)?)?;
    let n_seqs = read32(bns + 8)? as usize;
    let n_holes = read32(bns + 24)? as usize;
    fields.extend([bns + 16..bns + 24, bns + 28..bns + 48]);
    let holes = bns + BNS_STRUCT_BYTES;
    for h in 0..n_holes {
        let at = holes + h * AMB_STRUCT_BYTES;
        fields.push(at + 13..at + 16);
    }
    let anns = holes + n_holes * AMB_STRUCT_BYTES;
    for a in 0..n_seqs {
        let at = anns + a * ANN_STRUCT_BYTES;
        fields.push(at + 24..at + 40);
    }
    if fields.iter().any(|range| range.end > image.len()) {
        return None;
    }
    Some(fields)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive(text: &[u32]) -> Vec<usize> {
        let mut sa: Vec<usize> = (0..text.len()).collect();
        sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
        sa
    }

    #[test]
    fn induced_sorting_agrees_with_sorting_the_suffixes() {
        let mut random = Rand48::new(7);
        for length in [0usize, 1, 2, 3, 5, 17, 64, 200, 1000] {
            for _ in 0..20 {
                let mut text: Vec<u32> = (0..length)
                    .map(|_| (random.next() % 4) as u32 + 1)
                    .collect();
                text.push(0);
                assert_eq!(suffix_array(&text, 5), naive(&text));
            }
        }
        // Long runs, which is where the recursion is exercised.
        let mut text: Vec<u32> = std::iter::repeat_n([1, 2, 1, 2, 3], 40).flatten().collect();
        text.push(0);
        assert_eq!(suffix_array(&text, 5), naive(&text));
    }

    #[test]
    fn the_image_is_the_size_bwa_writes() {
        // Five repeats of eight bases, one contig named chr1: 1333 bytes in the golden.
        let fasta = format!(">chr1\n{}\n", "ACGTTGCA".repeat(5));
        assert_eq!(build_image(fasta.as_bytes()).unwrap().len(), 1333);
        let fasta = format!(">chr1\n{}\n", "ACGTTGCA".repeat(20));
        assert_eq!(build_image(fasta.as_bytes()).unwrap().len(), 1551);
    }
}
