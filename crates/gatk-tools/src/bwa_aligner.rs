//! `BwaSparkEngine.ReadAligner`: one partition of reads through BWA-MEM and back into records.
//!
//! GATK calls BWA through `gatk-bwamem-jni` and builds each record itself in
//! `BwaMemAlignmentUtils.applyAlignment`. Here the aligning is `bwa-mem4`, which answers with SAM,
//! so the alignment's own fields are read from its lines and what `applyAlignment` takes from the
//! input (the name, the bases as stored, the qualities, the read group) comes from the input.
//! Measured against the JNI on 6000 paired and 6205 single-end records, those fields agree
//! exactly, except that `bwa-mem4` adds an `MC` tag the JNI has no field for.
//!
//! The batch is the partition, as the JNI's is: `-K` is set past anything a partition can hold,
//! because a batch's insert-size statistics decide the pairs in it.

use std::path::Path;

use clap::Parser;
use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};

/// `bwa-mem4 index`, whose arguments its own binary parses as a subcommand.
#[derive(Parser)]
struct IndexCommand {
    #[command(flatten)]
    args: bwa_mem4::cmd_index::IndexArgs,
}

/// `bwa-mem4 mem`, the same way.
#[derive(Parser)]
struct MemCommand {
    #[command(flatten)]
    args: bwa_mem4::cmd_mem::MemArgs,
}

/// The contigs an image holds, in its order.
pub fn image_contigs(image: &[u8]) -> Result<Vec<String>, String> {
    let fasta = fasta_from_image(image)?;
    Ok(String::from_utf8_lossy(&fasta)
        .lines()
        .filter_map(|line| line.strip_prefix('>').map(str::to_string))
        .collect())
}

/// `FmIndex` and `BntSeq` read a prefix on disk, and the image is the only thing naming the contigs.
pub fn fasta_from_image(image: &[u8]) -> Result<Vec<u8>, String> {
    let bad = || "not a BWA-MEM index image".to_string();
    let read64 = |at: usize| -> Result<u64, String> {
        image
            .get(at..at + 8)
            .map(|b| u64::from_le_bytes(b.try_into().expect("eight bytes")))
            .ok_or_else(bad)
    };
    let read32 = |at: usize| -> Result<u32, String> {
        image
            .get(at..at + 4)
            .map(|b| u32::from_le_bytes(b.try_into().expect("four bytes")))
            .ok_or_else(bad)
    };
    let bwt_size = read64(56)? as usize;
    let n_sa = read64(1104)? as usize;
    let bns = crate::bwa_index_image::BWT_STRUCT_BYTES
        .checked_add(bwt_size.checked_mul(4).ok_or_else(bad)?)
        .and_then(|at| at.checked_add(n_sa.checked_mul(8)?))
        .ok_or_else(bad)?;
    let l_pac = read64(bns)? as usize;
    let n_seqs = read32(bns + 8)? as usize;
    let n_holes = read32(bns + 24)? as usize;
    let holes = bns + crate::bwa_index_image::BNS_STRUCT_BYTES;
    let anns = holes + n_holes * crate::bwa_index_image::AMB_STRUCT_BYTES;
    let mut text = anns + n_seqs * crate::bwa_index_image::ANN_STRUCT_BYTES;
    let mut names = Vec::with_capacity(n_seqs);
    for _ in 0..n_seqs {
        let end_of = |from: usize| {
            image
                .get(from..)
                .and_then(|rest| rest.iter().position(|&byte| byte == 0))
                .map(|length| from + length)
                .ok_or_else(bad)
        };
        let name_end = end_of(text)?;
        names.push(image[text..name_end].to_vec());
        text = end_of(name_end + 1)? + 1;
    }
    let pac = image.get(text..text + l_pac / 4 + 1).ok_or_else(bad)?;
    let mut bases: Vec<u8> = (0..l_pac)
        .map(|i| b"ACGT"[usize::from(pac[i >> 2] >> ((!i & 3) << 1) & 3)])
        .collect();
    // A hole goes back to being a hole, so that the index's own `lrand48` fills it again.
    for hole in 0..n_holes {
        let at = holes + hole * crate::bwa_index_image::AMB_STRUCT_BYTES;
        let offset = read64(at)? as usize;
        let length = read32(at + 8)? as usize;
        bases
            .get_mut(offset..offset + length)
            .ok_or_else(bad)?
            .fill(b'N');
    }
    let mut fasta = Vec::with_capacity(l_pac + n_seqs * 32);
    for (index, name) in names.iter().enumerate() {
        let at = anns + index * crate::bwa_index_image::ANN_STRUCT_BYTES;
        let offset = read64(at)? as usize;
        let length = read32(at + 8)? as usize;
        fasta.push(b'>');
        fasta.extend_from_slice(name);
        fasta.push(b'\n');
        fasta.extend_from_slice(bases.get(offset..offset + length).ok_or_else(bad)?);
        fasta.push(b'\n');
    }
    Ok(fasta)
}

/// `BwaMemAligner.alignSeqs`: what BWA says of each read, a list of lines per read. `options` are
/// the aligner's setters as `bwa mem` spells them.
pub fn alignments(
    reads: &[BamRecord],
    paired: bool,
    image: &[u8],
    header: &SamHeader,
    options: &[String],
    threads: usize,
) -> Result<Vec<Vec<BamRecord>>, String> {
    let mut answers: Vec<Vec<BamRecord>> = vec![Vec::new(); reads.len()];
    if reads.is_empty() {
        return Ok(answers);
    }
    let directory = tempfile::tempdir().map_err(|error| error.to_string())?;
    let sam = run_bwa(reads, paired, image, options, directory.path(), threads)?;
    let resolve = |name: &str| {
        header
            .sequences
            .iter()
            .position(|sequence| sequence.name == name)
            .map(|index| index as i32)
    };
    for line in sam.lines() {
        if line.starts_with('@') {
            continue;
        }
        let mut fields: Vec<&str> = line.split('\t').collect();
        // `applyAlignment` looks the contig up by its name, whatever the header holds.
        if let Some(contig) = fields
            .get(2)
            .filter(|contig| **contig != "*" && resolve(contig).is_none())
        {
            return Err(format!(
                "Reference index for '{contig}' not found in sequence dictionary."
            ));
        }
        // A read of no bases is answered with nothing between the tabs.
        for field in [9, 10] {
            if let Some(value) = fields.get_mut(field).filter(|value| value.is_empty()) {
                *value = "*";
            }
        }
        let record = htsjdk_bam::text_parse::parse_line(&fields.join("\t"), resolve)
            .map_err(|error| format!("{error:?}"))?;
        let index = record_index(&record, paired)?;
        answers
            .get_mut(index)
            .ok_or_else(|| format!("BWA answered for read {index} of {}", reads.len()))?
            .push(record);
    }
    Ok(answers)
}

/// `ReadAligner.apply`: every alignment of every read, in the reads' order.
pub fn align(
    reads: &[BamRecord],
    paired: bool,
    image: &[u8],
    header: &SamHeader,
    threads: usize,
) -> Result<Vec<BamRecord>, String> {
    if paired {
        if !reads.len().is_multiple_of(2) {
            return Err(
                "We're supposed to be aligning paired reads, but there are an odd number of them."
                    .to_string(),
            );
        }
        for pair in reads.chunks(2) {
            if pair[0].read_name != pair[1].read_name {
                return Err(format!(
                    "Read pair has varying template name: {} .vs {}",
                    pair[0].read_name, pair[1].read_name
                ));
            }
        }
    }
    let answers = alignments(reads, paired, image, header, &[], threads)?;
    let mut aligned = Vec::with_capacity(reads.len());
    for (read, lines) in reads.iter().zip(answers) {
        for alignment in lines {
            aligned.push(apply_alignment(read, alignment, true));
        }
    }
    Ok(aligned)
}

/// `PSBwaAligner.apply`: each read takes its primary line and nothing else, under the aligner's own
/// options, and with no `AS` or `XS` on one that stayed unmapped.
pub fn align_for_pathseq(
    reads: &[BamRecord],
    paired: bool,
    image: &[u8],
    header: &SamHeader,
    options: &[String],
    threads: usize,
) -> Result<Vec<BamRecord>, String> {
    let answers = alignments(reads, paired, image, header, options, threads)?;
    let mut aligned = Vec::with_capacity(reads.len());
    for (read, lines) in reads.iter().zip(answers) {
        let primary = lines
            .into_iter()
            .find(|line| line.flags & (0x100 | 0x800) == 0);
        aligned.push(match primary {
            Some(line) => apply_alignment(read, line, false),
            None => read.clone(),
        });
    }
    Ok(aligned)
}

/// `PSBwaFilter.apply`: the reads no line of which is a good enough match to the host. `keeps` is
/// `HostAlignmentReadFilter.test` over a line's CIGAR and `NM`.
pub fn host_filter(
    reads: Vec<BamRecord>,
    image: &[u8],
    header: &SamHeader,
    min_seed_length: i32,
    keeps: impl Fn(&BamRecord) -> bool,
    threads: usize,
) -> Result<Vec<BamRecord>, String> {
    let options = [
        "-h".to_string(),
        "0,0".to_string(),
        "-k".to_string(),
        min_seed_length.to_string(),
        "-T".to_string(),
        "0".to_string(),
    ];
    let answers = alignments(&reads, false, image, header, &options, threads)?;
    Ok(reads
        .into_iter()
        .zip(answers)
        .filter(|(_, lines)| {
            lines
                .iter()
                .all(|line| line.cigar.is_empty() || keeps(line))
        })
        .map(|(read, _)| read)
        .collect())
}

/// Which input read a line answers: the names handed to BWA are positions, since a real name could
/// lose a `/1` to the FASTQ reader. A pair's name is its first read's, and the flag says which.
fn record_index(record: &BamRecord, paired: bool) -> Result<usize, String> {
    let position: usize = record.read_name.parse().map_err(|_| {
        format!(
            "BWA answered for a read it was not asked about: {}",
            record.read_name
        )
    })?;
    Ok(if paired && record.flags & 0x80 != 0 {
        position + 1
    } else {
        position
    })
}

/// `applyAlignment` over one line: the bases and qualities are the input's, turned round for a
/// reverse strand and cut to what a supplementary line's hard clips leave.
fn apply_alignment(read: &BamRecord, alignment: BamRecord, generic: bool) -> BamRecord {
    const REVERSE: u16 = 0x10;
    const SECONDARY: u16 = 0x100;
    let mut bases = read.read_bases.clone();
    let mut qualities = read.base_qualities.clone();
    if alignment.flags & REVERSE != 0 && alignment.flags & SECONDARY == 0 {
        bases.reverse();
        for base in &mut bases {
            *base = match *base {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                other => other,
            };
        }
        qualities.reverse();
    }
    let (start, length) = (
        leading_hard_clip(&alignment.cigar),
        alignment.read_bases.len(),
    );
    if alignment.read_bases.len() != bases.len() {
        bases = bases[start..start + length].to_vec();
        if !qualities.is_empty() {
            qualities = qualities[start..start + length].to_vec();
        }
    }
    // `createSATags` builds the list BWA wrote, and the record takes it last.
    let mut tags = Tags::new();
    if alignment.cigar.is_empty() {
        if generic {
            tags.insert(Tag::new(b"AS"), TagValue::Int(0));
            tags.insert(Tag::new(b"XS"), TagValue::Int(0));
        }
    } else {
        for name in [b"NM", b"AS", b"XS", b"MD", b"XA"] {
            if let Some(value) = alignment.tags.get(Tag::new(name)) {
                tags.insert(Tag::new(name), value.clone());
            }
        }
    }
    if let Some(group) = read.tags.get(Tag::new(b"RG")) {
        tags.insert(Tag::new(b"RG"), group.clone());
    }
    // `BwaMemAlignmentUtils.makeSAMRecords` writes it, which the PathSeq aligner does not go through.
    if generic {
        if let Some(value) = alignment.tags.get(Tag::new(b"SA")) {
            tags.insert(Tag::new(b"SA"), value.clone());
        }
    }
    BamRecord {
        read_name: read.read_name.clone(),
        read_bases: bases,
        base_qualities: qualities,
        tags,
        ..alignment
    }
}

fn leading_hard_clip(cigar: &Cigar) -> usize {
    match cigar.elements.first() {
        Some(element) if element.op.to_char() == b'H' => element.length as usize,
        _ => 0,
    }
}

/// The reads as FASTQ, the image as an index, and `bwa mem` over them: its SAM.
fn run_bwa(
    reads: &[BamRecord],
    paired: bool,
    image: &[u8],
    options: &[String],
    directory: &Path,
    threads: usize,
) -> Result<String, String> {
    use std::io::Write;
    let io = |error: std::io::Error| error.to_string();
    let fasta = directory.join("reference.fa");
    std::fs::write(&fasta, fasta_from_image(image)?).map_err(io)?;
    let index = IndexCommand::try_parse_from(["index", &*fasta.to_string_lossy()])
        .map_err(|error| error.to_string())?;
    bwa_mem4::cmd_index::run(index.args).map_err(|error| error.to_string())?;
    let fastq = directory.join("reads.fq");
    let mut out = std::io::BufWriter::new(std::fs::File::create(&fastq).map_err(io)?);
    for (index, read) in reads.iter().enumerate() {
        // BWA reads no quality, so none needs to be right here.
        let name = if paired { index / 2 * 2 } else { index };
        writeln!(out, "@{name}").map_err(io)?;
        out.write_all(&read.read_bases).map_err(io)?;
        writeln!(out, "\n+").map_err(io)?;
        out.write_all(&vec![b'I'; read.read_bases.len()])
            .map_err(io)?;
        out.write_all(b"\n").map_err(io)?;
    }
    out.flush().map_err(io)?;
    drop(out);
    let sam = directory.join("out.sam");
    let threads = threads.max(1).to_string();
    let mut argv = vec![
        "mem",
        "-v",
        "0",
        "-t",
        threads.as_str(),
        "-K",
        "1000000000000",
        "-o",
    ];
    let sam_text = sam.to_string_lossy().into_owned();
    argv.push(&sam_text);
    if paired {
        argv.push("-p");
    }
    argv.extend(options.iter().map(String::as_str));
    let fasta_text = fasta.to_string_lossy().into_owned();
    let fastq_text = fastq.to_string_lossy().into_owned();
    argv.push(&fasta_text);
    argv.push(&fastq_text);
    let argv: Vec<String> = argv.into_iter().map(str::to_string).collect();
    let mem = MemCommand::try_parse_from(&argv).map_err(|error| error.to_string())?;
    bwa_mem4::cmd_mem::run(mem.args, &argv).map_err(|error| error.to_string())?;
    std::fs::read_to_string(&sam).map_err(io)
}

/// What `BwaMemAlignment` says of one line of a contig, with the JNI's fields as `alignSeqs` fills
/// them in. `ref_id` is -1 on the one line of a contig that did not align.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigAlignment {
    pub ref_id: i32,
    pub start: i64,
    pub reverse: bool,
    pub supplementary: bool,
    pub mapq: u32,
    pub cigar: String,
    pub mismatches: i32,
    pub score: i32,
    pub suboptimal: i32,
    pub md: String,
    pub xa: Option<String>,
}

/// `BwaMemAligner` after `setIntraCtgOptions` and `setZDropOption`, for `FindBreakpointEvidenceSpark`.
///
/// `setIntraCtgOptions` writes its fields into the C struct through a `ByteBuffer`, and libbwa fills
/// its scoring matrix only once, when the options are initialized. So a mismatch still costs what
/// the matrix says, 4, while `b` is 9 for what reads it directly, the MAPQ's identity and the
/// near-tie window. `bwa mem -B 9` refills the matrix and could not give this, which is why this is
/// the library and not the command.
pub struct ContigAligner {
    fm: bwa_index::FmIndex,
    bns: bwa_index::BntSeq,
    options: bwa_mem4_core::MemOpt,
    _directory: tempfile::TempDir,
}

impl ContigAligner {
    pub fn from_image(image: &[u8], z_dropoff: i32) -> Result<ContigAligner, String> {
        let directory = tempfile::tempdir().map_err(|error| error.to_string())?;
        let fasta = directory.path().join("reference.fa");
        std::fs::write(&fasta, fasta_from_image(image)?).map_err(|error| error.to_string())?;
        let index = IndexCommand::try_parse_from(["index", &*fasta.to_string_lossy()])
            .map_err(|error| error.to_string())?;
        bwa_mem4::cmd_index::run(index.args).map_err(|error| error.to_string())?;
        let fm = bwa_index::FmIndex::load(&fasta).map_err(|error| error.to_string())?;
        let bns = bwa_index::BntSeq::load(&fasta).map_err(|error| error.to_string())?;
        let options = bwa_mem4_core::MemOpt {
            o_del: 16,
            o_ins: 16,
            b: 9,
            pen_clip5: 5,
            pen_clip3: 5,
            zdrop: z_dropoff,
            ..bwa_mem4_core::MemOpt::default()
        };
        Ok(ContigAligner {
            fm,
            bns,
            options,
            _directory: directory,
        })
    }

    /// `alignSeqs` of one assembly's contigs. A contig's number in the call breaks a score tie in
    /// `mem_mark_primary_se`, so the call has to be one assembly's, as it is in the reference.
    pub fn align(&self, contigs: &[Vec<u8>]) -> Vec<Vec<ContigAlignment>> {
        let (fm, bns, options) = (&self.fm, &self.bns, &self.options);
        contigs
            .iter()
            .enumerate()
            .map(|(id, contig)| {
                let codes: Vec<u8> = contig
                    .iter()
                    .map(|&base| bwa_mem4_core::dna::nt4(base))
                    .collect();
                let length = codes.len() as i32;
                let found = bwa_mem4_mem::align_read(fm, bns, options, &codes);
                let mut regions = bwa_mem4_mem::mem_sort_dedup_patch(fm, options, &codes, found);
                bwa_mem4_mem::stamp_is_alt(bns, &mut regions);
                bwa_mem4_mem::mem_mark_primary_se(options, &mut regions, id as u64);
                let alternates =
                    bwa_mem4_mem::alt::mem_gen_alt(fm, bns, options, &regions, length, &codes);
                let mut lines: Vec<ContigAlignment> = Vec::new();
                for (region, xa) in regions.iter().zip(alternates) {
                    if region.score < options.t || region.secondary >= 0 {
                        continue;
                    }
                    let aln = bwa_mem4_mem::reg2aln(fm, bns, options, length, &codes, region);
                    let mut mapq = bwa_mem4_mem::mem_approx_mapq_se(options, region);
                    if let Some(first) = lines.first() {
                        if !region.is_alt && mapq > first.mapq {
                            mapq = first.mapq;
                        }
                    }
                    lines.push(ContigAlignment {
                        ref_id: aln.rid,
                        start: aln.pos,
                        reverse: aln.is_rev,
                        supplementary: !lines.is_empty(),
                        mapq,
                        cigar: bwa_mem4_mem::cigar_string(&aln.cigar),
                        mismatches: aln.nm,
                        score: aln.score,
                        suboptimal: aln.sub,
                        md: aln.md,
                        xa,
                    });
                }
                if lines.is_empty() {
                    lines.push(ContigAlignment {
                        ref_id: -1,
                        start: -1,
                        reverse: false,
                        supplementary: false,
                        mapq: 0,
                        cigar: String::new(),
                        mismatches: 0,
                        score: 0,
                        suboptimal: 0,
                        md: String::new(),
                        xa: None,
                    });
                }
                lines
            })
            .collect()
    }
}
