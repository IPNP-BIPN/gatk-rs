//! `SplitVcfs` one record at a time, for the runs nothing refuses.
//!
//! Picard never holds the file: each record is read, typed, written and indexed before the next.
//! The exact path (`exact`, beside it) reads the whole input, keeps both outputs as text and indexes
//! them afterwards, which costs several times the input in memory. This module is the same run
//! streamed, and it answers only for a run it can finish: anything that would refuse, or that it
//! cannot show writes the same bytes, gives `false` before an output exists, and the exact path
//! runs from the start and refuses (or writes) as it always did. That covers every refusal there
//! is, the order of refusals included, because the exact path is the one that decides them.
//!
//! # Why it writes the same bytes
//!
//! * The input is read with `next_block`'s checks on every block, as `decompress_all` reads it,
//!   and split into lines at `\n`, which no UTF-8 sequence contains, so the lossy decoding of
//!   each piece is the lossy decoding of the whole.
//! * Each piece is parsed as the header followed by its lines, so every record is parsed against
//!   the same header as in the whole file; only a refusal's line number would differ, and a
//!   refusal is the exact path's.
//! * Each record is encoded with the header's encoder, as `write_vcf` encodes it, after the
//!   header's own text.
//! * A text file's `.idx` is [`OnTheFlyIndex`], which the exact path's index is built by too.
//! * A block-compressed file's `.tbi` takes each feature at the writer's file pointer before its
//!   line. `BgzfWriter` deflates a block the moment it is full, so that pointer is the block map's
//!   virtual offset of the line's start, which is what `build_tabix` uses; the index is closed at
//!   the start of the terminator block, as there. A feature the creator or `check_ordering` would
//!   refuse stops the stream instead.

use std::fs::File;
use std::io::{BufWriter, Read, Write};

use gatk_barclay::Parser;
use gatk_tools::feature_codec::{feature_of_line, Codec};
use gatk_tools::index_feature_file::default_output;
use gatk_tools::split_vcfs::{record_type, RecordType};
use htsjdk_bgzf::read::BgzfReader;
use htsjdk_bgzf::BgzfWriter;
use htsjdk_tribble::index_write::Feature;
use htsjdk_tribble::tabix::{FeatureRef, TabixFormat, TabixIndexCreator};
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::reader::read_vcf;

use super::super::vcf_format_converter::output_type;
use super::block;
use crate::runners::{
    argument, dictionary_from_any, flag, header_dictionary, java_absolute_path, modified_millis,
    output_compression, OnTheFlyIndex,
};

/// How much text each parse takes: enough that the header parsed again with every piece is lost
/// in it, little enough that the records of one piece are a few megabytes.
const PIECE: usize = 1 << 22;

/// The run, streamed. `false` when the exact path has to answer instead; nothing is left behind.
pub(super) fn split(parser: &Parser, input: &str, snp_output: &str, indel_output: &str) -> bool {
    let create_index = flag(parser, "CREATE_INDEX");
    // The builder's shared tabix creator (#1450) is the exact path's alone.
    if snp_output == indel_output || (create_index && block(snp_output) && !block(indel_output)) {
        return false;
    }
    let outputs = [snp_output, indel_output];
    let temps = outputs.map(|output| format!("{output}.split-vcfs-part"));
    let finished = run(parser, input, outputs, &temps, create_index).is_some();
    if !finished {
        for temp in &temps {
            let _ = std::fs::remove_file(temp);
        }
    }
    finished
}

fn run(
    parser: &Parser,
    input: &str,
    outputs: [&str; 2],
    temps: &[String; 2],
    create_index: bool,
) -> Option<()> {
    output_type(outputs[0]).ok()?;
    output_type(outputs[1]).ok()?;
    let given: Option<Vec<(String, i32)>> = match argument(parser, "SEQUENCE_DICTIONARY") {
        None => None,
        Some(path) => Some(
            dictionary_from_any(&path)
                .ok()?
                .into_iter()
                .map(|record| (record.name.clone(), record.length))
                .collect(),
        ),
    };
    let strict = flag(parser, "STRICT");

    let mut lines = Lines::open(input)?;
    let header_text = lines.header()?;
    let header = read_vcf(&header_text).ok()?.header;
    let has_dictionary = header
        .lines
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }));
    if create_index && given.is_none() && !has_dictionary {
        return None;
    }
    let encoder = VcfEncoder::new(&header);
    let written_header = header.write();

    let (level, deflater) = output_compression(parser);
    let sink = |at: usize| {
        let compression = block(outputs[at]).then_some((level, deflater));
        Sink::create(&temps[at], compression, create_index)
    };
    let (mut snps, mut indels) = (sink(0)?, sink(1)?);
    snps.write(&written_header)?;
    indels.write(&written_header)?;

    let mut text = header_text.clone();
    let mut line = String::new();
    while let Some(piece) = lines.piece(PIECE)? {
        text.truncate(header_text.len());
        text.push_str(&piece);
        for record in &read_vcf(&text).ok()?.records {
            // The indel test runs first, though no record can answer both.
            let sink = match record_type(record) {
                RecordType::Indel => &mut indels,
                RecordType::Snp => &mut snps,
                _ if strict => return None,
                _ => continue,
            };
            line.clear();
            encoder.encode_into(record, &mut line).ok()?;
            line.push('\n');
            sink.write(&line)?;
        }
    }

    // `snpWriter.close()`, then `indelWriter.close()`.
    let dictionary = given.unwrap_or_else(|| header_dictionary(&written_header));
    snps.finish(outputs[0], level, deflater, &dictionary)?;
    indels.finish(outputs[1], level, deflater, &dictionary)
}

/// The input's text, a line-aligned piece at a time.
struct Lines {
    source: Source,
    buffer: Vec<u8>,
    ended: bool,
}

enum Source {
    Plain(File),
    Block(BgzfReader<File>),
}

impl Lines {
    /// `vcf_text`'s choice: block compressed when the file starts as a gzip member does.
    fn open(path: &str) -> Option<Lines> {
        let mut file = File::open(path).ok()?;
        let mut magic = [0u8; 2];
        let compressed = file.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b];
        let file = File::open(path).ok()?;
        Some(Lines {
            source: if compressed {
                Source::Block(BgzfReader::new(file))
            } else {
                Source::Plain(file)
            },
            buffer: Vec::new(),
            ended: false,
        })
    }

    /// More bytes, or `false` at the end. `None` where `vcf_text` would have refused.
    fn fill(&mut self) -> Option<bool> {
        match &mut self.source {
            Source::Plain(file) => {
                let start = self.buffer.len();
                self.buffer.resize(start + (1 << 20), 0);
                let read = file.read(&mut self.buffer[start..]).ok()?;
                self.buffer.truncate(start + read);
                Some(read > 0)
            }
            Source::Block(reader) => match reader.next_block() {
                Ok(Some(block)) => {
                    self.buffer.extend_from_slice(&block.data);
                    Some(true)
                }
                Ok(None) => Some(false),
                Err(_) => None,
            },
        }
    }

    /// Every line up to and including the first `#CHROM` one, when each line before it is a
    /// `##` line. Any other shape is the exact path's to read.
    fn header(&mut self) -> Option<String> {
        let mut checked = 0;
        loop {
            while let Some(end) = self.buffer[checked..].iter().position(|&b| b == b'\n') {
                let line = &self.buffer[checked..checked + end];
                if line.starts_with(b"#CHROM") {
                    let header: Vec<u8> = self.buffer.drain(..checked + end + 1).collect();
                    return Some(String::from_utf8_lossy(&header).into_owned());
                }
                if !line.starts_with(b"##") {
                    return None;
                }
                checked += end + 1;
            }
            if !self.fill()? {
                return None;
            }
        }
    }

    /// At least `want` bytes of whole lines when the input holds that many, the rest of it
    /// otherwise; `Some(None)` once it is all read.
    fn piece(&mut self, want: usize) -> Option<Option<String>> {
        while !self.ended && self.buffer.len() < want {
            if !self.fill()? {
                self.ended = true;
            }
        }
        let cut = if self.ended {
            self.buffer.len()
        } else {
            match self.buffer.iter().rposition(|&b| b == b'\n') {
                Some(at) => at + 1,
                // One line longer than a piece: read on until it ends.
                None => return self.piece(self.buffer.len() + want),
            }
        };
        if cut == 0 {
            return Some(None);
        }
        let piece: Vec<u8> = self.buffer.drain(..cut).collect();
        Some(Some(String::from_utf8_lossy(&piece).into_owned()))
    }
}

/// One output: a scratch file renamed into place once it is complete, and its index as it grows.
struct Sink {
    temp: String,
    body: Body,
    index: Index,
}

enum Body {
    Text(BufWriter<File>),
    Block(Box<BgzfWriter<BufWriter<File>>>),
}

enum Index {
    None,
    Tribble(OnTheFlyIndex),
    Tabix(TabixIndexCreator, Order),
}

impl Sink {
    /// `compression` is the run's level and deflater for a block-compressed file, `None` for text.
    fn create(
        temp: &str,
        compression: Option<(u32, htsjdk_bgzf::Deflater)>,
        create_index: bool,
    ) -> Option<Sink> {
        let file = BufWriter::with_capacity(1 << 20, File::create(temp).ok()?);
        // Serially: the index asks for the file pointer before every record, which waits for the
        // blocks in flight anyway, and one thread is what the reference uses.
        let (body, index) = if let Some((level, deflater)) = compression {
            (
                Body::Block(Box::new(
                    BgzfWriter::with_deflater(file, level, deflater).parallel(false),
                )),
                if create_index {
                    Index::Tabix(TabixIndexCreator::new(TabixFormat::VCF), Order::default())
                } else {
                    Index::None
                },
            )
        } else {
            (
                Body::Text(file),
                if create_index {
                    Index::Tribble(OnTheFlyIndex::new())
                } else {
                    Index::None
                },
            )
        };
        Some(Sink {
            temp: temp.to_string(),
            body,
            index,
        })
    }

    /// Text that ends in a newline: the header, or one record.
    fn write(&mut self, text: &str) -> Option<()> {
        match (&mut self.body, &mut self.index) {
            (Body::Text(file), index) => {
                if let Index::Tribble(index) = index {
                    for line in text.split_inclusive('\n') {
                        index.add_line(line);
                    }
                }
                file.write_all(text.as_bytes()).ok()
            }
            (Body::Block(writer), index) => {
                if let Index::Tabix(creator, order) = index {
                    for line in text.split_inclusive('\n') {
                        let Some(feature) = feature_of_line(line, Codec::Vcf) else {
                            continue;
                        };
                        order.check(&feature)?;
                        let at = writer.file_pointer().ok()? as i64;
                        let entry = FeatureRef {
                            contig: &feature.contig,
                            start: feature.start,
                            end: feature.end,
                            // Quoted only by the refusal `Order` stops before.
                            description: "",
                            sequence_length: 0,
                        };
                        creator.add_feature(entry, at).ok()?;
                    }
                }
                writer.write_all(text.as_bytes()).ok()
            }
        }
    }

    /// `close()`: the file renamed into place, then its index written beside it.
    fn finish(
        self,
        output: &str,
        level: u32,
        deflater: htsjdk_bgzf::Deflater,
        dictionary: &[(String, i32)],
    ) -> Option<()> {
        match (self.body, self.index) {
            (Body::Text(file), index) => {
                let file = file.into_inner().ok()?;
                let size = file.metadata().ok()?.len() as i64;
                drop(file);
                std::fs::rename(&self.temp, output).ok()?;
                if let Index::Tribble(index) = index {
                    let path = java_absolute_path(output);
                    let bytes = index.finish(dictionary, &path, size, modified_millis(output));
                    std::fs::write(default_output(output), bytes).ok()?;
                }
                Some(())
            }
            (Body::Block(writer), index) => {
                let file = writer.into_inner().ok()?.into_inner().ok()?;
                let size = file.metadata().ok()?.len();
                drop(file);
                let tabix = match index {
                    Index::Tabix(creator, _) => {
                        // Closed before the terminator block, which does not move the address.
                        let end = size.checked_sub(htsjdk_bgzf::EMPTY_GZIP_BLOCK.len() as u64)?;
                        let index = creator.finish((end << 16) as i64).ok()?;
                        let mut writer = BgzfWriter::with_deflater(Vec::new(), level, deflater);
                        writer.write_all(&index.write_body()).ok()?;
                        Some(writer.into_inner().ok()?)
                    }
                    _ => None,
                };
                std::fs::rename(&self.temp, output).ok()?;
                if let Some(tabix) = tabix {
                    std::fs::write(default_output(output), tabix).ok()?;
                }
                Some(())
            }
        }
    }
}

/// `check_ordering` and the creator's own contig test, as one: a start smaller than the one before
/// it on the same contig, or a contig that comes back, stops the stream before either refuses.
#[derive(Default)]
struct Order {
    previous: Option<(String, i32)>,
    visited: std::collections::HashSet<String>,
}

impl Order {
    fn check(&mut self, feature: &Feature) -> Option<()> {
        match &mut self.previous {
            Some((contig, start)) if *contig == feature.contig => {
                if feature.start < *start {
                    return None;
                }
                *start = feature.start;
            }
            _ => {
                if !self.visited.insert(feature.contig.clone()) {
                    return None;
                }
                self.previous = Some((feature.contig.clone(), feature.start));
            }
        }
        Some(())
    }
}
