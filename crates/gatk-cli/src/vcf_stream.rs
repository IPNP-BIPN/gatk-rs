//! A VCF's records one at a time, for a tool that never needs two of them at once.
//!
//! `htsjdk_vcf::reader::read_vcf` decodes the whole file before the first record is looked at,
//! which on a chromosome of gnomAD is five million contexts held at once. This is its record loop
//! (the `Eager` branch of `read_vcf_with`) as an iterator: the same counter, the same `#` lines
//! skipped without a count, and the same lazy genotype context over the line's own text.

use htsjdk_vcf::genotype_parse::{parse_genotypes, GenotypeContext};
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::header_parse::VcfVersion;
use htsjdk_vcf::reader::VcfFile;
use htsjdk_vcf::record_parse::{decode_line, RecordError};
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::VcfHeader;
use std::sync::Arc;

/// Where the header ends: just past the `#CHROM` line, or `None` when there is none.
pub fn header_end(text: &str) -> Option<usize> {
    let mut at = 0;
    for line in text.split_inclusive('\n') {
        at += line.len();
        if line.starts_with("#CHROM") {
            return Some(at);
        }
    }
    None
}

/// A file's lines as `str::lines` splits them, read as they are needed.
///
/// A block compressed file is UTF-8 or refused, as `decompress_all` and `String::from_utf8` are
/// over the whole of it; a plain one is read lossily, line by line, which is the same text because
/// a newline is never part of a multi-byte sequence.
pub struct LineReader {
    reader: Box<dyn std::io::BufRead>,
    strict: bool,
    bytes: Vec<u8>,
}

impl LineReader {
    pub fn new(reader: Box<dyn std::io::BufRead>, strict: bool) -> Self {
        LineReader {
            reader,
            strict,
            bytes: Vec::new(),
        }
    }

    /// The next line into `line`, without its `\n` or `\r\n`; `false` at the end.
    pub fn next_line(&mut self, line: &mut String) -> std::io::Result<bool> {
        self.bytes.clear();
        if self.reader.read_until(b'\n', &mut self.bytes)? == 0 {
            return Ok(false);
        }
        if self.bytes.last() == Some(&b'\n') {
            self.bytes.pop();
            if self.bytes.last() == Some(&b'\r') {
                self.bytes.pop();
            }
        }
        line.clear();
        if self.strict {
            let text = std::str::from_utf8(&self.bytes)
                .map_err(|error| std::io::Error::new(std::io::ErrorKind::InvalidData, error))?;
            line.push_str(text);
        } else {
            line.push_str(&String::from_utf8_lossy(&self.bytes));
        }
        Ok(true)
    }

    /// The lines up to and including `#CHROM`, or every line when there is none, and how many.
    pub fn read_header(&mut self) -> std::io::Result<(String, usize)> {
        let (mut header, mut count, mut line) = (String::new(), 0, String::new());
        while self.next_line(&mut line)? {
            header.push_str(&line);
            header.push('\n');
            count += 1;
            if line.starts_with("#CHROM") {
                break;
            }
        }
        Ok((header, count))
    }
}

/// Where the records come from: the rest of a text already in memory, or the rest of a file.
pub enum Body<'a> {
    Text(std::str::Lines<'a>),
    Stream(LineReader),
}

/// Why a record could not be read: the codec refused its line, or the file could not be read.
#[derive(Debug)]
pub enum Failure {
    Record(RecordError),
    Read(std::io::Error),
}

/// The records after a header `read_vcf` established from `header_lines` lines of text.
pub struct Records<'a> {
    body: Body<'a>,
    line: String,
    header: &'a VcfHeader,
    version: VcfVersion,
    line_number: usize,
    samples: Arc<[String]>,
    samples_were_already_sorted: bool,
}

impl<'a> Records<'a> {
    /// `body` is what follows the header, which was `header_lines` lines long.
    pub fn new(body: Body<'a>, file: &'a VcfFile, header_lines: usize) -> Self {
        let header = &file.header;
        Records {
            body,
            line: String::new(),
            header,
            version: file.codec_version,
            line_number: header_lines,
            samples: Arc::from(header.samples.clone()),
            samples_were_already_sorted: header.samples.windows(2).all(|pair| {
                htsjdk_vcf::genotype_type::compare_sample_names(&pair[0], &pair[1])
                    != std::cmp::Ordering::Greater
            }),
        }
    }
}

impl Iterator for Records<'_> {
    type Item = Result<VariantContext, Failure>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match &mut self.body {
                Body::Text(lines) => {
                    let next = lines.next()?;
                    self.line.clear();
                    self.line.push_str(next);
                }
                Body::Stream(reader) => match reader.next_line(&mut self.line) {
                    Ok(true) => {}
                    Ok(false) => return None,
                    Err(error) => return Some(Err(Failure::Read(error))),
                },
            }
            let line = self.line.as_str();
            let decoded = match decode_line(line, self.header, self.line_number, self.version) {
                Ok(Some(decoded)) => decoded,
                // A `#` line: no record and no increment of the counter.
                Ok(None) => continue,
                Err(error) => return Some(Err(Failure::Record(error))),
            };
            self.line_number += 1;
            let mut variant = decoded.variant;
            if let Some(block) = decoded.genotype_block {
                let context = GenotypeContext {
                    line,
                    header: self.header,
                    version: self.version,
                    contig: &variant.contig,
                    pos: variant.start,
                    line_number: self.line_number,
                };
                match parse_genotypes(&block, &variant.alleles, &context) {
                    Ok(genotypes) => {
                        variant.genotypes =
                            GenotypesContext::lazy(genotypes, block, self.samples.clone());
                        if !self.samples_were_already_sorted {
                            variant.genotypes.decode();
                        }
                    }
                    Err(error) => return Some(Err(Failure::Record(error))),
                }
            }
            return Some(Ok(variant));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The iterator yields what `read_vcf` holds, record for record, refusals included, whether
    /// it reads the text or a stream of it.
    #[test]
    fn records_are_read_vcfs() {
        let text = "##fileformat=VCFv4.2\n\
                    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tb\ta\r\n\
                    1\t10\t.\tA\tC\t.\tPASS\t.\tGT:GQ\t0/1:3\t1/1:4\n\
                    #a comment\n\
                    1\t20\t.\tG\tT\t.\t.\t.\tGT\t0/0\t./.\r\n\
                    1\tx\t.\tG\tT\t.\t.\t.\tGT\t0/0\t./.";
        let whole = htsjdk_vcf::reader::read_vcf(text).unwrap_err();
        let expected = match &whole.error {
            htsjdk_vcf::reader::ReadError::Record(error) => error.message().to_string(),
            other => panic!("{other:?}"),
        };
        let check = |streamed: Vec<Result<VariantContext, Failure>>| {
            assert_eq!(streamed.len(), 3);
            let records: Vec<VariantContext> = streamed[..2]
                .iter()
                .map(|record| record.as_ref().unwrap().clone())
                .collect();
            assert_eq!(records, whole.records);
            match &streamed[2] {
                Err(Failure::Record(error)) => assert_eq!(error.message(), expected),
                other => panic!("{other:?}"),
            }
        };

        let end = header_end(text).unwrap();
        let file = htsjdk_vcf::reader::read_vcf(&text[..end]).unwrap();
        let lines = text[..end].lines().count();
        check(Records::new(Body::Text(text[end..].lines()), &file, lines).collect());

        let mut reader = LineReader::new(Box::new(text.as_bytes()), true);
        let (header, count) = reader.read_header().unwrap();
        assert_eq!(count, lines);
        let file = htsjdk_vcf::reader::read_vcf(&header).unwrap();
        check(Records::new(Body::Stream(reader), &file, count).collect());
    }
}
