//! A VCF written as its records arrive (#1452): the file, the md5 and the index
//! [`super::write_variant_file`] writes, without the records or the text ever held whole.
//!
//! The bytes go to a file beside the output and are renamed onto it when the run finishes, so a
//! run that ends in a refusal before the old writer would have written anything still leaves the
//! output as it found it. The index is built as the lines go out: a `.idx` from each line's byte
//! offset, a `.tbi` from the virtual offset the block writer reports before each line, which is
//! what [`super::write_variant_file`] recovers by reading the file back.

use std::io::Write;

use gatk_tools::index_feature_file::{self, TabixStream};
use md5::{Digest, Md5};

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Failure, Thrown, PORT_FAILURE};

/// The file's bytes as they reach it, counted and, when asked, digested.
struct Counted {
    file: std::io::BufWriter<std::fs::File>,
    written: u64,
    md5: Option<Md5>,
}

impl Write for Counted {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let n = self.file.write(bytes)?;
        self.written += n as u64;
        if let Some(md5) = &mut self.md5 {
            md5.update(&bytes[..n]);
        }
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.file.flush()
    }
}

enum Body {
    Plain(Counted),
    Block(Box<htsjdk_bgzf::BgzfWriter<Counted>>),
}

enum Index {
    None,
    Tribble(super::OnTheFlyIndex),
    Tabix(Box<TabixStream>),
    /// `FeatureManager` knows no codec for the name: the refusal `build_tabix` gives.
    Refused(index_feature_file::Refusal),
}

pub(crate) struct VcfStream {
    output: String,
    /// Where the bytes go until [`VcfStream::finish`]; the output itself when it is not a
    /// regular file (`/dev/stdout`), which cannot be renamed onto.
    partial: String,
    body: Option<Body>,
    index: Index,
    dictionary: Vec<(String, i32)>,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
    /// The first write that failed, reported when the file is finished.
    failed: Option<std::io::Error>,
}

impl VcfStream {
    /// The output opened and its header written, with the md5 and index the arguments ask for.
    pub(crate) fn open(parser: &Parser, output: &str, header_text: &str) -> VcfStream {
        let create_index = super::flag(parser, "create-output-variant-index");
        let create_md5 = super::flag(parser, "create-output-variant-md5");
        let (level, deflater) = super::output_compression(parser);
        let special = std::fs::metadata(output).is_ok_and(|meta| !meta.is_file());
        let partial = if special {
            output.to_string()
        } else {
            format!("{output}.partial-{}", std::process::id())
        };
        let index = if !create_index {
            Index::None
        } else {
            match index_feature_file::index_kind(output) {
                index_feature_file::IndexKind::Tabix => {
                    let path = super::java_absolute_path(output);
                    match index_feature_file::codec_for(output) {
                        Some(codec) => Index::Tabix(Box::new(TabixStream::new(codec, &path))),
                        None => {
                            Index::Refused(index_feature_file::Refusal::NoSuitableCodecs { path })
                        }
                    }
                }
                _ => Index::Tribble(super::OnTheFlyIndex::new()),
            }
        };
        let mut stream = VcfStream {
            output: output.to_string(),
            partial,
            body: None,
            index,
            // The header's `##contig` lines, as `write_variant_file` reads them.
            dictionary: header_text
                .lines()
                .take_while(|line| line.starts_with('#'))
                .filter_map(|line| {
                    let body = line.strip_prefix("##contig=<")?.trim_end_matches('>');
                    let name = body.split(',').find_map(|f| f.strip_prefix("ID="))?;
                    let length = body
                        .split(',')
                        .find_map(|f| f.strip_prefix("length="))
                        .and_then(|value| value.parse().ok())
                        .unwrap_or(0);
                    Some((name.to_string(), length))
                })
                .collect(),
            level,
            deflater,
            failed: None,
        };
        match std::fs::File::create(&stream.partial) {
            Ok(file) => {
                let counted = Counted {
                    file: std::io::BufWriter::with_capacity(1 << 20, file),
                    written: 0,
                    md5: create_md5.then(Md5::new),
                };
                stream.body = Some(if output.ends_with(".gz") || output.ends_with(".bgz") {
                    Body::Block(Box::new(htsjdk_bgzf::BgzfWriter::with_deflater(
                        counted, level, deflater,
                    )))
                } else {
                    Body::Plain(counted)
                });
            }
            Err(error) => stream.failed = Some(error),
        }
        for line in header_text.split_inclusive('\n') {
            stream.write_line(line);
        }
        stream
    }

    /// One line of the file, its newline included.
    pub(crate) fn write_line(&mut self, line: &str) {
        if self.failed.is_some() {
            return;
        }
        let Some(body) = &mut self.body else { return };
        let written = match body {
            Body::Plain(counted) => {
                if let Index::Tribble(index) = &mut self.index {
                    index.add_line(line);
                }
                counted.write_all(line.as_bytes())
            }
            Body::Block(writer) => {
                if let Index::Tabix(index) = &mut self.index {
                    if let Some(feature) =
                        index_feature_file::line_feature(line, index_feature_file::Codec::Vcf)
                    {
                        match writer.file_pointer() {
                            Ok(offset) => index.add(feature, offset as i64),
                            Err(error) => {
                                self.failed = Some(error);
                                return;
                            }
                        }
                    }
                }
                writer.write_all(line.as_bytes())
            }
        };
        if let Err(error) = written {
            self.failed = Some(error);
        }
    }

    /// The run ends without the file: what was written goes.
    pub(crate) fn abandon(self) {
        drop(self.body);
        if self.partial != self.output {
            let _ = std::fs::remove_file(&self.partial);
        }
    }

    /// The file closed and put in place, then its md5, then its index: the order and the
    /// refusals of [`super::write_variant_file`].
    pub(crate) fn finish(mut self) -> Result<(), Thrown> {
        let output = self.output.clone();
        let could_not = |error: std::io::Error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
        };
        let (counted, end) = match self.body.take() {
            Some(Body::Plain(counted)) => (counted, None),
            Some(Body::Block(mut writer)) => {
                // The index ends where the last block ends, before the empty block that closes
                // the file.
                let end = writer
                    .flush()
                    .and_then(|()| writer.file_pointer())
                    .map_err(|error| self.failed.get_or_insert(error).kind());
                let counted = (*writer).into_inner();
                match (end, counted) {
                    (Ok(end), Ok(counted)) => (counted, Some(end as i64)),
                    (_, Err(error)) => {
                        self.failed.get_or_insert(error);
                        return Err(self.fail());
                    }
                    (Err(_), Ok(_)) => return Err(self.fail()),
                }
            }
            None => return Err(self.fail()),
        };
        if let Some(error) = self.failed.take() {
            self.abandon();
            return Err(could_not(error));
        }
        let Counted {
            mut file,
            written,
            md5,
        } = counted;
        if let Err(error) = file.flush() {
            self.abandon();
            return Err(could_not(error));
        }
        drop(file);
        if self.partial != output {
            std::fs::rename(&self.partial, &output).map_err(could_not)?;
        }
        if let Some(md5) = md5 {
            let digest = u128::from_be_bytes(md5.finalize().into());
            super::write_file(
                &format!("{output}.md5"),
                format!("{digest:032x}").as_bytes(),
            )?;
        }
        let companion = index_feature_file::default_output(&output);
        let index = match self.index {
            Index::None => return Ok(()),
            Index::Tribble(index) => index.finish(
                &self.dictionary,
                &super::java_absolute_path(&output),
                written as i64,
                super::modified_millis(&output),
            ),
            Index::Tabix(index) => index
                .finish(end.unwrap_or(0), self.deflater, self.level)
                .map_err(refused)?,
            Index::Refused(refusal) => return Err(refused(refusal)),
        };
        std::fs::write(&companion, index).map_err(|error| {
            Thrown::non_user(
                PORT_FAILURE,
                format!("could not write {companion}: {error}"),
            )
        })
    }

    fn fail(self) -> Thrown {
        let message = format!(
            "could not write {}: {}",
            self.output,
            self.failed
                .as_ref()
                .map_or_else(|| "no file".to_string(), ToString::to_string)
        );
        self.abandon();
        Thrown::non_user(PORT_FAILURE, message)
    }
}

/// An index refusal as the tool reports it.
fn refused(refusal: index_feature_file::Refusal) -> Thrown {
    Thrown {
        failure: Failure::User,
        exception: refusal.java_class(),
        message: Some(refusal.message()),
    }
}
