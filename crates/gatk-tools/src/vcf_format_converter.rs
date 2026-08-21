//! `VcfFormatConverter`, ported from `picard.vcf.VcfFormatConverter` (Picard 3.4.0).
//!
//! A file rewritten in the format its name asks for. The port covers the text formats; the binary
//! one is a codec of its own and the suite records it as a digest.
//!
//! # Two defaults that both say `true` and mean different things
//!
//! `REQUIRE_INDEX` is a declared argument and defaults to true, so a plain vcf with no index beside
//! it is refused before a record is read, and refused BY TRIBBLE: the message is the reader's and
//! names the file it looked for. `CREATE_INDEX` is set in the constructor rather than by the
//! parser, and it is the reason a file with no contig lines is refused, the dictionary being null
//! and the `PicardException` naming the indexing rather than the file.
//!
//! The two are independent: dropping the requirement does not stop the indexing from refusing a
//! file with no dictionary, and turning the indexing off accepts it.
//!
//! # A conversion is a rewrite, not a copy
//!
//! Every record goes through the decoder and the encoder, and the header is copied through
//! `new VCFHeader(header)` and emitted in the writer's own order. A header written out of order
//! comes back sorted whatever the format asked for.

use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::vcf_file::write_vcf;

/// What the run refuses, neither of which is raised by the tool's own code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConvertError {
    /// The reader wanted an index and found none. The path is the input's, as the message quotes
    /// it.
    MissingIndex { path: String },
    /// `CREATE_INDEX` with no dictionary to index against.
    NoSequenceDictionary,
    /// The reader or the writer refused.
    Vcf(String, String),
}

impl ConvertError {
    pub fn java_class(&self) -> &str {
        match self {
            ConvertError::MissingIndex { .. } => "htsjdk.tribble.TribbleException",
            ConvertError::NoSequenceDictionary => "picard.PicardException",
            ConvertError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            ConvertError::MissingIndex { path } => format!(
                "An index is required, but none found with file ending .idx, for input source: \
                 file://{path}"
            ),
            ConvertError::NoSequenceDictionary => {
                "A sequence dictionary must be available in the input file when creating indexed \
                 output."
                    .to_string()
            }
            ConvertError::Vcf(_, message) => message.clone(),
        }
    }
}

/// The arguments, with the defaults the tool starts from: both of them true, one from the parser
/// and one from the constructor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arguments {
    pub require_index: bool,
    pub create_index: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            require_index: true,
            create_index: true,
        }
    }
}

/// The input, with the two facts about it the run reads from the filesystem.
#[derive(Debug, Clone, Copy)]
pub struct Input<'a> {
    /// The path the reader's message quotes.
    pub path: &'a str,
    pub text: &'a str,
    /// Whether an index sits beside it.
    pub indexed: bool,
}

/// `doWork()` for a text output: the whole run, text in and text out.
pub fn convert(input: &Input, arguments: &Arguments) -> Result<String, ConvertError> {
    // `new VCFFileReader(INPUT, REQUIRE_INDEX)`, which looks for the index before anything else.
    if arguments.require_index && !input.indexed {
        return Err(ConvertError::MissingIndex {
            path: input.path.to_string(),
        });
    }
    let file = read_vcf(input.text).map_err(|failure| {
        ConvertError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }));
    if arguments.create_index && !has_dictionary {
        return Err(ConvertError::NoSequenceDictionary);
    }
    write_vcf(&file.header, &file.records).map_err(|error| {
        ConvertError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{error:?}"),
        )
    })
}
