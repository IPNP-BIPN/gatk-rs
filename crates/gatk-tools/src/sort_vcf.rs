//! `SortVcf`, ported from `picard.vcf.SortVcf` (Picard 3.4.0).
//!
//! One or more VCFs read into a sorting collection and written back out in dictionary order under a
//! header merged from all of them.
//!
//! # The order is the dictionary's, and ties keep their input order
//!
//! The comparator is contig index then position, and nothing else. Two records at the same locus
//! therefore keep the order they were added in, which across several inputs is the order the inputs
//! were given: the same pair of files handed over the other way round writes the same records in a
//! different order. The golden's `two-files` and `two-files-reversed` are that pair.
//!
//! A record on a contig the dictionary does not declare does NOT come back with a message. htsjdk's
//! comparator looks the contig up in a map and unboxes what it finds, so the run dies with a null
//! pointer; this port carries that as a refusal of its own rather than sorting the record somewhere
//! arbitrary.
//!
//! # The samples are checked by their sorted names
//!
//! ```java
//! if (!sampleList.equals(header.getSampleNamesInOrder())) { throw ... }
//! ```
//!
//! `getSampleNamesInOrder` is SORTED, so two files listing the same samples in different column
//! orders agree and are not refused. The output's columns are those sorted names, and each input's
//! genotypes are placed by name rather than by position.
//!
//! # The dictionary is the first input's
//!
//! The first file's contig lines become the run's, and any later file declaring a different
//! dictionary is a refusal. A file with no dictionary at all is refused unless one was supplied
//! separately.

use htsjdk_vcf::comparator::VariantContextComparator;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::vcf_file::write_vcf;

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SortError {
    /// An input with no contig lines and no dictionary supplied. The path is masked, as the golden
    /// masks it: the message names wherever the run happened.
    MissingDictionary,
    /// A later input whose dictionary is not the first's.
    DifferentDictionaries,
    /// A later input whose sorted sample names are not the first's.
    DifferentSamples,
    /// A record on a contig the dictionary does not declare, which htsjdk meets as a null pointer.
    UndeclaredContig,
    /// The reader, the merge or the writer refused.
    Vcf(String, String),
}

impl SortError {
    pub fn java_class(&self) -> &str {
        match self {
            SortError::MissingDictionary
            | SortError::DifferentDictionaries
            | SortError::DifferentSamples => "java.lang.IllegalArgumentException",
            SortError::UndeclaredContig => "java.lang.NullPointerException",
            SortError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            SortError::MissingDictionary => "Sequence dictionary was missing or empty for the \
                 VCF: <masked> Please add a sequence dictionary to this VCF or specify \
                 SEQUENCE_DICTIONARY."
                .to_string(),
            SortError::DifferentDictionaries => {
                "java.lang.AssertionError: SAM dictionaries are not the same".to_string()
            }
            SortError::DifferentSamples => {
                "Input file <masked> has sample names that don't match the other files.".to_string()
            }
            SortError::UndeclaredContig => "Cannot invoke \"java.lang.Integer.intValue()\" \
                 because the return value of \"java.util.Map.get(Object)\" is null"
                .to_string(),
            SortError::Vcf(_, message) => message.clone(),
        }
    }
}

/// The contig lines of a header, as the dictionary check compares them.
fn contig_lines(header: &VcfHeader) -> Vec<String> {
    header
        .lines
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .map(|line| line.render())
        .collect()
}

/// `getSampleNamesInOrder()`, which is sorted.
fn sample_names_in_order(header: &VcfHeader) -> Vec<String> {
    let mut names = header.samples.clone();
    names.sort();
    names
}

/// `doWork()`: the whole run, every input's text in and one text out.
pub fn sort(inputs: &[String]) -> Result<String, SortError> {
    let mut headers: Vec<VcfHeader> = Vec::new();
    let mut records: Vec<VariantContext> = Vec::new();
    let mut samples: Option<Vec<String>> = None;
    let mut dictionary: Option<Vec<String>> = None;

    for input in inputs {
        let file = read_vcf(input).map_err(|failure| {
            SortError::Vcf(failure.error.class().to_string(), failure.error.message())
        })?;
        let contigs = contig_lines(&file.header);
        if contigs.is_empty() {
            return Err(SortError::MissingDictionary);
        }
        match &dictionary {
            None => dictionary = Some(contigs),
            Some(first) if *first != contigs => return Err(SortError::DifferentDictionaries),
            Some(_) => {}
        }
        let names = sample_names_in_order(&file.header);
        match &samples {
            None => samples = Some(names),
            Some(first) if *first != names => return Err(SortError::DifferentSamples),
            Some(_) => {}
        }
        records.extend(file.records.iter().cloned());
        headers.push(file.header);
    }

    let sources: Vec<Source> = headers
        .iter()
        .map(|header| Source {
            header,
            version: None,
        })
        .collect();
    let (merged, _warnings) = smart_merge_headers(&sources, false)
        .map_err(|error| SortError::Vcf(error.class().to_string(), error.message().to_string()))?;

    let mut header = VcfHeader {
        lines: merged,
        samples: samples.unwrap_or_default(),
    };
    // `new VCFHeader(lines, sampleList)`, whose sample list is what the checks agreed on.
    header.samples.sort();

    // `from_header_lines` is given the contig lines alone: it counts its lookup against the list
    // it was handed, so a whole header would look like a header full of duplicates.
    let contigs: Vec<HeaderLine> = header
        .lines
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .cloned()
        .collect();
    let comparator = VariantContextComparator::from_header_lines(&contigs)
        .map_err(|error| SortError::Vcf(error.class().to_string(), error.message().to_string()))?;
    // Every contig is looked up once before the sort: htsjdk meets the missing one inside the
    // comparator, where it is a null pointer, and a sort whose comparator can fail has nowhere to
    // put the failure.
    for record in &records {
        if comparator.compare(record, record).is_err() {
            return Err(SortError::UndeclaredContig);
        }
    }
    // A stable sort, which is what keeps two records at one locus in the order they were read.
    records.sort_by(|left, right| {
        comparator
            .compare(left, right)
            .expect("every contig was checked")
            .cmp(&0)
    });

    write_vcf(&header, &records).map_err(|error| {
        SortError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{error:?}"),
        )
    })
}
