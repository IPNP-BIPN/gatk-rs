//! `GetSampleName`, ported from `org.broadinstitute.hellbender.tools.GetSampleName` (GATK 4.6.2.0).
//!
//! The smallest tool of the reporting-walker archetype: its `traverse` is empty and everything
//! happens in `onTraversalStart`, which reads the header's read groups and writes their sample names
//! to a file. It reads no records at all, so a BAM holding nothing but a header still produces its
//! sample.
//!
//! # The file has no trailing newline
//!
//! ```java
//! sampleNames.stream().map(...).collect(Collectors.joining("\n"))
//! ```
//!
//! `joining` puts a separator **between** names and nothing after the last, so a one-sample file is
//! seven bytes for `sample1` rather than eight. The length is the only way to see it, which is why
//! the suite compares both the text and the byte count.
//!
//! # Two refusals, of which only the second is reachable from a BAM
//!
//! The first guard asks whether the header or its read group list is null; htsjdk returns an
//! **empty list** for a header with no `@RG` lines, so that branch never fires and the message is
//! always "The given bam input has no sample names.".
//!
//! A read group with no `SM` is not the same thing as no sample: `getSample()` returns null,
//! `distinct()` keeps the null, and `Collectors.joining` writes the four letters `null`. The tool
//! finishes and the file is not empty.

use gatk_engine::reads::ReadsDataSource;
use htsjdk_bam::header::SamHeader;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK GetSampleName";

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GetSampleNameError {
    /// The header, or its read group list, was null. Unreachable through a BAM.
    NoHeaderOrReadGroups,
    /// No read groups at all, which is what an `@RG`-less header actually produces.
    NoSampleNames,
}

impl GetSampleNameError {
    /// The message `UserException.BadInput` carries, without its "Bad input: " prefix.
    pub fn message(&self) -> String {
        match self {
            GetSampleNameError::NoHeaderOrReadGroups => {
                "The given input bam has no header or no read groups.  Cannot determine a sample name."
                    .to_string()
            }
            GetSampleNameError::NoSampleNames => {
                "The given bam input has no sample names.".to_string()
            }
        }
    }
}

/// The distinct sample names, in read group order.
///
/// `distinct()` keeps the **first** occurrence of each, so the order is the header's and not
/// alphabetical, and a read group with no `SM` contributes a `None` that survives to the output.
pub fn sample_names(header: &SamHeader) -> Vec<Option<String>> {
    let mut seen: Vec<Option<String>> = Vec::new();
    for group in &header.read_groups {
        let sample = group.attributes.get("SM").map(|s| s.to_string());
        if !seen.contains(&sample) {
            seen.push(sample);
        }
    }
    seen
}

/// `GetSampleName.onTraversalStart`: the file's whole content.
///
/// Returns the text as written, with no trailing newline.
pub fn get_sample_name(
    source: &ReadsDataSource,
    url_encode: bool,
) -> Result<String, GetSampleNameError> {
    sample_name_text(source.header(), url_encode)
}

/// The same from a header alone, which is all this tool ever looks at.
pub fn sample_name_text(
    header: &SamHeader,
    url_encode: bool,
) -> Result<String, GetSampleNameError> {
    let samples = sample_names(header);
    if samples.is_empty() {
        return Err(GetSampleNameError::NoSampleNames);
    }
    Ok(samples
        .iter()
        .map(|sample| match sample {
            // `String.valueOf(null)` inside `Collectors.joining`.
            None => "null".to_string(),
            Some(name) if url_encode => url_encode_utf8(name),
            Some(name) => name.clone(),
        })
        .collect::<Vec<String>>()
        .join("\n"))
}

/// `IOUtils.urlEncode`, which is `java.net.URLEncoder.encode(text, UTF-8)`.
///
/// The form-encoding rules, not the path ones: a **space becomes `+`**, `*`, `-`, `.` and `_` pass
/// through beside the alphanumerics, and everything else becomes `%XX` of its UTF-8 bytes with
/// upper-case hexadecimal.
pub fn url_encode_utf8(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for byte in text.bytes() {
        match byte {
            b'a'..=b'z' | b'A'..=b'Z' | b'0'..=b'9' | b'.' | b'-' | b'*' | b'_' => {
                out.push(byte as char)
            }
            b' ' => out.push('+'),
            other => {
                out.push('%');
                out.push_str(&format!("{other:02X}"));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::ReadGroup;

    fn header(groups: &[(&str, Option<&str>)]) -> SamHeader {
        let mut header = SamHeader::default();
        for (id, sample) in groups {
            let mut group = ReadGroup::new(id);
            if let Some(sample) = sample {
                group.attributes.set("SM", sample);
            }
            header.read_groups.push(group);
        }
        header
    }

    #[test]
    fn one_sample_has_no_trailing_newline() {
        let text = sample_name_text(&header(&[("rg1", Some("sample1"))]), false).expect("a sample");
        assert_eq!(text, "sample1");
        assert_eq!(text.len(), 7);
    }

    #[test]
    fn the_order_is_the_headers_and_repeats_collapse() {
        let text = sample_name_text(
            &header(&[("rg1", Some("zebra")), ("rg2", Some("alpha"))]),
            false,
        )
        .expect("two samples");
        assert_eq!(text, "zebra\nalpha");
        assert_eq!(text.len(), 11);

        let repeated =
            sample_name_text(&header(&[("rg1", Some("s1")), ("rg2", Some("s1"))]), false)
                .expect("one sample");
        assert_eq!(repeated, "s1");
    }

    #[test]
    fn a_read_group_with_no_sample_writes_the_word_null() {
        let text = sample_name_text(&header(&[("rg1", None)]), false).expect("it finishes");
        assert_eq!(text, "null");
    }

    #[test]
    fn no_read_groups_is_the_second_refusal() {
        assert_eq!(
            sample_name_text(&header(&[]), false).unwrap_err(),
            GetSampleNameError::NoSampleNames
        );
    }

    #[test]
    fn a_space_becomes_a_plus_and_not_a_percent_twenty() {
        assert_eq!(
            url_encode_utf8("a sample/with+odd chars & more"),
            "a+sample%2Fwith%2Bodd+chars+%26+more"
        );
    }
}
