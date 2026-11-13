//! `RenameSampleInVcf.doWork`, over [`gatk_tools::rename_sample_in_vcf`].

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::rename_sample_in_vcf::{rename, RenameError};

use super::{assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let text = vcf_text(&input)?;
    let old = argument(parser, "OLD_SAMPLE_NAME");
    let written = rename(
        &text,
        &argument(parser, "NEW_SAMPLE_NAME").unwrap_or_default(),
        old.as_deref(),
    )
    .map_err(|error| match error {
        // `getGenotypeSamples().get(0)` of a file with no sample at all.
        RenameError::UnexpectedSample { contained }
            if contained.is_empty() && no_samples(&text) =>
        {
            thrown(
                "java.lang.IndexOutOfBoundsException",
                "Index 0 out of bounds for length 0".to_string(),
            )
        }
        other => thrown(other.java_class(), other.message()),
    })?;
    write_vcf(parser, &output, &written, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}

fn no_samples(text: &str) -> bool {
    text.lines()
        .find(|line| line.starts_with("#CHROM"))
        .is_none_or(|line| line.split('\t').count() <= 9)
}
