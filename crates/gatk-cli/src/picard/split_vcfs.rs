//! `SplitVcfs.doWork`, over [`gatk_tools::split_vcfs`].
//!
//! `SEQUENCE_DICTIONARY` only replaces the dictionary the two indexes are built against, which is
//! why a VCF with no contig lines is refused without it only when indexing.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_LIMITATION};
use gatk_tools::split_vcfs::{split, SplitError};

use super::vcf_format_converter::{bcf_limitation, is_bcf, output_type};
use super::{assert_readable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, dictionary_from_any, flag, write_variant_file};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let snp_output = argument(parser, "SNP_OUTPUT").unwrap_or_default();
    let indel_output = argument(parser, "INDEL_OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    if is_bcf(&input) {
        return Err(bcf_limitation());
    }
    let text = vcf_text(&input)?;
    let create_index = flag(parser, "CREATE_INDEX");
    // Read through `SamReaderFactory`, as a SAM header.
    let given: Option<Vec<(String, i32)>> = match argument(parser, "SEQUENCE_DICTIONARY") {
        None => None,
        Some(path) => Some(
            dictionary_from_any(&path)?
                .into_iter()
                .map(|record| (record.name.clone(), record.length))
                .collect(),
        ),
    };
    // The library's own dictionary check is the header's; a given dictionary answers it instead.
    let result = split(
        &text,
        flag(parser, "STRICT"),
        create_index && given.is_none(),
    );
    let refuse = |error: SplitError| thrown(error.java_class(), error.message());
    // The reader and the dictionary refuse before the two writers are built and refuse a name;
    // a record of the wrong type only once they are.
    let written = match result {
        Err(SplitError::UnexpectedType(kind)) => Err(SplitError::UnexpectedType(kind)),
        other => Ok(other.map_err(refuse)?),
    };
    output_type(&snp_output)?;
    output_type(&indel_output)?;
    let block = |path: &str| path.ends_with(".vcf.gz") || path.ends_with(".vcf.bgz");
    if create_index && block(&snp_output) && !block(&indel_output) {
        // The builder hands the TabixIndexCreator it made for the SNP file to the text INDEL
        // writer it builds next, so two files' positions reach one index.
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "An indexed .vcf.gz SNP_OUTPUT beside a text INDEL_OUTPUT shares one tabix index in htsjdk, which this port does not reproduce. This message is the port's own and not GATK's.",
        ));
    }
    let written = written.map_err(refuse)?;
    for (output, text) in [
        (&snp_output, &written.snps),
        (&indel_output, &written.indels),
    ] {
        if block(output) {
            // A tabix index takes its contigs in the order they come, whatever dictionary it has.
            write_vcf(parser, output, text, flag(parser, "CREATE_INDEX"))?;
        } else {
            write_variant_file(parser, output, text, create_index, false, given.clone())?;
        }
    }
    Ok(0)
}
