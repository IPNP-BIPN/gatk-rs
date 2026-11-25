//! `SplitVcfs.doWork`, over [`gatk_tools::split_vcfs`].
//!
//! `SEQUENCE_DICTIONARY` only replaces the dictionary the two indexes are built against, which is
//! why a VCF with no contig lines is refused without it only when indexing.
//!
//! Picard streams: each record goes to its writer, and its writer's on-the-fly index, before the
//! next is read, so an index refusing an earlier record comes before `STRICT` refuses a later one.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};
use gatk_tools::split_vcfs::{split_stream, RecordType, Split, SplitError};

use super::vcf_format_converter::output_type;
use super::{
    assert_readable, record_features, thrown, vcf_text, virtual_offsets, write_vcf, TabixChecks,
};
use crate::runners::{argument, dictionary_from_any, flag, output_compression, write_variant_file};

fn block(path: &str) -> bool {
    path.ends_with(".vcf.gz") || path.ends_with(".vcf.bgz")
}

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let snp_output = argument(parser, "SNP_OUTPUT").unwrap_or_default();
    let indel_output = argument(parser, "INDEL_OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
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
    let refuse = |error: SplitError| thrown(error.java_class(), error.message());
    // The library's own dictionary check is the header's; a given dictionary answers it instead.
    // The reader and the dictionary refuse before the two writers are built and refuse a name.
    let written = split_stream(
        &text,
        flag(parser, "STRICT"),
        create_index && given.is_none(),
    )
    .map_err(refuse)?;
    output_type(&snp_output)?;
    output_type(&indel_output)?;
    // The builder hands the TabixIndexCreator it made for a `.vcf.gz` SNP file to the text INDEL
    // writer it builds next, so two files' positions reach one index.
    let shared = create_index && block(&snp_output) && !block(&indel_output);
    if create_index {
        index_refusal(parser, &snp_output, &indel_output, &written)?;
    }
    if let Some(kind) = written.unexpected {
        return Err(refuse(SplitError::UnexpectedType(kind)));
    }
    for (output, text) in [
        (&snp_output, &written.snps),
        (&indel_output, &written.indels),
    ] {
        if block(output) {
            // A tabix index takes its contigs in the order they come, whatever dictionary it has.
            write_vcf(parser, output, text, create_index)?;
        } else {
            write_variant_file(
                parser,
                output,
                text,
                create_index && !shared,
                false,
                given.clone(),
            )?;
        }
    }
    if shared {
        // Only a run with no record at all closes a shared index twice without a refusal, and an
        // index of nothing is the same bytes whichever file it is written beside.
        let (from, to) = (format!("{snp_output}.tbi"), format!("{indel_output}.tbi"));
        std::fs::copy(&from, &to).map_err(|error| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {to}: {error}"))
        })?;
    }
    Ok(0)
}

/// The on-the-fly indexes' refusals, fed as the two writers feed them: each record at its own
/// file's position, in input order, then `close()` on the SNP writer and on the INDEL writer.
/// A Tribble index (a text file's own) refuses no order, so only the tabix creators are followed.
fn index_refusal(
    parser: &Parser,
    snp_output: &str,
    indel_output: &str,
    written: &Split,
) -> Result<(), Thrown> {
    let (snp_block, indel_block) = (block(snp_output), block(indel_output));
    // The creator each writer feeds: the SNP file's own, the INDEL file's own, or the SNP file's.
    let creators = [
        snp_block.then_some(0),
        if indel_block {
            Some(1)
        } else {
            snp_block.then_some(0)
        },
    ];
    let mut checks = [TabixChecks::default(), TabixChecks::default()];
    let positions = [
        file_positions(parser, &written.snps, snp_output)?,
        file_positions(parser, &written.indels, indel_output)?,
    ];
    let mut features = [
        record_features(&written.snps),
        record_features(&written.indels),
    ];
    for kind in &written.routed {
        let writer = usize::from(*kind == RecordType::Indel);
        let Some((at, contig, start, end)) = features[writer].next() else {
            break;
        };
        if let Some(creator) = creators[writer] {
            checks[creator].add(contig, start, end, (positions[writer].0)(at))?;
        }
    }
    // `snpWriter.close()`, then `indelWriter.close()`.
    for writer in 0..2 {
        if let Some(creator) = creators[writer] {
            checks[creator].finalize(positions[writer].1)?;
        }
    }
    Ok(())
}

type Positions = (Box<dyn Fn(usize) -> i64>, i64);

/// Where a writer stands at each offset of its file's text, and when it closes: the byte count of
/// a text file; a BCF's offset of the record written from that line; a block-compressed file's
/// BGZF virtual offset, which at close is the start of the empty block that ends it, because the
/// index is finalised before that block counts.
fn file_positions(parser: &Parser, text: &str, output: &str) -> Result<Positions, Thrown> {
    if crate::bcf::is_bcf(output) {
        let (starts, end) = crate::bcf::positions(text)?;
        let at: std::collections::HashMap<usize, i64> = record_features(text)
            .map(|(offset, ..)| offset)
            .zip(starts)
            .collect();
        return Ok((
            Box::new(move |offset| at.get(&offset).copied().unwrap_or(end)),
            end,
        ));
    }
    if !block(output) {
        return Ok((Box::new(|offset| offset as i64), text.len() as i64));
    }
    let (level, deflater) = output_compression(parser);
    let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
    std::io::Write::write_all(&mut writer, text.as_bytes())
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
    let compressed = writer
        .into_inner()
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error}")))?;
    let end = compressed
        .len()
        .saturating_sub(htsjdk_bgzf::EMPTY_GZIP_BLOCK.len());
    Ok((Box::new(virtual_offsets(&compressed)), (end as i64) << 16))
}
