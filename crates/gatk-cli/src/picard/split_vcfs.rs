//! `SplitVcfs.doWork`, over [`gatk_tools::split_vcfs`].
//!
//! `SEQUENCE_DICTIONARY` only replaces the dictionary the two indexes are built against, which is
//! why a VCF with no contig lines is refused without it only when indexing.
//!
//! Picard streams: each record goes to its writer, and its writer's on-the-fly index, before the
//! next is read, so an index refusing an earlier record comes before `STRICT` refuses a later one.
//! A run nothing refuses is streamed the same way ([`stream`]); the rest is decided here.

mod stream;

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};
use gatk_tools::split_vcfs::{split_stream, RecordType, Split, SplitError};

use super::vcf_format_converter::{bcf_limitation, is_bcf, output_type};
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
    if is_bcf(&input) {
        return Err(bcf_limitation());
    }
    if stream::split(parser, &input, &snp_output, &indel_output) {
        return Ok(0);
    }
    exact(parser, &input, &snp_output, &indel_output)
}

/// The run with the whole input in hand, which is what decides every refusal.
fn exact(
    parser: &Parser,
    input: &str,
    snp_output: &str,
    indel_output: &str,
) -> Result<i32, Thrown> {
    let text = vcf_text(input)?;
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
    output_type(snp_output)?;
    output_type(indel_output)?;
    // The builder hands the TabixIndexCreator it made for a `.vcf.gz` SNP file to the text INDEL
    // writer it builds next, so two files' positions reach one index.
    let shared = create_index && block(snp_output) && !block(indel_output);
    if create_index {
        index_refusal(parser, snp_output, indel_output, &written)?;
    }
    if let Some(kind) = written.unexpected {
        return Err(refuse(SplitError::UnexpectedType(kind)));
    }
    for (output, text) in [(snp_output, &written.snps), (indel_output, &written.indels)] {
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
        file_positions(parser, &written.snps, snp_block)?,
        file_positions(parser, &written.indels, indel_block)?,
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
/// a text file; a block-compressed file's BGZF virtual offset, which at close is the start of the
/// empty block that ends it, because the index is finalised before that block counts.
fn file_positions(
    parser: &Parser,
    text: &str,
    block_compressed: bool,
) -> Result<Positions, Thrown> {
    if !block_compressed {
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

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    /// An empty directory of the test's own.
    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("split-vcfs-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// SNPs, indels and a MIXED record over two contigs, one sample, long enough to fill several
    /// BGZF blocks so the tabix offsets cross block boundaries.
    fn input() -> String {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##contig=<ID=chr1,length=1000000>\n\
             ##contig=<ID=chr2,length=1000000>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
        );
        for contig in ["chr1", "chr2"] {
            for at in 1..4000 {
                let (reference, alternate) = match at % 5 {
                    0 => ("A", "AT"),
                    1 => ("AC", "A"),
                    2 => ("A", "C,AT"),
                    _ => ("G", "T"),
                };
                text.push_str(&format!(
                    "{contig}\t{}\trs{at}\t{reference}\t{alternate}\t50\tPASS\tDP={at}\tGT\t0/1\n",
                    at * 10
                ));
            }
        }
        text
    }

    /// Every file a run left in `dir`, by name, with an `.idx`'s timestamp zeroed: the two runs
    /// write their files at different moments.
    fn files(dir: &Path) -> Vec<(String, Vec<u8>)> {
        let mut found: Vec<(String, Vec<u8>)> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| {
                let path = entry.unwrap().path();
                let name = path.file_name().unwrap().to_string_lossy().into_owned();
                let mut bytes = std::fs::read(&path).unwrap();
                if name.ends_with(".idx") {
                    // magic, type, version, the indexed path and its NUL, then the size and the
                    // timestamp, eight bytes each.
                    let nul = 12 + bytes[12..].iter().position(|&b| b == 0).unwrap();
                    bytes[nul + 9..nul + 17].fill(0);
                }
                (name, bytes)
            })
            .collect();
        found.sort();
        found
    }

    /// The streamed run writes what the exact path writes, byte for byte, for every pairing of
    /// a text and a block-compressed output, indexed or not.
    #[test]
    fn the_stream_writes_what_the_exact_path_writes() {
        let scratch = scratch("same-bytes");
        let source = scratch.join("in.vcf");
        std::fs::write(&source, input()).unwrap();
        let compressed = scratch.join("in.vcf.gz");
        let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
        std::io::Write::write_all(&mut writer, input().as_bytes()).unwrap();
        std::fs::write(&compressed, writer.into_inner().unwrap()).unwrap();

        for input in [&source, &compressed] {
            for (snp, indel) in [("vcf.gz", "vcf.gz"), ("vcf", "vcf"), ("vcf", "vcf.gz")] {
                for index in ["true", "false"] {
                    let mut answers = Vec::new();
                    for streamed in [true, false] {
                        let out = scratch.join(format!("out-{snp}-{indel}-{index}"));
                        let _ = std::fs::remove_dir_all(&out);
                        std::fs::create_dir(&out).unwrap();
                        let (snp_output, indel_output) = (
                            out.join(format!("snp.{snp}"))
                                .to_string_lossy()
                                .into_owned(),
                            out.join(format!("indel.{indel}"))
                                .to_string_lossy()
                                .into_owned(),
                        );
                        let arguments: Vec<String> = [
                            "--INPUT",
                            &input.to_string_lossy(),
                            "--SNP_OUTPUT",
                            &snp_output,
                            "--INDEL_OUTPUT",
                            &indel_output,
                            "--STRICT",
                            "false",
                            "--CREATE_INDEX",
                            index,
                        ]
                        .map(String::from)
                        .to_vec();
                        let parser = super::super::parse("SplitVcfs", &arguments).unwrap();
                        let input = input.to_string_lossy();
                        if streamed {
                            assert!(super::stream::split(
                                &parser,
                                &input,
                                &snp_output,
                                &indel_output
                            ));
                        } else {
                            super::exact(&parser, &input, &snp_output, &indel_output).unwrap();
                        }
                        answers.push(files(&out));
                    }
                    assert_eq!(answers[0], answers[1], "{input:?} {snp} {indel} {index}");
                    assert!(answers[0].len() >= 2);
                }
            }
        }
    }

    /// A record STRICT refuses leaves the run to the exact path, and nothing of the stream behind.
    #[test]
    fn a_refusal_is_left_to_the_exact_path() {
        let scratch = scratch("refusal");
        let source = scratch.join("in.vcf");
        std::fs::write(&source, input()).unwrap();
        let snp_output = scratch.join("snp.vcf").to_string_lossy().into_owned();
        let indel_output = scratch.join("indel.vcf").to_string_lossy().into_owned();
        let input = source.to_string_lossy().into_owned();
        let arguments: Vec<String> = [
            "--INPUT",
            &input,
            "--SNP_OUTPUT",
            &snp_output,
            "--INDEL_OUTPUT",
            &indel_output,
        ]
        .map(String::from)
        .to_vec();
        let parser = super::super::parse("SplitVcfs", &arguments).unwrap();
        assert!(!super::stream::split(
            &parser,
            &input,
            &snp_output,
            &indel_output
        ));
        let left: Vec<String> = std::fs::read_dir(scratch)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        assert_eq!(left, ["in.vcf"]);
    }
}
