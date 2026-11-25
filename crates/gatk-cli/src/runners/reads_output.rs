//! Reads outputs written as `GATKTool.createSAMWriter` writes them: SAM, BAM or CRAM by the
//! output's name.
//!
//! The walkers here build their output as a finished BAM in memory, so this module takes that
//! BAM and writes what htsjdk 5.0.0's `SAMFileWriterFactory.makeWriter` would have written
//! instead: the BAM as it is, its records as SAM text for a `.sam` name, or as htsjdk's CRAM 3.1
//! for a `.cram` name. Re-reading the BAM costs a decode, and keeps every tool's own writer as it
//! is; the header and records the BAM carries are the ones `createSAMWriter` was handed.
//!
//! # What `createSAMWriter` checks, in its order
//!
//! * a `.cram` output without `-R` is `MISSING_REFERENCE_FOR_CRAM`, before any file exists;
//! * `ReadUtils.getCRAMVersion` reads `--output-cram-version` for EVERY output, so a value that
//!   is neither `3.0` nor `3.1` refuses a BAM run too;
//! * `createCommonSAMWriter` drops the index of a header that is not `SO:coordinate`.
//!
//! # Where the companions go
//!
//! * a BAM's index: `BAMFileWriter.createBamIndex` strips a trailing `.bam` and appends `.bai`,
//!   so `out.bam` is indexed by `out.bai` but `out.tmp` by `out.tmp.bai`;
//! * a CRAM's index: `IOUtil.addExtension(output, ".bai")`, `out.cram.bai`, built by
//!   `CRAMBAIIndexer` from the containers;
//! * a SAM text file has no index of any kind;
//! * the digest, for all three: `Md5CalculatingOutputStream` over the bytes written, in
//!   `<output>.md5`.

use super::*;

/// `MISSING_REFERENCE_FOR_CRAM`: `MissingReference(message, true)` appends where `-R` goes.
const MISSING_REFERENCE_FOR_CRAM: &str = "A reference file is required when using CRAM files.  \
     A reference is specified with the -R command line argument.";

/// The CRAM minor version `ReadUtils.getCRAMVersion` returns for `--output-cram-version`, or its
/// `BadInput` refusal. The default is `GATKConfig.outputCramVersion()`, which is `3.1`.
pub(super) fn cram_minor_version(parser: &Parser) -> Result<u8, Thrown> {
    let version = argument(parser, "output-cram-version").unwrap_or_else(|| "3.1".to_string());
    match version.trim() {
        "3.0" => Ok(0),
        "3.1" => Ok(1),
        _ => Err(Thrown::user(format!(
            "Bad input: Unsupported CRAM version '{version}'. Supported values are '3.0' and '3.1'."
        ))),
    }
}

/// The index `BAMFileWriter` writes beside a BAM named `output`.
pub(super) fn bam_index_path(output: &str) -> String {
    format!("{}.bai", output.strip_suffix(".bam").unwrap_or(output))
}

/// `Path.toUri().toString()`, the id `SAMFileWriterFactory` gives a CRAM it writes to a file.
fn file_id(output: &str) -> String {
    funcotator::java_file_uri(output)
}

/// The header and records of a finished BAM, which this port's writers have just encoded.
fn decode_bam(bam: &[u8]) -> Result<(SamHeader, Vec<BamRecord>), Thrown> {
    let failed = |error: String| Thrown::non_user(PORT_FAILURE, error);
    let plain = htsjdk_bgzf::read::decompress_all(bam).map_err(|e| failed(format!("{e:?}")))?;
    let reader =
        htsjdk_bam::reader::BamReader::new(&plain).map_err(|e| failed(format!("{e:?}")))?;
    let header = reader.header.text.clone();
    let records = reader
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| failed(format!("{e:?}")))?;
    Ok((header, records))
}

pub(crate) fn write_file(path: &str, bytes: &[u8]) -> Result<(), Thrown> {
    std::fs::write(path, bytes)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("could not write {path}: {error}")))
}

/// The writer options every CRAM here is written under: `minor_version` from
/// `--output-cram-version`, and the GZIP blocks deflated by `java.util.zip.Deflater`, which is what
/// htsjdk 5.0.0's `new DeflaterFactory()` returns under `-Dsamjdk.use_libdeflate=false`, the
/// property the pinned oracle runs every tool under (`default_java_props` in the manifest).
/// `--use-jdk-deflater` does not reach it: GATK's flag only swaps the BGZF factory. The SAM header
/// block is at `samjdk.compression_level`, which GATK's configuration sets to 2.
pub(crate) fn cram_options(minor_version: u8) -> htsjdk_cram::writer::WriterOptions {
    let level = std::env::var(gatk_tools::gatk_config::COMPRESSION_LEVEL).ok();
    htsjdk_cram::writer::WriterOptions {
        minor_version,
        gzip: htsjdk_cram::gzip::Deflater::Jdk,
        header_gzip_level: gatk_tools::gatk_config::compression_level(level.as_deref()) as i32,
    }
}

/// `header` and `records` as htsjdk's `CRAMFileWriter` writes them to `output` against the
/// reference `reference`, with the `.bai` `CRAMBAIIndexer` builds when `index` is set. A writer
/// that refuses leaves both files there and empty, as the streams it opened first do in Java.
///
/// With no reference, which only a Picard tool reaches, the writer has htsjdk's
/// `CRAMLazyReferenceSource`, which refuses the first mapped record in its own words.
pub(crate) fn write_cram(
    output: &str,
    reference: Option<&str>,
    header: &SamHeader,
    records: &[BamRecord],
    index: bool,
    minor_version: u8,
) -> Result<Vec<u8>, Thrown> {
    let index_path = format!("{output}.bai");
    let mut source: Box<dyn htsjdk_cram::reader::ReferenceSource> = match reference {
        Some(path) => Box::new(htsjdk_cram::reader::FastaReference::open(
            std::path::Path::new(path),
        )),
        None => Box::new(htsjdk_cram::reader::LazyReference),
    };
    let written = htsjdk_cram::writer::write_cram_with(
        header,
        records,
        &file_id(output),
        source.as_mut(),
        cram_options(minor_version),
        index,
    );
    let (cram, bai) = match written {
        Ok(written) => written,
        Err(error) => {
            if index {
                write_file(&index_path, &[])?;
            }
            write_file(output, &[])?;
            return Err(Thrown::non_user(error.java_exception(), error.message()));
        }
    };
    write_file(output, &cram)?;
    if let Some(bai) = bai {
        write_file(&index_path, &bai)?;
    }
    Ok(cram)
}

/// A finished BAM written to `output` as `createSAMWriter`'s writer would have written its header
/// and records, with the index and digest beside it.
pub(super) fn write_reads(
    parser: &Parser,
    output: &str,
    bam: &[u8],
    bai: Option<Vec<u8>>,
) -> Result<Option<String>, Thrown> {
    let is_cram = output.ends_with(".cram");
    let reference = argument(parser, "reference");
    if is_cram && reference.is_none() {
        return Err(Thrown::user(MISSING_REFERENCE_FOR_CRAM.to_string()));
    }
    let minor_version = cram_minor_version(parser)?;
    // htsjdk's writer builds an index only for a `SO:coordinate` file, whatever
    // `--create-output-bam-index` says: an unsorted or queryname file is written without one.
    let bai = bai.filter(|_| is_coordinate_sorted(bam));

    let written: std::borrow::Cow<[u8]> = if let (true, Some(reference)) = (is_cram, &reference) {
        let (header, records) = decode_bam(bam)?;
        write_cram(
            output,
            Some(reference),
            &header,
            &records,
            bai.is_some(),
            minor_version,
        )?
        .into()
    } else if output.ends_with(".sam") {
        let (header, records) = decode_bam(bam)?;
        let text = htsjdk_bam::sam_file::write_sam(&header, &records).ok_or_else(|| {
            Thrown::non_user(PORT_FAILURE, format!("could not write {output} as SAM"))
        })?;
        write_file(output, text.as_bytes())?;
        text.into_bytes().into()
    } else {
        write_file(output, bam)?;
        if let Some(bai) = &bai {
            write_file(&bam_index_path(output), bai)?;
        }
        bam.into()
    };
    if flag(parser, "create-output-bam-md5") {
        let digest = format!("{output}.md5");
        write_file(
            &digest,
            gatk_tools::gather_bam_files::md5_file(&written).as_bytes(),
        )?;
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bam_index_replaces_only_a_bam_extension() {
        assert_eq!(bam_index_path("out/x.bam"), "out/x.bai");
        assert_eq!(bam_index_path("out/x.tmp"), "out/x.tmp.bai");
        assert_eq!(bam_index_path("out/x"), "out/x.bai");
    }
}
