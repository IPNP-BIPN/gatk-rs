//! `FilterAlignmentArtifacts`: each call's supporting reads reassembled into unitigs, the unitigs
//! realigned with BWA-MEM, and the call filtered when they map elsewhere as well.
//!
//! Only the walker around that is ported. The reassembly is the Mutect2 read-threading assembler,
//! the realignment is BWA-MEM's C through JNI, and neither exists in this port; what does is the
//! rule that decides which reads go into them, [`gatk_tools::filter_alignment_artifacts`]'s
//! `supportsVariant`. So a call that SOME read supports is the port's limitation, refused by name,
//! and a call no read supports takes the reference's own path for it: no read reaches the
//! assembler, the pileup over none gives no unitig, BWA is never asked, and the call is written
//! with `UNITIGS` empty and `NALIGNS=0`.
//!
//! The order is the reference's:
//!
//! * **startup**: the reads and their dictionaries validated as a variant walker validates them,
//!   the reference required, and the BWA image required to be readable before the traversal;
//! * **the traversal**: the calls in input order, a filtered one skipped unless
//!   `--dont-skip-filtered-variants` is given, and each one's supporting reads looked for among
//!   the reads that cover its start;
//! * **the output**: the input's header with the tool's FILTER and three INFO lines, every call
//!   with its two attributes, and the writer's standard arguments.
use super::*;
use gatk_tools::filter_alignment_artifacts as artifacts;

/// A read as `supportsVariant` sees it, from the BAM record.
fn support_read(record: &htsjdk_bam::record::BamRecord) -> artifacts::Read {
    let text = record.cigar.to_text();
    let mut cigar = Vec::new();
    let mut length = 0i32;
    for character in text.chars() {
        if let Some(digit) = character.to_digit(10) {
            length = length * 10 + digit as i32;
        } else {
            cigar.push((character, length));
            length = 0;
        }
    }
    artifacts::Read {
        name: record.read_name.clone(),
        start: record.alignment_start,
        cigar,
        bases: record.read_bases.clone(),
    }
}

pub fn filter_alignment_artifacts(parser: &Parser) -> Outcome {
    let inputs = arguments(parser, "variant");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --variant is a GATK feature that this port does not carry yet. This \
             message is the port's own and not GATK's.",
        ));
    }
    let resolved = resolve_read_filters(parser, "FilterAlignmentArtifacts")?;
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;
    if arguments(parser, "input").is_empty() {
        return Err(Thrown::command_line(
            "Argument input was missing: Argument 'input' is required",
        ));
    }
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup_over(parser, input)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    if argument(parser, "reference").is_none() {
        return Err(Thrown::non_user(
            "java.lang.NullPointerException",
            "Cannot invoke \"org.broadinstitute.hellbender.engine.GATKPath.toPath()\" because \
             \"referencePathSpecifier\" is null",
        ));
    }
    // `new BwaMemIndex(image)` opens the image at startup.
    let image = argument(parser, "bwa-mem-index-image").ok_or_else(|| {
        Thrown::command_line(
            "Argument bwa-mem-index-image was missing: Argument 'bwa-mem-index-image' is required",
        )
    })?;
    if std::fs::metadata(&image)
        .map(|meta| !meta.is_file())
        .unwrap_or(true)
    {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            format!(
                "--bwa-mem-index-image {image} cannot be opened, and the message GATK's native \
                 loader prints for that is not one this port reproduces. This message is the \
                 port's own and not GATK's."
            ),
        ));
    }
    if argument(parser, "bam-output").is_some() {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--bam-output writes the haplotypes the Mutect2 assembler builds, which this port does \
             not carry yet. This message is the port's own and not GATK's.",
        ));
    }

    // The reads, every one that could cover a call, filtered as the tool's filters say.
    let mut reads: Vec<(String, artifacts::Read)> = Vec::new();
    for path in arguments(parser, "input") {
        let source =
            gatk_engine::reads::ReadsDataSource::open_unindexed(std::path::Path::new(&path))
                .map_err(|error| Thrown::user(format!("{error:?}")))?;
        let sequences = source.header().sequences.clone();
        let header = source.header().clone();
        let filter = read_filter(parser, &resolved, &header)?;
        for read in gatk_tools::read_walker::traverse(&source, &[], &filter)
            .map_err(|error| Thrown::user(format!("{error:?}")))?
        {
            let contig = sequences
                .get(read.reference_index as usize)
                .map(|sequence| sequence.name.clone())
                .unwrap_or_default();
            reads.push((contig, support_read(&read)));
        }
    }

    let mut file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: failure.error.class(),
        message: Some(failure.error.message()),
    })?;
    let traversed: Vec<htsjdk_vcf::variant::VariantContext> =
        variants_in_traversal(&file.records, intervals.as_deref(), &input)?
            .into_iter()
            .cloned()
            .collect();
    let skip_filtered = !flag(parser, "dont-skip-filtered-variants");
    let tolerance = number_or(parser, "indel-start-tolerance", 5);
    let mut written = Vec::new();
    for record in traversed {
        if skip_filtered && record.filters.as_ref().is_some_and(|f| !f.is_empty()) {
            continue;
        }
        let variant = artifacts::Variant {
            start: record.start as i32,
            reference: record.reference().base_string().into_bytes(),
            alternates: record
                .alternate_alleles()
                .iter()
                .map(|allele| allele.display_string().into_bytes())
                .collect(),
        };
        let supported = reads.iter().any(|(contig, read)| {
            *contig == record.contig && artifacts::supports_variant(read, &variant, tolerance)
        });
        if supported {
            return Err(Thrown::non_user(
                PORT_LIMITATION,
                format!(
                    "Reads support the call at {}:{}, and realigning them needs the Mutect2 \
                     assembler and BWA-MEM, which this port does not carry yet. This message is \
                     the port's own and not GATK's.",
                    record.contig, record.start
                ),
            ));
        }
        let mut out = record.clone();
        out.attributes
            .push(("NALIGNS".to_string(), htsjdk_vcf::variant::Value::Int(0)));
        out.attributes
            .push(("UNITIGS".to_string(), htsjdk_vcf::variant::Value::Missing));
        written.push(out);
    }

    // The header: the input's lines, the tool's own, and the defaults.
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
    let compound = |id: &str, number: Cardinality, description: &str| HeaderLine::Compound {
        key: "INFO".to_string(),
        id: id.to_string(),
        number,
        line_type: LineType::Integer,
        description: description.to_string(),
        extra: Vec::new(),
    };
    let added = [
        HeaderLine::Filter {
            id: "alignment".to_string(),
            description: "Alignment artifact".to_string(),
        },
        compound(
            "UNITIGS",
            Cardinality::Unbounded,
            "Sizes of reassembled unitigs",
        ),
        compound(
            "ALIGN_DIFF",
            Cardinality::Fixed(1),
            "Difference in alignment score between best and next-best alignment",
        ),
        compound(
            "NALIGNS",
            Cardinality::Fixed(1),
            "Number of joint alignments",
        ),
    ];
    let same = |a: &HeaderLine, b: &HeaderLine| a.render() == b.render();
    for line in added {
        if !file
            .header
            .lines
            .iter()
            .any(|existing| same(existing, &line))
        {
            file.header.lines.push(line);
        }
    }
    file.header.lines.extend(default_tool_vcf_header_lines(
        parser,
        "FilterAlignmentArtifacts",
    ));

    let keep = variant_output_filter(parser, intervals.as_deref())?;
    written.retain(|record| keep(record));
    apply_sites_only(parser, &mut file.header, &mut written);
    let rendered =
        htsjdk_vcf::vcf_file::write_vcf(&file.header, &written).map_err(|error| Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException",
            message: Some(format!("{error:?}")),
        })?;
    write_variant_output(parser, &output, &rendered)?;
    // `onTraversalSuccess` returns "SUCCESS", which `Main` prints.
    Ok(Some("SUCCESS".to_string()))
}
