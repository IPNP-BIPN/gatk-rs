//! `SelectVariants.createVCFHeaderLineList` and the header it hands the writer.
//!
//! The five existing `select-variants-*` suites all skip the `#` lines of the output, and for a
//! good reason: the default header carries a `##GATKCommandLine` line holding the wall-clock time
//! of the run, so a golden of the whole header would flake by construction. That is why this is
//! the part of the tool nothing measured, and it is not a small part -- the method merges, adds,
//! REPLACES and removes lines, and the replacement is the half a port gets wrong.
//!
//! Four rules, in the order the reference applies them.
//!
//! 1. The input's own header lines, through `VCFUtils.smartMergeHeaders` over a set of one.
//! 2. The tool's default lines: `##source=SelectVariants` and `##GATKCommandLine`, and BOTH of
//!    them are gone when `--add-output-vcf-command-line` is false, because
//!    `getDefaultToolVCFHeaderLines` returns an empty set rather than dropping one line.
//! 3. `--keep-original-ac` adds three lines and `--keep-original-dp` adds one, whether or not any
//!    record ends up carrying them.
//! 4. AN, AC, AF and then DP are REMOVED and re-added from htsjdk's standard definitions. This is
//!    not a merge: an input whose `##INFO=<ID=AC>` says `Description="Not the standard count"`
//!    comes out with htsjdk's description, and AF appears in a file that never had one. The order
//!    is `ChromosomeCounts.keyNames`, which is AN, AC, AF rather than the alphabetical reading.
//!
//! The drops run LAST, after the four replacements, so `--drop-info-annotation AC` removes the
//! standard line that had just been put there.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants` of
//! GATK 4.6.2.0.

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use htsjdk_vcf::standard_header_lines::standard_info_line;

/// `GATKVCFConstants`: the four keys `--keep-original-ac` and `--keep-original-dp` add.
pub const ORIGINAL_AC_KEY: &str = "AC_Orig";
pub const ORIGINAL_AF_KEY: &str = "AF_Orig";
pub const ORIGINAL_AN_KEY: &str = "AN_Orig";
pub const ORIGINAL_DP_KEY: &str = "DP_Orig";

/// `ChromosomeCounts.keyNames`, in its own order: allele NUMBER first, then count, then frequency.
pub const CHROMOSOME_COUNT_KEYS: [&str; 3] = ["AN", "AC", "AF"];

/// What the header construction reads off the command line.
pub struct HeaderArguments {
    pub keep_original_chr_counts: bool,
    pub keep_original_depth: bool,
    pub info_annotations_to_drop: Vec<String>,
    pub genotype_annotations_to_drop: Vec<String>,
    /// `--add-output-vcf-command-line`. False removes the `##source` line as well.
    pub add_output_vcf_command_line: bool,
    /// The `##GATKCommandLine` value, which the caller builds because it holds the run's own time.
    pub tool_command_line: Option<HeaderLine>,
    /// The samples the output declares, which are the SELECTED ones rather than the input's.
    pub samples: Vec<String>,
}

/// `GATKVCFHeaderLines.getInfoLine`, for the four keys this tool can add.
fn gatk_info_line(id: &str) -> HeaderLine {
    let (number, line_type, description) = match id {
        ORIGINAL_AC_KEY => (Cardinality::A, LineType::Integer, "Original AC"),
        ORIGINAL_AF_KEY => (Cardinality::A, LineType::Float, "Original AF"),
        ORIGINAL_AN_KEY => (Cardinality::Fixed(1), LineType::Integer, "Original AN"),
        ORIGINAL_DP_KEY => (Cardinality::Fixed(1), LineType::Integer, "Original DP"),
        _ => unreachable!("no GATK header line for {id}"),
    };
    HeaderLine::Compound {
        key: "INFO".to_string(),
        id: id.to_string(),
        number,
        line_type,
        description: description.to_string(),
        extra: Vec::new(),
    }
}

/// True when the line is an `##INFO` (or `##FORMAT`) line with this id.
fn is_compound(line: &HeaderLine, key: &str, id: &str) -> bool {
    matches!(
        line,
        HeaderLine::Compound { key: k, id: i, .. } if k == key && i == id
    )
}

/// `createVCFHeaderLineList`, then `new VCFHeader(lines, samples)`.
pub fn output_header(input: &VcfHeader, arguments: &HeaderArguments) -> VcfHeader {
    // One source, which is what `Collections.singletonMap` gives `smartMergeHeaders`. The merge
    // over a single header is not the identity: it is what repairs a line whose declaration
    // disagrees with the standard one for its id.
    let (mut lines, _warnings) = smart_merge_headers(
        &[Source {
            header: input,
            version: None,
        }],
        true,
    )
    .unwrap_or_else(|_| (input.lines.clone(), Vec::new()));

    if arguments.add_output_vcf_command_line {
        // `##source=<tool>` and the command line, and the pair is all-or-nothing.
        lines.push(HeaderLine::Unstructured {
            key: "source".to_string(),
            value: "SelectVariants".to_string(),
        });
        if let Some(command_line) = &arguments.tool_command_line {
            lines.push(command_line.clone());
        }
    }

    if arguments.keep_original_chr_counts {
        for key in [ORIGINAL_AC_KEY, ORIGINAL_AF_KEY, ORIGINAL_AN_KEY] {
            lines.push(gatk_info_line(key));
        }
    }
    if arguments.keep_original_depth {
        lines.push(gatk_info_line(ORIGINAL_DP_KEY));
    }

    // The replacement, which is a removal and an addition rather than a merge.
    for key in CHROMOSOME_COUNT_KEYS {
        lines.retain(|line| !is_compound(line, "INFO", key));
        if let Some(standard) = standard_info_line(key) {
            lines.push(standard);
        }
    }
    lines.retain(|line| !is_compound(line, "INFO", "DP"));
    if let Some(standard) = standard_info_line("DP") {
        lines.push(standard);
    }

    // Last, so that dropping AC removes the line the loop above had just added.
    lines.retain(|line| {
        !arguments
            .info_annotations_to_drop
            .iter()
            .any(|id| is_compound(line, "INFO", id))
    });
    lines.retain(|line| {
        !arguments
            .genotype_annotations_to_drop
            .iter()
            .any(|id| is_compound(line, "FORMAT", id))
    });

    VcfHeader {
        lines,
        samples: arguments.samples.clone(),
    }
}
