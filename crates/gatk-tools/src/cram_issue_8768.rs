//! `CRAMIssue8768Detector`, ported from the tool and `CRAMIssue8768Analyzer` (GATK 4.6.2.0).
//!
//! A CRAM walked container by container, looking for the containers that GATK issue 8768 would
//! have corrupted. Reading the containers is not ported; deciding which of them is a suspect, and
//! writing the report, are.
//!
//! # The suspect is the container AFTER the one that opens at position 1
//!
//! ```java
//! } else if (previousAlignmentContext.getReferenceContext().isMappedSingleRef() &&
//!         (previousAlignmentContext.getAlignmentStart() == 1)) {
//!     recordContainerStats(badContainersForContig, true, container, containerOrdinalForContig);
//! }
//! ```
//!
//! The test reads the PREVIOUS container's alignment start, so a contig whose first container
//! opens at 1 gets exactly one bad container, the second, whatever that one's own start is. The
//! container that opens at 1 is itself always reported as good, by the branch above.
//!
//! # The count printed beside the rate is the total base count
//!
//! ```java
//! return new Tuple<>(totalBases, misMatches/(double) totalBases);
//! ...
//! containerStats.a,   // mismatches
//! ```
//!
//! `.a` is `totalBases` and it is stored into a field named `misMatchCount`, so the report's
//! `Mismatch Rate/Count: 0.733333/30` pairs a rate with its own denominator. [`ContainerStats`]
//! keeps the reference's name for that field and says what it holds.
//!
//! # A bad contig is keyed by reference id, never by name
//!
//! The map key is `ReferenceContext.toString()`, which reads `SINGLE_REFERENCE: 0`. The text
//! report never resolves a contig name; only the TSV does, and only through the sequence
//! dictionary.
//!
//! # Four good containers per contig, and the counter starts at one
//!
//! `nGoodContainersReportedForContig = 1` in the new-context branch rather than 0, so the fifth
//! container of a contig is the first one dropped. `--verbose` is the only thing that shows it.
//!
//! # Multi-ref and unmapped containers are counted but never recorded
//!
//! `recordContainerStats` returns without adding a row unless the context is single-ref, so the
//! ordinal advances for a container that never appears, and the ordinals in the report skip. Those
//! containers still become the previous context, which is what flushes the contig before them.
//!
//! # The last contig is flushed by the EOF container, not by the loop
//!
//! Nothing after the loop emits `badContainersForContig`. It survives only because the EOF
//! container is read like any other and its context is UNMAPPED_UNPLACED, which is a change of
//! context for any mapped contig before it.
//!
//! # A foreign CRAM stops the analysis dead
//!
//! `isForeignCRAM` returns true from inside the container loop and `doAnalysis` RETURNS, so there
//! is no report body, no averages, and a return code of 0 as though nothing were wrong.

use gatk_engine::java_format::format_decimals;
use gatk_engine::tsv_table::{java_double_to_string, quote_if_needed};

/// `Slice.EMBEDDED_REFERENCE_ABSENT_CONTENT_ID`.
pub const EMBEDDED_REFERENCE_ABSENT_CONTENT_ID: i32 = -1;

/// `NUMBER_OF_GOOD_CONTAINERS_PER_CONTIG_TO_REPORT`.
pub const GOOD_CONTAINERS_PER_CONTIG: i32 = 4;

/// `CRAMIssue8768Detector.DEFAULT_MISMATCH_RATE_THRESHOLD`.
pub const DEFAULT_MISMATCH_RATE_THRESHOLD: f64 = 0.05;

/// `ReferenceContext`, reduced to what the report reads off it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefContext {
    Single(i32),
    Multiple,
    UnmappedUnplaced,
}

impl RefContext {
    /// `ReferenceContext.UNINITIALIZED_REFERENCE_ID`, which the good list starts at.
    pub const UNINITIALIZED_REFERENCE_ID: i32 = -3;

    /// `getReferenceContextID`. The reference throws for the two non-single cases; here they carry
    /// the ids the constants give them, because the good-container loop compares them by id.
    pub fn id(self) -> i32 {
        match self {
            RefContext::Single(id) => id,
            RefContext::Multiple => -2,
            RefContext::UnmappedUnplaced => -1,
        }
    }

    pub fn is_mapped_single_ref(self) -> bool {
        matches!(self, RefContext::Single(_))
    }

    /// `toString`, which is what the map key and the report line hold.
    pub fn text(self) -> String {
        match self {
            RefContext::Single(id) => format!("SINGLE_REFERENCE: {id}"),
            RefContext::Multiple => "MULTIPLE_REFERENCE".to_string(),
            RefContext::UnmappedUnplaced => "UNMAPPED_UNPLACED".to_string(),
        }
    }
}

/// One container, as the analyzer sees it: its alignment context, the shape checks read off it,
/// and the two numbers computed from its records.
#[derive(Debug, Clone, PartialEq)]
pub struct ContainerMeta {
    pub context: RefContext,
    pub start: i32,
    pub span: i32,
    pub slices: usize,
    /// `getCompressionHeader().isReferenceRequired()`. The EOF container carries no compression
    /// header at all, so this is `None` there, and the reference only ever asks for it when the
    /// context is single-ref, which the EOF container is not.
    pub reference_required: Option<bool>,
    pub embedded_reference: i32,
    /// The total base count of the container's records, which the report prints as a count.
    pub bases: i64,
    pub mismatches: i64,
    pub is_eof: bool,
}

/// `ContainerStats`, with the reference's field names.
#[derive(Debug, Clone, PartialEq)]
pub struct ContainerStats {
    pub container_ordinal: i32,
    pub is_bad: bool,
    pub context: RefContext,
    pub alignment_start: i32,
    pub alignment_span: i32,
    /// Named as the reference names it. It holds the TOTAL BASE COUNT: see the module header.
    pub mismatch_count: i64,
    pub mismatch_rate: f64,
}

impl ContainerStats {
    /// The report's line for one container, used for both lists.
    ///
    /// `AlignmentContext.toString` is `sequenceId=%s, start=%d, span=%d`, and the whole of it sits
    /// inside the parentheses.
    pub fn line(&self) -> String {
        format!(
            "  Ordinal: {} (sequenceId={}, start={}, span={}) Mismatch Rate/Count: {}/{}",
            self.container_ordinal,
            self.context.text(),
            self.alignment_start,
            self.alignment_span,
            format_decimals(self.mismatch_rate, 6),
            self.mismatch_count,
        )
    }
}

/// What one walk produced.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Analysis {
    /// The bad contigs in the order they were closed, which is a `LinkedHashMap`'s order. The key
    /// is `ReferenceContext.toString()`, not a contig name.
    pub bad_contigs: Vec<(String, Vec<ContainerStats>)>,
    pub good: Vec<ContainerStats>,
    /// The line `isForeignCRAM` emitted, if it stopped the walk.
    pub foreign: Option<String>,
}

/// `recordContigStats` refusing to close the same contig twice.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DuplicateBadContig {
    pub contig: String,
}

impl DuplicateBadContig {
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalStateException"
    }

    pub fn message(&self) -> String {
        format!(
            "Attempt to add a bad contig ({}) more than once",
            self.contig
        )
    }
}

/// `isForeignCRAM`, which emits its reason and stops the analysis.
///
/// The reference-less test is guarded by the context being single-ref, which is what keeps it away
/// from the EOF container's absent compression header.
pub fn foreign_message(meta: &ContainerMeta) -> Option<String> {
    if meta.slices > 1 {
        return Some(
            "Multi-slice container detected. This file was not written by GATK or Picard."
                .to_string(),
        );
    }
    if meta.context.is_mapped_single_ref()
        && !meta
            .reference_required
            .expect("a compression header on a mapped container")
    {
        return Some(
            "Reference-less container detected. This file was not written by GATK or Picard."
                .to_string(),
        );
    }
    if meta.slices > 0 && meta.embedded_reference != EMBEDDED_REFERENCE_ABSENT_CONTENT_ID {
        return Some(format!(
            "Embedded reference block (ID {}) detected. This file was not written by GATK or \
             Picard.",
            meta.embedded_reference
        ));
    }
    None
}

/// `doAnalysis`, over containers already read.
pub fn analyse(
    containers: &[ContainerMeta],
    verbose: bool,
) -> Result<Analysis, DuplicateBadContig> {
    let mut analysis = Analysis::default();
    let mut bad_for_contig: Vec<ContainerStats> = Vec::new();
    let mut ordinal = 0;
    let mut reported_good = 0;
    let mut previous: Option<(RefContext, i32)> = None;

    for meta in containers {
        ordinal += 1;
        if let Some(message) = foreign_message(meta) {
            analysis.foreign = Some(message);
            return Ok(analysis);
        }
        match previous {
            // The first container of the whole file cannot be bad.
            None => {
                record(&mut analysis.good, false, meta, ordinal);
                reported_good += 1;
            }
            Some((previous_context, previous_start)) => {
                if previous_context != meta.context {
                    if !bad_for_contig.is_empty() {
                        let key = previous_context.text();
                        if analysis.bad_contigs.iter().any(|(name, _)| name == &key) {
                            return Err(DuplicateBadContig { contig: key });
                        }
                        analysis
                            .bad_contigs
                            .push((key, std::mem::take(&mut bad_for_contig)));
                    }
                    ordinal = 1;
                    record(&mut analysis.good, false, meta, ordinal);
                    reported_good = 1;
                } else if previous_context.is_mapped_single_ref() && previous_start == 1 {
                    record(&mut bad_for_contig, true, meta, ordinal);
                } else if verbose || reported_good < GOOD_CONTAINERS_PER_CONTIG {
                    record(&mut analysis.good, false, meta, ordinal);
                    reported_good += 1;
                }
            }
        }
        previous = Some((meta.context, meta.start));
        if meta.is_eof {
            break;
        }
    }
    Ok(analysis)
}

/// `recordContainerStats`, which drops anything that is not single-ref.
fn record(into: &mut Vec<ContainerStats>, is_bad: bool, meta: &ContainerMeta, ordinal: i32) {
    if !meta.context.is_mapped_single_ref() {
        return;
    }
    into.push(ContainerStats {
        container_ordinal: ordinal,
        is_bad,
        context: meta.context,
        alignment_start: meta.start,
        alignment_span: meta.span,
        // `.a` of the tuple, which is the total base count.
        mismatch_count: meta.bases,
        mismatch_rate: meta.mismatches as f64 / meta.bases as f64,
    });
}

/// The three lines `analyzeCRAMHeader` writes before any container is read, which is why a foreign
/// CRAM still has them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CramHeaderInfo {
    pub file_name: String,
    pub version: String,
    pub id_base64: String,
}

/// The text file, what reached `System.out`, and what `doWork` returned.
#[derive(Debug, Clone, PartialEq)]
pub struct Report {
    pub text: String,
    pub stdout: String,
    pub code: i32,
}

/// `printTextResults`, with the header lines `doAnalysis` emitted before it.
pub fn report(
    header: &CramHeaderInfo,
    analysis: &Analysis,
    mismatch_rate_threshold: f64,
    echo_to_console: bool,
) -> Report {
    let mut text = String::new();
    let mut stdout = String::new();
    let emit = |line: &str, text: &mut String| {
        text.push_str(line);
        text.push('\n');
    };

    emit(&format!("CRAM File Name: {}", header.file_name), &mut text);
    emit(&format!("CRAM Version: {}", header.version), &mut text);
    emit(
        &format!("CRAM ID Contents: {}", header.id_base64),
        &mut text,
    );

    if let Some(message) = &analysis.foreign {
        emit(message, &mut text);
        return Report {
            text,
            stdout,
            code: 0,
        };
    }

    // Both banners begin with a newline, so the report carries a blank line before them, and the
    // corrupt one ends with a colon and a newline, so it carries one after.
    let (banner, code) = if analysis.bad_contigs.is_empty() {
        (
            "\n**********************NO CORRUPT CONTAINERS DETECTED**********************"
                .to_string(),
            0,
        )
    } else {
        (
            "\n**********************!!!!!Possible CORRUPT CONTAINERS DETECTED!!!!!\
             **********************:\n"
                .to_string(),
            1,
        )
    };
    emit(&banner, &mut text);
    // The banner reaches the console whether or not --echo-to-stdout was given.
    stdout.push_str(&banner);
    stdout.push('\n');

    // An average over no containers at all is 0.0/0, which prints as NaN rather than refusing.
    let total_good = analysis.good.len();
    let sum_good: f64 = analysis.good.iter().map(|stats| stats.mismatch_rate).sum();
    let average_good = sum_good / total_good as f64;
    let good_line = format!(
        "Average mismatch rate for presumed good containers: {}",
        format_decimals(average_good, 6)
    );
    emit(&good_line, &mut text);
    if echo_to_console {
        stdout.push_str(&good_line);
        stdout.push('\n');
    }

    if !analysis.bad_contigs.is_empty() {
        let total_bad: usize = analysis
            .bad_contigs
            .iter()
            .map(|(_, containers)| containers.len())
            .sum();
        let sum_bad: f64 = analysis
            .bad_contigs
            .iter()
            .map(|(_, containers)| {
                containers
                    .iter()
                    .map(|stats| stats.mismatch_rate)
                    .sum::<f64>()
            })
            .sum();
        let average_bad = sum_bad / total_bad as f64;
        let bad_line = format!(
            "Average mismatch rate for suspected bad containers: {}",
            format_decimals(average_bad, 6)
        );
        emit(&bad_line, &mut text);
        if echo_to_console {
            stdout.push_str(&bad_line);
            stdout.push('\n');
        }

        if average_bad > mismatch_rate_threshold {
            // The same kind of number twice in two formats: the measured rate in `%f`, the
            // threshold in `%1.2f`.
            let exceeded = format!(
                "The average base mismatch rate of {} for suspected bad containers exceeds the \
                 threshold rate of {}, and indicates this file may be corrupt.",
                format_decimals(average_bad, 6),
                format_decimals(mismatch_rate_threshold, 2),
            );
            emit(&exceeded, &mut text);
            if echo_to_console {
                stdout.push_str(&exceeded);
                stdout.push('\n');
            }
        }

        // The two headings are written to the file only, never echoed.
        emit("\nSuspected CORRUPT Containers:", &mut text);
        for (_, containers) in &analysis.bad_contigs {
            for stats in containers {
                let line = stats.line();
                emit(&line, &mut text);
                if echo_to_console {
                    stdout.push_str(&line);
                    stdout.push('\n');
                }
            }
        }
    }

    emit("\nPresumed GOOD Containers:", &mut text);
    let mut last_contig = RefContext::UNINITIALIZED_REFERENCE_ID;
    for stats in &analysis.good {
        if last_contig != RefContext::UNINITIALIZED_REFERENCE_ID
            && last_contig != stats.context.id()
        {
            emit("", &mut text);
            if echo_to_console {
                stdout.push('\n');
            }
        }
        last_contig = stats.context.id();
        let line = stats.line();
        emit(&line, &mut text);
        if echo_to_console {
            stdout.push_str(&line);
            stdout.push('\n');
        }
    }

    Report { text, stdout, code }
}

/// `printTSVResults`: the same rows, and the only place a contig NAME is ever resolved.
pub fn tsv(analysis: &Analysis, file_name: &str, dictionary: &[String]) -> String {
    let mut out = String::new();
    out.push_str(
        "file_name\tcontig_name\tcontainer_ordinal\tcontainer_is_bad\tmismatch_rate\t\
         alignment_start\talignment_span\n",
    );
    let row = |stats: &ContainerStats, out: &mut String| {
        let values = [
            file_name.to_string(),
            dictionary[stats.context.id() as usize].clone(),
            stats.container_ordinal.to_string(),
            if stats.is_bad { "1" } else { "0" }.to_string(),
            java_double_to_string(stats.mismatch_rate),
            stats.alignment_start.to_string(),
            stats.alignment_span.to_string(),
        ];
        let quoted: Vec<String> = values.iter().map(|value| quote_if_needed(value)).collect();
        out.push_str(&quoted.join("\t"));
        out.push('\n');
    };

    if analysis.bad_contigs.is_empty() {
        out.push_str("#No bad containers detected\n");
    } else {
        out.push_str("#Bad containers:\n");
        for (_, containers) in &analysis.bad_contigs {
            for stats in containers {
                row(stats, &mut out);
            }
        }
    }
    if analysis.good.is_empty() {
        out.push_str("#No good mapped containers detected\n");
    } else {
        out.push_str("#Good containers:\n");
        for stats in &analysis.good {
            row(stats, &mut out);
        }
    }
    out
}
