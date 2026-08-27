//! `StructuralVariantDiscoverer`: how the split alignments of an assembled contig become a
//! structural variant.
//!
//! A contig that aligns in two pieces is read as an adjacency, and which variant it is comes from
//! the SIGNATURE of the pair rather than from any argument.
//!
//! Reading the BAM and writing the VCF are not ported, and neither is the complex-event path. The
//! sort requirement, the read filters and the two simple signatures are.

/// One alignment of an assembled contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub contig_name: String,
    pub reference: String,
    pub start: i32,
    pub end: i32,
    pub reverse_strand: bool,
    pub supplementary: bool,
    pub secondary: bool,
    pub unmapped: bool,
    pub mapping_quality: i32,
}

/// The tool's own read filters, applied before any signature is read.
pub fn passes_read_filters(alignment: &Alignment) -> bool {
    !alignment.unmapped && !alignment.secondary
}

/// The sort order the input must be in.
///
/// The tool gathers a contig's alignments by walking CONSECUTIVE records of the same name, so a
/// coordinate-sorted file would scatter them and it refuses one outright.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Queryname,
    Coordinate,
    Unsorted,
}

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiscovererError {
    NotQuerynameSorted,
}

impl DiscovererError {
    pub fn message(&self) -> String {
        match self {
            DiscovererError::NotQuerynameSorted => {
                "This tool requires a queryname-sorted source of reads.".to_string()
            }
        }
    }
}

pub fn check_sort_order(order: SortOrder) -> Result<(), DiscovererError> {
    if order == SortOrder::Queryname {
        Ok(())
    } else {
        Err(DiscovererError::NotQuerynameSorted)
    }
}

/// What one contig's pair of alignments is read as.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Signature {
    /// The pieces leave a gap on the reference.
    Deletion { start: i32, end: i32 },
    /// The pieces overlap on the reference: the overlap is the repeat unit.
    TandemDuplication { start: i32, end: i32 },
    /// Anything else, including a pair whose strands differ, which produces no call at all.
    None,
}

/// The signature of one contig's alignments.
///
/// A contig with anything but two surviving alignments produces nothing. A pair whose strands
/// DIFFER produces nothing either: a strand flip alone is not an inversion signature, which is the
/// measurement rather than an omission.
pub fn signature(alignments: &[Alignment]) -> Signature {
    let kept: Vec<&Alignment> = alignments
        .iter()
        .filter(|alignment| passes_read_filters(alignment))
        .collect();
    if kept.len() != 2 {
        return Signature::None;
    }
    let (first, second) = (kept[0], kept[1]);
    if first.reference != second.reference {
        return Signature::None;
    }
    if first.reverse_strand != second.reverse_strand {
        return Signature::None;
    }
    let (left, right) = if first.start <= second.start {
        (first, second)
    } else {
        (second, first)
    };
    if right.start > left.end + 1 {
        // The gap between them, whose ends are the last base of the left piece and the base before
        // the first of the right one.
        Signature::Deletion {
            start: left.end,
            end: right.start - 1,
        }
    } else if right.start <= left.end {
        Signature::TandemDuplication {
            start: right.start - 1,
            end: left.end,
        }
    } else {
        Signature::None
    }
}

/// One called variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    /// `<DEL>` or `<DUP>`.
    pub alternate: String,
    pub id: String,
    /// `CTG_NAMES`, so each call says which contig made it.
    pub contig_names: Vec<String>,
}

/// The calls one queryname-sorted file produces, one record per adjacency.
pub fn discover(alignments: &[Alignment]) -> Vec<Variant> {
    let mut out = Vec::new();
    let mut index = 0;
    while index < alignments.len() {
        let name = &alignments[index].contig_name;
        let mut end = index;
        while end < alignments.len() && alignments[end].contig_name == *name {
            end += 1;
        }
        let group = &alignments[index..end];
        index = end;
        let reference = group[0].reference.clone();
        match signature(group) {
            Signature::Deletion { start, end } => out.push(Variant {
                contig: reference.clone(),
                start,
                end,
                alternate: "<DEL>".to_string(),
                id: format!("DEL_{reference}_{start}_{end}"),
                contig_names: vec![name.clone()],
            }),
            Signature::TandemDuplication { start, end } => out.push(Variant {
                contig: reference.clone(),
                start,
                end,
                alternate: "<DUP>".to_string(),
                id: format!(
                    "INS-DUPLICATION-TANDEM-EXPANSION_{reference}_{}_{end}",
                    start + 1
                ),
                contig_names: vec![name.clone()],
            }),
            Signature::None => {}
        }
    }
    out.sort_by_key(|variant| variant.start);
    out
}
