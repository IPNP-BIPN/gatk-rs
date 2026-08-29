//! Which read filter owns which controlled argument, and which of them are required.
//!
//! The trim in [`gatk_barclay`] is the rule: a controlled argument whose plugin nobody selected
//! and nobody set leaves the definition list before the required check, so a required argument of
//! an unselected read filter does not fire. The rule needs a table, and this is it.
//!
//! The table cannot be derived from the declarations. [`crate::tool_declarations`] records that an
//! argument is controlled by `GATKReadFilterPluginDescriptor`, which is the DESCRIPTOR and not the
//! filter that contributed it, and the names do not name their filters: `--library` is
//! `LibraryReadFilter`'s, `--platform-filter-name` is `PlatformReadFilter`'s,
//! `--black-listed-lanes` is `PlatformUnitReadFilter`'s. So it is measured, in the
//! `plugin-argument-ownership` golden, and transcribed here in the golden's order.
//!
//! Twelve of the twenty-eight are required, which is what makes the trim load-bearing: without it
//! a plain walker command line would be asked for `--keep-intervals`, `--library`, `--sample`,
//! `--read-name` and eight more that GATK never asks for.
//!
//! Ported from `org.broadinstitute.hellbender.engine.filters` through
//! `org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor`.

/// The long name of the argument that names read filters, which is what the trim asks about.
///
/// `isDependentArgumentAllowed` is the descriptor's, and every descriptor in GATK answers it by
/// looking at its own selector: a plugin is selected when this argument's values name it, or when
/// the tool handed the descriptor an instance of it as a default.
pub const SELECTOR: &str = "read-filter";

/// One controlled argument, and the filter that declared it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ownership {
    /// The filter class's simple name, which is what a refusal names and what `--read-filter`
    /// takes.
    pub owner: &'static str,
    /// The argument's long name.
    pub long_name: &'static str,
    /// Its short name, empty where the golden records a bare dash.
    pub short_name: &'static str,
    /// Whether the reference reports it as required.
    pub required: bool,
}

/// The table, in the golden's order, which is by owner and then by long name.
pub const OWNERSHIP: [Ownership; 28] = [
    o("AmbiguousBaseReadFilter", "ambig-filter-bases", "", false),
    o("AmbiguousBaseReadFilter", "ambig-filter-frac", "", false),
    o(
        "ExcessiveEndClippedReadFilter",
        "max-clipped-bases",
        "",
        false,
    ),
    o(
        "FlowBasedTPAttributeValidReadFilter",
        "read-filter-max-hmer",
        "",
        false,
    ),
    o("FragmentLengthReadFilter", "max-fragment-length", "", false),
    o("FragmentLengthReadFilter", "min-fragment-length", "", false),
    o("IntervalOverlapReadFilter", "keep-intervals", "", true),
    o(
        "JexlExpressionReadTagValueFilter",
        "read-filter-expression",
        "",
        true,
    ),
    o("LibraryReadFilter", "library", "library", true),
    o(
        "MappingQualityReadFilter",
        "maximum-mapping-quality",
        "",
        false,
    ),
    o(
        "MappingQualityReadFilter",
        "minimum-mapping-quality",
        "",
        false,
    ),
    o(
        "MateDistantReadFilter",
        "mate-too-distant-length",
        "",
        false,
    ),
    o(
        "OverclippedReadFilter",
        "dont-require-soft-clips-both-ends",
        "",
        false,
    ),
    o("OverclippedReadFilter", "filter-too-short", "", false),
    o("PlatformReadFilter", "platform-filter-name", "", true),
    o("PlatformUnitReadFilter", "black-listed-lanes", "", true),
    o(
        "ReadGroupBlackListReadFilter",
        "read-group-black-list",
        "",
        true,
    ),
    o("ReadGroupReadFilter", "keep-read-group", "", true),
    o("ReadLengthReadFilter", "max-read-length", "", true),
    o("ReadLengthReadFilter", "min-read-length", "", false),
    o("ReadNameReadFilter", "read-name", "", true),
    o("ReadStrandFilter", "keep-reverse-strand-only", "", true),
    o("ReadTagValueFilter", "read-filter-tag", "", true),
    o("ReadTagValueFilter", "read-filter-tag-comp", "", false),
    o("ReadTagValueFilter", "read-filter-tag-op", "", false),
    o("SampleReadFilter", "sample", "sample", true),
    o(
        "SoftClippedReadFilter",
        "max-soft-clipped-leading-trailing-ratio",
        "",
        false,
    ),
    o("SoftClippedReadFilter", "max-soft-clipped-ratio", "", false),
];

const fn o(
    owner: &'static str,
    long_name: &'static str,
    short_name: &'static str,
    required: bool,
) -> Ownership {
    Ownership {
        owner,
        long_name,
        short_name,
        required,
    }
}

/// The filter that declared an argument, by long name.
pub fn owner(long_name: &str) -> Option<&'static str> {
    OWNERSHIP
        .iter()
        .find(|entry| entry.long_name == long_name)
        .map(|entry| entry.owner)
}

/// The row for an argument, by long name.
pub fn ownership(long_name: &str) -> Option<&'static Ownership> {
    OWNERSHIP.iter().find(|entry| entry.long_name == long_name)
}
