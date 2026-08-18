//! Every `Mutect2Filter`'s identity and the header lines that describe them, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering` and `GATKVCFConstants`
//! (GATK 4.6.2.0).
//!
//! # Three error types, and they are not decoration
//!
//! `ErrorProbabilities` combines the filters of one type by their **maximum** and the types with
//! each other as **independent**. A filter's type therefore decides which other filters can mask it:
//! `SEQUENCING` holds the tumour-evidence filter alone, `NON_SOMATIC` holds contamination and
//! germline, and `ARTIFACT` holds the other sixteen.
//!
//! # Nine per allele and nine per site
//!
//! A per-site filter's one probability is copied to every alternate allele; a per-allele filter
//! answers a list. The two are the `Mutect2VariantFilter` and `Mutect2AlleleFilter` subclasses.
//!
//! # The header's list is not the engine's list
//!
//! `MUTECT_FILTER_NAMES` includes `PASS` and `FAIL`, which no filter reports, and every entry
//! becomes a `##FILTER` line whether or not the filter runs. `MUTECT_AS_FILTER_NAMES` is one entry
//! and becomes an `##INFO` line instead.
//!
//! # Which filters a run builds
//!
//! `buildFiltersList` guards six of them with `if (!MTFAC.mitochondria && !MTFAC.microbial)`, adds
//! one back under `if (MTFAC.microbial)`, and one more under an argument. Neither the engine nor its
//! stats file will say which were built: the stats file names only the filters that **fired**. See
//! [`build_filters_list`], which is measured against the `filter-list-by-mode` golden.

pub use crate::error_probabilities::ErrorType;

/// Whether a filter answers one probability per allele or one for the site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arity {
    PerAllele,
    PerSite,
}

impl Arity {
    pub fn name(self) -> &'static str {
        match self {
            Arity::PerAllele => "per-allele",
            Arity::PerSite => "per-site",
        }
    }
}

/// One filter's identity: its Java class, the name it reports, its error type and its arity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FilterIdentity {
    pub class: &'static str,
    pub filter_name: &'static str,
    pub error_type: ErrorType,
    pub arity: Arity,
}

/// Every filter `buildFiltersList` can construct, in the order it constructs them.
///
/// All eighteen names are distinct: `StrandArtifactFilter` reports `strand_bias` and
/// `StrictStrandBiasFilter` reports `strict_strand`, which are not the same name.
pub const FILTERS: [FilterIdentity; 18] = [
    FilterIdentity {
        class: "TumorEvidenceFilter",
        filter_name: "weak_evidence",
        error_type: ErrorType::Sequencing,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "BaseQualityFilter",
        filter_name: "base_qual",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "MappingQualityFilter",
        filter_name: "map_qual",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "DuplicatedAltReadFilter",
        filter_name: "duplicate",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "StrandArtifactFilter",
        filter_name: "strand_bias",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "ContaminationFilter",
        filter_name: "contamination",
        error_type: ErrorType::NonSomatic,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "StrictStrandBiasFilter",
        filter_name: "strict_strand",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "ReadPositionFilter",
        filter_name: "position",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "MinAlleleFractionFilter",
        filter_name: "low_allele_frac",
        error_type: ErrorType::Artifact,
        arity: Arity::PerAllele,
    },
    FilterIdentity {
        class: "NormalArtifactFilter",
        filter_name: "normal_artifact",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "NRatioFilter",
        filter_name: "n_ratio",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "PanelOfNormalsFilter",
        filter_name: "panel_of_normals",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "ClusteredEventsFilter",
        filter_name: "clustered_events",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "MultiallelicFilter",
        filter_name: "multiallelic",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "FragmentLengthFilter",
        filter_name: "fragment",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "PolymeraseSlippageFilter",
        filter_name: "slippage",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "FilteredHaplotypeFilter",
        filter_name: "haplotype",
        error_type: ErrorType::Artifact,
        arity: Arity::PerSite,
    },
    FilterIdentity {
        class: "GermlineFilter",
        filter_name: "germline",
        error_type: ErrorType::NonSomatic,
        arity: Arity::PerSite,
    },
];

/// `GATKVCFConstants.MUTECT_FILTER_NAMES`, in its own order, which is not the filters' order.
pub const MUTECT_FILTER_NAMES: [&str; 22] = [
    "PASS",
    "slippage",
    "panel_of_normals",
    "clustered_events",
    "weak_evidence",
    "germline",
    "multiallelic",
    "strand_bias",
    "normal_artifact",
    "base_qual",
    "map_qual",
    "fragment",
    "position",
    "contamination",
    "duplicate",
    "orientation",
    "haplotype",
    "strict_strand",
    "n_ratio",
    "low_allele_frac",
    "possible_numt",
    "FAIL",
];

/// `GATKVCFConstants.MUTECT_AS_FILTER_NAMES`, one entry, which becomes an `##INFO` line.
pub const MUTECT_AS_FILTER_NAMES: [&str; 1] = ["AS_FilterStatus"];

/// The `Description` of each `##FILTER` line, by filter name.
///
/// `normal_artifact`'s description is the string `artifact_in_normal`, which is another filter's
/// old name rather than a sentence: the header line describes itself with an identifier.
pub const FILTER_DESCRIPTIONS: [(&str, &str); 22] = [
    (
        "PASS",
        "Site contains at least one allele that passes filters",
    ),
    (
        "slippage",
        "Site filtered due to contraction of short tandem repeat region",
    ),
    ("panel_of_normals", "Blacklisted site in panel of normals"),
    ("clustered_events", "Clustered events observed in the tumor"),
    (
        "weak_evidence",
        "Mutation does not meet likelihood threshold",
    ),
    (
        "germline",
        "Evidence indicates this site is germline, not somatic",
    ),
    (
        "multiallelic",
        "Site filtered because too many alt alleles pass tumor LOD",
    ),
    (
        "strand_bias",
        "Evidence for alt allele comes from one read direction only",
    ),
    ("normal_artifact", "artifact_in_normal"),
    ("base_qual", "alt median base quality"),
    ("map_qual", "ref - alt median mapping quality"),
    ("fragment", "abs(ref - alt) median fragment length"),
    (
        "position",
        "median distance of alt variants from end of reads",
    ),
    ("contamination", "contamination"),
    (
        "duplicate",
        "evidence for alt allele is overrepresented by apparent duplicates",
    ),
    (
        "orientation",
        "orientation bias detected by the orientation bias mixture model",
    ),
    (
        "haplotype",
        "Variant near filtered variant on same haplotype.",
    ),
    (
        "strict_strand",
        "Evidence for alt allele is not represented in both directions",
    ),
    ("n_ratio", "Ratio of N to alt exceeds specified ratio"),
    (
        "low_allele_frac",
        "Allele fraction is below specified threshold",
    ),
    (
        "possible_numt",
        "Allele depth is below expected coverage of NuMT in autosome",
    ),
    (
        "FAIL",
        "Fail the site if all alleles fail but for different reasons.",
    ),
];

/// The two `##INFO` lines the tool declares beside the filters, rendered as `toString` does.
pub const INFO_LINES: [(&str, &str); 2] = [
    ("AS_FilterStatus", "INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=\"Filter status for each allele, as assessed by ApplyVQSR. Note that the VCF filter field will reflect the most lenient/sensitive status across all alleles.\">"),
    ("STRQ", "INFO=<ID=STRQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors\">"),
];

/// `FilterMutectCalls.FILTERING_STATUS_VCF_KEY`.
pub const FILTERING_STATUS_VCF_KEY: &str = "filtering_status";

/// The `##filtering_status` value Mutect2 writes, which the tool strips.
pub const MUTECT2_FILTERING_STATUS: &str =
    "Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove false positives.";

/// The value `FilterMutectCalls` writes over it, under the same key.
pub const FILTERED_FILTERING_STATUS: &str =
    "These calls have been filtered by FilterMutectCalls to \
     label false positives with a list of failed filters and true positives with PASS.";

/// The `##FILTER` line a name produces, as `VCFFilterHeaderLine.toString` renders it.
///
/// The rendering has **no leading `##`**: `toString` is the key and the value, and the writer adds
/// the hashes.
pub fn filter_line(name: &str) -> Option<String> {
    FILTER_DESCRIPTIONS
        .iter()
        .find(|(key, _)| *key == name)
        .map(|(_, description)| format!("FILTER=<ID={name},Description=\"{description}\">"))
}

/// `M2FiltersArgumentCollection.DEFAULT_MIN_MEDIAN_MAPPING_QUALITY`.
pub const DEFAULT_MIN_MEDIAN_MAPPING_QUALITY: i32 = 30;

/// `DEFAULT_MIN_MEDIAN_MAPPING_QUALITY_FOR_MICROBIAL`.
pub const DEFAULT_MIN_MEDIAN_MAPPING_QUALITY_FOR_MICROBIAL: i32 = 20;

/// The sentinel `minMedianMappingQuality` carries until it is asked for.
pub const UNSET_MIN_MEDIAN_MAPPING_QUALITY: i32 = -1;

/// `DEFAULT_LOG_SNV_PRIOR`, `MathUtils.log10ToLog(-6)`.
pub fn default_log_snv_prior() -> f64 {
    crate::allele_likelihoods::log10_to_log(-6.0)
}

/// `DEFAULT_LOG_INDEL_PRIOR`, `MathUtils.log10ToLog(-7)`.
pub fn default_log_indel_prior() -> f64 {
    crate::allele_likelihoods::log10_to_log(-7.0)
}

/// `DEFAULT_LOG_SNV_PRIOR_FOR_MITO`, `MathUtils.log10ToLog(-2.5)`.
pub fn default_log_snv_prior_for_mito() -> f64 {
    crate::allele_likelihoods::log10_to_log(-2.5)
}

/// `DEFAULT_LOG_INDEL_PRIOR_FOR_MITO`, `MathUtils.log10ToLog(-3.75)`.
pub fn default_log_indel_prior_for_mito() -> f64 {
    crate::allele_likelihoods::log10_to_log(-3.75)
}

/// `ReadOrientationFilter`, which is not in [`FILTERS`] because it is built only when a tar.gz of
/// artifact priors is supplied.
pub const READ_ORIENTATION_FILTER_CLASS: &str = "ReadOrientationFilter";

/// As much of `M2FiltersArgumentCollection` as decides which filters exist and how two of them are
/// configured.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterArguments {
    pub mitochondria: bool,
    pub microbial: bool,
    /// `minMedianMappingQuality`, [`UNSET_MIN_MEDIAN_MAPPING_QUALITY`] until asked for.
    pub min_median_mapping_quality: i32,
    pub log_snv_prior: f64,
    pub log_indel_prior: f64,
    /// Whether `readOrientationPriorTarGzs` is non-empty.
    pub read_orientation_priors: bool,
}

impl Default for FilterArguments {
    fn default() -> Self {
        Self {
            mitochondria: false,
            microbial: false,
            min_median_mapping_quality: UNSET_MIN_MEDIAN_MAPPING_QUALITY,
            log_snv_prior: default_log_snv_prior(),
            log_indel_prior: default_log_indel_prior(),
            read_orientation_priors: false,
        }
    }
}

impl FilterArguments {
    /// `getMinMedianMappingQuality()`, which **writes to the field it reads**.
    ///
    /// It takes `&mut self` for that reason. The memoisation is observable: a collection asked once
    /// while `microbial` is unset answers `30` and goes on answering `30` after the flag is set,
    /// where the same collection with the flag set first answers `20`.
    pub fn min_median_mapping_quality(&mut self) -> i32 {
        if self.min_median_mapping_quality == UNSET_MIN_MEDIAN_MAPPING_QUALITY {
            self.min_median_mapping_quality = if self.microbial {
                DEFAULT_MIN_MEDIAN_MAPPING_QUALITY_FOR_MICROBIAL
            } else {
                DEFAULT_MIN_MEDIAN_MAPPING_QUALITY
            };
        }
        self.min_median_mapping_quality
    }

    /// `getLogSnvPrior()`, which compares against the default **by value**.
    ///
    /// So passing the default explicitly under mitochondrial mode is indistinguishable from not
    /// passing it, and yields the mitochondrial prior.
    pub fn log_snv_prior(&self) -> f64 {
        if self.mitochondria && self.log_snv_prior == default_log_snv_prior() {
            default_log_snv_prior_for_mito()
        } else {
            self.log_snv_prior
        }
    }

    /// `getLogIndelPrior()`, the same shape.
    pub fn log_indel_prior(&self) -> f64 {
        if self.mitochondria && self.log_indel_prior == default_log_indel_prior() {
            default_log_indel_prior_for_mito()
        } else {
            self.log_indel_prior
        }
    }
}

/// `buildFiltersList`: the classes a run constructs, in construction order.
///
/// That order is the engine's: `ErrorProbabilities` collects the list into a `LinkedHashMap`, so it
/// is the iteration order everywhere downstream.
///
/// Eighteen by default, twelve in mitochondrial mode, thirteen in microbial. The two modes drop the
/// same six; microbial adds `PolymeraseSlippageFilter` back **at the end**, so that filter is
/// sixteenth by default and last in microbial mode.
///
/// `NormalArtifactFilter` is built in every mode, directly under a comment saying it does not apply
/// to mitochondria. The `add` has no guard.
///
/// The `ReadOrientationFilter` branch needs a tar.gz of artifact priors, which is a fixture rather
/// than an argument; no golden exercises it, and it is written here rather than refused because it
/// is a list membership with an unambiguous position, not an arithmetic path.
pub fn build_filters_list(arguments: &FilterArguments) -> Vec<&'static str> {
    // The first twelve, unconditionally, in the reference's order.
    let mut classes: Vec<&'static str> = FILTERS[..12].iter().map(|f| f.class).collect();

    if arguments.read_orientation_priors {
        classes.push(READ_ORIENTATION_FILTER_CLASS);
    }

    if !arguments.mitochondria && !arguments.microbial {
        classes.extend(FILTERS[12..].iter().map(|f| f.class));
    }

    if arguments.microbial {
        classes.push("PolymeraseSlippageFilter");
    }

    classes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_filter_name_is_distinct() {
        for (index, filter) in FILTERS.iter().enumerate() {
            assert!(
                !FILTERS[..index]
                    .iter()
                    .any(|other| other.filter_name == filter.filter_name),
                "{} repeats {}",
                filter.class,
                filter.filter_name
            );
        }
    }

    #[test]
    fn the_error_types_are_split_one_two_and_sixteen() {
        let sequencing = FILTERS
            .iter()
            .filter(|f| f.error_type == ErrorType::Sequencing)
            .count();
        let non_somatic = FILTERS
            .iter()
            .filter(|f| f.error_type == ErrorType::NonSomatic)
            .count();
        let artifact = FILTERS
            .iter()
            .filter(|f| f.error_type == ErrorType::Artifact)
            .count();
        assert_eq!((sequencing, non_somatic, artifact), (1, 2, 15));
    }

    #[test]
    fn the_header_names_are_not_the_filter_names() {
        // Two of the header's names belong to no filter at all.
        for name in ["PASS", "FAIL"] {
            assert!(!FILTERS.iter().any(|f| f.filter_name == name));
            assert!(MUTECT_FILTER_NAMES.contains(&name));
        }
        // And two filters the header names are built only under an argument this list cannot see.
        for name in ["orientation", "possible_numt"] {
            assert!(MUTECT_FILTER_NAMES.contains(&name));
            assert!(!FILTERS.iter().any(|f| f.filter_name == name));
        }
    }

    /// The three lists, and where the slippage filter sits in each.
    #[test]
    fn the_mode_decides_the_list_and_the_order() {
        let default = build_filters_list(&FilterArguments::default());
        assert_eq!(default.len(), 18);
        assert_eq!(
            default[15], "PolymeraseSlippageFilter",
            "sixteenth by default"
        );

        let mitochondria = build_filters_list(&FilterArguments {
            mitochondria: true,
            ..FilterArguments::default()
        });
        assert_eq!(mitochondria.len(), 12);
        assert!(!mitochondria.contains(&"PolymeraseSlippageFilter"));
        // The comment above the `add` says this one does not apply to mitochondria. It is here.
        assert!(mitochondria.contains(&"NormalArtifactFilter"));

        let microbial = build_filters_list(&FilterArguments {
            microbial: true,
            ..FilterArguments::default()
        });
        assert_eq!(microbial.len(), 13);
        assert_eq!(
            microbial[12], "PolymeraseSlippageFilter",
            "and last in microbial mode"
        );
    }

    /// The getter writes to the field it reads, so the order of calls decides the answer.
    #[test]
    fn the_mapping_quality_getter_remembers_its_first_answer() {
        let mut asked_first = FilterArguments::default();
        assert_eq!(asked_first.min_median_mapping_quality(), 30);
        asked_first.microbial = true;
        assert_eq!(
            asked_first.min_median_mapping_quality(),
            30,
            "the flag came too late"
        );

        let mut set_first = FilterArguments {
            microbial: true,
            ..FilterArguments::default()
        };
        assert_eq!(set_first.min_median_mapping_quality(), 20);
    }

    /// Passing the default explicitly is indistinguishable from not passing it.
    #[test]
    fn an_explicit_default_prior_is_the_mitochondrial_one() {
        let explicit = FilterArguments {
            mitochondria: true,
            log_snv_prior: default_log_snv_prior(),
            ..FilterArguments::default()
        };
        assert_eq!(explicit.log_snv_prior(), default_log_snv_prior_for_mito());

        let other = FilterArguments {
            mitochondria: true,
            log_snv_prior: -12.0,
            ..FilterArguments::default()
        };
        assert_eq!(other.log_snv_prior(), -12.0);
    }
}
