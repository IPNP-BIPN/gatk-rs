//! The annotations GATK discovers, and which of them a command line asks for.
//!
//! Ported from `org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor`
//! and `org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine` (GATK
//! 4.6.2.0), as far as a tool that never computes an annotation still reaches them. `CombineGVCFs` is such a tool:
//! it annotates nothing, and yet the annotations it was asked for decide two things in its output.
//! Their **descriptions** become header lines, and their **raw keys** are the INFO fields its merge
//! combines rather than takes the median of.
//!
//! # The table is what the reference instantiates
//!
//! [`CATALOGUE`] is every class `ClassFinder` finds under
//! `org.broadinstitute.hellbender.tools.walkers.annotator` that is neither abstract nor the base
//! interface, each instantiated with its no-argument constructor and asked for `getDescriptions()`,
//! `getRawDescriptions()`, `getRawKeyNames()` and `getKeyNames()`. The lines are kept as the
//! reference renders them and parsed on use, so a description is compared as the reference wrote it
//! rather than as this port would have rebuilt it.
//!
//! # A group is any interface that extends `Annotation`, found one class deep
//!
//! ```java
//! Collections.addAll(interfaces, annot.getClass().getInterfaces());
//! ```
//!
//! `getInterfaces()` is the class's **own** declaration, not its superclass's. So
//! `BaseQualityRankSumTest`, which declares `StandardAnnotation` and inherits `InfoFieldAnnotation`
//! from `RankSumTest`, belongs to one group, while `QualByDepth`, which declares both, belongs to
//! four. `-G InfoFieldAnnotation` is therefore a valid group that selects some info annotations
//! and not others, and `ClippingRankSumTest` belongs to no group at all.
//!
//! # The resolved list is sorted by simple name
//!
//! `getResolvedInstances` collects into a `TreeSet` ordered by `getSimpleName()`, so the order a
//! command line names annotations in is lost. It matters in one place: two reducible annotations
//! that share a raw key (`AS_FisherStrand` and `AS_StrandOddsRatio` both read `AS_SB_TABLE`) are
//! combined by whichever sorts first, and the second finds the key already consumed.

use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::header_lines::parse_meta_line;
use htsjdk_vcf::header_parse::VcfVersion;

/// Which of the four annotation interfaces a class implements, which is the list the engine files
/// it under.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Info,
    Genotype,
    JumboInfo,
    JumboGenotype,
}

/// One discovered annotation class.
#[derive(Debug)]
pub struct Entry {
    /// `getClass().getSimpleName()`, which is what `--annotation` names it by.
    pub name: &'static str,
    pub kind: Kind,
    /// The groups it is discovered in, breadth first over its own interfaces.
    pub groups: &'static [&'static str],
    /// `getRawKeyNames()` for a `ReducibleAnnotation`, and `None` for any other.
    pub raw_keys: Option<&'static [&'static str]>,
    /// `getKeyNames()`.
    pub key_names: &'static [&'static str],
    /// `getDescriptions()`, rendered without the leading `##`.
    pub descriptions: &'static [&'static str],
    /// `getRawDescriptions()`, rendered the same way, for a reducible annotation.
    pub raw_descriptions: &'static [&'static str],
}

impl Entry {
    /// The descriptions as header lines.
    pub fn description_lines(&self) -> Vec<HeaderLine> {
        parse(self.descriptions)
    }

    /// The raw descriptions as header lines.
    pub fn raw_description_lines(&self) -> Vec<HeaderLine> {
        parse(self.raw_descriptions)
    }

    /// `instanceof ReducibleAnnotation`.
    pub fn is_reducible(&self) -> bool {
        self.raw_keys.is_some()
    }
}

fn parse(lines: &[&str]) -> Vec<HeaderLine> {
    lines
        .iter()
        .map(|line| {
            parse_meta_line(&format!("##{line}"), VcfVersion::Vcf4_2, 0)
                .expect("a catalogue line the reference rendered")
        })
        .collect()
}

/// The entry for a simple class name.
pub fn entry(name: &str) -> Option<&'static Entry> {
    CATALOGUE.iter().find(|entry| entry.name == name)
}

/// Whether any discovered annotation names this group.
pub fn is_group(name: &str) -> bool {
    CATALOGUE.iter().any(|entry| entry.groups.contains(&name))
}

/// The annotations that declare an `@Argument` of their own, and so are the predecessors of
/// dependent arguments: `AllelePseudoDepth`'s Dirichlet prior arguments,
/// `--assembly-complexity-reference-mode`, the two `--denovo-*` thresholds and
/// `--allow-old-rms-mapping-quality-annotation-data`.
const WITH_ARGUMENTS: [&str; 4] = [
    "AllelePseudoDepth",
    "AssemblyComplexity",
    "PossibleDeNovo",
    "RMSMappingQuality",
];

/// `GATKAnnotationArgumentCollection`: what the command line said about annotations.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AnnotationArguments {
    /// `--annotation`, in command-line order.
    pub annotations: Vec<String>,
    /// `--annotation-group`.
    pub groups: Vec<String>,
    /// `--annotations-to-exclude`.
    pub excluded: Vec<String>,
    /// `--disable-tool-default-annotations`.
    pub disable_tool_defaults: bool,
    /// `--enable-all-annotations`.
    pub enable_all: bool,
}

/// What `validateAndResolvePlugins` refuses, each a `CommandLineException`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResolutionError {
    /// `CommandLineException.BadArgumentValue`, whose message Barclay prefixes.
    BadArgumentValue(String),
    /// A plain `CommandLineException`.
    CommandLine(String),
}

impl ResolutionError {
    /// The message as the reference reports it: `BadArgumentValue(String)` prefixes its own.
    pub fn message(&self) -> String {
        match self {
            ResolutionError::BadArgumentValue(message) => {
                format!("Illegal argument value: {message}")
            }
            ResolutionError::CommandLine(message) => message.clone(),
        }
    }
}

/// `Utils.getDuplicatedItems`, each name repeated on the command line, in the order its second
/// occurrence was seen. The reference collects them into a `HashSet`, so a line repeating two
/// names lists them in hash order instead; a single repeated name, which is every row measured,
/// reads the same either way.
fn duplicated(names: &[String]) -> Vec<String> {
    let mut seen: Vec<&String> = Vec::new();
    let mut duplicates: Vec<String> = Vec::new();
    for name in names {
        if seen.contains(&name) {
            if !duplicates.contains(name) {
                duplicates.push(name.clone());
            }
        } else {
            seen.push(name);
        }
    }
    duplicates
}

/// `validateAndResolvePlugins` and then `getResolvedInstances`, for a tool whose defaults are the
/// groups and annotations given.
///
/// The checks run in the reference's order, so a command line with two faults is refused for the
/// first of them. The result is sorted by simple name and has the excluded annotations removed.
pub fn resolve(
    arguments: &AnnotationArguments,
    default_groups: &[&str],
    default_annotations: &[&str],
) -> Result<Vec<&'static Entry>, ResolutionError> {
    let enabled_twice = duplicated(&arguments.annotations);
    if !enabled_twice.is_empty() {
        return Err(ResolutionError::BadArgumentValue(format!(
            "The annotation(s) are enabled more than once: {}",
            enabled_twice.join(", ")
        )));
    }
    let disabled_twice = duplicated(&arguments.excluded);
    if !disabled_twice.is_empty() {
        return Err(ResolutionError::BadArgumentValue(format!(
            "The annotation(s) are disabled more than once: {}",
            disabled_twice.join(", ")
        )));
    }
    let both: Vec<String> = arguments
        .annotations
        .iter()
        .filter(|name| arguments.excluded.contains(name))
        .cloned()
        .collect();
    if !both.is_empty() {
        return Err(ResolutionError::CommandLine(format!(
            "The annotation(s): {} are both enabled and disabled",
            both.join(", ")
        )));
    }
    for name in &arguments.excluded {
        if entry(name).is_none() {
            return Err(ResolutionError::BadArgumentValue(format!(
                "Disabled annotation ({name}) does not exist"
            )));
        }
    }
    // `requiredPredecessors`: Barclay asks `isDependentArgumentAllowed` of every argument an
    // annotation declares, given or not, so each annotation that HAS arguments and is selected
    // (by name, as a tool default, or through a group the command line or the tool enables) is
    // recorded as a predecessor. Disabling one then refuses the run, unless the tool names it
    // individually as a default, where the reference only warns (barclay#23).
    for name in &arguments.excluded {
        let selected = arguments.annotations.contains(name)
            || default_annotations.contains(&name.as_str())
            || arguments
                .groups
                .iter()
                .map(String::as_str)
                .chain(default_groups.iter().copied())
                .any(|group| entry(name).is_some_and(|found| found.groups.contains(&group)));
        if WITH_ARGUMENTS.contains(&name.as_str())
            && selected
            && !default_annotations.contains(&name.as_str())
        {
            return Err(ResolutionError::CommandLine(format!(
                "Values were supplied for ({name}) that is also disabled"
            )));
        }
    }
    for name in &arguments.annotations {
        if entry(name).is_none() && !default_annotations.contains(&name.as_str()) {
            return Err(ResolutionError::CommandLine(format!(
                "Unrecognized annotation name: {name}"
            )));
        }
    }
    for group in &arguments.groups {
        if !is_group(group) {
            return Err(ResolutionError::CommandLine(format!(
                "Unrecognized annotation group name: {group}"
            )));
        }
    }

    let mut names: Vec<&'static str> = Vec::new();
    let mut add = |entry: &'static Entry| {
        if !names.contains(&entry.name) {
            names.push(entry.name);
        }
    };
    if !arguments.disable_tool_defaults {
        for name in default_annotations {
            if let Some(found) = entry(name) {
                add(found);
            }
        }
        for group in default_groups {
            for found in CATALOGUE
                .iter()
                .filter(|entry| entry.groups.contains(group))
            {
                add(found);
            }
        }
    }
    for group in &arguments.groups {
        for found in CATALOGUE
            .iter()
            .filter(|entry| entry.groups.contains(&group.as_str()))
        {
            add(found);
        }
    }
    if arguments.enable_all {
        CATALOGUE.iter().for_each(&mut add);
    } else {
        for name in &arguments.annotations {
            if let Some(found) = entry(name) {
                add(found);
            }
        }
    }
    names.retain(|name| !arguments.excluded.iter().any(|excluded| excluded == name));
    // The `TreeSet` by simple name, which is `String.compareTo`: UTF-16 code units, which for these
    // ASCII names is byte order.
    names.sort_unstable();
    Ok(names
        .into_iter()
        .map(|name| entry(name).expect("a catalogued name"))
        .collect())
}

/// `VariantAnnotatorEngine.getVCFAnnotationDescriptions(useRaw)` over a resolved list, without the
/// overlap and expression lines the caller adds.
///
/// The info annotations come first, each reducible one contributing its raw lines when asked to
/// keep them and its finished lines unless `use_raw`; then the genotype annotations, then the two
/// jumbo kinds. Duplicates are the caller's to collapse, as they are the reference's `Set`'s.
pub fn descriptions(
    resolved: &[&'static Entry],
    use_raw: bool,
    keep_raw_combined: bool,
) -> Vec<HeaderLine> {
    let mut lines = Vec::new();
    for entry in resolved.iter().filter(|entry| entry.kind == Kind::Info) {
        if entry.is_reducible() {
            if use_raw || keep_raw_combined {
                lines.extend(entry.raw_description_lines());
            }
            if !use_raw {
                lines.extend(entry.description_lines());
            }
        } else {
            lines.extend(entry.description_lines());
        }
    }
    for kind in [Kind::Genotype, Kind::JumboInfo, Kind::JumboGenotype] {
        for entry in resolved.iter().filter(|entry| entry.kind == kind) {
            lines.extend(entry.description_lines());
        }
    }
    lines
}

/// Every annotation the reference discovers, sorted by simple name.
pub static CATALOGUE: &[Entry] = &[
    Entry {
        name: "AS_BaseQualityRankSumTest",
        kind: Kind::Info,
        groups: &["AS_StandardAnnotation"],
        raw_keys: Some(&["AS_RAW_BaseQRankSum"]),
        key_names: &["AS_BaseQRankSum"],
        descriptions: &["INFO=<ID=AS_BaseQRankSum,Number=A,Type=Float,Description=\"allele specific Z-score from Wilcoxon rank sum test of each Alt Vs. Ref base qualities\">"],
        raw_descriptions: &["INFO=<ID=AS_RAW_BaseQRankSum,Number=1,Type=String,Description=\"raw data for allele specific rank sum test of base qualities\">"],
    },
    Entry {
        name: "AS_FisherStrand",
        kind: Kind::Info,
        groups: &["AS_StandardAnnotation"],
        raw_keys: Some(&["AS_SB_TABLE"]),
        key_names: &["AS_FS"],
        descriptions: &["INFO=<ID=AS_FS,Number=A,Type=Float,Description=\"allele specific phred-scaled p-value using Fisher's exact test to detect strand bias of each alt allele\">"],
        raw_descriptions: &["INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.\">"],
    },
    Entry {
        name: "AS_InbreedingCoeff",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "AS_StandardAnnotation", "AlleleSpecificAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["AS_InbreedingCoeff"],
        descriptions: &["INFO=<ID=AS_InbreedingCoeff,Number=A,Type=Float,Description=\"Allele-specific inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "AS_MappingQualityRankSumTest",
        kind: Kind::Info,
        groups: &["AS_StandardAnnotation"],
        raw_keys: Some(&["AS_RAW_MQRankSum"]),
        key_names: &["AS_MQRankSum"],
        descriptions: &["INFO=<ID=AS_MQRankSum,Number=A,Type=Float,Description=\"Allele-specific Mapping Quality Rank Sum\">"],
        raw_descriptions: &["INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description=\"Allele-specfic raw data for Mapping Quality Rank Sum\">"],
    },
    Entry {
        name: "AS_QualByDepth",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "ReducibleAnnotation", "AS_StandardAnnotation", "AlleleSpecificAnnotation", "VariantAnnotation"],
        raw_keys: Some(&["AS_QUALapprox", "AS_QUAL", "AS_VarDP"]),
        key_names: &["AS_QD"],
        descriptions: &["INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Allele-specific Variant Confidence/Quality by Depth\">"],
        raw_descriptions: &["INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Allele-specific Variant Confidence/Quality by Depth\">"],
    },
    Entry {
        name: "AS_RMSMappingQuality",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "AS_StandardAnnotation", "ReducibleAnnotation", "AlleleSpecificAnnotation", "VariantAnnotation"],
        raw_keys: Some(&["AS_RAW_MQ"]),
        key_names: &["AS_MQ"],
        descriptions: &["INFO=<ID=AS_MQ,Number=A,Type=Float,Description=\"Allele-specific RMS Mapping Quality\">"],
        raw_descriptions: &["INFO=<ID=AS_RAW_MQ,Number=1,Type=String,Description=\"Allele-specfic raw data for RMS Mapping Quality\">"],
    },
    Entry {
        name: "AS_ReadPosRankSumTest",
        kind: Kind::Info,
        groups: &["AS_StandardAnnotation"],
        raw_keys: Some(&["AS_RAW_ReadPosRankSum"]),
        key_names: &["AS_ReadPosRankSum"],
        descriptions: &["INFO=<ID=AS_ReadPosRankSum,Number=A,Type=Float,Description=\"allele specific Z-score from Wilcoxon rank sum test of each Alt vs. Ref read position bias\">"],
        raw_descriptions: &["INFO=<ID=AS_RAW_ReadPosRankSum,Number=1,Type=String,Description=\"allele specific raw data for rank sum test of read position bias\">"],
    },
    Entry {
        name: "AS_StrandBiasMutectAnnotation",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardMutectAnnotation", "AlleleSpecificAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["AS_SB_TABLE"],
        descriptions: &["INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "AS_StrandOddsRatio",
        kind: Kind::Info,
        groups: &["AS_StandardAnnotation"],
        raw_keys: Some(&["AS_SB_TABLE"]),
        key_names: &["AS_SOR"],
        descriptions: &["INFO=<ID=AS_SOR,Number=A,Type=Float,Description=\"Allele specific strand Odds Ratio of 2x|Alts| contingency table to detect allele specific strand bias\">"],
        raw_descriptions: &["INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.\">"],
    },
    Entry {
        name: "AlleleFraction",
        kind: Kind::Genotype,
        groups: &["GenotypeAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["AF"],
        descriptions: &["FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "AllelePseudoDepth",
        kind: Kind::Genotype,
        groups: &["GenotypeAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["DD", "DF"],
        descriptions: &["FORMAT=<ID=DD,Number=R,Type=Float,Description=\"Allele depth based on Dirichlet posterior pseudo-counts\">", "FORMAT=<ID=DF,Number=R,Type=Float,Description=\"Allele Fraction based on Dirichlet posterior pseudo-counts\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "AssemblyComplexity",
        kind: Kind::JumboInfo,
        groups: &["JumboInfoAnnotation", "AlleleSpecificAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["HEC", "HAPCOMP", "HAPDOM"],
        descriptions: &["INFO=<ID=HEC,Number=.,Type=Integer,Description=\"Counts of support for haplotype groups excluding difference at the site in question.\">", "INFO=<ID=HAPCOMP,Number=A,Type=Integer,Description=\"Edit distances of each alt allele's most common supporting haplotype from closest germline haplotype, excluding differences at the site in question.\">", "INFO=<ID=HAPDOM,Number=A,Type=Float,Description=\"For each alt allele, fraction of read support that best fits the most-supported haplotype containing the allele\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "BaseQuality",
        kind: Kind::Info,
        groups: &["StandardMutectAnnotation"],
        raw_keys: None,
        key_names: &["MBQ"],
        descriptions: &["INFO=<ID=MBQ,Number=R,Type=Integer,Description=\"median base quality by allele\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "BaseQualityHistogram",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["BQHIST"],
        descriptions: &["INFO=<ID=BQHIST,Number=A,Type=Integer,Description=\"Base quality counts for each allele represented sparsely as alternating entries of qualities and counts for each allele.For example [10,1,0,20,0,1] means one ref base with quality 10 and one alt base with quality 20.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "BaseQualityRankSumTest",
        kind: Kind::Info,
        groups: &["StandardAnnotation"],
        raw_keys: None,
        key_names: &["BaseQRankSum"],
        descriptions: &["INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "ChromosomeCounts",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["AN", "AC", "AF"],
        descriptions: &["INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">", "INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">", "INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "ClippingRankSumTest",
        kind: Kind::Info,
        groups: &[],
        raw_keys: None,
        key_names: &["ClippingRankSum"],
        descriptions: &["INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "CountNs",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["NCount"],
        descriptions: &["INFO=<ID=NCount,Number=1,Type=Integer,Description=\"Count of N bases in the pileup\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "Coverage",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["DP"],
        descriptions: &["INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "CycleSkipStatus",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_CSS"],
        descriptions: &["INFO=<ID=X_CSS,Number=A,Type=String,Description=\"Flow: cycle skip status: cycle-skip, possible-cycle-skip, non-skip\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "DepthPerAlleleBySample",
        kind: Kind::Genotype,
        groups: &["GenotypeAnnotation", "StandardAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["AD"],
        descriptions: &["FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "DepthPerSampleHC",
        kind: Kind::Genotype,
        groups: &["GenotypeAnnotation", "StandardHCAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["DP"],
        descriptions: &["FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "ExcessHet",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["ExcessHet"],
        descriptions: &["INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "FisherStrand",
        kind: Kind::Info,
        groups: &["StandardAnnotation"],
        raw_keys: None,
        key_names: &["FS"],
        descriptions: &["INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "FragmentDepthPerAlleleBySample",
        kind: Kind::JumboGenotype,
        groups: &["JumboGenotypeAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["FAD"],
        descriptions: &["FORMAT=<ID=FAD,Number=R,Type=Integer,Description=\"Count of fragments supporting each allele.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "FragmentLength",
        kind: Kind::Info,
        groups: &["StandardMutectAnnotation"],
        raw_keys: None,
        key_names: &["MFRL"],
        descriptions: &["INFO=<ID=MFRL,Number=R,Type=Integer,Description=\"median fragment length by allele\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "GcContent",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_GCC"],
        descriptions: &["INFO=<ID=X_GCC,Number=1,Type=Float,Description=\"Flow: percentage of G or C in the window around hmer\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "GenotypeSummaries",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["NCC", "GQ_MEAN", "GQ_STDDEV"],
        descriptions: &["INFO=<ID=NCC,Number=1,Type=Integer,Description=\"Number of no-called samples\">", "INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description=\"Mean of all GQ values\">", "INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description=\"Standard deviation of all GQ values\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "HaplotypeFilteringAnnotation",
        kind: Kind::JumboInfo,
        groups: &["JumboInfoAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["ASSEMBLED_HAPS", "FILTERED_HAPS"],
        descriptions: &["INFO=<ID=ASSEMBLED_HAPS,Number=1,Type=Integer,Description=\"Haplotypes detected by the assembly region before haplotype filtering is applied\">", "INFO=<ID=FILTERED_HAPS,Number=1,Type=Integer,Description=\"Haplotypes filtered out by the haplotype filtering code\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "HmerIndelLength",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_HIL"],
        descriptions: &["INFO=<ID=X_HIL,Number=A,Type=Integer,Description=\"Flow: length of the hmer indel, if so\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "HmerIndelNuc",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_HIN"],
        descriptions: &["INFO=<ID=X_HIN,Number=A,Type=String,Description=\"Flow: nucleotide of the hmer indel, if so\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "HmerMotifs",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_LM", "X_RM"],
        descriptions: &["INFO=<ID=X_LM,Number=A,Type=String,Description=\"Flow: motif to the left of the indel\">", "INFO=<ID=X_RM,Number=A,Type=String,Description=\"Flow: motif to the right of the indel\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "InbreedingCoeff",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["InbreedingCoeff"],
        descriptions: &["INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "IndelClassify",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_IC"],
        descriptions: &["INFO=<ID=X_IC,Number=A,Type=String,Description=\"Flow: indel class: ins, del, NA\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "IndelLength",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["X_IL"],
        descriptions: &["INFO=<ID=X_IL,Number=A,Type=Integer,Description=\"Flow: length of indel\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "LikelihoodRankSumTest",
        kind: Kind::Info,
        groups: &[],
        raw_keys: None,
        key_names: &["LikelihoodRankSum"],
        descriptions: &["INFO=<ID=LikelihoodRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref haplotype likelihoods\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "MappingQuality",
        kind: Kind::Info,
        groups: &["StandardMutectAnnotation"],
        raw_keys: None,
        key_names: &["MMQ"],
        descriptions: &["INFO=<ID=MMQ,Number=R,Type=Integer,Description=\"median mapping quality by allele\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "MappingQualityRankSumTest",
        kind: Kind::Info,
        groups: &["StandardAnnotation"],
        raw_keys: None,
        key_names: &["MQRankSum"],
        descriptions: &["INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "MappingQualityZero",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["MQ0"],
        descriptions: &["INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "OrientationBiasReadCounts",
        kind: Kind::JumboGenotype,
        groups: &["JumboGenotypeAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["F1R2", "F2R1"],
        descriptions: &["FORMAT=<ID=F1R2,Number=R,Type=Integer,Description=\"Count of reads in F1R2 pair orientation supporting each allele\">", "FORMAT=<ID=F2R1,Number=R,Type=Integer,Description=\"Count of reads in F2R1 pair orientation supporting each allele\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "OriginalAlignment",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["OCM"],
        descriptions: &["INFO=<ID=OCM,Number=1,Type=Integer,Description=\"Number of alt reads whose original alignment doesn't match the current contig.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "PossibleDeNovo",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["hiConfDeNovo", "loConfDeNovo"],
        descriptions: &["INFO=<ID=hiConfDeNovo,Number=1,Type=String,Description=\"High confidence possible de novo mutation (GQ >= 20 for all trio members)=[comma-delimited list of child samples]\">", "INFO=<ID=loConfDeNovo,Number=1,Type=String,Description=\"Low confidence possible de novo mutation (GQ >= 10 for child, GQ > 0 for parents)=[comma-delimited list of child samples]\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "QualByDepth",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["QD"],
        descriptions: &["INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "RMSMappingQuality",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardAnnotation", "ReducibleAnnotation", "VariantAnnotation"],
        raw_keys: Some(&["RAW_MQandDP"]),
        key_names: &["MQ", "RAW_MQandDP"],
        descriptions: &["INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">"],
        raw_descriptions: &["INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description=\"Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.\">"],
    },
    Entry {
        name: "RawGtCount",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "ReducibleAnnotation", "VariantAnnotation"],
        raw_keys: Some(&["RAW_GT_COUNT"]),
        key_names: &["RAW_GT_COUNT"],
        descriptions: &["INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer,Description=\"Counts of genotypes w.r.t. the reference allele in the following order: 0/0, 0/*, */*, i.e. all alts lumped together; for use in calculating excess heterozygosity\">"],
        raw_descriptions: &["INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer,Description=\"Counts of genotypes w.r.t. the reference allele in the following order: 0/0, 0/*, */*, i.e. all alts lumped together; for use in calculating excess heterozygosity\">"],
    },
    Entry {
        name: "ReadPosRankSumTest",
        kind: Kind::Info,
        groups: &["StandardAnnotation"],
        raw_keys: None,
        key_names: &["ReadPosRankSum"],
        descriptions: &["INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "ReadPosition",
        kind: Kind::Info,
        groups: &["StandardMutectAnnotation"],
        raw_keys: None,
        key_names: &["MPOS"],
        descriptions: &["INFO=<ID=MPOS,Number=A,Type=Integer,Description=\"median distance from end of read\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "ReferenceBases",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["REF_BASES"],
        descriptions: &["INFO=<ID=REF_BASES,Number=1,Type=String,Description=\"local reference bases.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "SampleList",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["Samples"],
        descriptions: &["INFO=<ID=Samples,Number=.,Type=String,Description=\"List of polymorphic samples\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "StrandBiasBySample",
        kind: Kind::Genotype,
        groups: &["GenotypeAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["SB"],
        descriptions: &["FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "StrandOddsRatio",
        kind: Kind::Info,
        groups: &["StandardAnnotation"],
        raw_keys: None,
        key_names: &["SOR"],
        descriptions: &["INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "TandemRepeat",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "StandardMutectAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["STR", "RU", "RPA"],
        descriptions: &["INFO=<ID=STR,Number=0,Type=Flag,Description=\"Variant is a short tandem repeat\">", "INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">", "INFO=<ID=RPA,Number=R,Type=Integer,Description=\"Number of times tandem repeat unit is repeated, for each allele (including reference)\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "TransmittedSingleton",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["transmittedSingleton", "nonTransmittedSingleton"],
        descriptions: &["INFO=<ID=transmittedSingleton,Number=1,Type=String,Description=\"Possible transmitted singleton (site with AC=2 from parent and child). Parent ID is listed.\">", "INFO=<ID=nonTransmittedSingleton,Number=1,Type=String,Description=\"Possible non transmitted singleton (site with AC=1 in just one parent). Parent ID is listed.\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "UniqueAltReadCount",
        kind: Kind::Info,
        groups: &["InfoFieldAnnotation", "AlleleSpecificAnnotation", "VariantAnnotation"],
        raw_keys: None,
        key_names: &["AS_UNIQ_ALT_READ_COUNT"],
        descriptions: &["INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description=\"Number of reads with unique start and mate end positions for each alt at a variant site\">"],
        raw_descriptions: &[],
    },
    Entry {
        name: "VariantType",
        kind: Kind::Info,
        groups: &["StandardFlowBasedAnnotation"],
        raw_keys: None,
        key_names: &["VARIANT_TYPE"],
        descriptions: &["INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description=\"Flow: type of variant: SNP/NON-H-INDEL/H-INDEL\">"],
        raw_descriptions: &[],
    },
];
