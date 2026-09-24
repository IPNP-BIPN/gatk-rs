//! `ReferenceConfidenceVariantContextMerger`, ported from
//! `org.broadinstitute.hellbender.tools.walkers` (GATK 4.6.2.0).
//!
//! How the records of several GVCFs that meet at one position become one record. `CombineGVCFs`
//! calls it wherever a stopped record carries a real alternate; a site of reference blocks alone
//! takes the tool's own cheaper merge instead.
//!
//! # The genotypes are NOT remapped before they are called
//!
//! ```java
//! return merge(vcs, loc, refBase, removeNonRefSymbolicAllele, samplesAreUniquified, false);
//! ```
//!
//! `CombineGVCFs` calls the overload that passes `useRemappedAllelesForGenotyping == false`, the one
//! the reference's own TODO (#8317) wants retired. So under `--call-genotypes` a genotype is matched
//! against the merged alleles by its ORIGINAL alleles: a deletion `AT>A` merged against a longer
//! reference `ATG` is no longer in the list as `A`, and the call falls back to the reference.
//!
//! # A number that is not reducible is replaced by its median, parsed the way `toString` prints it
//!
//! ```java
//! if (!value.toString().contains(",")) { ... }
//! else { String[] valueArray = value.toString().split("\\[|" + AnnotationUtils.LIST_DELIMITER +  "|\\]"); ... }
//! ```
//!
//! A multi-valued INFO field is decoded as a `List`, and a `List`'s `toString` is `[a, b]` with a
//! SPACE after each comma. The split is on the bracket and the comma alone, so every element after
//! the first keeps a leading space. `Double.parseDouble` trims it and `Integer.parseInt` does not:
//! a `Float` list contributes every element, an `Integer` list contributes its FIRST element and
//! then throws, and the throw is caught and logged once. `MBQ=40,40,0` therefore takes part in the
//! median as the single value `40`.
//!
//! # The median is the upper one
//!
//! `Utils.getMedianValue` sorts and takes element `size / 2`, so two values give the larger, and a
//! value that is a `Double` in one input keeps its `Double` rendering: `0.0` becomes `0.00`, and any
//! negative value takes `formatVCFDouble`'s exponent branch, `-0.524` becoming `-5.240e-01`.

use gatk_annotation::catalogue::Entry;
use gatk_engine::genotype_index::{
    genotypes_in_canonical_order, index_of_first_genotype_with_allele, subsetted_pl_indices,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext, NO_LOG10_PERROR};

/// `GATKVariantContextUtils.SUM_GL_THRESH_NOCALL`.
pub const SUM_GL_THRESH_NOCALL: f64 = -0.1;

/// `GATKVariantContextUtils.DEFAULT_PLOIDY`.
pub const DEFAULT_PLOIDY: usize = 2;

/// `Allele.NON_REF_ALLELE`.
pub fn non_ref() -> Allele {
    Allele::from_str("<NON_REF>", false).expect("a symbolic allele")
}

/// `Allele.SPAN_DEL`.
pub fn span_del() -> Allele {
    Allele::from_str("*", false).expect("the spanning deletion")
}

/// `GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED`.
fn deprecated_span_del() -> Allele {
    Allele::from_str("<*:DEL>", false).expect("a symbolic allele")
}

/// `VariantContext.Type`, as `determineType` decides it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// `typeOfBiallelicVariant(ref, allele)`: a symbolic alternate is symbolic, one of the reference's
/// length is a SNP or an MNP, and anything else is an indel. `*` against a one-base reference is a
/// SNP by that rule, and against a longer one an indel.
fn biallelic_type(reference: &Allele, allele: &Allele) -> VariantType {
    if allele.is_symbolic() {
        VariantType::Symbolic
    } else if reference.len() == allele.len() {
        if allele.len() == 1 {
            VariantType::Snp
        } else {
            VariantType::Mnp
        }
    } else {
        VariantType::Indel
    }
}

/// `VariantContext.getType()`.
pub fn variant_type(vc: &VariantContext) -> VariantType {
    if vc.alleles.len() == 1 {
        return VariantType::NoVariation;
    }
    let reference = vc.reference();
    let mut decided: Option<VariantType> = None;
    for allele in vc.alternate_alleles() {
        let this = biallelic_type(reference, allele);
        match decided {
            None => decided = Some(this),
            Some(previous) if previous != this => return VariantType::Mixed,
            Some(_) => {}
        }
    }
    decided.unwrap_or(VariantType::NoVariation)
}

/// `Allele.getBases()` for a base allele.
fn bases_of(allele: &Allele) -> Vec<u8> {
    allele.display_string().into_bytes()
}

/// `GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef`.
pub fn is_unmixed_mnp_ignoring_non_ref(vc: &VariantContext) -> bool {
    let length = vc.reference().len();
    if length < 2 {
        return false;
    }
    let non_ref = non_ref();
    vc.alleles.iter().all(|allele| {
        if allele.is_symbolic() {
            *allele == non_ref
        } else {
            allele.len() == length
        }
    })
}

/// What the merge refuses, in the reference's own classes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeError {
    /// `IllegalStateException` from `determineReferenceAllele`, two references of one length that
    /// differ.
    InconsistentReferences {
        contig: String,
        start: i64,
        first: String,
        second: String,
    },
    /// `IllegalStateException("the wrong reference was selected")`, a record whose reference is
    /// longer than the merged one.
    WrongReference,
    /// `UserException` from `getIndexesOfRelevantAllelesForGVCF`: a record without `<NON_REF>`.
    MissingNonRef { position: i64 },
    /// An index the reference would have thrown `ArrayIndexOutOfBoundsException` on, a PL shorter
    /// than its ploidy and allele count say.
    IndexOutOfBounds { index: usize, length: usize },
    /// `ClassCastException` from sorting a median list that mixed `Integer` and `Double`.
    MixedMedian { key: String },
    /// A raw annotation value the reducible parser refuses.
    Annotation { class: String, message: String },
}

impl MergeError {
    /// The Java class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            MergeError::InconsistentReferences { .. } | MergeError::WrongReference => {
                "java.lang.IllegalStateException"
            }
            MergeError::MissingNonRef { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            MergeError::IndexOutOfBounds { .. } => "java.lang.ArrayIndexOutOfBoundsException",
            MergeError::MixedMedian { .. } => "java.lang.ClassCastException",
            MergeError::Annotation { .. } => "java.lang.NumberFormatException",
        }
    }

    /// The exception's message.
    pub fn message(&self) -> String {
        match self {
            MergeError::InconsistentReferences {
                contig,
                start,
                first,
                second,
            } => format!(
                "The provided variant file(s) have inconsistent references for the same position(s) \
                 at {contig}:{start}, {first} vs. {second}"
            ),
            MergeError::WrongReference => "the wrong reference was selected".to_string(),
            MergeError::MissingNonRef { position } => format!(
                "The list of input alleles must contain <NON_REF> as an allele but that is not the \
                 case at position {position}; please use the Haplotype Caller with gVCF output to \
                 generate appropriate records"
            ),
            MergeError::IndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            MergeError::MixedMedian { key } => format!(
                "the values of {key} mix java.lang.Integer and java.lang.Double, which do not compare"
            ),
            MergeError::Annotation { message, .. } => message.clone(),
        }
    }
}

/// `Allele.toString()`, which is what the inconsistent-reference message prints: the bases, and a
/// `*` after a reference.
fn allele_to_string(allele: &Allele) -> String {
    let mut text = allele.display_string();
    if allele.is_reference() {
        text.push('*');
    }
    text
}

/// A record and the alleles it is remapped to, `VCWithNewAlleles`.
struct Remapped<'a> {
    vc: &'a VariantContext,
    new_alleles: Vec<Allele>,
    is_spanning_event: bool,
}

impl Remapped<'_> {
    /// `filterAllelesForFinalSet`: the called, non-reference alleles other than `<NON_REF>`, and no
    /// symbolic one from a record that is itself symbolic.
    fn alleles_for_final_set(&self) -> impl Iterator<Item = &Allele> + '_ {
        let non_ref = non_ref();
        let symbolic_record = variant_type(self.vc) == VariantType::Symbolic;
        self.new_alleles.iter().filter(move |allele| {
            **allele != non_ref
                && !allele.is_reference()
                && !(allele.is_symbolic() && symbolic_record)
                && !allele.is_no_call()
        })
    }

    fn is_spanning_deletion(&self) -> bool {
        (self.is_spanning_event && variant_type(self.vc) == VariantType::Mixed)
            || self.vc.alleles.contains(&span_del())
            || self.vc.alleles.contains(&deprecated_span_del())
    }

    fn is_non_spanning_event(&self) -> bool {
        !self.is_spanning_event && variant_type(self.vc) == VariantType::Mixed
    }
}

/// `replaceWithNoCallsAndDels`: a spanning record's reference becomes a no-call, and each alternate
/// a `*` if it is shorter than the reference, `<NON_REF>` if it is `<NON_REF>`, and a no-call
/// otherwise. A somatic merge never writes `*`.
fn replace_with_no_calls_and_dels(vc: &VariantContext, somatic: bool) -> Vec<Allele> {
    let non_ref = non_ref();
    let mut result = vec![Allele::no_call()];
    for allele in vc.alternate_alleles() {
        if *allele == non_ref {
            result.push(allele.clone());
        } else if !somatic && allele.len() < vc.reference().len() {
            result.push(span_del());
        } else {
            result.push(Allele::no_call());
        }
    }
    result
}

/// `remapAlleles`: the merged reference, and each alternate extended by the reference bases the
/// record's own reference lacked.
pub fn remap_alleles(vc: &VariantContext, reference: &Allele) -> Result<Vec<Allele>, MergeError> {
    let reference_bases = bases_of(reference);
    let extra = reference_bases.len() as i64 - vc.reference().len() as i64;
    if extra < 0 {
        return Err(MergeError::WrongReference);
    }
    let extra = extra as usize;
    let span_del = span_del();
    let mut result = vec![reference.clone()];
    for allele in vc.alternate_alleles() {
        if allele.is_symbolic() || *allele == span_del || allele.is_no_call() {
            result.push(allele.clone());
        } else if extra > 0 {
            let mut bases = bases_of(allele);
            bases.extend_from_slice(&reference_bases[reference_bases.len() - extra..]);
            result.push(Allele::create(&bases, false).expect("the bases of an allele, extended"));
        } else {
            result.push(allele.clone());
        }
    }
    Ok(result)
}

/// `GATKVariantContextUtils.determineReferenceAllele(VCs, loc)`: the longest reference among the
/// records that start at the position, or `None` when none does.
fn determine_reference_allele(
    vcs: &[VariantContext],
    contig: &str,
    start: i64,
) -> Result<Option<Allele>, MergeError> {
    let mut reference: Option<Allele> = None;
    for vc in vcs {
        // `contextMatchesLoc`: the same contig and the same start.
        if vc.contig != contig || vc.start != start {
            continue;
        }
        let mine = vc.reference().clone();
        reference = Some(match reference {
            None => mine,
            Some(current) if current.len() < mine.len() => mine,
            Some(current) if mine.len() < current.len() => current,
            Some(current) if current != mine => {
                return Err(MergeError::InconsistentReferences {
                    contig: vc.contig.clone(),
                    start: vc.start,
                    first: allele_to_string(&current),
                    second: allele_to_string(&mine),
                })
            }
            Some(current) => current,
        });
    }
    Ok(reference)
}

/// `ReferenceConfidenceVariantContextMerger.getBestDepthValue`: `MIN_DP` when the genotype has one,
/// its `DP` otherwise.
fn best_depth_value(genotype: &Genotype) -> Result<i64, MergeError> {
    if let Some(value) = genotype.get("MIN_DP") {
        let text = value_to_string(value);
        return text.parse::<i64>().map_err(|_| MergeError::Annotation {
            class: "java.lang.NumberFormatException".to_string(),
            message: format!("For input string: \"{text}\""),
        });
    }
    Ok(genotype.dp.unwrap_or(0) as i64)
}

/// `calculateVCDepth`: the record's `DP` when it has one, the sum of its genotypes' best depths
/// otherwise.
fn vc_depth(vc: &VariantContext) -> Result<i64, MergeError> {
    if let Some((_, value)) = vc.attributes.iter().find(|(key, _)| key == "DP") {
        // `getAttributeAsInt(DP, 0)`: a missing value is the default, anything else is parsed.
        return Ok(match value {
            Value::Missing => 0,
            Value::Int(number) => *number,
            other => {
                let text = value_to_string(other);
                if text == "." {
                    0
                } else {
                    text.parse::<i64>().unwrap_or_else(|_| {
                        htsjdk_vcf::genotype_likelihoods::parse_java_double(&text)
                            .map(|d| d as i64)
                            .unwrap_or(0)
                    })
                }
            }
        });
    }
    let mut depth = 0;
    for genotype in vc.genotypes.iter() {
        depth += best_depth_value(genotype)?;
    }
    Ok(depth)
}

/// `Object.toString()` of a decoded attribute: a string is itself, a list is `[a, b]`.
pub fn value_to_string(value: &Value) -> String {
    match value {
        Value::Missing => ".".to_string(),
        Value::Int(number) => number.to_string(),
        Value::Double(number) => number.to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Str(text) => text.clone(),
        Value::List(items) => format!(
            "[{}]",
            items
                .iter()
                .map(value_to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// `VariantContext.getAttributeAsList(key)`: a list is its elements, anything else one element.
fn attribute_as_list(value: &Value) -> Vec<String> {
    match value {
        Value::List(items) => items.iter().map(value_to_string).collect(),
        other => vec![value_to_string(other)],
    }
}

/// A value the median sorts, boxed as the header declares it.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Numeric {
    Integer(i32),
    Double(f64),
}

/// `Double.compareTo`: `-0.0` before `0.0`, and every `NaN` equal to every other and after
/// everything else.
fn compare_doubles(left: f64, right: f64) -> std::cmp::Ordering {
    match (left.is_nan(), right.is_nan()) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) => left.total_cmp(&right),
    }
}

/// The type the header declares for an INFO key, `None` when it declares none.
fn info_type(header: &VcfHeader, key: &str) -> Option<LineType> {
    header.lines.iter().find_map(|line| match line {
        HeaderLine::Compound {
            key: kind,
            id,
            line_type,
            ..
        } if kind == "INFO" && id == key => Some(*line_type),
        _ => None,
    })
}

/// A FORMAT line's count type.
fn format_cardinality(header: &VcfHeader, key: &str) -> Option<Cardinality> {
    header.lines.iter().find_map(|line| match line {
        HeaderLine::Compound {
            key: kind,
            id,
            number,
            ..
        } if kind == "FORMAT" && id == key => Some(*number),
        _ => None,
    })
}

/// An INFO line's count type.
fn info_cardinality(header: &VcfHeader, key: &str) -> Option<Cardinality> {
    header.lines.iter().find_map(|line| match line {
        HeaderLine::Compound {
            key: kind,
            id,
            number,
            ..
        } if kind == "INFO" && id == key => Some(*number),
        _ => None,
    })
}

/// `Integer.parseInt`: an optional sign and decimal digits, nothing else, not even a space.
fn java_parse_int(text: &str) -> Option<i32> {
    let digits = text
        .strip_prefix('-')
        .or_else(|| text.strip_prefix('+'))
        .unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    text.parse::<i32>().ok()
}

/// `parseNumericInfoAttributeValue`: the header's type decides the box, and a key the header does
/// not declare is a `Double` when it holds a `.`.
fn parse_numeric(header: &VcfHeader, key: &str, text: &str) -> Option<Numeric> {
    match info_type(header, key) {
        None => {
            if text.contains('.') {
                htsjdk_vcf::genotype_likelihoods::parse_java_double(text).map(Numeric::Double)
            } else {
                java_parse_int(text).map(Numeric::Integer)
            }
        }
        Some(LineType::Integer) => java_parse_int(text).map(Numeric::Integer),
        Some(LineType::Float) => {
            htsjdk_vcf::genotype_likelihoods::parse_java_double(text).map(Numeric::Double)
        }
        Some(_) => None,
    }
}

/// `SOMATIC_INFO_ANNOTATIONS_TO_MOVE`.
const SOMATIC_INFO_TO_MOVE: [&str; 1] = ["TLOD"];
/// `SOMATIC_INFO_ANNOTATIONS_TO_DROP`.
const SOMATIC_INFO_TO_DROP: [&str; 1] = ["POPAF"];
/// `Mutect2FilteringEngine.STANDARD_MUTECT_INFO_FIELDS_FOR_FILTERING`.
pub const MUTECT_FILTERING_INFO: [&str; 4] = ["MMQ", "MBQ", "MPOS", "MFRL"];
/// `SOMATIC_FORMAT_ANNOTATIONS_TO_KEEP`.
const SOMATIC_FORMAT_TO_KEEP: [&str; 4] = ["OCM", "PGT", "PID", "PS"];

/// One key's collected values, which are either raw strings for a reducible annotation or numbers.
enum Collected {
    /// The allele list each raw string was written against, and the string.
    Reducible(Vec<(Vec<Allele>, String)>),
    Numbers(Vec<Numeric>),
}

/// The merger, configured the way `CombineGVCFs.onTraversalStart` builds it.
pub struct Merger<'a> {
    /// The inputs' merged header, whose INFO and FORMAT lines type the values.
    pub header: &'a VcfHeader,
    /// The resolved annotations, in the engine's order. Only the reducible info ones are read.
    pub annotations: Vec<&'static Entry>,
    /// `--input-is-somatic`.
    pub somatic: bool,
    /// `--drop-somatic-filtering-annotations`.
    pub drop_somatic_filtering_annotations: bool,
    /// `--call-genotypes`.
    pub call_genotypes: bool,
}

/// What `merge` hands back beside the record: the warnings the reference logs once.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct MergeNotes {
    pub invalid_annotation: bool,
}

impl Merger<'_> {
    /// `annotatorEngine.isRequestedReducibleRawKey(key)`.
    fn is_reducible_raw_key(&self, key: &str) -> bool {
        self.annotations.iter().any(|entry| {
            entry.kind == gatk_annotation::catalogue::Kind::Info
                && entry.raw_keys.is_some_and(|keys| keys.contains(&key))
        })
    }

    /// `merge(vcs, loc, refBase, removeNonRefSymbolicAllele, samplesAreUniquified)`.
    ///
    /// `None` is the reference's `null`, which only a site with no record starting at it and no
    /// reference base can produce.
    pub fn merge(
        &self,
        vcs: &[VariantContext],
        contig: &str,
        start: i64,
        reference_base: Option<u8>,
        remove_non_ref: bool,
        notes: &mut MergeNotes,
    ) -> Result<Option<VariantContext>, MergeError> {
        let reference = match determine_reference_allele(vcs, contig, start)? {
            Some(found) => found,
            None => match reference_base {
                Some(base) => {
                    Allele::create(&[base], true).map_err(|error| MergeError::Annotation {
                        class: "java.lang.IllegalArgumentException".to_string(),
                        message: error.to_string(),
                    })?
                }
                None => return Ok(None),
            },
        };

        let mut remapped: Vec<Remapped> = Vec::with_capacity(vcs.len());
        let mut aggregated_filters: Vec<String> = Vec::new();
        let mut saw_pass_sample = false;
        for vc in vcs {
            let is_spanning_event = start != vc.start;
            let new_alleles = if is_spanning_event {
                replace_with_no_calls_and_dels(vc, self.somatic)
            } else {
                remap_alleles(vc, &reference)?
            };
            remapped.push(Remapped {
                vc,
                new_alleles,
                is_spanning_event,
            });
            if self.somatic && vc.filters_were_applied() {
                if vc.is_filtered() {
                    for filter in vc.filters.iter().flatten() {
                        if !aggregated_filters.contains(filter) {
                            aggregated_filters.push(filter.clone());
                        }
                    }
                } else {
                    saw_pass_sample = true;
                }
            }
        }

        let alleles = collect_target_alleles(&remapped, &reference, remove_non_ref);

        let mut ids: Vec<String> = Vec::new();
        let mut depth: i64 = 0;
        let mut collected: Vec<(String, Collected)> = Vec::new();
        let mut genotypes: Vec<Genotype> = Vec::new();
        for pair in &remapped {
            genotypes.extend(self.merge_genotypes(pair.vc, &pair.new_alleles, &alleles)?);
            depth += vc_depth(pair.vc)?;
            if start != pair.vc.start {
                continue;
            }
            if pair.vc.id != "." && !ids.contains(&pair.vc.id) {
                ids.push(pair.vc.id.clone());
            }
            self.add_attributes(pair, &mut collected, notes);
        }

        let mut attributes = self.merge_attributes(depth, &alleles, collected)?;
        attributes.sort_by(|a, b| a.0.cmp(&b.0));

        // `computeEndFromAlleles(nonSymbolicAlleles(allelesList), start, start)`, whose first
        // non-symbolic allele is always the reference.
        let stop = start + reference.len().max(1) as i64 - 1;
        let mut merged = VariantContext::new(contig, start, alleles);
        merged.stop = stop;
        merged.id = if ids.is_empty() {
            ".".to_string()
        } else {
            ids.join(",")
        };
        merged.log10_p_error = NO_LOG10_PERROR;
        merged.attributes = attributes;
        merged.genotypes = genotypes.into();
        merged.filters = None;
        if self.somatic {
            if aggregated_filters.is_empty() || saw_pass_sample {
                merged.filters = Some(Vec::new());
            } else {
                merged.filters = Some(aggregated_filters);
            }
        }
        Ok(Some(merged))
    }

    /// `addReferenceConfidenceAttributes`: each INFO field of a record that starts here, collected
    /// either as raw data for a requested reducible annotation or as numbers for the median.
    fn add_attributes(
        &self,
        pair: &Remapped,
        collected: &mut Vec<(String, Collected)>,
        notes: &mut MergeNotes,
    ) {
        for (key, value) in &pair.vc.attributes {
            if SOMATIC_INFO_TO_MOVE.contains(&key.as_str())
                || SOMATIC_INFO_TO_DROP.contains(&key.as_str())
                || MUTECT_FILTERING_INFO.contains(&key.as_str())
            {
                continue;
            }
            let slot = match collected.iter().position(|(existing, _)| existing == key) {
                Some(index) => index,
                None => {
                    collected.push((
                        key.clone(),
                        if self.is_reducible_raw_key(key) {
                            Collected::Reducible(Vec::new())
                        } else {
                            Collected::Numbers(Vec::new())
                        },
                    ));
                    collected.len() - 1
                }
            };
            match &mut collected[slot].1 {
                Collected::Reducible(values) => {
                    values.push((pair.new_alleles.clone(), attribute_as_list(value).join(",")));
                }
                Collected::Numbers(values) => {
                    let text = value_to_string(value);
                    if !text.contains(',') {
                        match parse_numeric(self.header, key, &text) {
                            Some(number) => values.push(number),
                            None => notes.invalid_annotation = true,
                        }
                    } else {
                        // `split("\\[|,|\\]")`, which drops trailing empty strings.
                        let mut pieces: Vec<&str> = text.split(['[', ',', ']']).collect();
                        while pieces.last().is_some_and(|piece| piece.is_empty()) {
                            pieces.pop();
                        }
                        for piece in pieces.into_iter().filter(|piece| !piece.is_empty()) {
                            match parse_numeric(self.header, key, piece) {
                                Some(number) => values.push(number),
                                None => {
                                    notes.invalid_annotation = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// `mergeAttributes`: the reducible annotations combined, every other number replaced by its
    /// median, `DP` the summed depth, and the stale keys removed.
    fn merge_attributes(
        &self,
        depth: i64,
        alleles: &[Allele],
        mut collected: Vec<(String, Collected)>,
    ) -> Result<Vec<(String, Value)>, MergeError> {
        let mut attributes: Vec<(String, Value)> = Vec::new();
        let put = |attributes: &mut Vec<(String, Value)>, key: &str, value: Value| {
            if let Some(slot) = attributes.iter_mut().find(|(existing, _)| existing == key) {
                slot.1 = value;
            } else {
                attributes.push((key.to_string(), value));
            }
        };

        // `combineAnnotations`, over the requested info annotations in the engine's order.
        for entry in &self.annotations {
            if entry.kind != gatk_annotation::catalogue::Kind::Info {
                continue;
            }
            let Some(raw_keys) = entry.raw_keys else {
                continue;
            };
            for raw_key in raw_keys {
                let Some(position) = collected.iter().position(|(key, _)| key == raw_key) else {
                    continue;
                };
                let raw: Vec<(Vec<Allele>, String)> = match &collected[position].1 {
                    Collected::Reducible(values) => values.clone(),
                    Collected::Numbers(_) => Vec::new(),
                };
                for (key, value) in combine_raw_data(entry.name, alleles, &raw)? {
                    put(&mut attributes, &key, value);
                }
                collected.retain(|(key, _)| !raw_keys.contains(&key.as_str()));
            }
        }

        for (key, values) in collected {
            let Collected::Numbers(mut values) = values else {
                continue;
            };
            if values.is_empty() {
                continue;
            }
            let chosen = if values.len() == 1 {
                values[0]
            } else {
                let integers = values
                    .iter()
                    .filter(|value| matches!(value, Numeric::Integer(_)))
                    .count();
                if integers != 0 && integers != values.len() {
                    return Err(MergeError::MixedMedian { key });
                }
                values.sort_by(|a, b| match (a, b) {
                    (Numeric::Integer(x), Numeric::Integer(y)) => x.cmp(y),
                    (Numeric::Double(x), Numeric::Double(y)) => compare_doubles(*x, *y),
                    _ => std::cmp::Ordering::Equal,
                });
                values[values.len() / 2]
            };
            put(
                &mut attributes,
                &key,
                match chosen {
                    Numeric::Integer(number) => Value::Int(number as i64),
                    Numeric::Double(number) => Value::Double(number),
                },
            );
        }

        if depth > 0 {
            put(&mut attributes, "DP", Value::Str(depth.to_string()));
        }
        // `removeStaleAttributesAfterMerge`.
        for stale in ["AC", "AF", "AN", "MLEAC", "MLEAF", "END", "ECNT", "ECNTH"] {
            attributes.retain(|(key, _)| key != stale);
        }
        Ok(attributes)
    }

    /// `mergeRefConfidenceGenotypes`: each genotype's likelihoods and depths moved onto the merged
    /// alleles, and its call made again.
    fn merge_genotypes(
        &self,
        vc: &VariantContext,
        remapped: &[Allele],
        targets: &[Allele],
    ) -> Result<Vec<Genotype>, MergeError> {
        let mut merged = Vec::with_capacity(vc.genotypes.len());
        for original in vc.genotypes.iter() {
            let ploidy = original.ploidy();
            let mut genotype = original.clone();
            if !self.somatic {
                if original.pl.is_some() || original.ad.is_some() {
                    let relevant = indexes_of_relevant_alleles_for_gvcf(
                        remapped, targets, vc.start, original,
                    )?;
                    if let Some(pl) = &original.pl {
                        let map = subsetted_pl_indices(ploidy, &relevant).map_err(|_| {
                            MergeError::IndexOutOfBounds {
                                index: relevant.len(),
                                length: pl.len(),
                            }
                        })?;
                        let mut pls = Vec::with_capacity(map.len());
                        for index in map {
                            pls.push(*pl.get(index).ok_or(MergeError::IndexOutOfBounds {
                                index,
                                length: pl.len(),
                            })?);
                        }
                        genotype.pl = Some(pls);
                    }
                    if let Some(ad) = &original.ad {
                        genotype.ad = Some(remap_list(ad, &relevant, 0, 0));
                    }
                }
                if exclude_from_annotations(original) {
                    genotype.alleles = vec![Allele::no_call(); ploidy];
                }
            } else {
                genotype.extended.clear();
                if let Some(dp) = original.dp {
                    genotype.dp = Some(dp);
                }
                for key in SOMATIC_FORMAT_TO_KEEP {
                    if let Some(value) = original.get(key) {
                        set_extended(&mut genotype, key, value.clone());
                    }
                }
                let relevant =
                    indexes_of_relevant_alleles_for_gvcf(remapped, targets, vc.start, original)?;
                if let Some(ad) = &original.ad {
                    genotype.ad = Some(remap_list(ad, &relevant, 0, 0));
                } else if let Some(dp) = original.dp {
                    let mut ad = vec![0; targets.len()];
                    ad[0] = dp;
                    genotype.ad = Some(ad);
                }
                if let Some(af) = original.get("AF") {
                    let fractions = attribute_to_doubles(af);
                    let remapped_fractions = remap_list(&fractions, &relevant, 1, 0.0);
                    set_extended(
                        &mut genotype,
                        "AF",
                        Value::List(remapped_fractions.into_iter().map(Value::Double).collect()),
                    );
                } else if (is_hom_ref(original) || is_no_call(original))
                    && vc.alternate_alleles().len() == 1
                {
                    set_extended(
                        &mut genotype,
                        "AF",
                        Value::List(vec![Value::Double(0.0); targets.len() - 1]),
                    );
                }
                for key in SOMATIC_INFO_TO_MOVE {
                    self.move_somatic_attribute(vc, &relevant, original, &mut genotype, key);
                }
                if !self.drop_somatic_filtering_annotations {
                    for key in MUTECT_FILTERING_INFO {
                        self.move_somatic_attribute(vc, &relevant, original, &mut genotype, key);
                    }
                }
                if vc.filters_were_applied() && vc.genotypes.len() == 1 && !is_hom_ref(original) {
                    let filters = vc.filters.clone().unwrap_or_default();
                    if !filters.is_empty() {
                        genotype.filters = Some(filters.join(";"));
                    }
                }
            }

            let method = if self.call_genotypes && should_be_called(original) {
                Method::BestMatchToOriginal
            } else {
                Method::SetToNoCall
            };
            make_genotype_call(
                ploidy,
                &mut genotype,
                method,
                original.pl.as_deref(),
                targets,
                original,
            );
            merged.push(genotype);
        }
        Ok(merged)
    }

    /// `setPerSampleSomaticAttributes`: a genotype's own value remapped when it has one, the
    /// record's INFO value remapped and moved to the genotype when it does not.
    fn move_somatic_attribute(
        &self,
        vc: &VariantContext,
        relevant: &[usize],
        original: &Genotype,
        genotype: &mut Genotype,
        key: &str,
    ) {
        let (source, cardinality) = if let Some(value) = original.get(key) {
            (value.clone(), format_cardinality(self.header, key))
        } else if let Some((_, value)) = vc.attributes.iter().find(|(k, _)| k == key) {
            (value.clone(), info_cardinality(self.header, key))
        } else {
            return;
        };
        let list = attribute_to_list(&source);
        let values = match cardinality {
            Some(Cardinality::A) => remap_optional_list(&list, relevant, 1),
            Some(Cardinality::R) => remap_optional_list(&list, relevant, 0),
            _ => list.into_iter().map(Some).collect(),
        };
        set_extended(
            genotype,
            key,
            Value::List(
                values
                    .into_iter()
                    .map(|value| value.map(Value::Str).unwrap_or(Value::Missing))
                    .collect(),
            ),
        );
    }
}

/// `GenotypeBuilder.attribute(key, value)`: replaces a value already set under the key.
fn set_extended(genotype: &mut Genotype, key: &str, value: Value) {
    if let Some(slot) = genotype.extended.iter_mut().find(|(k, _)| k == key) {
        slot.1 = value;
    } else {
        genotype.extended.push((key.to_string(), value));
    }
}

/// `VariantContextGetters.attributeToList`: a string is split on commas, a list is its elements.
fn attribute_to_list(value: &Value) -> Vec<String> {
    match value {
        Value::List(items) => items.iter().map(value_to_string).collect(),
        Value::Str(text) => text.split(',').map(str::to_string).collect(),
        other => vec![value_to_string(other)],
    }
}

/// `getAttributeAsDoubleArray(g, key, () -> new double[]{0.0}, 0.0)` over a string attribute.
fn attribute_to_doubles(value: &Value) -> Vec<f64> {
    let text = value_to_string(value);
    let cleaned: String = text
        .trim()
        .chars()
        .filter(|c| *c != '[' && *c != ']')
        .collect();
    cleaned
        .split(',')
        .map(|piece| {
            if piece == "." {
                0.0
            } else {
                htsjdk_vcf::genotype_likelihoods::parse_java_double(piece).unwrap_or(0.0)
            }
        })
        .collect()
}

/// `AlleleSubsettingUtils.remapList`: `offset` 0 for a per-allele list, 1 for a per-alternate one,
/// and `filler` wherever the source allele has no entry.
fn remap_list<T: Clone>(original: &[T], relevant: &[usize], offset: usize, filler: T) -> Vec<T> {
    (offset..relevant.len())
        .map(|i| {
            let old = relevant[i];
            if old >= original.len() + offset {
                filler.clone()
            } else {
                original[old - offset].clone()
            }
        })
        .collect()
}

/// [`remap_list`] with the null filler `generateAnnotationValueVector` passes.
fn remap_optional_list(
    original: &[String],
    relevant: &[usize],
    offset: usize,
) -> Vec<Option<String>> {
    let wrapped: Vec<Option<String>> = original.iter().cloned().map(Some).collect();
    remap_list(&wrapped, relevant, offset, None)
}

/// `AlleleSubsettingUtils.getIndexesOfRelevantAllelesForGVCF`, with `doSomaticMerge` false as both
/// of the merger's calls pass it.
fn indexes_of_relevant_alleles_for_gvcf(
    remapped: &[Allele],
    targets: &[Allele],
    position: i64,
    genotype: &Genotype,
) -> Result<Vec<usize>, MergeError> {
    let non_ref = non_ref();
    let Some(index_of_non_ref) = remapped.iter().position(|allele| *allele == non_ref) else {
        return Err(MergeError::MissingNonRef { position });
    };
    let span_del = span_del();
    let mut mapping = vec![0usize; targets.len()];
    for i in 1..targets.len() {
        if targets[i] == span_del && genotype.pl.is_some() {
            let occurrences = remapped
                .iter()
                .filter(|allele| **allele == span_del)
                .count();
            if occurrences > 1 {
                let best = index_of_best_del(
                    remapped,
                    genotype.pl.as_deref().unwrap_or(&[]),
                    genotype.ploidy(),
                )?;
                mapping[i] = best.unwrap_or(index_of_non_ref);
                continue;
            }
        }
        mapping[i] = remapped
            .iter()
            .position(|allele| *allele == targets[i])
            .unwrap_or(index_of_non_ref);
    }
    Ok(mapping)
}

/// `indexOfBestDel`: the `*` whose homozygous genotype has the smallest PL.
fn index_of_best_del(
    alleles: &[Allele],
    pls: &[i32],
    ploidy: usize,
) -> Result<Option<usize>, MergeError> {
    let span_del = span_del();
    let mut best: Option<usize> = None;
    let mut best_pl = i32::MAX;
    for (i, allele) in alleles.iter().enumerate() {
        if *allele != span_del {
            continue;
        }
        let hom_alt = index_of_first_genotype_with_allele(ploidy, i + 1).unwrap_or(0) as i64 - 1;
        let pl = usize::try_from(hom_alt)
            .ok()
            .and_then(|index| pls.get(index).copied())
            .ok_or(MergeError::IndexOutOfBounds {
                index: hom_alt.max(0) as usize,
                length: pls.len(),
            })?;
        if pl < best_pl {
            best = Some(i);
            best_pl = pl;
        }
    }
    Ok(best)
}

/// `collectTargetAlleles`: the merged reference, every record's called alternates in order, and
/// `*` and `<NON_REF>` at the end when they are wanted and not already there.
fn collect_target_alleles(
    remapped: &[Remapped],
    reference: &Allele,
    remove_non_ref: bool,
) -> Vec<Allele> {
    let mut alleles: Vec<Allele> = vec![reference.clone()];
    for pair in remapped {
        for allele in pair.alleles_for_final_set() {
            if !alleles.contains(allele) {
                alleles.push(allele.clone());
            }
        }
    }
    let saw_spanning_deletion = remapped.iter().any(Remapped::is_spanning_deletion);
    let saw_non_spanning_event = remapped.iter().any(Remapped::is_non_spanning_event);
    if saw_spanning_deletion && (saw_non_spanning_event || !remove_non_ref) {
        let span_del = span_del();
        if !alleles.contains(&span_del) {
            alleles.push(span_del);
        }
    }
    if !remove_non_ref {
        let non_ref = non_ref();
        if !alleles.contains(&non_ref) {
            alleles.push(non_ref);
        }
    }
    alleles
}

/// `Genotype.isHomRef()`.
pub fn is_hom_ref(genotype: &Genotype) -> bool {
    genotype.is_hom_ref()
}

/// `Genotype.isNoCall()`.
pub fn is_no_call(genotype: &Genotype) -> bool {
    genotype.is_no_call()
}

/// `GenotypeGVCFsEngine.excludeFromAnnotations`: a hom-ref or no-call with no depth and a GQ of 0.
fn exclude_from_annotations(genotype: &Genotype) -> bool {
    (is_hom_ref(genotype) || is_no_call(genotype))
        && genotype.dp.is_none_or(|dp| dp == 0)
        && genotype.gq == Some(0)
}

/// `Genotype.isNonInformative()`: no PL, or every PL zero.
fn is_non_informative(genotype: &Genotype) -> bool {
    genotype
        .pl
        .as_ref()
        .is_none_or(|pl| pl.iter().all(|value| *value == 0))
}

/// `GenotypeUtils.shouldBeCalled`.
fn should_be_called(genotype: &Genotype) -> bool {
    !is_non_informative(genotype) || genotype.gq.is_some_and(|gq| gq > 0)
}

/// `GenotypeAssignmentMethod`, the three values `CombineGVCFs` reaches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    SetToNoCall,
    BestMatchToOriginal,
    PreferPls,
}

/// `GenotypeLikelihoods.fromPLs(pl).getAsVector()`: each PL over minus ten.
pub(crate) fn likelihoods_of(pls: &[i32]) -> Vec<f64> {
    pls.iter().map(|pl| *pl as f64 / -10.0).collect()
}

/// `MathUtils.maxElementIndex`: the first index of the maximum.
pub(crate) fn max_element_index(values: &[f64]) -> usize {
    let mut best = 0;
    for (index, value) in values.iter().enumerate() {
        if *value > values[best] {
            best = index;
        }
    }
    best
}

/// `GenotypeLikelihoods.getGQLog10FromLikelihoods`.
pub(crate) fn gq_log10_from_likelihoods(chosen: usize, likelihoods: &[f64]) -> f64 {
    let mut other = f64::NEG_INFINITY;
    for (index, value) in likelihoods.iter().enumerate() {
        if index != chosen && *value >= other {
            other = *value;
        }
    }
    let quality = likelihoods[chosen] - other;
    if quality < 0.0 {
        let max = likelihoods
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let total: f64 = likelihoods
            .iter()
            .map(|value| 10f64.powf(value - max))
            .sum();
        let normalised = 10f64.powf(likelihoods[chosen] - max) / total;
        (1.0 - normalised).log10()
    } else {
        -quality
    }
}

/// `Math.round`, which is `floor(x + 0.5)`.
pub(crate) fn java_round(value: f64) -> i64 {
    (value + 0.5).floor() as i64
}

/// `bestMatchToOriginalGT`: each original allele the new list still holds, or a no-call, or the
/// reference in place of anything else.
pub(crate) fn best_match_to_original(targets: &[Allele], original: &[Allele]) -> Vec<Allele> {
    original
        .iter()
        .map(|allele| {
            if targets.contains(allele) || allele.is_no_call() {
                allele.clone()
            } else {
                targets[0].clone()
            }
        })
        .collect()
}

/// `GATKVariantContextUtils.makeGenotypeCall`, for the three methods this tool reaches.
pub fn make_genotype_call(
    ploidy: usize,
    genotype: &mut Genotype,
    method: Method,
    pls: Option<&[i32]>,
    targets: &[Allele],
    original: &Genotype,
) {
    if method != Method::SetToNoCall
        && (is_hom_ref(original) || is_no_call(original))
        && original.gq == Some(0)
    {
        genotype.alleles = vec![Allele::no_call(); ploidy];
        if original.dp == Some(0) {
            genotype.pl = None;
            genotype.dp = None;
            genotype.ad = None;
            genotype.gq = None;
            genotype.extended.clear();
            return;
        }
    }
    match method {
        Method::SetToNoCall => genotype.alleles = vec![Allele::no_call(); ploidy],
        Method::BestMatchToOriginal => {
            let uninformative = original.gq == Some(0)
                && original.pl.as_ref().is_none_or(|pl| pl.first() == Some(&0));
            if uninformative {
                genotype.alleles = vec![Allele::no_call(); ploidy];
            } else {
                genotype.alleles = best_match_to_original(targets, &original.alleles);
            }
        }
        Method::PreferPls => {
            let likelihoods = pls.map(likelihoods_of);
            let informative = likelihoods
                .as_ref()
                .is_some_and(|gls| gls.iter().sum::<f64>() < SUM_GL_THRESH_NOCALL);
            match likelihoods {
                Some(gls) if informative => {
                    let best = max_element_index(&gls);
                    let called: Vec<Allele> = genotypes_in_canonical_order(ploidy, targets.len())
                        .into_iter()
                        .nth(best)
                        .unwrap_or_default()
                        .into_iter()
                        .map(|index| targets[index].clone())
                        .collect();
                    let gq = gq_log10_from_likelihoods(best, &gls);
                    if called.contains(&non_ref()) {
                        genotype.alleles = vec![targets[0].clone(); ploidy];
                        genotype.pl = Some(vec![0; gls.len()]);
                        genotype.gq = Some(0);
                    } else if best == 0 && gq > SUM_GL_THRESH_NOCALL {
                        genotype.alleles = vec![Allele::no_call(); ploidy];
                    } else {
                        genotype.alleles = called;
                    }
                    if targets.len() > 1 {
                        // `gb.log10PError(gq)`, which caps the GQ at 99.
                        genotype.gq = if gq == NO_LOG10_PERROR {
                            None
                        } else {
                            Some(crate::genotyping_engine::gq_of_log10(gq))
                        };
                    }
                }
                _ => genotype.alleles = best_match_to_original(targets, &original.alleles),
            }
        }
    }
}

/// `combineRawData` for one reducible annotation, keyed by its simple name.
///
/// Each raw string is parsed against the alleles ITS record was remapped to, and combined into
/// slots keyed by the merged alleles, so an allele one input has and another lacks is summed only
/// where it exists.
fn combine_raw_data(
    name: &str,
    alleles: &[Allele],
    raw: &[(Vec<Allele>, String)],
) -> Result<Vec<(String, Value)>, MergeError> {
    use gatk_annotation::{
        allele_specific_rank_sum as rank_sum, allele_specific_site_statistics as site,
        allele_specific_strand_bias as strand_bias, mapping_quality, raw_gt_count,
    };
    let annotation_error = |class: &str, message: String| MergeError::Annotation {
        class: class.to_string(),
        message,
    };
    match name {
        "RMSMappingQuality" => {
            let mut tuples = Vec::new();
            for (_, text) in raw {
                tuples.push(
                    mapping_quality::parse_raw_data_string(text).map_err(|error| {
                        annotation_error("java.lang.NumberFormatException", format!("{error:?}"))
                    })?,
                );
            }
            let Some((sum, depth)) = mapping_quality::RmsMappingQuality::combine_raw_data(&tuples)
            else {
                return Ok(Vec::new());
            };
            Ok(vec![(
                "RAW_MQandDP".to_string(),
                Value::Str(mapping_quality::raw_annotation_string(sum, depth)),
            )])
        }
        "RawGtCount" => {
            let texts: Vec<String> = raw.iter().map(|(_, text)| text.clone()).collect();
            let combined = raw_gt_count::combine_raw_data(&texts)
                .map_err(|error| annotation_error(error.class(), error.message()))?;
            Ok(vec![("RAW_GT_COUNT".to_string(), Value::Str(combined))])
        }
        "AS_RMSMappingQuality" => {
            let mut combined: Vec<(Allele, Option<f64>)> = alleles
                .iter()
                .map(|allele| (allele.clone(), None))
                .collect();
            for (source_alleles, text) in raw {
                let parsed = site::as_rms_parse_raw(source_alleles, text).map_err(|error| {
                    annotation_error("java.lang.NumberFormatException", format!("{error:?}"))
                })?;
                for slot in combined.iter_mut() {
                    // `toAdd.getAttribute(currentAllele)`: the LAST entry for an allele the source
                    // lists twice, which is what a `HashMap` put leaves.
                    let add = parsed
                        .iter()
                        .rev()
                        .find(|(allele, _)| *allele == slot.0)
                        .and_then(|(_, value)| *value);
                    if let Some(add) = add {
                        slot.1 = Some(slot.1.map_or(add, |current| current + add));
                    }
                }
            }
            Ok(vec![(
                "AS_RAW_MQ".to_string(),
                Value::Str(site::as_rms_raw_string(alleles, &combined)),
            )])
        }
        "AS_BaseQualityRankSumTest" | "AS_MappingQualityRankSumTest" | "AS_ReadPosRankSumTest" => {
            let key = match name {
                "AS_BaseQualityRankSumTest" => "AS_RAW_BaseQRankSum",
                "AS_MappingQualityRankSumTest" => "AS_RAW_MQRankSum",
                _ => "AS_RAW_ReadPosRankSum",
            };
            let mut combined: Vec<(Allele, gatk_engine::histogram::Histogram)> = alleles
                .iter()
                .map(|allele| (allele.clone(), gatk_engine::histogram::Histogram::new()))
                .collect();
            for (source_alleles, text) in raw {
                let parsed =
                    rank_sum::parse_raw_data_string(source_alleles, text).map_err(|error| {
                        annotation_error("java.lang.IllegalStateException", format!("{error:?}"))
                    })?;
                for slot in combined.iter_mut() {
                    if let Some((_, histogram)) =
                        parsed.iter().rev().find(|(allele, _)| *allele == slot.0)
                    {
                        slot.1.add_histogram(histogram).map_err(|error| {
                            annotation_error(
                                "java.lang.IllegalStateException",
                                format!("{error:?}"),
                            )
                        })?;
                    }
                }
            }
            let text =
                rank_sum::make_combined_annotation_string(alleles, &combined).map_err(|error| {
                    annotation_error("java.lang.IllegalStateException", format!("{error:?}"))
                })?;
            Ok(vec![(key.to_string(), Value::Str(text))])
        }
        "AS_FisherStrand" | "AS_StrandOddsRatio" => {
            let mut combined: Vec<(Allele, Option<Vec<i32>>)> = alleles
                .iter()
                .map(|allele| (allele.clone(), None))
                .collect();
            for (source_alleles, text) in raw {
                let parsed =
                    strand_bias::parse_raw_data_string(source_alleles, text).map_err(|error| {
                        annotation_error("java.lang.IllegalStateException", format!("{error:?}"))
                    })?;
                for slot in combined.iter_mut() {
                    let Some((_, counts)) =
                        parsed.iter().rev().find(|(allele, _)| *allele == slot.0)
                    else {
                        continue;
                    };
                    match &mut slot.1 {
                        Some(existing) if !existing.is_empty() => {
                            if counts.len() >= 2 {
                                existing[0] += counts[0];
                                existing[1] += counts[1];
                            }
                        }
                        _ => {
                            slot.1 = Some(if counts.len() >= 2 {
                                vec![counts[0], counts[1]]
                            } else {
                                Vec::new()
                            });
                        }
                    }
                }
            }
            Ok(vec![(
                "AS_SB_TABLE".to_string(),
                Value::Str(strand_bias::make_raw_annotation_string(&combined)),
            )])
        }
        "AS_QualByDepth" => {
            let mut combined: Vec<(Allele, Option<i64>)> = alleles
                .iter()
                .map(|allele| (allele.clone(), None))
                .collect();
            for (source_alleles, text) in raw {
                let mut tokens: Vec<&str> = text.split('|').collect();
                while tokens.last().is_some_and(|token| token.is_empty()) {
                    tokens.pop();
                }
                let mut parsed: Vec<(Allele, Option<i64>)> = Vec::new();
                for (index, token) in tokens.iter().enumerate() {
                    let allele = source_alleles.get(index).ok_or_else(|| {
                        annotation_error(
                            "java.lang.IndexOutOfBoundsException",
                            format!(
                                "Index {index} out of bounds for length {}",
                                source_alleles.len()
                            ),
                        )
                    })?;
                    let value = if token.is_empty() || *token == "." {
                        None
                    } else {
                        Some(java_parse_int(token).ok_or_else(|| {
                            annotation_error(
                                "java.lang.NumberFormatException",
                                format!("For input string: \"{token}\""),
                            )
                        })? as i64)
                    };
                    parsed.push((allele.clone(), value));
                }
                for slot in combined.iter_mut() {
                    let add = parsed
                        .iter()
                        .rev()
                        .find(|(allele, _)| *allele == slot.0)
                        .and_then(|(_, value)| *value);
                    if let Some(add) = add {
                        slot.1 = Some((slot.1.unwrap_or(0) as i32).wrapping_add(add as i32) as i64);
                    }
                }
            }
            let text = combined
                .iter()
                .map(|(_, value)| value.unwrap_or(0).to_string())
                .collect::<Vec<_>>()
                .join("|");
            Ok(vec![("AS_QUALapprox".to_string(), Value::Str(text))])
        }
        _ => Ok(Vec::new()),
    }
}
