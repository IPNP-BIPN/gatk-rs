//! The driving variants of a `MultiVariantWalker`: every `--variant` opened, cross-checked and
//! merged into one header and one stream, by [`gatk_tools::multi_variant`].
//!
//! What is here is the half that touches files: reading the inputs, their Tribble indexes for a
//! header with no `##contig` line, and turning the engine's refusals into what the dispatcher
//! reports. Every multi-input tool opens its inputs through [`open`] and walks them through
//! [`DrivingVariants::merged`].

use gatk_barclay::{Parser, Value};
use gatk_tools::main_entry::{Failure, Thrown};
use gatk_tools::multi_variant::{self as engine, FeatureInput, Merged, Refusal};
use htsjdk_bam::header::SamHeader;
use htsjdk_vcf::variant::VariantContext;

/// One opened `--variant`.
pub struct Driving {
    pub input: FeatureInput,
    /// The path as given on the command line.
    pub path: String,
    pub codec: gatk_tools::feature_codec::Codec,
    /// The file's text, decompressed if it was block compressed.
    pub text: String,
    pub file: htsjdk_vcf::reader::VcfFile,
}

/// `MultiVariantDataSource` once constructed.
pub struct DrivingVariants {
    pub inputs: Vec<Driving>,
    pub merged: Merged,
}

/// An engine refusal as the dispatcher reports it.
pub fn thrown(refusal: Refusal) -> Thrown {
    let exception: &'static str = match &refusal {
        Refusal::BadInput(_) => "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        Refusal::Incompatible(_) => {
            "org.broadinstitute.hellbender.exceptions.UserException$IncompatibleSequenceDictionaries"
        }
        Refusal::NoDictionary => "org.broadinstitute.hellbender.exceptions.UserException",
        Refusal::DuplicateName { .. } => "java.lang.IllegalStateException",
        Refusal::Other { class, .. } => super::java_class_name(class),
    };
    Thrown {
        failure: if refusal.is_user() {
            Failure::User
        } else {
            Failure::Other
        },
        exception,
        message: Some(refusal.message()),
    }
}

/// `--variant`'s values as `FeatureInput`s, tags included.
pub fn feature_inputs(parser: &Parser) -> Vec<FeatureInput> {
    let one = |value: &Value| -> Option<FeatureInput> {
        let (raw, tag) = match value {
            Value::Tagged { value, tag, .. } => (value.clone(), tag.clone()),
            Value::Str(text) => (text.clone(), None),
            _ => return None,
        };
        Some(FeatureInput {
            absolute: super::absolute_path(&raw),
            tag,
            raw,
        })
    };
    parser
        .definitions()
        .iter()
        .find(|definition| definition.long_name() == "variant")
        .map(|definition| match &definition.value {
            Value::List(values) => values.iter().filter_map(one).collect(),
            other => one(other).into_iter().collect(),
        })
        .unwrap_or_default()
}

/// `initializeDrivingVariants` and the data source's constructor: the inputs checked for
/// repeats, opened, their dictionaries cross-checked unless `cross_validate` is off, and merged.
pub fn open(inputs: Vec<FeatureInput>, cross_validate: bool) -> Result<DrivingVariants, Thrown> {
    engine::require_unique(&inputs).map_err(thrown)?;
    open_unchecked(inputs, cross_validate)
}

/// [`open`] without the repeat check, which a tool that fills the driving inputs itself instead of
/// through `--variant` (VariantEval) never reaches.
pub fn open_unchecked(
    inputs: Vec<FeatureInput>,
    cross_validate: bool,
) -> Result<DrivingVariants, Thrown> {
    let mut opened = Vec::with_capacity(inputs.len());
    for input in inputs {
        let path = input.raw.clone();
        let (codec, text) = super::open_feature_input(&path)?;
        let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })?;
        opened.push(Driving {
            input,
            path,
            codec,
            text,
            file,
        });
    }
    let versions: Vec<Option<String>> = opened
        .iter()
        .map(|driving| {
            driving
                .file
                .header_version
                .map(|version| version.version_string().to_string())
        })
        .collect();
    let sources: Vec<engine::Source> = opened
        .iter()
        .zip(&versions)
        .map(|(driving, version)| engine::Source {
            name: driving.input.name(),
            header: &driving.file.header,
            version: version.as_deref(),
            index_contigs: super::index_dictionary(&driving.path).map(|dictionary| {
                dictionary
                    .sequences
                    .into_iter()
                    .map(|sequence| sequence.name)
                    .collect()
            }),
        })
        .collect();
    if cross_validate {
        let dictionaries: Vec<(String, Option<Vec<engine::Contig>>)> = opened
            .iter()
            .zip(&sources)
            .map(|(driving, source)| (driving.input.name().to_string(), source.dictionary()))
            .collect();
        engine::cross_validate(&dictionaries).map_err(thrown)?;
    }
    let merged = engine::merge(&sources).map_err(|refusal| match refusal {
        Refusal::DuplicateName {
            name,
            first,
            second,
        } => Thrown::non_user(
            "java.lang.IllegalStateException",
            format!(
                "Duplicate key {name} (attempted merging values {} and {})",
                super::java_vcf_header_string(&opened[first].text),
                super::java_vcf_header_string(&opened[second].text)
            ),
        ),
        other => thrown(other),
    })?;
    Ok(DrivingVariants {
        inputs: opened,
        merged,
    })
}

/// [`open`] over the command line's `--variant`, with the cross-check the tool's
/// `doDictionaryCrossValidation` and `--disable-sequence-dictionary-validation` leave on.
pub fn open_from(parser: &Parser, tool_cross_validates: bool) -> Result<DrivingVariants, Thrown> {
    let inputs = feature_inputs(parser);
    if inputs.is_empty() {
        return Err(Thrown::command_line(
            "Argument variant was missing: Argument 'variant' is required",
        ));
    }
    let cross_validate =
        tool_cross_validates && !super::flag(parser, "disable-sequence-dictionary-validation");
    open(inputs, cross_validate)
}

impl DrivingVariants {
    /// The merged dictionary as a SAM header, for the validation every variant walker runs
    /// against the reference and the reads. Empty when no input carries a `##contig` line, which
    /// is what that validation reads as "derive one from the index".
    pub fn dictionary(&self) -> SamHeader {
        let any_own = self.inputs.iter().any(|driving| {
            driving
                .file
                .header
                .lines
                .iter()
                .any(|line| matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }))
        });
        let mut header = SamHeader::default();
        if any_own {
            for contig in self.merged.dictionary.iter().flatten() {
                header
                    .sequences
                    .push(htsjdk_bam::header::SequenceRecord::new(
                        &contig.name,
                        contig.length as i32,
                    ));
            }
        }
        header
    }

    /// A variant walker's startup validation over the merged inputs: the first input stands for
    /// the driving file and the merged dictionary for its dictionary.
    pub(super) fn startup(&self, parser: &Parser) -> Result<super::VariantWalkerStart, Thrown> {
        let first = &self.inputs[0];
        super::variant_walker_validation(
            parser,
            first.path.clone(),
            first.text.clone(),
            first.codec,
            self.dictionary(),
        )
    }

    /// The merged header, samples included, as a writer starts from it.
    pub fn header(&self) -> htsjdk_vcf::header::VcfHeader {
        htsjdk_vcf::header::VcfHeader {
            lines: self.header_lines(),
            samples: self.merged.samples.clone(),
        }
    }

    /// The merged header's lines as a writer emits them.
    pub fn header_lines(&self) -> Vec<htsjdk_vcf::header::HeaderLine> {
        engine::written_lines(&self.merged.lines)
    }

    /// Each input's records the traversal reaches, input by input: each restricted to `-L`, as
    /// each data source is. `-L` queries every input by interval, so each needs an index before a
    /// record is read.
    pub fn reached(
        &self,
        intervals: Option<&[gatk_engine::interval::SimpleInterval]>,
    ) -> Result<Vec<Vec<&VariantContext>>, Thrown> {
        let spans = gatk_engine::variant_source::intervals_for_traversal(intervals);
        if spans.is_some() {
            if let Some(driving) = self
                .inputs
                .iter()
                .find(|driving| !super::has_feature_index(&driving.path))
            {
                return Err(Thrown::user(
                    gatk_tools::count_variants::CountVariantsError::IntervalsWithoutRandomAccess {
                        path: driving.path.clone(),
                    }
                    .message(),
                ));
            }
        }
        Ok(self
            .inputs
            .iter()
            .map(|driving| {
                let loci: Vec<Spanned> = driving.file.records.iter().map(Spanned).collect();
                gatk_engine::variant_source::traverse(&loci, spans)
                    .into_iter()
                    .map(|spanned| spanned.0)
                    .collect()
            })
            .collect())
    }

    /// Every record the traversal reaches, as `(input, record)`, in the merge's order: the
    /// [`reached`](Self::reached) records of every input, merged. A refusal from the merge
    /// arrives after the records handed out before it.
    #[allow(clippy::type_complexity)]
    pub fn merged(
        &self,
        intervals: Option<&[gatk_engine::interval::SimpleInterval]>,
    ) -> Result<(Vec<(usize, &VariantContext)>, Option<Thrown>), Thrown> {
        let reached = self.reached(intervals)?;
        let contigs = self.merged.contigs();
        let keys: Vec<Vec<engine::Key>> = reached
            .iter()
            .map(|records| {
                records
                    .iter()
                    .map(|record| engine::key_of(&contigs, &record.contig, record.start))
                    .collect()
            })
            .collect();
        let merge = engine::merging_order(&keys);
        let order = merge
            .order
            .into_iter()
            .map(|(input, index)| (input, reached[input][index]))
            .collect();
        Ok((order, merge.failure.map(thrown)))
    }

    /// `MultiVariantWalker.traverse`'s wrapper around a failure inside `apply`.
    pub fn apply_failure(&self, input: usize, record: &VariantContext) -> Thrown {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            engine::apply_failure_message(
                &record.contig,
                record.start,
                &super::java_variant_context_string(record, self.inputs[input].input.name()),
            ),
        )
    }
}

/// A record as the per-input traversal sees it: its span, `END` included.
struct Spanned<'a>(&'a VariantContext);

impl gatk_engine::variant_source::Located for Spanned<'_> {
    fn contig(&self) -> &str {
        &self.0.contig
    }
    fn start(&self) -> i32 {
        self.0.start as i32
    }
    fn stop(&self) -> i32 {
        self.0.stop as i32
    }
}
