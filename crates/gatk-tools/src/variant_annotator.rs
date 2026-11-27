//! `VariantAnnotator`: how one VCF is annotated from another.
//!
//! A resource file is tagged with a name, an expression names one of its fields, and the value
//! lands on any record at the same position. What is ported is which records are annotated, under
//! what key, and what happens when the two files disagree about the alleles.
//!
//! Reading the VCFs and the read-based annotations are not ported. The expression machinery is.
//!
//! The reference is `VariantAnnotatorEngine` (GATK 4.7.0.0): `VAExpression` parses an
//! expression, `annotateExpressions` copies the values, and `VariantOverlapAnnotator` adds the
//! `--comp` flags. The model here is narrower than that Java in the places each item notes, and
//! the conformance suite `tests/variant_annotator.rs` measures it against its golden.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing", Java's nullable reference checked by the compiler. A
//!   `?` after an `Option` returns `None` from the function at once when there is nothing.
//! - `let Some(x) = y else { continue; };` takes the value out of an `Option`, or skips to the
//!   next loop iteration when there is none. Java has no single construct for it.
//! - `match` must cover every case; `ID_FIELD if ...` is a case with an extra condition.
//! - `impl Fn(&str) -> Arity` is a parameter that accepts any function or lambda from a field
//!   name to its arity, as a Java `Function<String, Arity>` would.

/// One record, reduced to what an expression reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub contig: String,
    pub position: i32,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    pub filters: Vec<String>,
    pub attributes: Vec<(String, String)>,
}

impl Record {
    /// The first INFO value under `key`, or `None` when the record does not carry it.
    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    /// The alternates two records share, which is what a per-allele value has to be mapped through.
    ///
    /// This is equality of the whole ALT list, in order. The Java instead splits both records
    /// into minimal-representation biallelics (`getMinRepresentationBiallelics`) and matches
    /// them one by one.
    pub fn shares_alternates(&self, other: &Record) -> bool {
        self.alternates == other.alternates
    }
}

/// How many values a field carries, which is what decides whether it can cross to a different
/// alternate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arity {
    /// `Number=1` or any fixed count: the value means the same whatever the alternate is.
    Scalar,
    /// `Number=A`: one value per alternate, so it cannot be carried to a different one.
    PerAllele,
}

/// One `--expression`, which is a resource TAG and a field name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Expression {
    pub tag: String,
    pub field: String,
}

impl Expression {
    /// `tag.field`, which is both how the expression is written and the key it produces.
    ///
    /// `split_once('.')` cuts at the first dot. `VAExpression` cuts at the last:
    ///
    /// ```java
    /// final int indexOfDot = fullExpression.lastIndexOf(".");
    /// ```
    ///
    /// so the two agree only while neither the tag nor the field contains a dot. No dot at all is
    /// `None` here and `UserException.BadInput` in the Java.
    pub fn parse(text: &str) -> Option<Expression> {
        let (tag, field) = text.split_once('.')?;
        Some(Expression {
            tag: tag.to_string(),
            field: field.to_string(),
        })
    }

    /// The annotation key, which is the expression as written (`expression.fullName`).
    pub fn key(&self) -> String {
        format!("{}.{}", self.tag, self.field)
    }
}

/// The three fields an expression may name that are not INFO attributes.
pub const ID_FIELD: &str = "ID";
/// The ALT column, special-cased by `annotateExpressions` as `ID` is.
pub const ALT_FIELD: &str = "ALT";
/// The FILTER column, special-cased by `annotateExpressions` as `ID` is.
pub const FILTER_FIELD: &str = "FILTER";

/// The VCF's own text for a field that has no value.
pub const MISSING: &str = ".";

/// The value one expression takes off a resource record.
///
/// `ID`, `ALT` and `FILTER` are read off the record itself; anything else is an INFO attribute, and
/// a field the resource does not carry yields nothing at all rather than a refusal.
///
/// A record whose ID is the VCF's missing marker yields NOTHING rather than an annotation whose
/// value is a dot: the tool asks whether the record HAS an id, not what its text is.
///
/// ```java
/// if (expression.fieldName.equals("ID")) {
///     if (expressionVC.hasID()) {
///         attributes.put(expression.fullName, expressionVC.getID());
///     }
/// } else if (expression.fieldName.equals("ALT")) {
///     attributes.put(expression.fullName, expressionVC.getAlternateAllele(0).getDisplayString());
/// ```
///
/// The Java's `ALT` is the first alternate only, and its `FILTER` is
/// `getFilters().stream().collect(Collectors.joining(","))`, or `"PASS"` for an unfiltered
/// record. This port joins every alternate with `,` and the filters with `;`.
pub fn value_of(record: &Record, field: &str) -> Option<String> {
    match field {
        ID_FIELD if record.id == MISSING => None,
        ID_FIELD => Some(record.id.clone()),
        ALT_FIELD => Some(record.alternates.join(",")),
        FILTER_FIELD => Some(record.filters.join(";")),
        other => record.attribute(other).map(str::to_string),
    }
}

/// Whether one expression's value reaches one input record.
///
/// A PER-ALLELE value cannot cross to a different alternate and is withheld whatever the arguments
/// say. A SCALAR one crosses freely unless `--resource-allele-concordance` is given, and so do the
/// three record fields, whose meaning does not depend on an alternate either.
///
/// In `annotateExpressions` the concordance test sits only in the INFO branch, gated by
/// `if ((useAltAlleles || useRefAndAltAlleles) || expressionAlleleConcordance)`; the `ID`,
/// `ALT` and `FILTER` branches before it take the value with no allele test.
pub fn annotates(
    input: &Record,
    resource: &Record,
    field: &str,
    arity: Arity,
    allele_concordance: bool,
) -> bool {
    let concordant = input.shares_alternates(resource);
    match field {
        ID_FIELD | ALT_FIELD | FILTER_FIELD => concordant || !allele_concordance,
        _ => match arity {
            Arity::PerAllele => concordant,
            Arity::Scalar => concordant || !allele_concordance,
        },
    }
}

/// One annotation the tool adds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    /// `tag.field` for an expression, or the tag alone for `--comp`.
    pub key: String,
    /// `None` for a `--comp` flag, which is a bare key.
    pub value: Option<String>,
}

/// The annotations one input record collects.
///
/// Every annotation the input already carried is kept: these are added to them, never instead of
/// them. A record with no resource at its position collects nothing rather than an empty value.
///
/// `resources` pairs each resource record with its tag. The Java takes, per expression,
/// `features.getValues(expression.binding, vc.getStart())` and then its first element; the
/// `find` below is that first record with the expression's tag at the input's position.
pub fn annotate(
    input: &Record,
    resources: &[(String, &Record)],
    expressions: &[Expression],
    arity: impl Fn(&str) -> Arity,
    allele_concordance: bool,
) -> Vec<Annotation> {
    let mut out = Vec::new();
    for expression in expressions {
        let Some((_, resource)) = resources
            .iter()
            .find(|(tag, record)| *tag == expression.tag && record.position == input.position)
        else {
            continue;
        };
        if !annotates(
            input,
            resource,
            &expression.field,
            arity(&expression.field),
            allele_concordance,
        ) {
            continue;
        }
        let Some(value) = value_of(resource, &expression.field) else {
            continue;
        };
        out.push(Annotation {
            key: expression.key(),
            value: Some(value),
        });
    }
    out
}

/// `--comp`, which adds a BARE FLAG under the tag alone rather than a key and a value, and which
/// is subject to the same allele test as a scalar expression.
///
/// `VariantOverlapAnnotator.annotateOverlap` puts `attribute(attributeKey, true)`, which the VCF
/// writer prints as a flag. Its own allele test matches any one alternate of the input against
/// the comparison's; the whole-list equality here is narrower.
pub fn comparison_annotation(
    input: &Record,
    tag: &str,
    comparison: Option<&Record>,
) -> Option<Annotation> {
    let comparison = comparison?;
    if comparison.position != input.position || !input.shares_alternates(comparison) {
        return None;
    }
    Some(Annotation {
        key: tag.to_string(),
        value: None,
    })
}

/// The one thing the expression machinery refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnotatorError {
    /// A tag that names no resource. An unknown FIELD is silent; an unknown FILE is not.
    UnknownResource { expression: String },
}

impl AnnotatorError {
    /// The `UserException.BadInput` text from the `VAExpression` constructor, `Bad input: ` being
    /// the prefix that class adds.
    pub fn message(&self) -> String {
        match self {
            AnnotatorError::UnknownResource { expression } => format!(
                "Bad input: The requested expression '{expression}' is invalid, could not find vcf \
                 input file"
            ),
        }
    }
}

/// `addExpressions`' own check, which runs before a record is read.
///
/// The check itself is in the `VAExpression` constructor that `addExpressions` calls once per
/// expression, so the first unknown tag in iteration order is the one refused.
pub fn check_expressions(
    expressions: &[Expression],
    tags: &[String],
) -> Result<(), AnnotatorError> {
    for expression in expressions {
        if !tags.contains(&expression.tag) {
            return Err(AnnotatorError::UnknownResource {
                expression: expression.key(),
            });
        }
    }
    Ok(())
}
