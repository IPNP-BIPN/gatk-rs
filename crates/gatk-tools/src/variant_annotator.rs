//! `VariantAnnotator`: how one VCF is annotated from another.
//!
//! A resource file is tagged with a name, an expression names one of its fields, and the value
//! lands on any record at the same position. What is ported is which records are annotated, under
//! what key, and what happens when the two files disagree about the alleles.
//!
//! Reading the VCFs and the read-based annotations are not ported. The expression machinery is.

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
    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    /// The alternates two records share, which is what a per-allele value has to be mapped through.
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
    pub fn parse(text: &str) -> Option<Expression> {
        let (tag, field) = text.split_once('.')?;
        Some(Expression {
            tag: tag.to_string(),
            field: field.to_string(),
        })
    }

    pub fn key(&self) -> String {
        format!("{}.{}", self.tag, self.field)
    }
}

/// The three fields an expression may name that are not INFO attributes.
pub const ID_FIELD: &str = "ID";
pub const ALT_FIELD: &str = "ALT";
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
    pub key: String,
    /// `None` for a `--comp` flag, which is a bare key.
    pub value: Option<String>,
}

/// The annotations one input record collects.
///
/// Every annotation the input already carried is kept: these are added to them, never instead of
/// them. A record with no resource at its position collects nothing rather than an empty value.
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
