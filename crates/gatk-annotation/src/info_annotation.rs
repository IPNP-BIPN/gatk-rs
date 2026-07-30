//! The annotation interface, ported from `InfoFieldAnnotation` and its parent `VariantAnnotation`.
//!
//! ```java
//! public interface InfoFieldAnnotation extends VariantAnnotation {
//!     Map<String, Object> annotate(ReferenceContext ref, VariantContext vc,
//!                                  AlleleLikelihoods<GATKRead, Allele> likelihoods);
//!     List<String> getKeyNames();
//! }
//! ```
//!
//! Two shapes of that signature are load-bearing.
//!
//! `Map<String, Object>` means the **Java type** of each value travels with it. `Coverage` puts a
//! `String`, `CountNs` puts a `Long`, `ChromosomeCounts` puts an `Integer` for one alternate allele
//! and an `ArrayList` for two or more. A port that returned `Map<String, String>` would agree on
//! every ordinary record and disagree wherever the encoder treats a list differently from a scalar,
//! so [`AnnotationValue`] keeps the distinction.
//!
//! And an empty map is not "all zeroes": it means the keys are **absent** from the record. Every
//! annotation here has at least one guard reaching that branch.
//!
//! # The likelihoods argument
//!
//! `null` likelihoods are a normal input, not an error: three of the annotations here answer an
//! empty map for them, and each one tests something slightly different. `Coverage` tests
//! `likelihoods == null || likelihoods.evidenceCount() == 0`; `MappingQualityZero` tests
//! `!vc.isVariant() || likelihoods == null` and **not** the evidence count, so an empty matrix
//! makes it write a zero where `Coverage` writes nothing at all; `CountNs` tests only the null.

use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;

/// A value an annotation puts into the INFO map, keeping the Java type it was boxed as.
///
/// The variants are the ones the ported annotations actually produce. They are distinguished
/// because the encoder is: an `Integer` and a one-element `ArrayList` are written the same way but
/// reached by different code, and a `String` that happens to hold digits is not an `Integer` to a
/// consumer that fetched it.
#[derive(Debug, Clone, PartialEq)]
pub enum AnnotationValue {
    /// A boxed `Integer`.
    Int(i32),
    /// A boxed `Long`. `CountNs` produces one, and it is not an `Integer`.
    Long(i64),
    /// A boxed `Double`.
    Double(f64),
    /// A `String`, including the ones built by `String.format("%d", ...)`, which are strings and
    /// not numbers however they look.
    Str(String),
    /// A boxed `Boolean`. `TandemRepeat` puts `true` under a `Flag` key.
    Flag(bool),
    /// An `ArrayList`, which is what an annotation produces once there is more than one alternate
    /// allele to report on.
    List(Vec<AnnotationValue>),
}

impl AnnotationValue {
    /// The Java class name of the boxed value, which is what an oracle dump can report and a port
    /// can therefore be held to.
    pub fn java_class(&self) -> &'static str {
        match self {
            AnnotationValue::Int(_) => "java.lang.Integer",
            AnnotationValue::Long(_) => "java.lang.Long",
            AnnotationValue::Double(_) => "java.lang.Double",
            AnnotationValue::Str(_) => "java.lang.String",
            AnnotationValue::Flag(_) => "java.lang.Boolean",
            AnnotationValue::List(_) => "java.util.ArrayList",
        }
    }

    /// `String.valueOf(value)`, which is how a dump renders it and how several annotations compose
    /// their own output.
    ///
    /// `None` for a `Double`, and for any list containing one. `Double.toString` is its own
    /// algorithm, not a format string: it prints the shortest decimal that round-trips, with
    /// Java's own rules for when to switch to `E` notation and for the trailing `.0`. It is not
    /// ported yet and none of the annotations here need it, since the encoder renders `AF` rather
    /// than this method. Producing a plausible-looking rendering would be inventing a golden.
    pub fn to_java_string(&self) -> Option<String> {
        match self {
            AnnotationValue::Int(value) => Some(value.to_string()),
            AnnotationValue::Long(value) => Some(value.to_string()),
            AnnotationValue::Double(_) => None,
            AnnotationValue::Str(value) => Some(value.clone()),
            AnnotationValue::Flag(value) => Some(value.to_string()),
            // `AbstractCollection.toString`: square brackets and ", " between elements.
            AnnotationValue::List(values) => {
                let inner: Option<Vec<String>> =
                    values.iter().map(|v| v.to_java_string()).collect();
                inner.map(|parts| format!("[{}]", parts.join(", ")))
            }
        }
    }
}

/// `InfoFieldAnnotation`.
///
/// The returned vector is the map's contents in insertion order. The encoder sorts, so the order
/// here is not observable in a file; it is kept so that a divergence can be traced to the branch
/// that produced it.
pub trait InfoFieldAnnotation {
    /// `getKeyNames()`: the keys this annotation may write, in declaration order.
    fn key_names(&self) -> Vec<&'static str>;

    /// `annotate(ref, vc, likelihoods)`.
    ///
    /// An empty result means the keys are absent from the record, which is a different statement
    /// from writing zero. `None` likelihoods are the reference's `null`, which several annotations
    /// are given in practice and each one guards differently.
    fn annotate(
        &self,
        reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)>;
}
