//! `VariantContext.getType()` named, from the alleles' display strings (htsjdk 5.0.0;
//! `determineType` and `typeOfBiallelicVariant` are unchanged since 4.2.0 but for formatting).
//!
//! The callers hold a record as its VCF text (the REF string and the ALT strings), not as an
//! htsjdk `VariantContext`, so the rule is applied to the strings directly.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` is a borrowed, read-only string; `&[String]` is a borrowed, read-only list of owned
//!   strings. The function reads them and the caller keeps them.
//! - `&'static str` is a string literal that lives for the whole program, so the function can
//!   return one of its constant names without allocating.
//! - `Option<T>` is "a value or nothing" (Java's nullable reference, checked by the compiler).

/// `VariantContext.getType()`, by the same rule [`crate::remove_nearby_indels`] ports: length
/// decides, and two alternates that disagree are MIXED.
///
/// **What**: the name of the `VariantContext.Type` constant htsjdk would report for a record
/// with this REF and these ALTs.
///
/// **How**: each alternate is typed against the reference as `typeOfBiallelicVariant` does,
/// and the first disagreement returns MIXED, as `determinePolymorphicType` does:
///
/// ```java
/// if (allele.isSymbolic()) return Type.SYMBOLIC;
///
/// if (ref.length() == allele.length()) {
///     if (allele.length() == 1) return Type.SNP;
///     else return Type.MNP;
/// }
/// ```
///
/// An ALT column of `.` is a record with only its reference allele, which `determineType`
/// reports as `NO_VARIATION`. The symbolic test covers `<...>` and breakpoint (`[`, `]`)
/// alleles; it does not cover a single breakend such as `A.`, which htsjdk's
/// `wouldBeSingleBreakend` also counts as symbolic.
pub fn type_name(reference: &str, alternates: &[String]) -> &'static str {
    if alternates.is_empty() || (alternates.len() == 1 && alternates[0] == ".") {
        return "NO_VARIATION";
    }
    // The type of the alternates seen so far; `None` until the first one, as Java's `type == null`.
    let mut kind: Option<&'static str> = None;
    for alternate in alternates {
        let this =
            if alternate.starts_with('<') || alternate.contains('[') || alternate.contains(']') {
                "SYMBOLIC"
            } else if alternate.len() == reference.len() {
                if reference.len() == 1 {
                    "SNP"
                } else {
                    "MNP"
                }
            } else {
                "INDEL"
            };
        // `match` must cover every case. The middle arm has a guard (`if seen == this`): the
        // same type again changes nothing; any other type is the reference's early `MIXED`.
        match kind {
            None => kind = Some(this),
            Some(seen) if seen == this => {}
            Some(_) => return "MIXED",
        }
    }
    kind.unwrap_or("NO_VARIATION")
}
