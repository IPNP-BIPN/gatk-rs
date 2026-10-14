//! `VariantContext.getType()` named, from the alleles' display strings (htsjdk 4.2.0).

/// `VariantContext.getType()`, by the same rule [`crate::remove_nearby_indels`] ports: length
/// decides, and two alternates that disagree are MIXED.
pub fn type_name(reference: &str, alternates: &[String]) -> &'static str {
    if alternates.is_empty() || (alternates.len() == 1 && alternates[0] == ".") {
        return "NO_VARIATION";
    }
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
        match kind {
            None => kind = Some(this),
            Some(seen) if seen == this => {}
            Some(_) => return "MIXED",
        }
    }
    kind.unwrap_or("NO_VARIATION")
}
