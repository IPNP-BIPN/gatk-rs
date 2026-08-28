//! Conformance for `Funcotator` against GATK 4.6.2.0, compared as the strings a GENCODE
//! funcotation is built out of.
//!
//! Golden from `tools/readfilter-conformance/FuncotatorDump.java`, whose GENCODE source is
//! FuncotateSegments' fixture reused whole. Walking the transcripts is not ported, so what the
//! port is asked for is the arithmetic and the formatting: given the transcript positions the
//! golden itself carries, it must produce the same four change strings.
//!
//! # What this suite is for
//!
//!  * **every alternate getting its own bracketed funcotation, joined by a comma**;
//!  * **the field order being the one the header line declares, and the header ending on an empty
//!    field**;
//!  * **the codon change writing the changed bases in upper case and the rest in lower**;
//!  * **the protein change putting the position in the middle for one amino acid**;
//!  * **the cDNA change and the genomic change being two different strings for one variant**;
//!  * **an intronic variant being named off its closest exon**;
//!  * **a variant outside every gene carrying the literals `Unknown` and `no_transcript`**;
//!  * **and the genetic code answering the letters the golden's protein changes carry.**

use gatk_corpus as corpus;
use gatk_tools::funcotator::{
    aligned_end_position, aligned_position, closest_exon_index,
    coding_sequence_change_string_for_xnp, codon_change_string_for_onp,
    eukaryotic_amino_acid_by_codon, extract_funcotator_keys_from_header_description,
    genome_change_string_for_xnp, intronic_cdna_string, is_position_in_frame,
    protein_change_end_position, protein_change_position, protein_change_string_for_onp,
    render_sanitized_funcotation_for_vcf, sanitize_funcotation_field_for_maf,
    sanitize_funcotation_field_for_vcf, transcript_id_without_version_number, CODON_LENGTH,
    NO_TRANSCRIPT, UNKNOWN_GENE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/funcotator.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The payload of the one line whose kind and label are given, the label being `=`-terminated.
fn payload(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("{kind} {label} is in the golden"));
    unescape(&line[prefix.len()..])
}

/// The `##INFO=<ID=FUNCOTATION` line of one case.
fn header(text: &str, label: &str) -> String {
    let prefix = format!("header\t{label}\t");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("header {label} is in the golden"));
    unescape(&line[prefix.len()..])
}

/// The declared funcotation field names, taken out of the header's description.
fn field_names(text: &str, label: &str) -> Vec<String> {
    let line = header(text, label);
    let description = line
        .split_once("Description=\"")
        .expect("a description")
        .1
        .rsplit_once('"')
        .expect("a closing quote")
        .0;
    extract_funcotator_keys_from_header_description(description)
}

/// One annotated record, as its position, its alleles and its funcotations.
struct Record {
    position: i32,
    reference: String,
    alternates: Vec<String>,
    funcotations: Vec<Vec<String>>,
}

fn records(text: &str, label: &str) -> Vec<Record> {
    payload(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let info = columns[7];
            let value = info
                .split(';')
                .find_map(|part| part.strip_prefix("FUNCOTATION="))
                .expect("a FUNCOTATION field");
            Record {
                position: columns[1].parse().expect("a position"),
                reference: columns[3].to_string(),
                alternates: columns[4].split(',').map(str::to_string).collect(),
                funcotations: value
                    .split("],[")
                    .map(|one| {
                        one.trim_start_matches('[')
                            .trim_end_matches(']')
                            .split('|')
                            .map(str::to_string)
                            .collect()
                    })
                    .collect(),
            }
        })
        .collect()
}

fn get<'a>(names: &[String], funcotation: &'a [String], field: &str) -> &'a str {
    let index = names
        .iter()
        .position(|name| name == &format!("Gencode_1_{field}"))
        .unwrap_or_else(|| panic!("{field} is a declared field"));
    &funcotation[index]
}

/// The header declares twenty-two fields and ends on an empty one, which is why every funcotation
/// ends on a trailing `|`.
#[test]
fn the_header_declares_the_field_order() {
    let text = golden();
    let names = field_names(&text, "annotated");
    assert_eq!(
        names.first().map(String::as_str),
        Some("Gencode_1_hugoSymbol")
    );
    assert_eq!(
        names.last().map(String::as_str),
        Some("Gencode_1_otherTranscripts")
    );
    assert_eq!(names.len(), 22);
    // Every record's funcotation has exactly one value per declared field, the last being empty.
    for record in records(&text, "annotated") {
        for funcotation in &record.funcotations {
            assert_eq!(funcotation.len(), names.len(), "at {}", record.position);
            assert_eq!(funcotation.last().map(String::as_str), Some(""));
        }
    }
}

/// Every alternate gets its own bracketed funcotation, so the two-alternate site carries two.
#[test]
fn every_alternate_gets_its_own_funcotation() {
    let text = golden();
    for record in records(&text, "annotated") {
        assert_eq!(
            record.funcotations.len(),
            record.alternates.len(),
            "at {}",
            record.position
        );
    }
    let two = records(&text, "annotated")
        .into_iter()
        .find(|record| record.alternates.len() == 2)
        .expect("a two-alternate site");
    assert_eq!(two.position, 1060);
}

/// The genomic change is rebuilt from the record's own position and alleles.
#[test]
fn the_genomic_change_is_the_position_and_the_alleles() {
    let text = golden();
    let names = field_names(&text, "annotated");
    for record in records(&text, "annotated") {
        for (alternate, funcotation) in record.alternates.iter().zip(&record.funcotations) {
            assert_eq!(
                genome_change_string_for_xnp("chr1", record.position, &record.reference, alternate),
                get(&names, funcotation, "genomeChange"),
                "at {}",
                record.position
            );
        }
    }
}

/// The codon change is rebuilt from the aligned codon the golden's own transcript position names,
/// with the changed base in upper case and the rest in lower.
#[test]
fn the_codon_change_upper_cases_the_changed_bases() {
    let text = golden();
    let names = field_names(&text, "annotated");
    for record in records(&text, "annotated") {
        for (alternate, funcotation) in record.alternates.iter().zip(&record.funcotations) {
            let codon_change = get(&names, funcotation, "codonChange");
            if codon_change.is_empty() {
                continue;
            }
            let transcript_position: i32 = get(&names, funcotation, "transcriptPos")
                .parse()
                .expect("a transcript position");
            let start = aligned_position(transcript_position);
            let stop = aligned_end_position(transcript_position);
            // The reference codon is the one the golden itself prints, read back in upper case.
            let (window, bases) = codon_change.split_once(')').expect("a codon window");
            assert_eq!(
                window,
                format!("c.({start}-{stop}"),
                "at {}",
                record.position
            );
            let (reference_codon, alternate_codon) = bases.split_once('>').expect("two codons");
            assert_eq!(
                codon_change_string_for_onp(
                    &reference_codon.to_lowercase(),
                    &alternate_codon.to_lowercase(),
                    start,
                    stop
                ),
                codon_change,
                "at {}",
                record.position
            );
            // And the changed base is the alternate allele, in the codon's own frame.
            let offset = (transcript_position - start) as usize;
            assert_eq!(
                alternate_codon.to_uppercase().as_bytes()[offset],
                alternate.as_bytes()[0],
                "at {}",
                record.position
            );
        }
    }
}

/// The protein change puts the position in the middle, and its letters are the genetic code's
/// answer for the two codons the codon change prints.
#[test]
fn the_protein_change_follows_the_genetic_code() {
    let text = golden();
    let names = field_names(&text, "annotated");
    for record in records(&text, "annotated") {
        for funcotation in &record.funcotations {
            let protein_change = get(&names, funcotation, "proteinChange");
            let codon_change = get(&names, funcotation, "codonChange");
            if protein_change.is_empty() || codon_change.is_empty() {
                continue;
            }
            let transcript_position: i32 = get(&names, funcotation, "transcriptPos")
                .parse()
                .expect("a transcript position");
            let start = aligned_position(transcript_position);
            let bases = codon_change.split_once(')').expect("a codon window").1;
            let (reference_codon, alternate_codon) = bases.split_once('>').expect("two codons");
            let reference_amino_acid =
                eukaryotic_amino_acid_by_codon(reference_codon).expect("a reference amino acid");
            let alternate_amino_acid =
                eukaryotic_amino_acid_by_codon(alternate_codon).expect("an alternate amino acid");
            let position = protein_change_position(start);
            assert_eq!(
                protein_change_string_for_onp(
                    reference_amino_acid,
                    alternate_amino_acid,
                    position,
                    position
                ),
                protein_change,
                "at {}",
                record.position
            );
        }
    }
}

/// The cDNA change is the transcript position and the two alleles, which is a different string
/// from the genomic change even though both name one base.
#[test]
fn the_cdna_change_is_the_transcript_position() {
    let text = golden();
    let names = field_names(&text, "annotated");
    for record in records(&text, "annotated") {
        for (alternate, funcotation) in record.alternates.iter().zip(&record.funcotations) {
            let cdna_change = get(&names, funcotation, "cDnaChange");
            let classification = get(&names, funcotation, "variantClassification");
            if classification == "INTRON" || classification == "IGR" {
                continue;
            }
            let transcript_position: i32 = get(&names, funcotation, "transcriptPos")
                .parse()
                .expect("a transcript position");
            assert_eq!(
                coding_sequence_change_string_for_xnp(
                    transcript_position,
                    &record.reference,
                    alternate
                ),
                cdna_change,
                "at {}",
                record.position
            );
            assert_ne!(cdna_change, get(&names, funcotation, "genomeChange"));
        }
    }
}

/// An intronic variant is named off its closest exon: the fixture's one exon runs to 1200, and
/// the variant at 1500 is three hundred bases past its end.
#[test]
fn an_intronic_variant_is_named_off_its_closest_exon() {
    let text = golden();
    let names = field_names(&text, "annotated");
    let intron = records(&text, "annotated")
        .into_iter()
        .find(|record| record.position == 1500)
        .expect("the intronic record");
    let funcotation = &intron.funcotations[0];
    assert_eq!(get(&names, funcotation, "variantClassification"), "INTRON");
    let exons = [(1000, 1200)];
    assert_eq!(closest_exon_index(1500, &exons), Some(0));
    assert_eq!(
        intronic_cdna_string(1500, &exons, &intron.reference, &intron.alternates[0]),
        get(&names, funcotation, "cDnaChange")
    );
    assert_eq!(get(&names, funcotation, "cDnaChange"), "c.e1+300T>A");
    // Its codon and protein changes are empty: there is no codon to change.
    assert_eq!(get(&names, funcotation, "codonChange"), "");
    assert_eq!(get(&names, funcotation, "proteinChange"), "");
}

/// A variant outside every gene still gets a funcotation, whose gene and transcript are literals.
#[test]
fn a_variant_outside_every_gene_carries_the_literals() {
    let text = golden();
    let names = field_names(&text, "annotated");
    let outside = records(&text, "annotated")
        .into_iter()
        .find(|record| record.position == 5000)
        .expect("the intergenic record");
    let funcotation = &outside.funcotations[0];
    assert_eq!(get(&names, funcotation, "variantClassification"), "IGR");
    assert_eq!(get(&names, funcotation, "hugoSymbol"), UNKNOWN_GENE);
    assert_eq!(
        get(&names, funcotation, "annotationTranscript"),
        NO_TRANSCRIPT
    );
}

/// Removing the filtered variants takes the file from five records to four, and the record it
/// drops is the one the unfiltered run annotated with a filter of its own.
#[test]
fn removing_the_filtered_variants_drops_a_record() {
    let text = golden();
    let all = records(&text, "annotated");
    let kept = records(&text, "remove-filtered");
    assert_eq!(all.len(), 5);
    assert_eq!(kept.len(), 4);
    let dropped: Vec<i32> = all
        .iter()
        .map(|record| record.position)
        .filter(|position| !kept.iter().any(|record| record.position == *position))
        .collect();
    assert_eq!(dropped, vec![1070]);
}

/// The transcript the golden names carries a version, which the reference strips off by name.
#[test]
fn the_transcript_version_is_stripped_by_name() {
    let text = golden();
    let names = field_names(&text, "annotated");
    let record = &records(&text, "annotated")[0];
    let transcript = get(&names, &record.funcotations[0], "annotationTranscript");
    assert_eq!(transcript, "ENST00000000001.1");
    assert_eq!(
        transcript_id_without_version_number(transcript),
        "ENST00000000001"
    );
    // Only a trailing dot and digits go: a name that ends on letters is left alone.
    assert_eq!(
        transcript_id_without_version_number(NO_TRANSCRIPT),
        NO_TRANSCRIPT
    );
    assert_eq!(transcript_id_without_version_number("ENST1.2.3"), "ENST1.2");
}

/// A reference version with no data-source directory is refused by name.
#[test]
fn a_reference_version_with_no_directory_is_refused() {
    let text = golden();
    let line = text
        .lines()
        .find(|line| line.starts_with("error\twrong-ref-version\t"))
        .expect("the refusal");
    assert!(
        line.ends_with("Could not find any data sources for given reference: hg38"),
        "{line}"
    );
}

/// The sanitiser encodes eight characters for a VCF and two for a MAF, so a field carrying a pipe
/// survives the join it would otherwise break.
#[test]
fn the_sanitiser_encodes_what_would_break_the_format() {
    assert_eq!(sanitize_funcotation_field_for_vcf("a|b"), "a_%7C_b");
    assert_eq!(
        sanitize_funcotation_field_for_vcf("a,b;c=d"),
        "a_%2C_b_%3B_c_%3D_d"
    );
    assert_eq!(sanitize_funcotation_field_for_vcf("a b"), "a_%20_b");
    assert_eq!(sanitize_funcotation_field_for_vcf("a\tb"), "a_%09_b");
    assert_eq!(sanitize_funcotation_field_for_vcf("a\nb"), "a_%0A_b");
    assert_eq!(sanitize_funcotation_field_for_vcf("a#b"), "a_%23_b");
    // The MAF one leaves the pipe and the comma alone.
    assert_eq!(sanitize_funcotation_field_for_maf("a|b,c"), "a|b,c");
    assert_eq!(
        sanitize_funcotation_field_for_maf("a\tb\nc"),
        "a_%09_b_%0A_c"
    );
}

/// The renderer joins the INCLUDED fields in the funcotation's own order, and renders nothing at
/// all when nothing is included.
#[test]
fn the_renderer_follows_the_funcotations_own_order() {
    let fields = vec![
        ("a".to_string(), "one".to_string()),
        ("b".to_string(), "two".to_string()),
        ("c".to_string(), "three".to_string()),
    ];
    let included = vec!["c".to_string(), "a".to_string()];
    assert_eq!(
        render_sanitized_funcotation_for_vcf(&fields, &included),
        "one|three"
    );
    assert_eq!(render_sanitized_funcotation_for_vcf(&fields, &[]), "");
}

/// The codon arithmetic, which everything above rests on.
#[test]
fn the_codon_arithmetic_is_one_based() {
    assert_eq!(CODON_LENGTH, 3);
    assert_eq!(aligned_position(1), 1);
    assert_eq!(aligned_position(51), 49);
    assert_eq!(aligned_position(61), 61);
    assert_eq!(aligned_end_position(51), 51);
    assert_eq!(aligned_end_position(61), 63);
    assert!(is_position_in_frame(61));
    assert!(!is_position_in_frame(51));
    assert_eq!(protein_change_position(49), 17);
    assert_eq!(protein_change_position(61), 21);
    // A change one codon long ends where it starts.
    assert_eq!(protein_change_end_position(21, 3), 21);
    assert_eq!(protein_change_end_position(21, 6), 22);
    // A position at or before the start folds the other way, three at a time: nought, minus one
    // and minus two all align to minus two, and minus three to minus five.
    assert_eq!(aligned_position(0), -2);
    assert_eq!(aligned_position(-1), -2);
    assert_eq!(aligned_position(-2), -2);
    assert_eq!(aligned_position(-3), -5);
}
