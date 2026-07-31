//! Conformance for the eight flow-based annotations, against the oracle.
//!
//! Golden from `tools/annotation-conformance/FlowAnnotationDump.java`.
//!
//! ```text
//! flow  VariantType      hmer-insertion  VARIANT_TYPE=h-indel
//! flow  HmerIndelLength  hmer-insertion  X_HIL=[8]        the run, not the indel
//! flow  HmerMotifs       hmer-insertion  X_LM=[CGTAA];X_RM=[GGGGT]
//! flow  CycleSkipStatus  snp             (absent)
//! ```
//!
//! `CYCLESKIP_STATUS` is **empty on every row**: it is the one annotation that declares an actual
//! flow order is required, and without a `--flow-order` argument the base class sets
//! `generateAnnotation = false` and the whole map comes back empty. So the annotation exists, is
//! registered, and produces nothing at all unless the tool was told the machine's flow order.

use std::io::Read;

use gatk_annotation::flow::{self, Window};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

const BASES: &str = "ACGTACGTACAAAAAGGGGTTTTCCCCACGTACGTACGT";
const WINDOW_START: i64 = 90;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/flow.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::from_str(bases, is_ref).expect("an allele")
}

/// The sites the dump measured, by label.
fn site(label: &str) -> VariantContext {
    let (start, reference, alts): (i64, &str, Vec<&str>) = match label {
        "snp" => (100, "A", vec!["C"]),
        "insertion" => (100, "A", vec!["AC"]),
        "deletion" => (100, "AC", vec!["A"]),
        "hmer-insertion" => (99, "A", vec!["AA"]),
        "hmer-deletion" => (99, "AA", vec!["A"]),
        "mixed" => (100, "A", vec!["C", "AC"]),
        "spanning-deletion" => (100, "A", vec!["C", "*"]),
        "non-ref" => (100, "A", vec!["C", "<NON_REF>"]),
        "mnp" => (100, "AC", vec!["GT"]),
        "near-window-start" => (91, "A", vec!["C"]),
        "near-window-end" => (127, "A", vec!["C"]),
        "long-insertion" => (105, "A", vec!["ACGTACGT"]),
        other => panic!("{other} has no fixture"),
    };
    let mut alleles = vec![allele(reference, true)];
    for alt in alts {
        alleles.push(allele(alt, false));
    }
    let mut vc = VariantContext::new("chr1", start, alleles);
    vc.stop = start + reference.len() as i64 - 1;
    vc
}

/// `List.toString` for the lists these annotations produce, with `null` for an absent entry.
fn list<T: std::fmt::Display>(values: &[Option<T>]) -> String {
    let rendered: Vec<String> = values
        .iter()
        .map(|v| match v {
            Some(value) => value.to_string(),
            None => "null".to_string(),
        })
        .collect();
    format!("[{}]", rendered.join(", "))
}

fn plain_list<T: std::fmt::Display>(values: &[T]) -> String {
    let rendered: Vec<String> = values.iter().map(|v| v.to_string()).collect();
    format!("[{}]", rendered.join(", "))
}

/// `Float.toString`, in the shapes a GC content can take.
fn java_float(value: f32) -> String {
    if value == value.trunc() && value.abs() < 1e7 {
        return format!("{value:.1}");
    }
    format!("{value}")
}

#[test]
fn every_flow_key_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("key\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let bases = fields.next().expect("the bases");
        let order = fields.next().expect("a flow order");
        let expected = fields.next().unwrap_or("");
        let ours = match flow::base_array_to_key(bases.as_bytes(), order) {
            Some(key) => key
                .iter()
                .map(|k| k.to_string())
                .collect::<Vec<_>>()
                .join(","),
            None => "E".to_string(),
        };
        assert_eq!(ours, expected, "key of {bases:?} under {order}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no key rows");
    println!("{count} flow keys identical");
}

#[test]
fn every_annotation_answers_as_the_reference_answers() {
    let text = golden();
    let mut count = 0;
    let mut cycle_skip_rows = 0;
    let mut cycle_skip_throws = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("flow\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let name = fields.next().expect("an annotation");
        let label = fields.next().expect("a label");
        let expected = fields.next().unwrap_or("");
        let vc = site(label);
        let window = Window {
            start: WINDOW_START,
            bases: BASES.as_bytes(),
        };
        let computed = flow::annotate(&vc, &window, flow::DEFAULT_FLOW_ORDER);
        let ours = match name {
            "VariantType" => format!(
                "VARIANT_TYPE={}",
                computed.variant_type.clone().unwrap_or_default()
            ),
            "IndelClassify" => format!("X_IC={}", plain_list(&computed.indel)),
            "IndelLength" => format!("X_IL={}", list(&computed.indel_length)),
            "HmerIndelLength" => format!("X_HIL={}", plain_list(&computed.hmer_indel_length)),
            "HmerIndelNuc" => format!("X_HIN={}", list(&computed.hmer_indel_nuc)),
            "HmerMotifs" => {
                // The two motifs are put into the attribute map independently, so a left motif
                // that ran off the window leaves the right one in place rather than dropping
                // both. The golden's near-window-start row is that.
                let mut parts = Vec::new();
                if let Some(left) = &computed.left_motif {
                    parts.push(format!("X_LM={}", plain_list(left)));
                }
                if let Some(right) = &computed.right_motif {
                    parts.push(format!("X_RM={}", list(right)));
                }
                parts.join(";")
            }
            "GcContent" => match computed.gc_content {
                Some(gc) => format!("X_GCC={}", java_float(gc)),
                None => String::new(),
            },
            "CycleSkipStatus" => {
                cycle_skip_rows += 1;
                match &computed.cycle_skip {
                    // The annotation declares that an actual flow order is required, and no
                    // argument supplied one, so a successful computation still produces nothing.
                    Ok(_) => String::new(),
                    // But the motif lists are dereferenced before that check, so a variant whose
                    // left motif ran off the window throws first.
                    Err(gatk_annotation::flow::CycleSkipError::MissingMotifList) => {
                        cycle_skip_throws += 1;
                        "E:java.lang.NullPointerException".to_string()
                    }
                    // A null motif entry becomes the literal "null" in the concatenation, whose
                    // characters are not bases, so the flow key's period guard trips instead.
                    Err(gatk_annotation::flow::CycleSkipError::KeyFromNullMotif) => {
                        cycle_skip_throws += 1;
                        "E:org.broadinstitute.hellbender.exceptions.GATKException".to_string()
                    }
                }
            }
            other => panic!("unknown annotation {other}"),
        };
        assert_eq!(ours, expected, "{name} on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no flow rows");
    assert!(cycle_skip_rows > 0, "the golden lost its cycle-skip rows");
    println!(
        "{count} flow answers identical, {cycle_skip_rows} cycle-skip rows of which \
         {cycle_skip_throws} throw in the reference"
    );
}
