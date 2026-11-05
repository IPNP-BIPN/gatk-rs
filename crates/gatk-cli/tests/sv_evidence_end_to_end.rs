//! FindBreakpointEvidenceSpark and StructuralVariationDiscoveryPipelineSpark, replayed against
//! GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SvEvidenceDump.java`, which draws a reference and the
//! reads of a donor genome from seeded generators, and lists everything each run left under its
//! output directory: every debug file, the FASTQs and GFAs, the contigs, the VCFs, and for a run
//! that failed whatever it had written by then.
//!
//! Every path in the golden is relative, so the test changes into a scratch directory and replays
//! each `argv` row as it is. An interval file's duration and a VCF's date are masked as the dump
//! masks them.

use gatk_corpus as corpus;

mod sv_replay;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/sv_evidence.txt.gz";

#[test]
fn every_run_replays_as_the_reference_ran_it() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
        return;
    }
    let text = corpus::read_golden(&golden);
    let root = std::env::temp_dir().join("gatk-cli-sv-evidence");
    let _ = std::fs::remove_dir_all(&root);
    sv_replay::replay(&text, &root);
}
