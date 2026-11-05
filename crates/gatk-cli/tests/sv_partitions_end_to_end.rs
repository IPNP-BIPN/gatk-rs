//! FindBreakpointEvidenceSpark and StructuralVariationDiscoveryPipelineSpark on BAMs Disq cuts into
//! several partitions, and under dictionaries that make KmerCleaner's table large, replayed against
//! GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SvPartitionsDump.java`, replayed as
//! `sv_evidence_end_to_end` replays its own.

use gatk_corpus as corpus;

mod sv_replay;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "../gatk-tools/tests/data/sv_partitions.txt.gz";

#[test]
fn every_run_replays_as_the_reference_ran_it() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    if !golden.exists() {
        // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
        return;
    }
    let text = corpus::read_golden(&golden);
    let root = std::env::temp_dir().join("gatk-cli-sv-partitions");
    let _ = std::fs::remove_dir_all(&root);
    sv_replay::replay(&text, &root);
}
