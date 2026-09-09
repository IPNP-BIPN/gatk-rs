//! Which of two out-of-range arguments the parser reports, against the reference's own answers.
//!
//! Measured in the pinned container on `ReadAnonymizer` and `CallableLoci`, in both command-line
//! orders, and recorded in #1128. The order is a `java.util.HashMap`'s over the option specs, so a
//! test that pinned it to declaration order would pass on a port that is wrong.

fn reported(tool: &str, argv: &[&str]) -> String {
    let args: Vec<String> = argv.iter().map(|a| (*a).to_string()).collect();
    match gatk_cli::parse_for(tool, &args) {
        Ok(_) => "no refusal".to_string(),
        Err(message) => message,
    }
}

#[test]
fn the_first_out_of_range_argument_is_the_one_the_hash_map_reaches_first() {
    // `ref-base-quality` is declared at index 1 and `max-variants-per-shard` at index 20, so
    // declaration order would report the first of the two. The reference reports the second.
    for argv in [
        vec![
            "--input",
            "/x.bam",
            "--output",
            "/o.bam",
            "--reference",
            "/x.fasta",
            "--ref-base-quality",
            "61",
            "--max-variants-per-shard",
            "-1",
        ],
        vec![
            "--input",
            "/x.bam",
            "--output",
            "/o.bam",
            "--reference",
            "/x.fasta",
            "--max-variants-per-shard",
            "-1",
            "--ref-base-quality",
            "61",
        ],
    ] {
        assert_eq!(
            reported("ReadAnonymizer", &argv),
            "Argument max-variants-per-shard has a bad value: -1. minimum allowed value 0"
        );
    }
}

#[test]
fn one_bad_argument_is_still_its_own_message() {
    assert_eq!(
        reported(
            "ReadAnonymizer",
            &[
                "--input",
                "/x.bam",
                "--output",
                "/o.bam",
                "--reference",
                "/x.fasta",
                "--ref-base-quality",
                "61",
            ]
        ),
        "Argument ref-base-quality has a bad value: 61. allowed range [0, 60]."
    );
}
