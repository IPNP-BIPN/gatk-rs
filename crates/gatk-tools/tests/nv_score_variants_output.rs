//! The port's scores and VCFs against what `tools/nvscorevariants/NVScoreVariantsDump.py` wrote.
//!
//! The weights and the fixture are not in the repository. Set `NVSCORE_OUTPUT` to
//! `<dump> <fixture directory> <weights directory>`. Without it this passes having compared nothing,
//! and says so. A difference on any line of any case fails it.

use gatk_tools::nv_score_variants::{run, Arguments, Stop};
use std::collections::BTreeMap;

fn file_lines(path: &str) -> Option<Vec<String>> {
    let text = std::fs::read_to_string(path).ok()?;
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    lines.pop();
    Some(lines)
}

#[test]
fn scores_and_vcfs_against_the_reference() {
    let Ok(request) = std::env::var("NVSCORE_OUTPUT") else {
        eprintln!("NVSCORE_OUTPUT is not set: nothing compared");
        return;
    };
    let parts: Vec<&str> = request.split(' ').collect();
    let (dump, fixture, weights) = (parts[0], parts[1], parts[2]);

    let mut reference: BTreeMap<(String, String), BTreeMap<String, Vec<String>>> = BTreeMap::new();
    let text = std::fs::read_to_string(dump).unwrap();
    for line in text.lines() {
        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        if fields.len() < 4 {
            continue;
        }
        reference
            .entry((fields[0].to_string(), fields[1].to_string()))
            .or_default()
            .entry(fields[2].to_string())
            .or_default()
            .push(fields[3].to_string());
    }
    assert!(!reference.is_empty(), "{dump} holds no case");

    let work = tempfile::tempdir().unwrap();
    let directory = work.path().to_string_lossy().into_owned();
    let mut failed = 0;
    let mut score_rows = 0;
    for ((name, tag), kinds) in &reference {
        let output = format!("{directory}/{name}.{tag}.vcf");
        let scores = format!("{directory}/{name}.{tag}.tmp");
        let read_tensor = tag == "2d";
        let result = run(&Arguments {
            output: output.clone(),
            variant: format!("{fixture}/{name}.vcf"),
            reference: format!("{fixture}/reference.fasta"),
            input: read_tensor.then(|| format!("{fixture}/reads.bam")),
            read_tensor,
            batch_size: 64,
            tmp_file: Some(scores.clone()),
            accelerator: "auto".to_string(),
            weights_directory: weights.to_string(),
        });
        let returned = match &result {
            Ok(()) => "0".to_string(),
            Err(Stop::Python(why)) => {
                println!("{name} {tag}: stopped: {why}");
                "1".to_string()
            }
            Err(Stop::Limitation(why)) => panic!("{name} {tag}: {why}"),
        };
        let mut compare = |kind: &str, mine: Vec<String>| {
            let want = kinds.get(kind).cloned().unwrap_or_default();
            let differing = mine.iter().zip(&want).filter(|(a, b)| a != b).count()
                + mine.len().abs_diff(want.len());
            if differing != 0 {
                failed += 1;
                println!(
                    "{name} {tag} {kind}: {differing} lines differ of {}",
                    want.len()
                );
                let first = mine.iter().zip(&want).find(|(a, b)| a != b);
                if let Some((a, b)) = first {
                    println!("  port      {a}");
                    println!("  reference {b}");
                } else {
                    println!("  port {} lines, reference {}", mine.len(), want.len());
                }
            }
        };
        compare("returned", vec![returned]);
        let absent = || vec!["absent".to_string()];
        let mine_scores = file_lines(&scores).unwrap_or_else(absent);
        if kinds.get("scores").is_some_and(|lines| lines != &absent()) {
            score_rows += mine_scores.len();
        }
        compare("scores", mine_scores);
        compare("vcf", file_lines(&output).unwrap_or_else(absent));
    }
    println!(
        "{} cases, {score_rows} score rows, {failed} differing parts",
        reference.len()
    );
    assert_eq!(failed, 0);
}
