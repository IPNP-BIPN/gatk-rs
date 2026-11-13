//! #1021's first measurement: how much of NVScoreVariants is inference, against the whole tool.
//!
//!     cargo run --release -p gatk-tools --example nvscore_profile -- <fixture> <weights> [repeats]
//!
//! `<fixture>` is what `tools/nvscorevariants/make_fixture.py` writes, `<weights>` what
//! `extract_weights.py` writes. Each mode is timed as the tool runs it (`run`, every core), then
//! split into encoding (FASTA, VCF and BAM reading included) and inference at one thread and at
//! every core, best of `repeats`.

use std::time::Instant;

use gatk_tools::nv_score_variants::model::Trace;
use gatk_tools::nv_score_variants::weights::{self, Weights};
use gatk_tools::nv_score_variants::{encode_for_comparison, run, score_batch, Arguments};

fn best<T>(repeats: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut out = None;
    let mut fastest = f64::INFINITY;
    for _ in 0..repeats {
        let start = Instant::now();
        let value = f();
        fastest = fastest.min(start.elapsed().as_secs_f64());
        out = Some(value);
    }
    (fastest, out.unwrap())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (fixture, weights_directory) = (&args[1], &args[2]);
    let repeats = args.get(3).map_or(3, |r| r.parse().unwrap());
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    let work = std::env::temp_dir();
    println!("| mode | variants | whole tool, {cores} threads | encoding | inference, 1 thread | inference, {cores} threads | inference share of the tool |");
    println!("|---|---:|---:|---:|---:|---:|---:|");
    for (tag, read_tensor, pinned) in [("1D", false, weights::ONE_D), ("2D", true, weights::TWO_D)]
    {
        let arguments = Arguments {
            output: work
                .join(format!("nvscore_profile.{tag}.vcf"))
                .to_string_lossy()
                .into_owned(),
            variant: format!("{fixture}/variants.vcf"),
            reference: format!("{fixture}/reference.fasta"),
            input: read_tensor.then(|| format!("{fixture}/reads.bam")),
            read_tensor,
            batch_size: 64,
            tmp_file: Some(
                work.join(format!("nvscore_profile.{tag}.tmp"))
                    .to_string_lossy()
                    .into_owned(),
            ),
            accelerator: "auto".to_string(),
            weights_directory: weights_directory.clone(),
        };
        let weights = Weights::load(weights_directory, pinned).unwrap();
        let (whole, result) = best(repeats, || run(&arguments));
        result.unwrap();
        let (encoding, batches) = best(repeats, || encode_for_comparison(&arguments).unwrap());
        let variants: usize = batches.iter().map(|(records, _)| records.len()).sum();
        let infer = |threads| {
            best(repeats, || {
                for (_, batch) in &batches {
                    score_batch(&weights, read_tensor, batch, threads, &mut Trace::default());
                }
            })
            .0
        };
        let (one, all) = (infer(1), infer(cores));
        println!(
            "| {tag} | {variants} | {whole:.3} s | {encoding:.3} s | {one:.3} s | {all:.3} s | {:.0}% |",
            100.0 * all / whole
        );
    }
}
