//! Milestone GPU.3: the speedup of each device kernel, measured beside its byte-equality.
//!
//! `GATK_RS_GPU=1 cargo run --release -p gatk-engine --features gpu --example gpu_speed`
//!
//! Prints one line per workload: the CPU path as the engine runs it (one thread), the same work
//! spread over every core, the device, and whether all three gave the same bits.

use std::time::Instant;

use gatk_engine::gpu;
use gatk_engine::pair_hmm::read_likelihood_given_haplotype_log10;

struct Random(u64);
impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn bases(&mut self, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }
    fn quals(&mut self, n: usize, low: u64, span: u64) -> Vec<u8> {
        (0..n).map(|_| (low + self.next() % span) as u8).collect()
    }
}

type Read = (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>);

fn cpu(haplotypes: &[Vec<u8>], reads: &[Read], threads: usize) -> Vec<f64> {
    let row = |r: &Read| -> Vec<f64> {
        haplotypes
            .iter()
            .map(|h| read_likelihood_given_haplotype_log10(h, &r.0, &r.1, &r.2, &r.3, &r.4))
            .collect()
    };
    if threads == 1 {
        return reads.iter().flat_map(row).collect();
    }
    let chunk = reads.len().div_ceil(threads);
    std::thread::scope(|s| {
        let parts: Vec<_> = reads
            .chunks(chunk)
            .map(|part| s.spawn(move || part.iter().flat_map(row).collect::<Vec<_>>()))
            .collect();
        parts.into_iter().flat_map(|p| p.join().unwrap()).collect()
    })
}

fn seconds(start: Instant) -> f64 {
    start.elapsed().as_secs_f64()
}

fn main() {
    let Some(device) = gpu::device_name() else {
        eprintln!("no device path: build with --features gpu and set GATK_RS_GPU=1");
        std::process::exit(2);
    };
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    println!("device: {device}; cores: {cores}");
    println!("workload | pairs | cells | cpu 1 thread s | cpu {cores} threads s | device s | x vs 1 thread | x vs {cores} threads | same bits");
    // An active region's shape: a few haplotypes of a few hundred bases, reads of 150.
    for (haplotype_count, haplotype_length, read_count, read_length) in [
        (4, 300, 100, 150),
        (16, 400, 500, 150),
        (32, 500, 2000, 150),
    ] {
        let mut random = Random(haplotype_count as u64);
        let haplotypes: Vec<Vec<u8>> = (0..haplotype_count)
            .map(|_| random.bases(haplotype_length))
            .collect();
        let reads: Vec<Read> = (0..read_count)
            .map(|_| {
                (
                    random.bases(read_length),
                    random.quals(read_length, 10, 30),
                    random.quals(read_length, 40, 5),
                    random.quals(read_length, 40, 5),
                    vec![10; read_length],
                )
            })
            .collect();
        let inputs: Vec<gpu::ReadInputs> = reads
            .iter()
            .map(|r| gpu::ReadInputs {
                bases: &r.0,
                qualities: &r.1,
                insertion: &r.2,
                deletion: &r.3,
                gap_continuation: &r.4,
            })
            .collect();
        // One warm-up so the kernel's compilation is not charged to the first workload.
        gpu::pair_hmm_log10(&haplotypes[..1], &inputs[..1]);

        let t = Instant::now();
        let one = cpu(&haplotypes, &reads, 1);
        let one_s = seconds(t);
        let t = Instant::now();
        let many = cpu(&haplotypes, &reads, cores);
        let many_s = seconds(t);
        let t = Instant::now();
        let dev = gpu::pair_hmm_log10(&haplotypes, &inputs).expect("device path");
        let dev_s = seconds(t);
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let same = bits(&one) == bits(&many) && bits(&one) == bits(&dev);
        println!(
            "{haplotype_count}x{haplotype_length} hap, {read_count}x{read_length} read | {} | {} | {one_s:.3} | {many_s:.3} | {dev_s:.3} | {:.2} | {:.2} | {same}",
            haplotype_count * read_count,
            haplotype_count * read_count * haplotype_length * read_length,
            one_s / dev_s,
            many_s / dev_s,
        );
    }
}
