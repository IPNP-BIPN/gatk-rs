//! The port against a dump of the reference's own layers, when one is there.
//!
//! `tools/nvscorevariants/dump_layers.py` writes the dump in the reference's image; the weights are
//! not in the repository. Set `NVSCORE_LAYERS` to `<reference|read_tensor> <vcf> <fasta> <bam or -> <batch size> <dump> <weights directory>`
//! to get a line per layer: how many of its values differ, the largest distance in units in the last
//! place, and the largest distance outright. Without the variable this passes having compared
//! nothing, and says so.
//!
//! It asserts what the reference keeps stable and the port has been measured to match: the inputs
//! and the batch norm everywhere, and the convolutions under a oneDNN dispatch where they were
//! measured exact. `NVSCORE_DISPATCH` carries the reference's `ONEDNN_VERBOSE` lines; under any other
//! dispatch the convolutions are printed and not asserted. Everything from the first `SGEMM` on is
//! printed only: the reference moves those bits itself with the thread count and the batch
//! (`docs/nvscorevariants.md`).
//!
//! For the 1D model it also recomputes each convolution of the first rows from the reference's own
//! previous layer under several candidate orders, and prints how many values each misses: which
//! order a dispatch the port does not model yet takes.

use gatk_tools::nv_score_variants::kernels::{fused, relu, ConvWeights};
use gatk_tools::nv_score_variants::model::Trace;
use gatk_tools::nv_score_variants::weights::{self, Weights};
use gatk_tools::nv_score_variants::{encode_for_comparison, score_batch, Arguments};
use std::collections::HashMap;
use std::io::Read;

struct Entry {
    values: Vec<f32>,
    integers: Vec<i64>,
}

fn take_u32(file: &mut impl Read) -> Option<u32> {
    let mut bytes = [0u8; 4];
    file.read_exact(&mut bytes).ok()?;
    Some(u32::from_le_bytes(bytes))
}

fn take_u64(file: &mut impl Read) -> Option<u64> {
    let mut bytes = [0u8; 8];
    file.read_exact(&mut bytes).ok()?;
    Some(u64::from_le_bytes(bytes))
}

/// The next record of a dump: its name and what it holds.
fn next_entry(file: &mut impl Read) -> Option<(String, Entry)> {
    let length = take_u32(file)?;
    let mut name = vec![0u8; length as usize];
    file.read_exact(&mut name).ok()?;
    let kind = take_u32(file)?;
    let dimensions = take_u32(file)?;
    let mut count = 1usize;
    for _ in 0..dimensions {
        count *= take_u64(file)? as usize;
    }
    let mut entry = Entry {
        values: Vec::new(),
        integers: Vec::new(),
    };
    if kind == 0 {
        let mut bytes = vec![0u8; count * 4];
        file.read_exact(&mut bytes).ok()?;
        entry.values = bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| f32::from_le_bytes(*b))
            .collect();
    } else {
        let mut bytes = vec![0u8; count * 8];
        file.read_exact(&mut bytes).ok()?;
        entry.integers = bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|b| i64::from_le_bytes(*b))
            .collect();
    }
    Some((String::from_utf8(name).ok()?, entry))
}

fn ordinal(value: f32) -> i64 {
    let bits = i64::from(value.to_bits());
    if bits & 0x8000_0000 != 0 {
        0x8000_0000 - bits
    } else {
        bits
    }
}

/// `(values that differ, largest ulp, largest distance)`, with -0 and +0 two different things.
fn distance(a: &[f32], b: &[f32]) -> (usize, i64, f64) {
    assert_eq!(a.len(), b.len());
    let mut result = (0usize, 0i64, 0f64);
    for (x, y) in a.iter().zip(b) {
        if x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()) {
            continue;
        }
        result.0 += 1;
        result.1 = result.1.max((ordinal(*x) - ordinal(*y)).abs());
        result.2 = result.2.max((f64::from(*x) - f64::from(*y)).abs());
    }
    result
}

#[test]
fn layers_against_the_reference() {
    let Ok(request) = std::env::var("NVSCORE_LAYERS") else {
        eprintln!("NVSCORE_LAYERS is not set: nothing compared");
        return;
    };
    let parts: Vec<&str> = request.split(' ').collect();
    let read_tensor = parts[0] == "read_tensor";
    let arguments = Arguments {
        output: String::new(),
        variant: parts[1].to_string(),
        reference: parts[2].to_string(),
        input: (parts[3] != "-").then(|| parts[3].to_string()),
        read_tensor,
        batch_size: parts[4].parse().unwrap(),
        tmp_file: None,
        accelerator: "cpu".to_string(),
        weights_directory: parts[6].to_string(),
    };
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let pinned = if read_tensor {
        weights::TWO_D
    } else {
        weights::ONE_D
    };
    let model = Weights::load(&arguments.weights_directory, pinned).unwrap();
    let batches = encode_for_comparison(&arguments).unwrap();
    let mut file = std::io::BufReader::new(std::fs::File::open(parts[5]).unwrap());

    let mut totals: Vec<(String, usize, usize, i64, f64)> = Vec::new();
    let mut add = |name: &str, total: usize, d: (usize, i64, f64)| {
        if let Some(row) = totals.iter_mut().find(|row| row.0 == name) {
            row.1 += total;
            row.2 += d.0;
            row.3 = row.3.max(d.1);
            row.4 = row.4.max(d.2);
        } else {
            totals.push((name.to_string(), total, d.0, d.1, d.2));
        }
    };

    // A batch's dump is in the order the reference made it, inputs first and `score` last.
    let mut pending: HashMap<String, Entry> = HashMap::new();
    let mut current = 0usize;
    let mut computed: Option<Trace> = None;
    let mut printed_rows = (0usize, 0usize);
    let mut first_batch: HashMap<String, Vec<f32>> = HashMap::new();
    while let Some((name, entry)) = next_entry(&mut file) {
        let (batch_text, layer) = name.split_once('/').unwrap();
        let index: usize = batch_text.parse().unwrap();
        if index != current {
            current = index;
            computed = None;
            pending.clear();
        }
        let (records, batch) = &batches[index];
        match layer {
            "input.type" => assert_eq!(entry.integers, batch.types, "types of {index}"),
            "input.best_practices" => {
                let d = distance(&entry.values, &batch.best_practices);
                add(layer, entry.values.len(), d);
                assert_eq!(d.0, 0, "annotations of batch {index}");
            }
            "input.reference" | "input.read_tensor" => {
                let d = distance(&entry.values, &batch.tensor);
                add(layer, entry.values.len(), d);
                if d.0 != 0 {
                    report_tensor(records, &entry.values, &batch.tensor, read_tensor);
                }
                assert_eq!(d.0, 0, "tensor of batch {index}");
            }
            _ => {
                let trace = computed.get_or_insert_with(|| {
                    let mut trace = Trace {
                        keep: true,
                        ..Trace::default()
                    };
                    score_batch(&model, read_tensor, batch, threads, &mut trace);
                    trace
                });
                let port_name = match layer {
                    "logits" if read_tensor => "dense_2.3",
                    "logits" => "dense_3",
                    other => other,
                };
                if index == 0 && layer.starts_with("conv.") {
                    first_batch.insert(layer.to_string(), entry.values.clone());
                }
                if let Some(mine) = trace.layers.get(port_name) {
                    add(layer, mine.len(), distance(mine, &entry.values));
                    // By batch too: the last batch of a file is the short one, and the short one is
                    // where the known differences are.
                    add(
                        &format!("{index}/{layer}"),
                        mine.len(),
                        distance(mine, &entry.values),
                    );
                    if layer == "score" {
                        for (a, b) in mine.iter().zip(&entry.values) {
                            printed_rows.0 += 1;
                            if format!("{a:.3}") != format!("{b:.3}") {
                                printed_rows.1 += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    for (name, total, differ, ulp, absolute) in &totals {
        println!("{name}\t{differ}/{total}\tmax ulp {ulp}\tmax abs {absolute:e}");
    }
    println!(
        "score at %.3f\t{} of {} rows differ",
        printed_rows.1, printed_rows.0
    );

    let dispatch = std::env::var("NVSCORE_DISPATCH").unwrap_or_default();
    // Measured exact: 1D under `jit:avx2` (Rosetta), 2D under `brgconv:avx2` (Rosetta) and under
    // AVX-512 (`jit:avx512_core` reported by the 1D run, Xeon 6975P, run 44799995145).
    let convolutions_exact = if read_tensor {
        dispatch.contains("brgconv:avx2") || dispatch.contains("avx512_core")
    } else {
        dispatch.contains("jit:avx2,")
    };
    println!(
        "dispatch\t{}\tconvolutions {}",
        dispatch.lines().next().unwrap_or("not given"),
        if convolutions_exact {
            "asserted"
        } else {
            "printed only"
        }
    );
    for (name, _, differ, _, _) in &totals {
        let always = name == "normalize" || (read_tensor && name == "dense_1.0");
        let convolution = name.starts_with("conv.");
        if always || (convolution && convolutions_exact) {
            assert_eq!(*differ, 0, "{name} is measured exact under this dispatch");
        }
    }
    if !read_tensor {
        probe_1d_orders(&model, &batches[0].1.tensor, &first_batch);
    }
}

/// One candidate order of a oneDNN convolution: where the bias goes, the input channel block, the
/// chunk of channels summed from zero and added through memory, and whether the channels of a
/// block are taken inside or outside the kernel positions.
struct Order {
    name: &'static str,
    bias_first: bool,
    block: usize,
    chunk: usize,
    channels_outside: bool,
}

const ORDERS: [Order; 10] = [
    Order {
        name: "bias first, block 8 (jit:avx2)",
        bias_first: true,
        block: 8,
        chunk: usize::MAX,
        channels_outside: false,
    },
    Order {
        name: "bias first, block 16",
        bias_first: true,
        block: 16,
        chunk: usize::MAX,
        channels_outside: false,
    },
    Order {
        name: "bias first, block 16, channels outside",
        bias_first: true,
        block: 16,
        chunk: usize::MAX,
        channels_outside: true,
    },
    Order {
        name: "bias last, block 8",
        bias_first: false,
        block: 8,
        chunk: usize::MAX,
        channels_outside: false,
    },
    Order {
        name: "bias last, block 16",
        bias_first: false,
        block: 16,
        chunk: usize::MAX,
        channels_outside: false,
    },
    Order {
        name: "bias last, block 16, channels outside",
        bias_first: false,
        block: 16,
        chunk: usize::MAX,
        channels_outside: true,
    },
    Order {
        name: "bias last, block 16, chunks of 16",
        bias_first: false,
        block: 16,
        chunk: 16,
        channels_outside: false,
    },
    Order {
        name: "bias last, block 16, chunks of 32",
        bias_first: false,
        block: 16,
        chunk: 32,
        channels_outside: false,
    },
    Order {
        name: "bias last, block 16, chunks of 64",
        bias_first: false,
        block: 16,
        chunk: 64,
        channels_outside: false,
    },
    Order {
        name: "bias last, block 16, chunks of 128",
        bias_first: false,
        block: 16,
        chunk: 128,
        channels_outside: false,
    },
];

/// `rows` images of `src`, `(rows, inputs, 1, width)`, convolved in `order`.
fn convolve_in(
    order: &Order,
    src: &[f32],
    width: usize,
    conv: &ConvWeights,
    rows: usize,
) -> Vec<f32> {
    let out_width = width + 1 - conv.width;
    let mut dst = vec![0f32; rows * conv.outputs * out_width];
    for b in 0..rows {
        let image = &src[b * conv.inputs * width..];
        for o in 0..conv.outputs {
            let weights = &conv.weight[o * conv.inputs * conv.width..];
            for x in 0..out_width {
                let term = |sum: f32, c: usize, k: usize| {
                    fused(image[c * width + x + k], weights[c * conv.width + k], sum)
                };
                let mut stored = 0f32;
                let mut chain = if order.bias_first { conv.bias[o] } else { 0.0 };
                let mut start = 0;
                while start < conv.inputs {
                    let end = start.saturating_add(order.chunk).min(conv.inputs);
                    let mut first = start;
                    while first < end {
                        let last = (first + order.block).min(end);
                        if order.channels_outside {
                            for c in first..last {
                                for k in 0..conv.width {
                                    chain = term(chain, c, k);
                                }
                            }
                        } else {
                            for k in 0..conv.width {
                                for c in first..last {
                                    chain = term(chain, c, k);
                                }
                            }
                        }
                        first = last;
                    }
                    if !order.bias_first {
                        stored = if start == 0 {
                            chain + conv.bias[o]
                        } else {
                            chain + stored
                        };
                        chain = 0.0;
                    }
                    start = end;
                }
                dst[(b * conv.outputs + o) * out_width + x] =
                    if order.bias_first { chain } else { stored };
            }
        }
    }
    dst
}

/// Each 1D convolution of the first rows of the first batch, from the reference's own input to it.
fn probe_1d_orders(model: &Weights, tensor: &[f32], reference: &HashMap<String, Vec<f32>>) {
    const ROWS: usize = 4;
    let batch = tensor.len() / (128 * 4);
    let mut input = vec![0f32; ROWS * 4 * 128];
    for b in 0..ROWS.min(batch) {
        for w in 0..128 {
            for c in 0..4 {
                input[(b * 4 + c) * 128 + w] = tensor[(b * 128 + w) * 4 + c];
            }
        }
    }
    let mut width = 128;
    let mut previous: Option<&str> = None;
    for layer in ["conv.0", "conv.2", "conv.5"] {
        let Some(want) = reference.get(layer) else {
            return;
        };
        let dimensions = model.shape(&format!("{layer}.weight"));
        let conv = ConvWeights {
            outputs: dimensions[0],
            inputs: dimensions[1],
            height: 1,
            width: *dimensions.last().unwrap(),
            weight: model.get(&format!("{layer}.weight")),
            bias: model.get(&format!("{layer}.bias")),
            bias_last: false,
        };
        if let Some(name) = previous {
            let before = &reference[name];
            input = before[..ROWS * conv.inputs * width].to_vec();
            relu(&mut input);
        }
        let out_width = width + 1 - conv.width;
        let want = &want[..ROWS * conv.outputs * out_width];
        for order in &ORDERS {
            let got = convolve_in(order, &input, width, &conv, ROWS);
            let missed = got
                .iter()
                .zip(want)
                .filter(|(a, b)| a.to_bits() != b.to_bits())
                .count();
            println!("probe\t{layer}\t{}\t{missed}/{}", order.name, want.len());
        }
        width = out_width;
        previous = Some(layer);
    }
}

fn report_tensor(
    records: &[gatk_tools::nv_score_variants::vcf::Record],
    want: &[f32],
    got: &[f32],
    read: bool,
) {
    let per_row = want.len() / records.len();
    let channels = if read { 15 } else { 4 };
    for (row, record) in records.iter().enumerate() {
        let (w, g) = (
            &want[row * per_row..(row + 1) * per_row],
            &got[row * per_row..(row + 1) * per_row],
        );
        if let Some(at) = w
            .iter()
            .zip(g)
            .position(|(a, b)| a.to_bits() != b.to_bits())
        {
            let cell = at / channels;
            eprintln!(
                "{}:{} {} {:?} first at cell {} ({} {}) channel {}: want {} got {}",
                record.contig,
                record.pos,
                record.reference,
                record.alternates,
                cell,
                cell / 128,
                cell % 128,
                at % channels,
                w[at],
                g[at]
            );
        }
    }
}
