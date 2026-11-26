//! The Metal host side of [`super`]: compiles `kernels.metal` once per thread, packs the inputs
//! and reads the bits back.

use std::io::Write;

use metal::objc::rc::autoreleasepool;
use metal::objc::runtime::Object;
use metal::objc::{msg_send, sel, sel_impl};
use metal::{
    Buffer, CommandBufferRef, CommandQueue, CompileOptions, ComputePipelineState, Device,
    MTLCommandBufferStatus, MTLResourceOptions, MTLSize,
};

use super::{ConvShape, ReadInputs};

const SOURCE: &str = include_str!("kernels.metal");
/// The scratch one dispatch may hold. A batch past it is split rather than refused.
const SCRATCH_BYTES: usize = 256 << 20;

struct Context {
    device: Device,
    queue: CommandQueue,
    f64_ops: ComputePipelineState,
    pair_hmm: ComputePipelineState,
    f32_ops: ComputePipelineState,
    conv2d_fused: ComputePipelineState,
}

thread_local! {
    // Metal's objects are not `Send` in these bindings, so each thread compiles its own.
    static CONTEXT: Option<Context> = autoreleasepool(build);
}

fn build() -> Option<Context> {
    let device = Device::system_default()?;
    let options = CompileOptions::new();
    // Fast math would license reassociation and contraction. The kernels hold no floating point
    // for it to touch, and tools/audit/gpu_kernels.py keeps both facts true.
    options.set_fast_math_enabled(false);
    let library = device
        .new_library_with_source(SOURCE, &options)
        .map_err(|e| writeln!(std::io::stderr(), "kernels.metal does not compile: {e}"))
        .ok()?;
    let pipeline = |name: &str| {
        let function = library.get_function(name, None).ok()?;
        device
            .new_compute_pipeline_state_with_function(&function)
            .ok()
    };
    Some(Context {
        f64_ops: pipeline("f64_ops")?,
        pair_hmm: pipeline("pair_hmm")?,
        f32_ops: pipeline("f32_ops")?,
        conv2d_fused: pipeline("conv2d_fused")?,
        queue: device.new_command_queue(),
        device,
    })
}

fn with_context<T>(run: impl FnOnce(&Context) -> Option<T>) -> Option<T> {
    CONTEXT.with(|context| context.as_ref().and_then(|c| autoreleasepool(|| run(c))))
}

/// Why a command buffer ended in error, in Metal's words.
// `msg_send!` tests a `cargo-clippy` feature this crate does not declare.
#[allow(unexpected_cfgs)]
fn failure(command: &CommandBufferRef) -> String {
    // SAFETY: `error` and `localizedDescription` are NSError and NSString accessors that return
    // nil or an object the autorelease pool owns; `UTF8String` is valid while that object is.
    unsafe {
        let error: *mut Object = msg_send![command, error];
        if error.is_null() {
            return "no error given".to_string();
        }
        let description: *mut Object = msg_send![error, localizedDescription];
        let text: *const std::os::raw::c_char = msg_send![description, UTF8String];
        std::ffi::CStr::from_ptr(text)
            .to_string_lossy()
            .into_owned()
    }
}

pub fn device_name() -> Option<String> {
    with_context(|c| Some(c.device.name().to_string()))
}

fn buffer<T: Copy>(device: &Device, data: &[T]) -> Buffer {
    // Metal refuses an empty buffer; one zeroed word stands in, and no thread reads it.
    let bytes = std::mem::size_of_val(data);
    if bytes == 0 {
        return device.new_buffer(8, MTLResourceOptions::StorageModeShared);
    }
    device.new_buffer_with_data(
        data.as_ptr().cast(),
        bytes as u64,
        MTLResourceOptions::StorageModeShared,
    )
}

fn read_back<T: Copy>(buffer: &Buffer, count: usize) -> Vec<T> {
    // SAFETY: the buffer is shared storage of at least `count` values of `T`, and the command that
    // wrote it has completed.
    unsafe { std::slice::from_raw_parts(buffer.contents() as *const T, count).to_vec() }
}

fn run(
    context: &Context,
    pipeline: &ComputePipelineState,
    buffers: &[&Buffer],
    words: &[u32],
    threads: usize,
) -> Option<()> {
    let command = context.queue.new_command_buffer();
    let encoder = command.new_compute_command_encoder();
    encoder.set_compute_pipeline_state(pipeline);
    for (index, b) in buffers.iter().enumerate() {
        encoder.set_buffer(index as u64, Some(b), 0);
    }
    for (offset, word) in words.iter().enumerate() {
        encoder.set_bytes(
            (buffers.len() + offset) as u64,
            4,
            (word as *const u32).cast(),
        );
    }
    let width = pipeline.thread_execution_width();
    encoder.dispatch_threads(
        MTLSize::new(threads as u64, 1, 1),
        MTLSize::new(width.min(threads as u64).max(1), 1, 1),
    );
    encoder.end_encoding();
    command.commit();
    command.wait_until_completed();
    // A command that fails leaves its output as it was allocated, zeros, which would read back as
    // likelihoods of zero. The caller then answers `None` and the CPU path runs instead.
    if command.status() == MTLCommandBufferStatus::Completed {
        return Some(());
    }
    // Straight to stderr rather than `eprintln!`, which libtest captures in a passing test, so
    // the `gpu` job can count these.
    let _ = writeln!(
        std::io::stderr(),
        "gatk-rs: a Metal command ended {:?} ({}); this call runs on the CPU path",
        command.status(),
        failure(command)
    );
    None
}

pub fn f64_ops(operation: u32, a: &[u64], b: &[u64]) -> Option<Vec<u64>> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    with_context(|c| {
        let out = buffer(&c.device, &vec![0u64; n]);
        if n > 0 {
            let (a, b) = (buffer(&c.device, a), buffer(&c.device, b));
            run(c, &c.f64_ops, &[&a, &b, &out], &[operation, n as u32], n)?;
        }
        Some(read_back(&out, n))
    })
}

pub fn pair_hmm_log10(haplotypes: &[Vec<u8>], reads: &[ReadInputs]) -> Option<Vec<f64>> {
    let p = super::PairHmmBatch::pack(haplotypes, reads);
    let per_dispatch = (SCRATCH_BYTES / 8 / p.stride).max(1);
    with_context(|c| {
        let d = &c.device;
        let inputs = [
            buffer(d, &p.read_bases),
            buffer(d, &p.read_spans),
            buffer(d, &p.read_rows),
            buffer(d, &p.haplotype_bases),
            buffer(d, &p.haplotype_spans),
            buffer(d, &p.haplotype_initial),
        ];
        let mut totals = Vec::with_capacity(p.jobs.len());
        for chunk in p.jobs.chunks(per_dispatch) {
            let jobs = buffer(d, chunk);
            let scratch = d.new_buffer(
                (chunk.len() * p.stride * 8) as u64,
                MTLResourceOptions::StorageModePrivate,
            );
            let out = buffer(d, &vec![0u64; chunk.len()]);
            let mut bound: Vec<&Buffer> = inputs.iter().collect();
            bound.extend([&jobs, &scratch, &out]);
            run(
                c,
                &c.pair_hmm,
                &bound,
                &[chunk.len() as u32, p.stride as u32],
                chunk.len(),
            )?;
            totals.extend(read_back::<u64>(&out, chunk.len()));
        }
        Some(super::PairHmmBatch::log10(totals))
    })
}

pub fn f32_ops(operation: u32, a: &[u32], b: &[u32], c: &[u32]) -> Option<Vec<u32>> {
    assert!(a.len() == b.len() && b.len() == c.len());
    let n = a.len();
    with_context(|context| {
        let out = buffer(&context.device, &vec![0u32; n]);
        if n > 0 {
            let d = &context.device;
            let (a, b, c) = (buffer(d, a), buffer(d, b), buffer(d, c));
            run(
                context,
                &context.f32_ops,
                &[&a, &b, &c, &out],
                &[operation, n as u32],
                n,
            )?;
        }
        Some(read_back(&out, n))
    })
}

/// The outputs one dispatch computes, so no command buffer runs long enough for the watchdog.
const CONV_OUTPUTS_PER_DISPATCH: usize = 1 << 18;

pub fn conv2d_fused(
    src: &[f32],
    weight: &[f32],
    bias: &[f32],
    shape: &ConvShape,
) -> Option<Vec<f32>> {
    let total = super::conv_outputs(shape);
    with_context(|c| {
        let d = &c.device;
        let (src, weight, bias) = (buffer(d, src), buffer(d, weight), buffer(d, bias));
        let mut dst = Vec::with_capacity(total);
        for first in (0..total).step_by(CONV_OUTPUTS_PER_DISPATCH) {
            let count = CONV_OUTPUTS_PER_DISPATCH.min(total - first);
            let (words, count_word) = super::conv_words(shape, first, count);
            let words = buffer(d, &words);
            let out = buffer(d, &vec![0u32; count]);
            run(
                c,
                &c.conv2d_fused,
                &[&src, &weight, &bias, &out, &words],
                &[count_word],
                count,
            )?;
            dst.extend(read_back::<f32>(&out, count));
        }
        Some(dst)
    })
}
