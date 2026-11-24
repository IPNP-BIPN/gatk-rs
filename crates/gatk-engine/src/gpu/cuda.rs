//! The CUDA host side of [`super`], on Linux: compiles `kernels.cu` (generated from
//! `kernels.metal`, so the device arithmetic is the Metal backend's line for line) with NVRTC, once
//! per process, packs the inputs and reads the bits back. The driver and NVRTC are loaded at run
//! time, so a host without them answers `None` and the CPU path runs.

use std::io::Write;
use std::sync::{Arc, OnceLock};

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchConfig, PushKernelArg,
    ValidAsZeroBits,
};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};

use super::{ConvShape, ReadInputs};

const SOURCE: &str = include_str!("kernels.cu");
/// The scratch one launch may hold. A batch past it is split rather than refused.
const SCRATCH_BYTES: usize = 256 << 20;
/// Threads per block: each thread holds a whole software binary64, so fewer than the maximum.
const BLOCK: u32 = 128;

struct Context {
    context: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    f64_ops: CudaFunction,
    pair_hmm: CudaFunction,
    f32_ops: CudaFunction,
    conv2d_fused: CudaFunction,
}

static CONTEXT: OnceLock<Option<Context>> = OnceLock::new();

fn build() -> Option<Context> {
    let context = CudaContext::new(0).ok()?;
    // No contraction and no fast math. The kernels hold no floating point for either to touch,
    // and tools/audit/gpu_kernels.py keeps both facts true.
    let options = CompileOptions {
        fmad: Some(false),
        use_fast_math: Some(false),
        ..Default::default()
    };
    let ptx = compile_ptx_with_opts(SOURCE, options)
        .map_err(|e| writeln!(std::io::stderr(), "kernels.cu does not compile: {e:?}"))
        .ok()?;
    let module = context.load_module(ptx).ok()?;
    Some(Context {
        f64_ops: module.load_function("f64_ops").ok()?,
        pair_hmm: module.load_function("pair_hmm").ok()?,
        f32_ops: module.load_function("f32_ops").ok()?,
        conv2d_fused: module.load_function("conv2d_fused").ok()?,
        stream: context.default_stream(),
        context,
    })
}

fn with_context<T>(
    run: impl FnOnce(&Context) -> Result<T, cudarc::driver::DriverError>,
) -> Option<T> {
    // cudarc panics when it cannot load libcuda or libnvrtc; that host has no device path.
    let context = CONTEXT.get_or_init(|| std::panic::catch_unwind(build).ok().flatten());
    let result = run(context.as_ref()?);
    // A failed launch would leave the output as allocated, zeros, which would read back as
    // likelihoods of zero. The caller answers `None` and the CPU path runs instead.
    result
        .map_err(|e| {
            let _ = writeln!(
                std::io::stderr(),
                "gatk-rs: a CUDA command ended in error ({e:?}); this call runs on the CPU path"
            );
        })
        .ok()
}

pub fn device_name() -> Option<String> {
    with_context(|c| c.context.name())
}

fn upload<T: DeviceRepr + ValidAsZeroBits>(
    c: &Context,
    data: &[T],
) -> Result<CudaSlice<T>, cudarc::driver::DriverError> {
    if data.is_empty() {
        return zeros(c, 0);
    }
    c.stream.clone_htod(data)
}

fn zeros<T: DeviceRepr + ValidAsZeroBits>(
    c: &Context,
    n: usize,
) -> Result<CudaSlice<T>, cudarc::driver::DriverError> {
    // CUDA refuses an empty allocation; one word stands in, and no thread reads it.
    c.stream.alloc_zeros(n.max(1))
}

fn grid(threads: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: ((threads as u32).div_ceil(BLOCK), 1, 1),
        block_dim: (BLOCK, 1, 1),
        shared_mem_bytes: 0,
    }
}

fn read_back<T: DeviceRepr>(
    c: &Context,
    out: &CudaSlice<T>,
    n: usize,
) -> Result<Vec<T>, cudarc::driver::DriverError> {
    let mut values = c.stream.clone_dtoh(out)?;
    values.truncate(n);
    Ok(values)
}

pub fn f64_ops(operation: u32, a: &[u64], b: &[u64]) -> Option<Vec<u64>> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    with_context(|c| {
        let mut out = zeros::<u64>(c, n)?;
        if n > 0 {
            let (a, b) = (upload(c, a)?, upload(c, b)?);
            let count = n as u32;
            let mut launch = c.stream.launch_builder(&c.f64_ops);
            launch
                .arg(&a)
                .arg(&b)
                .arg(&mut out)
                .arg(&operation)
                .arg(&count);
            // SAFETY: the arguments are the kernel's, in order and of its types.
            unsafe { launch.launch(grid(n)) }?;
        }
        read_back(c, &out, n)
    })
}

pub fn pair_hmm_log10(haplotypes: &[Vec<u8>], reads: &[ReadInputs]) -> Option<Vec<f64>> {
    let p = super::PairHmmBatch::pack(haplotypes, reads);
    let per_dispatch = (SCRATCH_BYTES / 8 / p.stride).max(1);
    let flat = |spans: &[[u32; 2]]| spans.iter().flatten().copied().collect::<Vec<u32>>();
    with_context(|c| {
        let read_bases = upload(c, &p.read_bases)?;
        let read_spans = upload(c, &flat(&p.read_spans))?;
        let read_rows = upload(c, &p.read_rows)?;
        let haplotype_bases = upload(c, &p.haplotype_bases)?;
        let haplotype_spans = upload(c, &flat(&p.haplotype_spans))?;
        let haplotype_initial = upload(c, &p.haplotype_initial)?;
        let mut totals = Vec::with_capacity(p.jobs.len());
        for chunk in p.jobs.chunks(per_dispatch) {
            let jobs = upload(c, &flat(chunk))?;
            let mut scratch = zeros::<u64>(c, chunk.len() * p.stride)?;
            let mut out = zeros::<u64>(c, chunk.len())?;
            let (count, stride) = (chunk.len() as u32, p.stride as u32);
            let mut launch = c.stream.launch_builder(&c.pair_hmm);
            launch
                .arg(&read_bases)
                .arg(&read_spans)
                .arg(&read_rows)
                .arg(&haplotype_bases)
                .arg(&haplotype_spans)
                .arg(&haplotype_initial)
                .arg(&jobs)
                .arg(&mut scratch)
                .arg(&mut out)
                .arg(&count)
                .arg(&stride);
            // SAFETY: the arguments are the kernel's, in order; `uint2` and `PairJob` are pairs
            // of `u32`, which is how the spans and jobs are flattened.
            unsafe { launch.launch(grid(chunk.len())) }?;
            totals.extend(read_back(c, &out, chunk.len())?);
        }
        Ok(super::PairHmmBatch::log10(totals))
    })
}

pub fn f32_ops(operation: u32, a: &[u32], b: &[u32], c: &[u32]) -> Option<Vec<u32>> {
    assert!(a.len() == b.len() && b.len() == c.len());
    let n = a.len();
    with_context(|context| {
        let mut out = zeros::<u32>(context, n)?;
        if n > 0 {
            let (a, b, c) = (
                upload(context, a)?,
                upload(context, b)?,
                upload(context, c)?,
            );
            let count = n as u32;
            let mut launch = context.stream.launch_builder(&context.f32_ops);
            launch
                .arg(&a)
                .arg(&b)
                .arg(&c)
                .arg(&mut out)
                .arg(&operation)
                .arg(&count);
            // SAFETY: the arguments are the kernel's, in order and of its types.
            unsafe { launch.launch(grid(n)) }?;
        }
        read_back(context, &out, n)
    })
}

/// The outputs one launch computes, as on the Metal backend.
const CONV_OUTPUTS_PER_DISPATCH: usize = 1 << 18;

pub fn conv2d_fused(
    src: &[f32],
    weight: &[f32],
    bias: &[f32],
    shape: &ConvShape,
) -> Option<Vec<f32>> {
    let total = super::conv_outputs(shape);
    let bits = |v: &[f32]| v.iter().map(|x| x.to_bits()).collect::<Vec<u32>>();
    with_context(|c| {
        let (src, weight, bias) = (
            upload(c, &bits(src))?,
            upload(c, &bits(weight))?,
            upload(c, &bits(bias))?,
        );
        let mut dst = Vec::with_capacity(total);
        for first in (0..total).step_by(CONV_OUTPUTS_PER_DISPATCH) {
            let count = CONV_OUTPUTS_PER_DISPATCH.min(total - first);
            let (words, count_word) = super::conv_words(shape, first, count);
            let words = upload(c, &words)?;
            let mut out = zeros::<u32>(c, count)?;
            let mut launch = c.stream.launch_builder(&c.conv2d_fused);
            launch
                .arg(&src)
                .arg(&weight)
                .arg(&bias)
                .arg(&mut out)
                .arg(&words)
                .arg(&count_word);
            // SAFETY: the arguments are the kernel's, in order; `ConvShape` is eleven `u32`s.
            unsafe { launch.launch(grid(count)) }?;
            dst.extend(read_back(c, &out, count)?.into_iter().map(f32::from_bits));
        }
        Ok(dst)
    })
}
