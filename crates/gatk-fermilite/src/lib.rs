//! fermi-lite, the assembler `FindBreakpointEvidenceSpark` calls through JNI, byte for byte.
//!
//! The library GATK 4.6.2.0 bundles answers `FermiLiteAssembler.getFermiLiteVersion()` with
//! 562638112b7992f4d7a243ccd9ab2b3bd08c6b08, r53, and that is the commit ported here. HEAD (r55)
//! differs only in returning nothing instead of dying when no read is left to assemble.
//!
//! [`assemble`] is `FermiLiteAssembler.createAssembly`: the buffer the Java side writes, the loop in
//! the JNI's C that reads it, `fml_assemble`, and the marshalling of what comes back.
//!
//! # The bundled library keeps its assertions, and an assertion kills the JVM
//!
//! The JNI Makefile passes `-O2 -Wall -fPIC` and no `-DNDEBUG`, so the 40-odd `assert`s call
//! `__assert_fail`. The one anybody will meet: when every read is filtered out, or there are no reads,
//! r53 reaches `mr_insert_multi` with an empty string, and `assert(len > 0 && ...)` aborts the
//! process. r55 fixed exactly that. Here each of those asserts is [`FmlError::Abort`], and it is on
//! the caller to decide what a dead JVM means for a byte-identical tool.
//!
//! # A `char` table index that is negative for some bytes
//!
//! `bfc_seq_conv` looks a base up with `seq_nt6_table[(int)s[i]]`, on a signed `char`, where
//! `worker_count` and `max_streak` cast to `uint8_t` first. Bytes from 0x80 read up to 128 bytes
//! before the table in `.data`. In `libfml.Linux.so` that is, going down from the table: 24 zero
//! bytes (0xE8 to 0xFF), the `__dso_handle` pointer (0xE0 to 0xE7), 16 zero bytes of padding (0xD0
//! to 0xDF) and the tail of `.got.plt` (0x80 to 0xCF). Zero gives `b = ob = (0 - 1) & 7 = 7`, not
//! 4, which is reproduced, and matters, as an uncorrectable base that is not "N" until it has been
//! reverse complemented twice. A pointer's bytes vary with ASLR and with lazy binding, and may even
//! decode as a base, so those bytes are refused as [`FmlError::Undefined`]. The low byte of
//! `__dso_handle` is the exception: page-aligned load address plus 0x226500, always 0.
//!
//! # Quals are read through `strdup`, and written back
//!
//! The JNI takes `l_seq` from the bases' `strlen` and `strdup`s the quals on their own. A NUL in the
//! quals leaves `qual` shorter than `l_seq`. `worker_count` and `bfc_seq_conv` then read past the
//! allocation and `bfc_ec1` writes past it, one byte per base, which is a heap overflow whose result
//! depends on what `malloc` put after it. That is refused as [`FmlError::Undefined`]. Note that a
//! read's qualities are not raw Phred scores here: a Phred 0 passed raw would be this NUL.
//!
//! A NUL in the bases is different and entirely defined: the JNI advances `readData` by the bases'
//! `strlen + 1` for both strings, so the parse goes out of step with the reads and carries on,
//! taking quals for bases. Java's zero-filled buffer makes that deterministic, so it is reproduced;
//! only running off the end of the buffer is refused. Quals longer than bases overflow Java's
//! buffer unless a shorter one elsewhere made room: [`FmlError::JavaBufferOverflow`].
//!
//! # Floats are x86-64's
//!
//! Three things rest on it. `(int)` of NaN is `INT_MIN` (`cvttsd2si`), and NaN is ordinary: `kcov`
//! is 0/0 when no k-mer reaches `min_cnt`, which is any assembly at low depth, and it goes through
//! `(int)(.1 * kcov + .499)` then a clamp. `min_elen` is `(int)(tot_len / n * 2.5 + .499)` with n
//! 0 as well, before the abort. And the two bubble poppers differ: `mag_vh_pop_simple` narrows
//! `n_diff` and `r_diff` to `float` and then compares `r_diff` with the `double` 0.1, while
//! `mag_v_pop_open` keeps both in `double`. Every such `float` here is one, and `.8` and `.7` are
//! fields of type `float` widened at the comparison, so 0.8 is 0.800000011920929.
//!
//! # What turned out not to matter
//!
//! `drand48` is imported, but only the `ks_shuffle_*` and `ks_sample_*` that `KSORT_INIT` expands to
//! call it, and nothing calls those. `qsort` is not imported at all, no `libm` function either
//! (`M_LN2` is a constant). `fm_verbose` is never set by the JNI, so nothing is global across
//! assemblies and there is no state to carry. No hash is iterated in an order that shows:
//! `bfc_ch_hist` only sums. `n_threads` stays at 1 in the JNI path, GATK never calls `setNThreads`,
//! and it would not be deterministic if it did: `bfc_kmer_bufclear` drops k-mers still buffered when
//! a worker finishes. So the threaded paths of `kt_for` and `mr_insert_multi` are not ported.
//!
//! The uninitialised memory that exists is not read before it is written: the `ec` and `absent` of
//! `ecbase_t`, `dummy` and `H1` in `ksw_i16`, `ok.info` in `unitig.c`. Two reads just outside an
//! array do nothing: `bc[c] -= l` with `c` still -1 and `l` 0 in `rle_insert_cached`, and `H0[-1]`
//! with an empty target in `ksw_i16`.

// Conditions stay as they are written in C, to be read side by side with it.
#![allow(clippy::explicit_auto_deref)]
#![allow(clippy::explicit_counter_loop)]
#![allow(clippy::if_same_then_else)]
#![allow(clippy::int_plus_one)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::needless_late_init)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::too_many_arguments)]

/// An `assert` that is compiled into the bundled library. Failure aborts the JVM there;
/// [`assemble`] turns it into [`FmlError::Abort`].
macro_rules! c_assert {
    ($cond:expr, $what:expr) => {
        if !$cond {
            std::panic::resume_unwind(Box::new($crate::CAbort($what)));
        }
    };
}

mod bfc;
mod bubble;
mod htab;
mod khash;
mod ksort;
mod ksw;
mod mag;
mod misc;
mod rld;
mod rope;
mod unitig;

use std::io::{self, Write};

/// `FML_COMMIT` of the bundled library.
pub const FML_COMMIT: &str = "562638112b7992f4d7a243ccd9ab2b3bd08c6b08";

/// `MAG_F_AGGRESSIVE`: pop variant bubbles (not default).
pub const MAG_F_AGGRESSIVE: i32 = 0x20;
/// `MAG_F_POPOPEN`: aggressive tip trimming (default).
pub const MAG_F_POPOPEN: i32 = 0x40;
/// `MAG_F_NO_SIMPL`: skip bubble simplification (default).
pub const MAG_F_NO_SIMPL: i32 = 0x80;

/// `magopt_t`, as `mag_init_opt` and `fml_opt_init` leave it. `trim_len` and `trim_depth` have no
/// field: `trim_len` is 0, for which `mag_g_trim_open` returns at once.
#[derive(Clone, Debug)]
pub(crate) struct MagOpt {
    pub(crate) flag: i32,
    pub(crate) min_ovlp: i32,
    pub(crate) min_elen: i32,
    pub(crate) min_ensr: i32,
    pub(crate) min_insr: i32,
    pub(crate) max_bdist: i32,
    pub(crate) max_bdiff: i32,
    pub(crate) max_bvtx: i32,
    pub(crate) min_merge_len: i32,
    pub(crate) min_dratio1: f32,
    pub(crate) max_bcov: f32,
    pub(crate) max_bfrac: f32,
}

/// `fml_opt_t`. Everything but the cleaning flag stays at what `fml_opt_init` says, because that is
/// the one option GATK changes; `n_threads` is 1 and there is no way to make it otherwise.
#[derive(Clone, Debug)]
pub struct FmlOpt {
    pub(crate) ec_k: i32,
    pub(crate) min_cnt: i32,
    pub(crate) max_cnt: i32,
    pub(crate) min_asm_ovlp: i32,
    pub(crate) min_merge_len: i32,
    pub(crate) mag_opt: MagOpt,
}

impl Default for FmlOpt {
    /// `fml_opt_init`, through `mag_init_opt`; `new FermiLiteAssembler()`.
    fn default() -> Self {
        FmlOpt {
            ec_k: 0,
            min_cnt: 4,
            max_cnt: 8,
            min_asm_ovlp: 33,
            min_merge_len: 0,
            mag_opt: MagOpt {
                flag: MAG_F_NO_SIMPL | MAG_F_POPOPEN,
                min_ovlp: 0,
                min_elen: 300,
                min_ensr: 4,
                min_insr: 3,
                max_bdist: 512,
                max_bdiff: 50,
                max_bvtx: 64,
                min_merge_len: 0,
                min_dratio1: 0.7,
                max_bcov: 10.,
                max_bfrac: 0.15,
            },
        }
    }
}

impl FmlOpt {
    /// `FermiLiteAssembler.setCleaningFlag`: `mag_opt.flag`, replaced, not or-ed in.
    pub fn with_cleaning_flag(mut self, flag: i32) -> Self {
        self.mag_opt.flag = flag;
        self
    }
    /// `FermiLiteAssembler.getCleaningFlag`.
    pub fn cleaning_flag(&self) -> i32 {
        self.mag_opt.flag
    }
}

/// `fml_ovlp_t` exactly as packed, the two ints `interpretAssemblyData` reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FmlOvlp {
    /// `len:31, from:1`: the overlap length, with the sign bit set when the overlap is at this
    /// contig's 3' end.
    pub raw_len: i32,
    /// `id:31, to:1`: the target contig, with the sign bit set when the overlap is at the target's
    /// 3' end.
    pub raw_id: i32,
}

impl FmlOvlp {
    /// `Connection.getTarget`, as an index into the contigs.
    pub fn target(&self) -> usize {
        (self.raw_id & i32::MAX) as usize
    }
    /// `Connection.getOverlapLen`.
    pub fn overlap_len(&self) -> i32 {
        self.raw_len & i32::MAX
    }
    /// `Connection.isRC`: the target is a predecessor, upstream of this contig's 5' end.
    pub fn is_rc(&self) -> bool {
        self.raw_len >= 0
    }
    /// `Connection.isTargetRC`.
    pub fn is_target_rc(&self) -> bool {
        self.raw_id < 0
    }
}

/// `fml_utg_t`, a `FermiLiteAssembly.Contig`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    /// `ACGT`.
    pub seq: Vec<u8>,
    /// Per-base coverage as `getPerBaseCoverage` hands it over: the count plus 33, at most 126.
    pub cov: Vec<u8>,
    /// `nsr`, `getNSupportingReads`.
    pub nsr: i32,
    /// The 5' connections then the 3' ones, in the order the JNI copies them.
    pub ovlp: Vec<FmlOvlp>,
}

/// What `createAssembly` does instead of returning.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FmlError {
    /// `assert` failed in the library, and the JVM aborts: the text says which one.
    Abort(&'static str),
    /// The bytes are such that C's result depends on the heap or on load addresses.
    Undefined(String),
    /// `java.nio.BufferOverflowException` in `makeReadData`: the quals outgrow the buffer sized by
    /// the bases.
    JavaBufferOverflow,
}

impl std::fmt::Display for FmlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FmlError::Abort(what) => write!(f, "libfml assertion fails, killing the JVM: {what}"),
            FmlError::Undefined(what) => write!(f, "undefined behaviour in libfml: {what}"),
            FmlError::JavaBufferOverflow => write!(f, "java.nio.BufferOverflowException"),
        }
    }
}

impl std::error::Error for FmlError {}

/// The payload of a failed [`c_assert!`].
pub(crate) struct CAbort(pub(crate) &'static str);

/// `bseq1_t`. `seq` and `qual` have the length `l_seq` and no terminator.
#[derive(Default)]
pub(crate) struct Bseq {
    pub(crate) l_seq: i32,
    pub(crate) seq: Vec<u8>,
    pub(crate) qual: Vec<u8>,
}

/// `seq_nt6_table`.
pub(crate) const SEQ_NT6_TABLE: [u8; 256] = {
    let mut t = [5u8; 256];
    t[b'A' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'G' as usize] = 3;
    t[b'T' as usize] = 4;
    t[b'a' as usize] = 1;
    t[b'c' as usize] = 2;
    t[b'g' as usize] = 3;
    t[b't' as usize] = 4;
    t
};

/// `seq_nt6_table[(int)c]` with `c` a signed `char`: what `libfml.Linux.so` has in `.data` under
/// the table, where it is the same every run. [`read_buffer`] has refused the rest.
#[inline]
pub(crate) fn seq_nt6_signed(c: u8) -> u8 {
    if c < 0x80 {
        SEQ_NT6_TABLE[c as usize]
    } else {
        0
    }
}

/// `seq_reverse`.
pub(crate) fn seq_reverse(s: &mut [u8]) {
    s.reverse();
}

/// `seq_revcomp6`.
pub(crate) fn seq_revcomp6(s: &mut [u8]) {
    s.reverse();
    for c in s.iter_mut() {
        if (1..=4).contains(c) {
            *c = 5 - *c;
        }
    }
}

/// x86-64 `cvttsd2si`, `(int)` of a `double`: out of range and NaN give `INT_MIN`, where Rust's
/// `as` would saturate.
#[inline]
pub(crate) fn cvttsd2si(x: f64) -> i32 {
    if x.is_nan() || x >= 2147483648.0 || x <= -2147483649.0 {
        i32::MIN
    } else {
        x as i32
    }
}

/// `strlen` from `pos`, refusing to leave `buf`, whose tail C would read as somebody else's heap.
fn strlen_in(buf: &[u8], pos: usize) -> Result<usize, FmlError> {
    buf.get(pos..)
        .and_then(|rest| rest.iter().position(|&b| b == 0))
        .ok_or_else(|| FmlError::Undefined("the JNI reads past the read buffer".to_string()))
}

/// `FermiLiteAssembler.makeReadData`: `allocateDirect` zero-fills, and its capacity counts the
/// bases only.
fn make_read_data(reads: &[(&[u8], &[u8])]) -> Result<Vec<u8>, FmlError> {
    let mut capacity = 4usize;
    for (bases, _) in reads {
        capacity += 2 * (bases.len() + 1);
    }
    let mut buf = vec![0u8; capacity];
    buf[0..4].copy_from_slice(&(reads.len() as i32).to_le_bytes());
    let mut at = 4;
    for (bases, quals) in reads {
        for s in [bases, quals] {
            let end = at + s.len() + 1;
            if end > capacity {
                return Err(FmlError::JavaBufferOverflow);
            }
            buf[at..end - 1].copy_from_slice(s);
            at = end;
        }
    }
    Ok(buf)
}

/// The first loop of `createAssemblyData`, and the refusals of what C does not define.
fn read_buffer(buf: &[u8], n: usize) -> Result<Vec<Bseq>, FmlError> {
    let mut seqs = Vec::with_capacity(n);
    let mut at = 4;
    for i in 0..n {
        let l = strlen_in(buf, at)?;
        let seq = buf[at..at + l].to_vec();
        let skip = l + 1;
        at += skip;
        let ql = strlen_in(buf, at)?;
        let mut qual = buf[at..at + ql].to_vec();
        at += skip;
        if ql < l {
            return Err(FmlError::Undefined(format!(
                "read {i} has {l} bases but {ql} quals before a NUL: bfc reads and writes the heap after them"
            )));
        }
        qual.truncate(l);
        if let Some(c) = seq
            .iter()
            .find(|&&c| (0x80..=0xCF).contains(&c) || (0xE1..=0xE7).contains(&c))
        {
            return Err(FmlError::Undefined(format!(
                "read {i} has the base byte {c:#x}, which indexes a pointer below seq_nt6_table"
            )));
        }
        seqs.push(Bseq {
            l_seq: l as i32,
            seq,
            qual,
        });
    }
    Ok(seqs)
}

/// `FermiLiteAssembler.createAssembly(reads)`, each read its `getBases()` and `getQuals()`.
pub fn assemble(opt: &FmlOpt, reads: &[(&[u8], &[u8])]) -> Result<Vec<Contig>, FmlError> {
    let buf = make_read_data(reads)?;
    let seqs = read_buffer(&buf, reads.len())?;
    drop(buf);
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        misc::fml_assemble(opt, seqs)
    })) {
        Ok(contigs) => Ok(contigs),
        Err(payload) => match payload.downcast::<CAbort>() {
            Ok(a) => Err(FmlError::Abort(a.0)),
            Err(other) => std::panic::resume_unwind(other),
        },
    }
}

/// One assembly as `tools/fermilite-conformance/FermiLiteDump.java` prints it.
pub fn write_dump(
    out: &mut dyn Write,
    index: u64,
    name: &str,
    flag: i32,
    n_reads: usize,
    contigs: &[Contig],
) -> io::Result<()> {
    writeln!(
        out,
        "A {index} {name} flag={flag} reads={n_reads} contigs={}",
        contigs.len()
    )?;
    for (i, c) in contigs.iter().enumerate() {
        write!(out, "C {i} {} {} ", c.seq.len(), c.nsr)?;
        out.write_all(&c.seq)?;
        out.write_all(b" ")?;
        out.write_all(&c.cov)?;
        out.write_all(b"\n")?;
        for o in &c.ovlp {
            writeln!(
                out,
                "L {} {} {} {}",
                o.target(),
                o.overlap_len(),
                u8::from(o.is_rc()),
                u8::from(o.is_target_rc()),
            )?;
        }
    }
    Ok(())
}

/// The line in place of an assembly the JVM did not survive.
pub fn write_dump_abort(
    out: &mut dyn Write,
    index: u64,
    name: &str,
    flag: i32,
    n_reads: usize,
) -> io::Result<()> {
    writeln!(out, "X {index} {name} flag={flag} reads={n_reads} ABORT")
}
