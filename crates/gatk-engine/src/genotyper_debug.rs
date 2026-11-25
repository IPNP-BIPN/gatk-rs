//! `--debug-genotyper-output`, ported from `HaplotypeCallerGenotypingDebugger` (GATK 4.7.0.0;
//! `HaplotypeCallerGenotypingDebugger.java` is unchanged since 4.6.2.0).
//!
//! The reference keeps one static `PrintStream` that HaplotypeCaller's engine opens and every
//! layer prints to when it is set: the trimmer, the likelihood engines, read disqualification, the
//! genotyping engine and the DRAGEN models. Here it is a buffer per thread, enabled for a run and
//! taken when the run ends, so that a layer can print without being handed a writer.
//!
//! ```java
//! private static PrintStream genotyperDebugOutStream;
//! ```
//!
//! One difference that no output shows: the reference's `close()` closes the stream but leaves the
//! field set, so its `isEnabled()` stays true afterwards, where [`take`] turns printing off.
//!
//! # Reading this file without reading Rust
//!
//! - `thread_local!` declares a static with one copy per thread: Java's `ThreadLocal`, without
//!   the `get`/`set` calls. `OUT.with(|out| ...)` runs the closure (Java's lambda) on this thread's
//!   copy.
//! - `RefCell<Option<String>>` is a cell whose contents may be changed through a shared
//!   reference, checked at run time (`borrow` to read, `borrow_mut` to write). `Option<String>` is
//!   "a buffer or nothing", and nothing is the reference's `null` stream.
//! - `Option::take` moves the value out and leaves `None` behind, in one step.

use std::cell::RefCell;

thread_local! {
    // `None` until `initialize`: the reference's field before `initialize` has run.
    static OUT: RefCell<Option<String>> = const { RefCell::new(None) };
}

/// `initialize(location)`: printing starts.
///
/// The reference opens the file here; the port collects the text and the caller writes it, so
/// the location is not this function's business.
pub fn initialize() {
    OUT.with(|out| *out.borrow_mut() = Some(String::new()));
}

/// `isEnabled()`: `genotyperDebugOutStream != null`.
pub fn is_enabled() -> bool {
    OUT.with(|out| out.borrow().is_some())
}

/// `println(text)`, a no-op when not enabled.
///
/// The `'\n'` is `PrintStream.println`'s line separator on the Linux the goldens are taken on.
pub fn println(text: &str) {
    OUT.with(|out| {
        // `if let Some(buffer) = ...`: run the block only when a buffer is there, binding it.
        if let Some(buffer) = out.borrow_mut().as_mut() {
            buffer.push_str(text);
            buffer.push('\n');
        }
    });
}

/// `close()`: what was printed, and printing stops.
pub fn take() -> Option<String> {
    OUT.with(|out| out.borrow_mut().take())
}
