//! `--debug-genotyper-output`, ported from `HaplotypeCallerGenotypingDebugger` (GATK 4.6.2.0).
//!
//! The reference keeps one static `PrintStream` that HaplotypeCaller's engine opens and every
//! layer prints to when it is set: the trimmer, the likelihood engines, read disqualification, the
//! genotyping engine and the DRAGEN models. Here it is a buffer per thread, enabled for a run and
//! taken when the run ends, so that a layer can print without being handed a writer.

use std::cell::RefCell;

thread_local! {
    static OUT: RefCell<Option<String>> = const { RefCell::new(None) };
}

/// `initialize(location)`: printing starts.
pub fn initialize() {
    OUT.with(|out| *out.borrow_mut() = Some(String::new()));
}

/// `isEnabled()`.
pub fn is_enabled() -> bool {
    OUT.with(|out| out.borrow().is_some())
}

/// `println(text)`, a no-op when not enabled.
pub fn println(text: &str) {
    OUT.with(|out| {
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
