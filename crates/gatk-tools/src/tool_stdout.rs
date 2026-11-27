//! What a tool prints to standard output while it runs, as `System.out.println` does.
//!
//! The reference's debug switches (`--print-pileupcalling-status` among them) print straight to
//! standard output in the middle of a traversal. A runner returns its result at the end, so the
//! lines are kept here, per thread, and the command line puts them before anything the tool
//! returns.
//!
//! # Reading this file without reading Rust
//!
//! - `thread_local!` declares a static with one independent copy per thread, as Java's
//!   `ThreadLocal`; `OUT.with(|out| ...)` runs the closure on this thread's copy.
//! - `RefCell` lets a value that is not declared mutable be borrowed for writing at run time
//!   (`borrow_mut`), and panics if two such borrows overlap.

use std::cell::RefCell;

// `const { ... }` initialises each thread's copy without a lazy first-use check.
thread_local! {
    static OUT: RefCell<String> = const { RefCell::new(String::new()) };
}

/// `System.out.println(line)`.
pub fn println(line: &str) {
    OUT.with(|out| {
        let mut out = out.borrow_mut();
        out.push_str(line);
        out.push('\n');
    });
}

/// Everything printed since the last take, which empties it.
pub fn take() -> String {
    // `std::mem::take` moves the buffer out and leaves an empty `String` in its place, so the
    // text is returned without a copy.
    OUT.with(|out| std::mem::take(&mut *out.borrow_mut()))
}
