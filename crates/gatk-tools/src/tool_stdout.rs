//! What a tool prints to standard output while it runs, as `System.out.println` does.
//!
//! The reference's debug switches (`--print-pileupcalling-status` among them) print straight to
//! standard output in the middle of a traversal. A runner returns its result at the end, so the
//! lines are kept here, per thread, and the command line puts them before anything the tool
//! returns.

use std::cell::RefCell;

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
    OUT.with(|out| std::mem::take(&mut *out.borrow_mut()))
}
