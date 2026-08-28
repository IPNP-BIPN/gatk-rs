//! `gatk-rs <Tool> <args>`.
//!
//! `mainEntry` ends in `System.exit`, and so does this: the status is the whole of what a caller
//! downstream of it can see.

use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let outcome = gatk_cli::run(&args);
    print!("{}", outcome.stdout);
    let _ = std::io::stdout().flush();
    eprint!("{}", outcome.stderr);
    let _ = std::io::stderr().flush();
    std::process::exit(outcome.status);
}
