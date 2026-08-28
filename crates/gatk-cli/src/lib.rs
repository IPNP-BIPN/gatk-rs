//! The dispatcher: what `gatk <Tool> <args>` does before a tool ever sees an argument.
//!
//! `Main.mainEntry` is three steps and an exit status: route the command line, run whatever the
//! route named, and turn whatever came back into a status. The routing and the statuses are ported
//! in [`gatk_tools::main_entry`]; what is here is the entry point that puts them together, writes
//! to the two streams the reference writes to, and returns the status the reference would have
//! exited with.
//!
//! # What this does not do yet
//!
//! Two things, both of them measured and neither of them portable from what is measured:
//!
//!  * **the main usage listing** is three hundred and seventy-three lines of tool names and their
//!    one-line summaries, and the summaries are not in any golden: what is reproduced here is the
//!    stream it goes to, the status that follows it, and its first line;
//!  * **a tool's own arguments** need a Barclay definition each, and the declarations golden
//!    carries names rather than types, so no command line can be handed to the parser yet.
//!
//! Both are deliberate boundaries rather than omissions, and the test beside this file states
//! them as such: it compares what the port claims against the golden, and says nothing about the
//! lines the port does not produce.
//!
//! Ported from `org.broadinstitute.hellbender.Main`.

pub mod definitions;

use gatk_tools::main_entry::{self, Failure, Route, Stream};

/// The reference's own pins, which `printVersionInfo` reads off the jar's manifest.
pub const TOOLKIT_NAME: &str = "The Genome Analysis Toolkit (GATK)";
pub const TOOLKIT_VERSION: &str = "4.6.2.0";
pub const HTSJDK_VERSION: &str = "4.2.0";
pub const PICARD_VERSION: &str = "3.4.0";

/// The first line of the main usage, which Barclay writes in colour.
///
/// The escapes are the reference's: bold red for the label, green for the program name. They are
/// bytes of the output like any other, which is why they are written out rather than described.
pub const USAGE_FIRST_LINE: &str =
    "\u{1b}[1m\u{1b}[31mUSAGE:  \u{1b}[32m<program name>\u{1b}[1m\u{1b}[31m [-h]";

/// What one run wrote and what it returned, which is what `mainEntry` turns into an exit status.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Outcome {
    pub stdout: String,
    pub stderr: String,
    pub status: i32,
}

/// A tool the port can actually run, and what it does with the arguments after its name.
///
/// The list is empty of walkers on purpose: a tool becomes reachable here the moment its name
/// resolves, and runnable only once its arguments have a declaration to parse against. The gap
/// between the two is Milestone C's, and it is visible rather than papered over.
pub type Runner = fn(&[String]) -> Result<Option<String>, (Failure, String)>;

/// `Main.instanceMain`, with the two streams returned rather than written.
pub fn run(args: &[String]) -> Outcome {
    let mut stdout = String::new();
    let mut stderr = String::new();
    match main_entry::route(args) {
        Route::Usage(stream) => {
            // `handleResult(null)` prints nothing, and the run is a success: asking for help is
            // not an error, which is the one path that reaches the runner with no program.
            let target = match stream {
                Stream::Stdout => &mut stdout,
                Stream::Stderr => &mut stderr,
            };
            target.push_str(&usage());
            Outcome {
                stdout,
                stderr,
                status: 0,
            }
        }
        Route::Version => {
            let (first, rest) = main_entry::version_lines(
                TOOLKIT_NAME,
                TOOLKIT_VERSION,
                HTSJDK_VERSION,
                PICARD_VERSION,
            );
            // The first line goes to `System.out` whatever stream the method was handed, and the
            // two under it go to the stream. Here both are stdout, which is what the command line
            // asks for, and the split is only visible when a caller hands it something else.
            stdout.push_str(&first);
            stdout.push_str(&rest);
            Outcome {
                stdout,
                stderr,
                status: 0,
            }
        }
        Route::Unknown { message } => {
            // The refusal is the usage AND the message: the usage goes to stderr first, and the
            // UserException follows it through `handleUserException`.
            stderr.push_str(&usage());
            stderr.push_str(&main_entry::user_exception_report(&message));
            Outcome {
                stdout,
                stderr,
                status: main_entry::exit_status(Failure::User),
            }
        }
        Route::Tool { name } => {
            // The parse comes before the run, and it only happens for a tool whose whole argument
            // surface converts: a parser missing a third of its arguments would refuse a command
            // line the reference accepts, which is a worse answer than refusing to parse at all.
            if let Some(error) = parse_failure(&name, main_entry::tool_arguments(args)) {
                stderr.push_str(&main_entry::user_exception_report(&error));
                return Outcome {
                    stdout,
                    stderr,
                    status: main_entry::exit_status(Failure::CommandLine),
                };
            }
            match runner(&name) {
                None => {
                    stderr.push_str(&main_entry::user_exception_report(&not_ported(&name)));
                    Outcome {
                        stdout,
                        stderr,
                        status: main_entry::exit_status(Failure::User),
                    }
                }
                Some(runner) => match runner(main_entry::tool_arguments(args)) {
                    Ok(result) => {
                        if let Some(printed) = main_entry::tool_returned(result.as_deref()) {
                            stdout.push_str(&printed);
                            stdout.push('\n');
                        }
                        Outcome {
                            stdout,
                            stderr,
                            status: 0,
                        }
                    }
                    Err((failure, message)) => {
                        stderr.push_str(&main_entry::user_exception_report(&message));
                        Outcome {
                            stdout,
                            stderr,
                            status: main_entry::exit_status(failure),
                        }
                    }
                },
            }
        }
    }
}

/// Whether this port can hand a tool's command line to the ported Barclay parser at all.
///
/// A tool is parseable when every argument it declares converts to a definition. None of the seven
/// tools whose declarations are carried is, today: each of them declares a `GATKPath`, and the
/// conversion a `GATKPath` goes through is not measured. See [`definitions::UNCONVERTIBLE_CLASSES`].
pub fn parseable(tool: &str) -> bool {
    gatk_tools::tool_declarations::declarations(tool)
        .map(|list| definitions::missing(list).is_empty())
        .unwrap_or(false)
}

/// The refusal the ported parser makes of a command line, if the tool is parseable at all.
pub fn parse_failure(tool: &str, args: &[String]) -> Option<String> {
    if !parseable(tool) {
        return None;
    }
    let list = gatk_tools::tool_declarations::declarations(tool)?;
    let mut parser = gatk_barclay::Parser::new(definitions::definitions(list));
    let borrowed: Vec<&str> = args.iter().map(String::as_str).collect();
    parser
        .parse_arguments(&borrowed)
        .err()
        .map(|error| error.message)
}

/// The tools this port can run, which is none of them yet.
///
/// A name that resolves and has no runner is not a name that does not resolve: the reference would
/// have run it, and saying so is the honest answer. What it is NOT is a refusal the reference
/// makes, which is why [`not_ported`] says whose refusal it is.
pub fn runner(_name: &str) -> Option<Runner> {
    None
}

/// The port's own refusal, which the reference has no equivalent of.
///
/// It is worded so that nobody mistakes it for GATK's: a tool that GATK runs and this does not is
/// a gap in the port, and a golden will never contain this string.
pub fn not_ported(name: &str) -> String {
    format!("{name} is a GATK tool that this port does not carry yet. This message is the port's own and not GATK's.")
}

/// The main usage, of which the first line is the reference's and the rest is not yet.
///
/// The reference prints three hundred and seventy-three lines here: a header, then every tool it
/// found on the class path under its program group, each with the one-line summary its annotation
/// carries. The names are ported in [`gatk_tools::main_catalogue`]; the summaries are not measured,
/// so the listing is the names alone and the golden is not asked to agree with it.
pub fn usage() -> String {
    let mut text = String::from(USAGE_FIRST_LINE);
    text.push('\n');
    for name in gatk_tools::main_catalogue::CATALOGUE {
        text.push_str("    ");
        text.push_str(name);
        text.push('\n');
    }
    text
}
