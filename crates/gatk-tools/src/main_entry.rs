//! `Main`'s entry: which stream each path writes to, what it returns, and what status an
//! exception carries.
//!
//! The class path search and the Barclay parse are not ported here; [`crate::main_catalogue`]
//! holds the names a token resolves against and [`crate::main_dispatch`] holds what happens when
//! one does not. What is ported is the routing between them, the five exit statuses, and the two
//! handlers `mainEntry` calls on its way out.
//!
//! Ported from `org.broadinstitute.hellbender.Main`.

use crate::main_catalogue;
use crate::main_dispatch;

/// A `CommandLineException`, which is a parse that failed.
pub const COMMANDLINE_EXCEPTION_EXIT_VALUE: i32 = 1;
/// A `UserException`, which is everything the user could have written differently.
pub const USER_EXCEPTION_EXIT_VALUE: i32 = 2;
/// Anything else, which is a bug rather than a refusal.
pub const ANY_OTHER_EXCEPTION_EXIT_VALUE: i32 = 3;
/// A Picard tool that returned non-zero, which arrives wrapped.
pub const PICARD_TOOL_EXCEPTION: i32 = 4;
/// An `OutOfMemoryError`, whose status is the shell's own convention for a fatal signal.
pub const OUT_OF_MEMORY_EXIT_VALUE: i32 = 137;

/// The property `handleUserException` names when it declines to print a stack trace.
pub const STACK_TRACE_ON_USER_EXCEPTION_PROPERTY: &str = "GATK_STACKTRACE_ON_USER_EXCEPTION";

/// The banner `printDecoratedExceptionMessage` puts around a message: seventy-one asterisks.
pub const BANNER: &str = "***********************************************************************";

/// The prefix `handleUserException` hands the decorator.
pub const USER_ERROR_PREFIX: &str = "A USER ERROR has occurred: ";

/// Which stream a path writes its usage to, which is the difference between asking for help and
/// being refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stream {
    Stdout,
    Stderr,
}

/// What `extractCommandLineProgram` decided, before anything has been printed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Route {
    /// The usage, and the run returns null: the one path that reaches `runCommandLineProgram`
    /// with no program.
    Usage(Stream),
    /// The version, printed instead of running anything.
    Version,
    /// A tool, whose arguments are everything after its name.
    Tool { name: String },
    /// The usage on stderr and then a `UserException` carrying this message.
    Unknown { message: String },
}

/// Which exception a path ends in, which is what the status is read off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Failure {
    CommandLine,
    User,
    PicardNonZero,
    OutOfMemory,
    Other,
}

/// `mainEntry`'s catch blocks, in the order they are written.
pub fn exit_status(failure: Failure) -> i32 {
    match failure {
        Failure::CommandLine => COMMANDLINE_EXCEPTION_EXIT_VALUE,
        Failure::User => USER_EXCEPTION_EXIT_VALUE,
        Failure::PicardNonZero => PICARD_TOOL_EXCEPTION,
        Failure::OutOfMemory => OUT_OF_MEMORY_EXIT_VALUE,
        Failure::Other => ANY_OTHER_EXCEPTION_EXIT_VALUE,
    }
}

/// `extractCommandLineProgram`'s three questions, in its own order.
///
/// The help test looks at the FIRST argument alone, so a tool name followed by `-h` resolves the
/// tool and lets it answer for itself. The version test looks at EVERY argument, so the same tool
/// name followed by `--version` never runs. And a name that resolves to nothing takes the usage
/// with it: the refusal is the usage on stderr and the message both.
pub fn route(args: &[String]) -> Route {
    if args.is_empty() || args[0] == "-h" || args[0] == "--help" {
        return Route::Usage(Stream::Stdout);
    }
    if args
        .iter()
        .any(|arg| arg == "-version" || arg == "--version")
    {
        return Route::Version;
    }
    if main_catalogue::resolves(&args[0]) {
        return Route::Tool {
            name: args[0].clone(),
        };
    }
    let catalogue: Vec<String> = main_catalogue::CATALOGUE
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    Route::Unknown {
        message: main_dispatch::unknown_command_message(&catalogue, &args[0])
            .expect("a name that resolves is not this path"),
    }
}

/// `runCommandLineProgram`: the first argument is the tool's name and never reaches the tool.
pub fn tool_arguments(args: &[String]) -> &[String] {
    if args.is_empty() {
        args
    } else {
        &args[1..]
    }
}

/// `handleResult`, which is silent for null and prints the value under a line of its own
/// otherwise.
pub fn tool_returned(result: Option<&str>) -> Option<String> {
    result.map(|value| format!("Tool returned:\n{value}"))
}

/// `printDecoratedExceptionMessage`, whose banner is printed twice with a blank line inside each
/// side of the message.
pub fn decorated_exception_message(prefix: &str, message: &str) -> String {
    format!("{BANNER}\n\n{prefix}{message}\n\n{BANNER}\n")
}

/// `handleUserException`: the decorated message, and then the property that would have printed a
/// stack trace.
pub fn user_exception_report(message: &str) -> String {
    format!(
        "{}Set the system property {STACK_TRACE_ON_USER_EXCEPTION_PROPERTY} \
         (--java-options '-D{STACK_TRACE_ON_USER_EXCEPTION_PROPERTY}=true') to print the stack \
         trace.\n",
        decorated_exception_message(USER_ERROR_PREFIX, message)
    )
}

/// `printVersionInfo`, which is split across two streams.
///
/// The first line goes to `System.out` whatever stream the method is handed, and only the two
/// lines under it go to the stream, so printing the version anywhere but stdout tears it in half.
/// The three versions are the reference's own pins.
pub fn version_lines(toolkit: &str, version: &str, htsjdk: &str, picard: &str) -> (String, String) {
    (
        format!("{toolkit} v{version}\n"),
        format!("HTSJDK Version: {htsjdk}\nPicard Version: {picard}\n"),
    )
}
