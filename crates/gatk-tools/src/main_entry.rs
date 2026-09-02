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

/// Which of the two handlers `mainEntry` gives a throwable to.
///
/// `handleUserException` catches `CommandLineException` and `UserException`; everything else,
/// `OutOfMemoryError` included, falls to `handleNonUserException`. The two write different things,
/// so the answer is not a detail of the status.
pub fn is_user_exception(failure: Failure) -> bool {
    matches!(failure, Failure::CommandLine | Failure::User)
}

/// What a run threw: which handler it reaches, the class the reference would name, and the message.
///
/// The class is carried because `handleNonUserException` PRINTS it. `handleUserException` does not,
/// so a user exception's class is here for the same reason its status is: it is what the reference
/// threw, and a port that only kept the message could not tell the two handlers apart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Thrown {
    pub failure: Failure,
    /// The exception's BINARY name, which is what `Throwable.toString` prints: a nested class
    /// carries a `$` where its source spells a dot.
    pub exception: &'static str,
    /// `getLocalizedMessage`, which is `None` where the reference's is null. An empty message is
    /// not the same thing and is `Some("")`.
    pub message: Option<String>,
}

/// `CommandLineException`, which is a parse that failed. Barclay's, not GATK's: `Main` imports it
/// from `org.broadinstitute.barclay.argparser`, and it is no `UserException`.
pub const COMMANDLINE_EXCEPTION: &str = "org.broadinstitute.barclay.argparser.CommandLineException";
/// `UserException`, which is everything the user could have written differently.
pub const USER_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.UserException";
/// A capability the reference has and this port does not, which is no exception of the reference's
/// at all. The name is deliberately not a Java one: `::` cannot appear in a binary class name, so
/// a reader and a diff both see at once that this line is the port's and not GATK's.
pub const PORT_LIMITATION: &str = "gatk_rs::PortLimitation";
/// The port's own plumbing failing where no golden says what the reference does, which is a third
/// thing again: not a refusal, and not a feature that was never carried.
pub const PORT_FAILURE: &str = "gatk_rs::PortFailure";

impl Thrown {
    /// A `CommandLineException`: the usage, then the decorated message, then status one.
    pub fn command_line(message: impl Into<String>) -> Self {
        Self {
            failure: Failure::CommandLine,
            exception: COMMANDLINE_EXCEPTION,
            message: Some(message.into()),
        }
    }

    /// A `UserException`: the decorated message, then status two.
    pub fn user(message: impl Into<String>) -> Self {
        Self {
            failure: Failure::User,
            exception: USER_EXCEPTION,
            message: Some(message.into()),
        }
    }

    /// Anything else, which is a bug rather than a refusal: the class, the message, status three.
    pub fn non_user(exception: &'static str, message: impl Into<String>) -> Self {
        Self {
            failure: Failure::Other,
            exception,
            message: Some(message.into()),
        }
    }

    /// Whether `handleUserException` is the handler this reaches.
    pub fn is_user(&self) -> bool {
        is_user_exception(self.failure)
    }

    /// The status `mainEntry` exits with.
    pub fn status(&self) -> i32 {
        exit_status(self.failure)
    }

    /// What the handler writes to stderr, which is not the same text for the two of them.
    pub fn report(&self) -> String {
        if self.is_user() {
            user_exception_report(self.message.as_deref().unwrap_or_default())
        } else {
            non_user_exception_report(self.exception, self.message.as_deref())
        }
    }
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

/// `handleNonUserException`, whose whole body is `printStackTrace`.
///
/// The first line is `Throwable.toString()`: the class's binary name and, where there is one, `: `
/// and the message. There is no banner, no `A USER ERROR has occurred:` prefix and no notice about
/// a system property, which is what separates this path from the other one and what a port with a
/// single banner gets wrong on every non-user refusal (`main-non-user`).
///
/// The frames under that line are the reference's own stack and this port has no equivalent to
/// print. That is a boundary rather than an omission: the suite records the first line, which is
/// what a row-by-row comparison reads, and that everything after it was a frame.
pub fn non_user_exception_report(exception: &str, message: Option<&str>) -> String {
    match message {
        Some(message) => format!("{exception}: {message}\n"),
        None => format!("{exception}\n"),
    }
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
