//! `SplitCRAM`, ported from `org.broadinstitute.hellbender.tools.SplitCRAM` (GATK 4.6.2.0).
//!
//! A CRAM cut into shards without a single record being decoded. The tool walks containers, reads
//! the record count each one declares in its header, and starts a new output once the count
//! reaches a threshold. The bytes it copies are htsjdk's; what is the tool's own is where the cuts
//! fall and what the outputs are called, which is all that is here.
//!
//! # The threshold is a minimum
//!
//! ```java
//! while (cramContainerIterator.hasNext() && (records < shardRecords)) {
//!     final Container container = cramContainerIterator.next();
//!     container.write(cramHeader.getCRAMVersion(), os);
//!     records += container.getContainerHeader().getNumberOfRecords();
//! }
//! ```
//!
//! The test comes before the container is read and the count is added after, so a shard overshoots
//! by up to one whole container. The test being strict, a threshold exactly the size of a container
//! still gives one container per shard: the first container takes the running count to the
//! threshold and the loop stops. One above it takes a second container.
//!
//! # `--shard-max-output-count` does not limit anything above one
//!
//! ```java
//! while (cramContainerIterator.hasNext()) {
//!     int shardOuputCount = 0;
//!     ...
//!     shardOuputCount++;
//!     if ( shardMaxOutputCount != 0 && shardOuputCount >= shardMaxOutputCount ) { break; }
//! }
//! ```
//!
//! The counter is declared inside the outer loop, so every shard resets it and it is always one
//! when tested. Only the value 1 ever stops the run: 2 and 3 leave every shard standing, which is
//! what the golden shows over a five-container input. [`plan`] keeps the reset where the reference
//! put it rather than repairing it.
//!
//! # The name is a format string, checked before anything is read
//!
//! `onStartup` refuses a template that `%[0-9]*d` cannot find anywhere in it, and does so with a
//! bare `IllegalArgumentException` rather than a `UserException`. So `%04d` passes, and a width
//! carrying a flag such as `%-4d` does not.

/// `SplitCRAM.DEFAULT_SHARD_RECORDS`.
pub const DEFAULT_SHARD_RECORDS: i64 = 10_000_000;

/// The default `--output`, relative to the working directory.
pub const DEFAULT_TEMPLATE: &str = "output_%04d.cram";

/// What the run refuses, which it does before opening the input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitError {
    TemplateMissingFormatter { template: String },
}

impl SplitError {
    pub fn java_class(&self) -> &str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> String {
        match self {
            SplitError::TemplateMissingFormatter { template } => {
                format!("output template missing a %d enumerator formatter: {template}")
            }
        }
    }
}

/// One `%[0-9]*d` found in a template: where it is, and the width it asks for.
struct Formatter {
    start: usize,
    end: usize,
    width: usize,
    zero_padded: bool,
}

/// `SplitCRAM.numeratorFormat.matcher(template).find()`, which is a search rather than a match, so
/// the formatter can be anywhere in the name.
fn find_formatter(template: &str) -> Option<Formatter> {
    let bytes = template.as_bytes();
    for start in 0..bytes.len() {
        if bytes[start] != b'%' {
            continue;
        }
        let mut end = start + 1;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
        if end < bytes.len() && bytes[end] == b'd' {
            let digits = &template[start + 1..end];
            return Some(Formatter {
                start,
                end: end + 1,
                width: digits.parse().unwrap_or(0),
                zero_padded: digits.starts_with('0'),
            });
        }
    }
    None
}

/// `onStartup`'s check.
pub fn accepts_template(template: &str) -> bool {
    find_formatter(template).is_some()
}

/// `String.format(cramOutputTemplate, shard)` for the one conversion the check allows.
pub fn format_name(template: &str, shard: i32) -> Option<String> {
    let formatter = find_formatter(template)?;
    let digits = shard.to_string();
    let padded = if digits.len() >= formatter.width {
        digits
    } else {
        let fill = if formatter.zero_padded { '0' } else { ' ' };
        let mut padded = String::new();
        for _ in digits.len()..formatter.width {
            padded.push(fill);
        }
        padded.push_str(&digits);
        padded
    };
    Some(format!(
        "{}{padded}{}",
        &template[..formatter.start],
        &template[formatter.end..]
    ))
}

/// One output: its name, and the record count of every container that went into it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Shard {
    pub name: String,
    pub containers: Vec<i32>,
}

impl Shard {
    /// What the inner loop counted, which is what decided the cut after it.
    pub fn records(&self) -> i64 {
        self.containers.iter().map(|count| *count as i64).sum()
    }
}

/// `doWork`: the containers of one CRAM grouped into shards.
///
/// `containers` is the record count each container declares, in file order. `max_output_count` is
/// `--shard-max-output-count`, where 0 is off and, through the reset above, every value but 1 is
/// off as well.
pub fn plan(
    containers: &[i32],
    shard_records: i64,
    max_output_count: i32,
    template: &str,
) -> Result<Vec<Shard>, SplitError> {
    if !accepts_template(template) {
        return Err(SplitError::TemplateMissingFormatter {
            template: template.to_string(),
        });
    }
    let mut shards = Vec::new();
    let mut shard = 0;
    let mut index = 0;
    // An empty CRAM never enters the outer loop, so it produces no shard rather than an empty one.
    while index < containers.len() {
        let name = format_name(template, shard).expect("the template was checked");
        shard += 1;
        // Declared here in the reference, and therefore reset for every shard.
        let mut shard_output_count = 0;
        let mut group = Vec::new();
        let mut records: i64 = 0;
        while index < containers.len() && records < shard_records {
            group.push(containers[index]);
            records += containers[index] as i64;
            index += 1;
        }
        shards.push(Shard {
            name,
            containers: group,
        });
        shard_output_count += 1;
        if max_output_count != 0 && shard_output_count >= max_output_count {
            break;
        }
    }
    Ok(shards)
}
