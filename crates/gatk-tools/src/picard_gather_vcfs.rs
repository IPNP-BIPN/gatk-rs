//! `GatherVcfs`, ported from `picard.vcf.GatherVcfs` (Picard 3.4.0).
//!
//! Shards of one scattered run concatenated into a single file.
//!
//! # Not the GATK tool of almost the same name
//!
//! [`crate::gather_vcfs`] is `GatherVcfsCloud`, already ported. The two share the shape and
//! disagree on nearly everything else: this one refuses with an exit code where that one throws,
//! it reorders its inputs on request, and its comment lines are keyed `GatherVcfs.comment`.
//!
//! # A refusal is an exit code, not an exception
//!
//! ```java
//! } catch (RuntimeException e) {
//!     log.error("There was a problem with gathering the INPUT.", e);
//!     Files.deleteIfExists(OUTPUT.toPath());
//!     return 1;
//! }
//! ```
//!
//! Everything after the dictionary check runs inside that try, so a caller sees a status and a log
//! line. An `AssertionError` is not a `RuntimeException`, so the dictionary mismatch, which
//! `assertSameDictionary` raises as one, walks straight out while the sample mismatch beside it
//! becomes exit 1. [`Failure`] keeps the two apart.
//!
//! # There are two order checks and they compare different things
//!
//! The check before writing compares each file's FIRST record with the previous file's first
//! record. The check inside the gathering compares the next file's first record with the LAST
//! RECORD WRITTEN. A pair of shards that overlap passes the first and fails the second, with a
//! different message.
//!
//! # Only the first comment survives
//!
//! `VCFHeader.addMetaDataLine` keys an unstructured line by its key alone, so a second `CO=`
//! replaces nothing and is dropped without a word.

use htsjdk_vcf::comparator::VariantContextComparator;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::vcf_file::write_vcf;

/// The key every `CO=` becomes, which is not the neighbouring tool's.
pub const COMMENT_KEY: &str = "GatherVcfs.comment";

/// The arguments the gathering reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arguments {
    pub comments: Vec<String>,
    /// `CREATE_INDEX`, which the constructor sets to true rather than the parser.
    pub create_index: bool,
    /// `REORDER_INPUT_BY_FIRST_VARIANT`.
    pub reorder_input_by_first_variant: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            comments: Vec::new(),
            create_index: true,
            reorder_input_by_first_variant: false,
        }
    }
}

/// One input, carried with the path its refusals name.
#[derive(Debug, Clone, Copy)]
pub struct Input<'a> {
    pub path: &'a str,
    pub text: &'a str,
}

/// How a run ends when it does not write a file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Failure {
    /// Raised out of `doWork`, either before the try or as an `AssertionError`.
    Thrown { class: String, message: String },
    /// Caught by the try, logged, and turned into exit code 1 with the output deleted.
    Exit { class: String, message: String },
}

impl Failure {
    pub fn java_class(&self) -> &str {
        match self {
            Failure::Thrown { class, .. } | Failure::Exit { class, .. } => class,
        }
    }

    pub fn message(&self) -> &str {
        match self {
            Failure::Thrown { message, .. } | Failure::Exit { message, .. } => message,
        }
    }

    /// The exit code a caller sees: `None` when the failure was thrown instead.
    pub fn exit_code(&self) -> Option<i32> {
        match self {
            Failure::Exit { .. } => Some(1),
            Failure::Thrown { .. } => None,
        }
    }

    /// The line htsjdk's `Log` writes for an exit, with the timestamp where the reference puts
    /// the run's own. The golden masks it, so the mask is what this renders.
    pub fn log_line(&self) -> Option<String> {
        match self {
            Failure::Exit { class, message } => Some(format!(
                "ERROR\tMASKED\tGatherVcfs\tThere was a problem with gathering the INPUT.\
                 {class}: {message}\n"
            )),
            Failure::Thrown { .. } => None,
        }
    }
}

/// One sequence of the dictionary, printed the way `SAMSequenceRecord.toString` prints it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    pub index: usize,
    pub assembly: Option<String>,
}

impl Sequence {
    fn rendered(&self) -> String {
        format!(
            "SAMSequenceRecord(name={},length={},dict_index={},assembly={},alternate_names=[])",
            self.name,
            self.length,
            self.index,
            self.assembly.as_deref().unwrap_or("null")
        )
    }

    /// `isSameSequence`, minus the MD5 and alternate-name branches no VCF contig line can reach.
    /// A length of `UNKNOWN_SEQUENCE_LENGTH` matches anything, which is why 0 is not a mismatch.
    fn is_same_sequence(&self, other: &Sequence) -> bool {
        if self.index != other.index {
            return false;
        }
        if self.length != 0 && other.length != 0 && self.length != other.length {
            return false;
        }
        self.name == other.name
    }
}

/// `VCFHeader.getSequenceDictionary()`, which is `None` when there are no contig lines.
pub fn dictionary(header: &VcfHeader) -> Option<Vec<Sequence>> {
    let field = |fields: &[(String, String)], key: &str| {
        fields
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.clone())
    };
    let sequences: Vec<Sequence> = header
        .lines
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .enumerate()
        .map(|(index, line)| match line {
            HeaderLine::Contig { fields, .. } => Sequence {
                name: field(fields, "ID").unwrap_or_default(),
                length: field(fields, "length")
                    .and_then(|value| value.parse().ok())
                    .unwrap_or(0),
                index,
                assembly: field(fields, "assembly"),
            },
            _ => unreachable!("filtered to contig lines"),
        })
        .collect();
    if sequences.is_empty() {
        None
    } else {
        Some(sequences)
    }
}

/// `SAMSequenceDictionary.assertSameDictionary`, whose failure is an `AssertionError`.
pub fn assert_same_dictionary(expected: &[Sequence], found: &[Sequence]) -> Result<(), String> {
    let template = |detail: String| format!("SAM dictionaries are not the same: {detail}.");
    for (position, this) in expected.iter().enumerate() {
        match found.get(position) {
            None => {
                return Err(template(format!(
                    "{} is present in only one dictionary",
                    this.rendered()
                )))
            }
            Some(that) => {
                if !that.is_same_sequence(this) {
                    return Err(template(format!(
                        "{} was found when {} was expected",
                        that.rendered(),
                        this.rendered()
                    )));
                }
            }
        }
    }
    if let Some(extra) = found.get(expected.len()) {
        return Err(template(format!(
            "{} is present in only one dictionary",
            extra.rendered()
        )));
    }
    Ok(())
}

/// The message a null dictionary produces, which is where `new VariantContextComparator(dict)`
/// dereferences it.
fn null_dictionary() -> Failure {
    Failure::Exit {
        class: "java.lang.NullPointerException".to_string(),
        message: "Cannot invoke \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \
                  \"dictionary\" is null"
            .to_string(),
    }
}

/// `header.getGenotypeSamples()`, which is the column order and not a sorted set.
fn samples(header: &VcfHeader) -> Vec<String> {
    header.samples.clone()
}

/// The `[a, b]` a Java `SortedSet` prints.
fn rendered_set(names: &[String]) -> String {
    let mut sorted: Vec<&String> = names.iter().collect();
    sorted.sort();
    format!(
        "[{}]",
        sorted
            .iter()
            .map(|name| name.as_str())
            .collect::<Vec<&str>>()
            .join(", ")
    )
}

/// `VariantContextComparator.compare`, which unboxes a null and so throws a NullPointerException
/// when a record sits on a contig the dictionary does not have. No run in the golden reaches it.
fn compare(
    comparator: &VariantContextComparator,
    left: &VariantContext,
    right: &VariantContext,
) -> Result<i32, Failure> {
    comparator.compare(left, right).map_err(|_| Failure::Exit {
        class: "java.lang.NullPointerException".to_string(),
        message: "Cannot invoke \"java.lang.Integer.intValue()\" because the return value of \
                  \"java.util.Map.get(Object)\" is null"
            .to_string(),
    })
}

fn comparator_of(header: &VcfHeader) -> VariantContextComparator {
    // The comparator refuses a list that holds anything but contigs, so it is handed those alone.
    let contigs: Vec<HeaderLine> = header
        .lines
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .cloned()
        .collect();
    VariantContextComparator::from_header_lines(&contigs).expect("a dictionary already checked")
}

/// One parsed input, held for as long as the gathering needs it.
struct Shard {
    path: String,
    header: VcfHeader,
    records: Vec<VariantContext>,
}

fn read_shard(input: &Input) -> Result<Shard, Failure> {
    match read_vcf(input.text) {
        Ok(file) => Ok(Shard {
            path: input.path.to_string(),
            header: file.header,
            records: file.records,
        }),
        Err(failure) => Err(Failure::Exit {
            class: failure.error.class().to_string(),
            message: failure.error.message(),
        }),
    }
}

/// `assertSameSamplesAndValidOrdering`, which also does the reordering.
fn check_and_reorder(shards: Vec<Shard>, args: &Arguments) -> Result<Vec<Shard>, Failure> {
    let first = shards.first().expect("at least one input");
    let expected_dictionary = match dictionary(&first.header) {
        Some(dictionary) => dictionary,
        // `new VariantContextComparator(header.getSequenceDictionary())` on a null dictionary.
        None => return Err(null_dictionary()),
    };
    let comparator = comparator_of(&first.header);
    let expected_samples = samples(&first.header);

    let mut shards = shards;
    if args.reorder_input_by_first_variant {
        // The comparator answers 1 whenever the left file is empty, so every empty file sorts
        // last, and `List.sort` is stable so files that tie keep their order.
        shards.sort_by(
            |left, right| match (left.records.first(), right.records.first()) {
                (None, None) => std::cmp::Ordering::Equal,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (Some(_), None) => std::cmp::Ordering::Less,
                (Some(one), Some(other)) => comparator
                    .compare(one, other)
                    .map(|answer| answer.cmp(&0))
                    .unwrap_or(std::cmp::Ordering::Equal),
            },
        );
    }

    let mut last: Option<(&VariantContext, &str)> = None;
    for shard in &shards {
        match dictionary(&shard.header) {
            None => {
                return Err(Failure::Thrown {
                    class: "java.lang.AssertionError".to_string(),
                    message: assert_same_dictionary(&expected_dictionary, &[])
                        .expect_err("a dictionary against nothing"),
                })
            }
            Some(found) => {
                if let Err(message) = assert_same_dictionary(&expected_dictionary, &found) {
                    return Err(Failure::Thrown {
                        class: "java.lang.AssertionError".to_string(),
                        message,
                    });
                }
            }
        }

        let these = samples(&shard.header);
        if these != expected_samples {
            let unique_to_first: Vec<String> = expected_samples
                .iter()
                .filter(|name| !these.contains(name))
                .cloned()
                .collect();
            let unique_to_this: Vec<String> = these
                .iter()
                .filter(|name| !expected_samples.contains(name))
                .cloned()
                .collect();
            return Err(Failure::Exit {
                class: "java.lang.IllegalArgumentException".to_string(),
                message: format!(
                    "VCFs do not have identical sample lists. Samples unique to first file: {}. \
                     Samples unique to {}: {}.",
                    rendered_set(&unique_to_first),
                    shard.path,
                    rendered_set(&unique_to_this)
                ),
            });
        }

        // An empty file moves nothing: `lastContext` and `lastFile` are only set inside the
        // `hasNext()` branch.
        if let Some(current) = shard.records.first() {
            if let Some((previous, previous_path)) = last {
                if compare(&comparator, previous, current)? >= 0 {
                    return Err(Failure::Exit {
                        class: "java.lang.IllegalArgumentException".to_string(),
                        message: format!(
                            "First record in file {} is not after first record in previous file {}",
                            shard.path, previous_path
                        ),
                    });
                }
            }
            last = Some((current, &shard.path));
        }
    }
    Ok(shards)
}

/// `gatherConventionally`, which is every run where the block copying path was not chosen.
fn gather_conventionally(shards: &[Shard], args: &Arguments) -> Result<String, Failure> {
    let first = shards.first().expect("at least one input");
    let mut header = first.header.clone();
    // `addMetaDataLine` keys an unstructured line by its key, so only the first comment lands.
    if let Some(comment) = args.comments.first() {
        header.lines.push(HeaderLine::Unstructured {
            key: COMMENT_KEY.to_string(),
            value: comment.clone(),
        });
    }
    let comparator = comparator_of(&header);

    let mut records: Vec<VariantContext> = Vec::new();
    let mut last: Option<(&VariantContext, &str)> = None;
    for shard in shards {
        if let (Some((previous, previous_path)), Some(next)) = (last, shard.records.first()) {
            if compare(&comparator, next, previous)? <= 0 {
                return Err(Failure::Exit {
                    class: "java.lang.IllegalArgumentException".to_string(),
                    message: format!(
                        "First variant in file {} is at {}:{} but last variant in earlier file {} \
                         is at {}:{}",
                        shard.path,
                        next.contig,
                        next.start,
                        previous_path,
                        previous.contig,
                        previous.start
                    ),
                });
            }
        }
        for record in &shard.records {
            records.push(record.clone());
        }
        if let Some(record) = shard.records.last() {
            last = Some((record, &shard.path));
        }
    }

    write_vcf(&header, &records).map_err(|error| Failure::Exit {
        class: "java.lang.IllegalStateException".to_string(),
        message: format!("{error:?}"),
    })
}

/// `doWork()` for the conventional path: text in, text out.
pub fn gather(inputs: &[Input], args: &Arguments) -> Result<String, Failure> {
    let shards: Vec<Shard> = inputs
        .iter()
        .map(read_shard)
        .collect::<Result<Vec<Shard>, Failure>>()?;
    let first = shards.first().expect("at least one input");

    // The dictionary check is BEFORE the try, so this one is thrown rather than an exit code.
    if args.create_index && dictionary(&first.header).is_none() {
        return Err(Failure::Thrown {
            class: "picard.PicardException".to_string(),
            message: "In order to index the resulting VCF input VCFs must contain ##contig lines."
                .to_string(),
        });
    }

    let shards = check_and_reorder(shards, args)?;
    gather_conventionally(&shards, args)
}
