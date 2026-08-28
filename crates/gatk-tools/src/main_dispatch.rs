//! `Main`'s tool resolution and its refusals, which is what `gatk <Tool> <args>` does before
//! Barclay ever sees an argument.
//!
//! Finding the tool classes on the class path is not ported: there is no class path here. What is
//! ported is what happens once a name has failed to resolve, which is the whole of the refusal.
//!
//! Ported from `org.broadinstitute.hellbender.Main`,
//! `org.broadinstitute.hellbender.cmdline.DeprecatedToolsRegistry` and
//! `htsjdk.samtools.util.StringUtil`.

/// `Main.HELP_SIMILARITY_FLOOR`, which the distance must be STRICTLY under.
pub const HELP_SIMILARITY_FLOOR: i32 = 7;
/// `Main.MINIMUM_SUBSTRING_LENGTH`, from which a substring anywhere in a name scores zero.
pub const MINIMUM_SUBSTRING_LENGTH: usize = 5;
/// The weights `getSuggestedAlternateCommand` hands to the Levenshtein distance, in its order.
pub const SWAP_COST: i32 = 0;
pub const SUBSTITUTION_COST: i32 = 2;
pub const INSERTION_COST: i32 = 1;
pub const DELETION_COST: i32 = 4;

/// `DeprecatedToolsRegistry.deprecatedTools`: the tool, the version it went out in, and what to do
/// instead.
pub const DEPRECATED_TOOLS: &[(&str, &str, &str)] = &[
    ("IndelRealigner", "4.0.0.0", "Please use GATK3 to run this tool"),
    (
        "RealignerTargetCreator",
        "4.0.0.0",
        "Please use GATK3 to run this tool",
    ),
    (
        "CNNScoreVariants",
        "4.6.1.0",
        "Please use the replacement tool NVScoreVariants instead, which produces virtually identical results",
    ),
    (
        "CNNVariantTrain",
        "4.6.1.0",
        "Please use a version of GATK prior to 4.6.1.0 to run this tool, or wait for the forthcoming Pytorch-based training tool for NVScoreVariants to be released",
    ),
    (
        "CNNVariantWriteTensors",
        "4.6.1.0",
        "Please use a version of GATK prior to 4.6.1.0 to run this tool, or wait for the forthcoming Pytorch-based training tool for NVScoreVariants to be released",
    ),
];

/// `getSuggestedAlternateCommand`, on a name that DOES resolve.
pub const COMMAND_MATCHES_PREFIX: &str = "Command matches: ";

/// `DeprecatedToolsRegistry.getToolDeprecationInfo`, which answers nothing for a live tool.
pub fn tool_deprecation_info(tool_name: &str) -> Option<String> {
    DEPRECATED_TOOLS
        .iter()
        .find(|(name, _, _)| *name == tool_name)
        .map(|(name, version, message)| {
            format!("{name} is no longer included in GATK as of version {version}. {message}")
        })
}

/// `StringUtil.levenshteinDistance`, whose four costs are separate and whose swap is a
/// transposition of two adjacent characters.
///
/// The rows are indexed over the SECOND string, so an insertion walks along it and a deletion
/// walks along the first. With the dispatcher's weights that means dropping a character from the
/// COMMAND costs four while adding one costs one.
///
/// The first row is seeded for indices 0 to `second.len() - 1` and the LAST cell is left at zero,
/// which the reference does too. An empty first string therefore answers zero whatever the second
/// one is, rather than the second's length in insertions.
pub fn levenshtein_distance(
    first: &str,
    second: &str,
    swap: i32,
    substitution: i32,
    insertion: i32,
    deletion: i32,
) -> i32 {
    let a = first.as_bytes();
    let b = second.as_bytes();
    let mut row0 = vec![0i32; b.len() + 1];
    let mut row1 = vec![0i32; b.len() + 1];
    let mut row2 = vec![0i32; b.len() + 1];
    for (j, cell) in row1.iter_mut().enumerate().take(b.len()) {
        *cell = j as i32 * insertion;
    }
    for i in 0..a.len() {
        row2[0] = (i as i32 + 1) * deletion;
        for j in 0..b.len() {
            row2[j + 1] = row1[j];
            if a[i] != b[j] {
                row2[j + 1] += substitution;
            }
            if i > 0
                && j > 0
                && a[i - 1] == b[j]
                && a[i] == b[j - 1]
                && row2[j + 1] > row0[j - 1] + swap
            {
                row2[j + 1] = row0[j - 1] + swap;
            }
            if row2[j + 1] > row1[j + 1] + deletion {
                row2[j + 1] = row1[j + 1] + deletion;
            }
            if row2[j + 1] > row2[j] + insertion {
                row2[j + 1] = row2[j] + insertion;
            }
        }
        std::mem::swap(&mut row0, &mut row1);
        std::mem::swap(&mut row1, &mut row2);
    }
    row1[b.len()]
}

/// One tool's score against a command, which is zero for a prefix or a long enough substring and
/// the weighted distance otherwise.
///
/// A name the command matches EXACTLY is not scored at all: the reference throws there, the search
/// being reached only once resolution has failed.
pub fn distance(command: &str, name: &str) -> Result<i32, String> {
    if name == command {
        return Err(format!("{COMMAND_MATCHES_PREFIX}{command}"));
    }
    if name.starts_with(command)
        || (MINIMUM_SUBSTRING_LENGTH <= command.len() && name.contains(command))
    {
        return Ok(0);
    }
    Ok(levenshtein_distance(
        command,
        name,
        SWAP_COST,
        SUBSTITUTION_COST,
        INSERTION_COST,
        DELETION_COST,
    ))
}

/// `getSuggestedAlternateCommand`, whose message always opens on the same line.
///
/// The suggestions are appended with eight spaces before each and NO separator after, so two of
/// them run together on one line. When every tool scores zero the best distance is bumped over
/// the floor instead of every tool being listed, which is what the empty command and a one-tool
/// catalogue the command prefixes both reach.
pub fn suggested_alternate_command(classes: &[String], command: &str) -> Result<String, String> {
    let mut distances = Vec::with_capacity(classes.len());
    let mut best_distance = i32::MAX;
    let mut best_count = 0usize;
    for name in classes {
        let d = distance(command, name)?;
        distances.push(d);
        if d < best_distance {
            best_distance = d;
            best_count = 1;
        } else if d == best_distance {
            best_count += 1;
        }
    }
    if best_distance == 0 && best_count == classes.len() {
        best_distance = HELP_SIMILARITY_FLOOR + 1;
    }
    let mut message = format!("'{command}' is not a valid command.\n");
    if best_distance < HELP_SIMILARITY_FLOOR {
        message.push_str(if best_count < 2 {
            "Did you mean this?\n"
        } else {
            "Did you mean one of these?\n"
        });
        for (name, d) in classes.iter().zip(&distances) {
            if *d == best_distance {
                message.push_str(&format!("        {name}"));
            }
        }
    }
    Ok(message)
}

/// `getUnknownCommandMessage`: the deprecation notice short-circuits the search entirely.
pub fn unknown_command_message(classes: &[String], command: &str) -> Result<String, String> {
    match tool_deprecation_info(command) {
        Some(message) => Ok(message),
        None => suggested_alternate_command(classes, command),
    }
}
