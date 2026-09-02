//! Conformance for `GATKConfig`'s system properties against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/GatkConfigDump.java`. This is not configuration a
//! user tunes: it is a table of defaults that reaches htsjdk, and one row of it decides the bytes
//! of every block-compressed file GATK writes.
//!
//! # What this suite is for
//!
//!  * **`samjdk.compression_level` defaulting to TWO**, where htsjdk's own default is five, and
//!    two is the one level pair Intel's GKL routes through igzip rather than zlib;
//!  * **every key being a `@SystemProperty`**, which is measured rather than assumed: a key added
//!    without the annotation would be read by GATK and never reach htsjdk;
//!  * **the property name being `@Key` and not the method name**;
//!  * **and injection leaving an already-set property alone**, so a `-D` on the command line wins.
//!
//! While the suite is `golden-pending` the dump is named by `GATK_CONFIG_DUMP`.

use gatk_tools::gatk_config::{compression_level, default, resolve, COMPRESSION_LEVEL, DEFAULTS};

fn rows<'a>(dump: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    dump.lines()
        .filter_map(|line| line.strip_prefix(&format!("{kind}\t")))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

#[test]
fn the_table_is_the_reference_one() {
    let dump = match std::env::var("GATK_CONFIG_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by GATK_CONFIG_DUMP"),
        Err(_) => {
            println!(
                "skipped: the gatk-config golden is still pending. Run the suite and point \
                 GATK_CONFIG_DUMP at tools/conformance/pending/gatk-config.GatkConfigDump.txt"
            );
            return;
        }
    };

    // The declaration: key, default, and whether it is injected at all.
    let declared = rows(&dump, "key");
    assert_eq!(
        declared.len(),
        DEFAULTS.len(),
        "the table has a row this port does not"
    );
    for (row, (key, value)) in declared.iter().zip(DEFAULTS) {
        assert_eq!(row[0], *key);
        assert_eq!(row[1], *value, "{key}");
        // Every key GATKConfig declares carries `@SystemProperty`. The port's table would be a
        // table of things htsjdk never sees if that stopped being true, so it is asserted.
        assert_eq!(row[2], "system", "{key} is no longer a system property");
    }

    // The effect: what `System.getProperties` holds once `Main` has injected them.
    for row in rows(&dump, "injected") {
        assert_eq!(
            resolve(row[0], None),
            Some(if row[1] == "<unset>" { "" } else { row[1] }),
            "{}",
            row[0]
        );
    }

    // And the precedence: an already-set property is not overwritten, so a `-D` wins.
    let precedence = rows(&dump, "precedence");
    assert_eq!(precedence.len(), 1);
    assert_eq!(precedence[0][0], COMPRESSION_LEVEL);
    assert_eq!(
        resolve(COMPRESSION_LEVEL, Some(precedence[0][1])),
        Some(precedence[0][1])
    );
}

/// The row this table exists for.
#[test]
fn the_compression_level_is_two_and_not_htsjdks_five() {
    assert_eq!(default(COMPRESSION_LEVEL), Some("2"));
    assert_eq!(compression_level(None), 2);
    // A `-D` on the command line wins, which is the only way a run writes at htsjdk's own level.
    assert_eq!(compression_level(Some("5")), 5);
}
