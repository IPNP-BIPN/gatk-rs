//! `GATKConfig`: the system properties `Main` installs before any tool runs.
//!
//! This is not configuration in the sense of something a user tunes. It is a table of defaults
//! that reaches **htsjdk**, and one of its rows decides the bytes of every block-compressed file
//! GATK writes:
//!
//! ```text
//! @SystemProperty
//! @Key("samjdk.compression_level")
//! @DefaultValue("2")
//! ```
//!
//! htsjdk's own default is **five**. Two is the one level pair Intel's GKL routes through ISA-L
//! igzip rather than through zlib, so the same eighty-one bytes deflate to 49 under a real `gatk`
//! invocation and to 43 under htsjdk's default. A covering-array row over `IndexFeatureFile` read
//! as a divergence for exactly that reason while its index body was already byte for byte (#1032).
//!
//! # Every key is a system property, and that is measured rather than assumed
//!
//! All twelve keys the interface declares carry `@SystemProperty`, so the table below is both the
//! configuration and the set of properties. The `gatk-config` suite dumps the annotation and the
//! effect separately, because a key added without the annotation would be read by GATK and never
//! seen by htsjdk, and nothing in the values would say so.
//!
//! # A property already set is not overwritten
//!
//! `injectSystemPropertiesFromConfig` leaves an existing value alone, so
//! `-Dsamjdk.compression_level=5` on the command line wins over the default. [`resolve`] is that
//! rule.
//!
//! Ported from `org.broadinstitute.hellbender.utils.config.GATKConfig` and
//! `org.broadinstitute.hellbender.utils.config.ConfigFactory.injectSystemPropertiesFromConfig`.
//!
//! GATK 4.7.0.0. `ConfigFactory.java` is unchanged since 4.6.2.0; `GATKConfig.java` gains one
//! method, `outputCramVersion()` with `@DefaultValue("3.1")`, which carries neither `@Key` nor
//! `@SystemProperty` and so adds no row here.
//!
//! # "Key" means a method with `@Key`
//!
//! The interface declares more settings than the twelve here (`codec_packages`,
//! `cloudPrefetchBuffer`, `createOutputBamIndex` and others). Those carry no `@Key` and no
//! `@SystemProperty`, are read only inside GATK, and never reach htsjdk. The golden's dumper
//! (`tools/readfilter-conformance/GatkConfigDump.java`) keeps only methods with `@Key`, and
//! every one of those is a system property.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing", Java's nullable reference with the null check
//!   enforced by the compiler: `Some(v)` holds a value, `None` is the absent case;
//! - `&str` is a borrowed, read-only string; `&'static str` is one baked into the binary;
//! - `.iter().find(..).map(..)` is Java's `stream().filter(..).findFirst().map(..)`;
//! - a function's last expression is its return value.

/// `samjdk.compression_level`, whose default is not htsjdk's.
pub const COMPRESSION_LEVEL: &str = "samjdk.compression_level";

/// The `@Key` methods `GATKConfig` declares, as (`@Key` value, `@DefaultValue`), sorted by key.
///
/// The key is the `@Key` string, not the method name: `samjdk.compression_level` is the
/// property, `samjdk_compression_level()` the method.
///
/// Sorted because that is the order the golden carries: `getDeclaredMethods` has no defined order,
/// and a table that depended on a JVM's reflection order would be a table about the JVM.
pub const DEFAULTS: &[(&str, &str)] = &[
    ("gatk_stacktrace_on_user_exception", "false"),
    ("samjdk.compression_level", "2"),
    ("samjdk.use_async_io_read_samtools", "false"),
    ("samjdk.use_async_io_write_samtools", "true"),
    ("samjdk.use_async_io_write_tribble", "false"),
    ("spark.driver.extraJavaOptions", ""),
    ("spark.driver.maxResultSize", "0"),
    ("spark.driver.userClassPathFirst", "true"),
    ("spark.executor.extraJavaOptions", ""),
    ("spark.executor.memoryOverhead", "600"),
    ("spark.io.compression.codec", "lzf"),
    ("spark.kryoserializer.buffer.max", "512m"),
];

/// The default for one key, or `None` for a key the config does not declare.
pub fn default(key: &str) -> Option<&'static str> {
    DEFAULTS
        .iter()
        .find(|(name, _)| *name == key)
        .map(|(_, value)| *value)
}

/// `injectSystemPropertiesFromConfig`, for one key: what is already set wins.
///
/// `set` is what the process already holds for this key, which is what a `-D` on the command line
/// put there.
///
/// The `'a` and `where 'static: 'a` only say the result may borrow either from `set` or from the
/// table compiled into the binary; they change nothing at run time.
///
/// ```java
/// if ( systemProperties.containsKey(entry.getKey()) ) {
///     logger.debug("System property already exists.  Not overriding: " + entry.getKey());
///     continue;
/// }
/// ```
pub fn resolve<'a>(key: &str, set: Option<&'a str>) -> Option<&'a str>
where
    'static: 'a,
{
    match set {
        Some(value) => Some(value),
        None => default(key),
    }
}

/// The BGZF compression level a tool run writes at, which is [`COMPRESSION_LEVEL`] as a number.
///
/// A value that is not a number is not this port's problem to interpret: htsjdk parses it with
/// `Integer.parseInt` and throws, and nothing here can throw, so the config's own default stands
/// in and the caller is none the wiser. That is a boundary and it is stated: no golden covers a
/// malformed property, because `GATKConfig` cannot produce one.
pub fn compression_level(set: Option<&str>) -> u32 {
    // `and_then` keeps going only on a value; `.parse().ok()` turns a parse failure into `None`;
    // `unwrap_or(2)` catches only the unparsable case, since an unset one already resolved to
    // the table's "2".
    resolve(COMPRESSION_LEVEL, set)
        .and_then(|value| value.parse().ok())
        .unwrap_or(2)
}
