//! `SetNmAndUqTags`: Picard's `@Deprecated` subclass of `SetNmMdAndUqTags`, which overrides
//! nothing.
//!
//! ```java
//! public class SetNmAndUqTags extends SetNmMdAndUqTags {
//! ```
//!
//! `pub(super) use` re-exports the parent's `run` under this module's name, so the dispatch table
//! can name either tool and reach the same function.

pub(super) use super::set_nm_md_and_uq_tags::run;
