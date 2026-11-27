//! `HtsgetReader`: what the tool decides before and after it talks to the server.
//!
//! The classes this is ported from are `HtsgetRequestBuilder`, `HtsgetResponse` and
//! `HtsgetErrorResponse`, as Jackson reads the last two. The transport is the runner's.
//!
//! Four things in them are not what they look like:
//!
//! * **a `.bam` id refuses every `--format`, `BAM` included.** The check's parentheses bind the
//!   format's test to `.sam` alone;
//! * **the scheme is always `https`**, whatever `--url` says, and a `--url` with no scheme is a
//!   path, so it reaches no host;
//! * **every failure to read the JSON is an `IOException`** and reports the request, not the JSON:
//!   an unknown property, a root other than `htsget`, a constant out of its enum;
//! * **a block's status is never looked at**: its body is written whatever the server answered.
//!
//! The reference is `org.broadinstitute.hellbender.tools.HtsgetReader` and the
//! `org.broadinstitute.hellbender.tools.htsgetreader` classes in GATK 4.7.0.0. The JSON is read
//! by a mapper built as
//!
//! ```java
//! .enable(DeserializationFeature.UNWRAP_ROOT_VALUE)
//! .configure(MapperFeature.ACCEPT_CASE_INSENSITIVE_PROPERTIES, true)
//! ```
//!
//! which is why every property name below is matched lowercased, and why the object must sit
//! under a single `htsget` key. Jackson, Jersey and the JDK's `java.net.URI` and `Base64` have no
//! source in this tree: what this file reproduces of them is their observable behaviour, as the
//! goldens record it.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" and `Result<T, E>` "a value or an error"; `?` after an
//!   expression hands the nothing or the error straight back to the caller;
//! - `&str` is a borrowed, read-only view of a string; `text[a..b]` is the substring between two
//!   BYTE offsets, and `find` returns a byte offset;
//! - `if let Some(x) = y { .. }` runs the block only when `y` holds a value, and `let Some(x) = y
//!   else { .. }` takes the value or the early exit;
//! - `match` must cover every case (the compiler checks), `_` is "anything else", and `x if
//!   cond =>` is an arm that applies only when `cond` holds;
//! - `|x| ...` is a lambda; `.iter().any(..)`, `.filter(..)`, `.map(..)`, `.collect()` are
//!   Java's stream operations.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_hash::hash_set_order;
use serde_json::Value;

/// `HtsgetRequestField`'s constants, in the order an `EnumSet` walks them, which is declaration
/// order: `--fields` given in any order is sent in this one.
pub const FIELDS: [&str; 11] = [
    "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL",
];

/// `SimpleInterval`, as far as a request asks of one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RequestInterval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// `HtsgetRequestBuilder`.
///
/// `endpoint` is `--url` as typed and `id` is `--id`; `format` and `data_class` are the enum
/// constants' names. `tags` and `notags` are kept as given; [`Request::validate`] and
/// [`Request::to_uri`] put them in the `HashSet` order the Java's sets have.
#[derive(Debug, Clone, Default)]
pub struct Request {
    pub endpoint: String,
    pub id: String,
    pub format: Option<String>,
    pub data_class: Option<String>,
    pub interval: Option<RequestInterval>,
    pub fields: Vec<String>,
    pub tags: Vec<String>,
    pub notags: Vec<String>,
}

/// A `HashSet<String>` of these strings, in its own order.
///
/// `unwrap_or_else` falls back to the order given when `hash_set_order` cannot compute one.
fn java_set(values: &[String]) -> Vec<String> {
    hash_set_order(values).unwrap_or_else(|_| values.to_vec())
}

impl Request {
    /// `validateRequest`: the message of the `UserException`, when there is one.
    ///
    /// ```java
    /// if ((this.id.endsWith(FileExtensions.BAM) || this.id.endsWith(FileExtensions.SAM) && (
    ///     this.format != HtsgetFormat.BAM && this.format != HtsgetFormat.CRAM))
    /// ```
    ///
    /// `&&` binds tighter than `||`, so the format test applies to `.sam` alone and a `.bam` id
    /// fails whatever the format: the module header's first point.
    pub fn validate(&self) -> Result<(), String> {
        if self.data_class.as_deref() == Some("header")
            && (self.interval.is_some()
                || !self.fields.is_empty()
                || !self.tags.is_empty()
                || !self.notags.is_empty())
        {
            return Err(
                "Invalid request: no query parameters except `format` may be specified when class=header"
                    .to_string(),
            );
        }
        if let Some(format) = &self.format {
            let name = format.as_str();
            let reads = self.id.ends_with(".bam")
                || (self.id.ends_with(".sam") && name != "BAM" && name != "CRAM");
            let variants = [".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"]
                .iter()
                .any(|extension| self.id.ends_with(extension))
                && name != "VCF"
                && name != "BCF";
            if reads || variants {
                return Err(format!(
                    "Specified format: {format} is incompatible with id's file extension"
                ));
            }
        }
        // `this.tags.stream().filter(getNoTags()::contains).collect(Collectors.joining(", "))`:
        // the overlap listed in the tag set's iteration order.
        let notags = java_set(&self.notags);
        let overlap: Vec<String> = java_set(&self.tags)
            .into_iter()
            .filter(|tag| notags.contains(tag))
            .collect();
        if !overlap.is_empty() {
            return Err(format!(
                "Invalid request: tags and notags overlap in the following fields: {}",
                overlap.join(", ")
            ));
        }
        Ok(())
    }

    /// `toURI`, past the validation: Jersey's `UriBuilder` over the endpoint. What it throws is an
    /// exception class and its message, since neither is a `UserException`.
    ///
    /// ```java
    /// final UriBuilder builder = UriBuilder.fromUri(this.endpoint)
    ///     .scheme("https")
    ///     .path(this.id);
    /// ```
    ///
    /// then one `queryParam` per argument given, in the order `format`, `class`, `referenceName`,
    /// `start`, `end`, `fields`, `tags`, `notags`, appended after any query the endpoint had.
    /// Whether the call to `validateRequest` that opens `toURI` has run is the caller's business.
    pub fn to_uri(&self) -> Result<String, UriFailure> {
        let endpoint = split_uri(&self.endpoint);
        if uri_scheme(&self.endpoint).is_some()
            && endpoint.authority.is_none()
            && !endpoint.path.starts_with('/')
        {
            return Err(UriFailure {
                class: "java.lang.IllegalArgumentException",
                message: "Schema specific part is opaque.".to_string(),
            });
        }
        check_braces(&self.id)?;
        let mut uri = String::from("https:");
        if let Some(authority) = endpoint.authority {
            uri.push_str("//");
            uri.push_str(authority);
        }
        let mut path = endpoint.path.to_string();
        if !self.id.is_empty() {
            // With no host and no path nothing separates the scheme from the id.
            if !path.ends_with('/')
                && !self.id.starts_with('/')
                && (endpoint.authority.is_some() || !path.is_empty())
            {
                path.push('/');
            }
            path.push_str(&encode(&self.id, Component::Path));
        }
        uri.push_str(&path);

        let mut query = endpoint.query.map(str::to_string);
        for value in self.tags.iter().chain(&self.notags) {
            check_braces(value)?;
        }
        // A closure that appends `name=value` to `query`. It is `mut` because it writes to the
        // `query` it captures; `query.take()` moves the current value out, leaving `None`.
        let mut parameter = |name: &str, value: &str| {
            let pair = format!("{name}={}", encode(value, Component::QueryParameter));
            query = Some(match query.take() {
                Some(earlier) if !earlier.is_empty() => format!("{earlier}&{pair}"),
                _ => pair,
            });
        };
        if let Some(format) = &self.format {
            parameter("format", format);
        }
        if let Some(data_class) = &self.data_class {
            parameter("class", data_class);
        }
        if let Some(interval) = &self.interval {
            parameter("referenceName", &interval.contig);
            // `getGA4GHStart` is zero based, `getGA4GHEnd` is the end as it stands.
            parameter("start", &(interval.start - 1).to_string());
            parameter("end", &interval.end.to_string());
        }
        if !self.fields.is_empty() {
            let chosen: Vec<&str> = FIELDS
                .iter()
                .copied()
                .filter(|field| self.fields.iter().any(|name| name == field))
                .collect();
            parameter("fields", &chosen.join(","));
        }
        if !self.tags.is_empty() {
            parameter("tags", &java_set(&self.tags).join(","));
        }
        if !self.notags.is_empty() {
            parameter("notags", &java_set(&self.notags).join(","));
        }
        if let Some(query) = query {
            uri.push('?');
            uri.push_str(&query);
        }
        if let Some(fragment) = endpoint.fragment {
            uri.push('#');
            uri.push_str(fragment);
        }
        unresolved_template(&uri)?;
        Ok(uri)
    }
}

/// What `UriBuilder` threw: the exception's class and message, since neither is a
/// `UserException`.
#[derive(Debug, PartialEq, Eq)]
pub struct UriFailure {
    pub class: &'static str,
    pub message: String,
}

/// A `{` with no `}` after it, which the template parser refuses as soon as it is handed the text.
fn check_braces(text: &str) -> Result<(), UriFailure> {
    let mut rest = text;
    while let Some(open) = rest.find('{') {
        match rest[open..].find('}') {
            Some(close) => rest = &rest[open + close + 1..],
            None => {
                return Err(UriFailure {
                    class: "java.lang.IllegalArgumentException",
                    message: format!(
                        "Invalid syntax in the template \"{text}\". Check if a path parameter is terminated with a \"}}\"."
                    ),
                })
            }
        }
    }
    Ok(())
}

/// `build()` with no values: a template names a variable nothing gives, and a lone `}` reaches
/// `new URI`, which says where it is.
fn unresolved_template(uri: &str) -> Result<(), UriFailure> {
    if let Some(open) = uri.find('{') {
        let close = uri[open..].find('}').unwrap_or(uri.len() - open - 1);
        let name = &uri[open + 1..open + close];
        return Err(UriFailure {
            class: "java.lang.IllegalArgumentException",
            message: format!("The template variable '{name}' has no value"),
        });
    }
    if let Some(at) = uri.find('}') {
        let part = if uri[..at].contains('?') {
            "query"
        } else {
            "path"
        };
        return Err(UriFailure {
            class: "javax.ws.rs.core.UriBuilderException",
            message: format!(
                "java.net.URISyntaxException: Illegal character in {part} at index {at}: {uri}"
            ),
        });
    }
    Ok(())
}

/// The five parts RFC 3986's appendix B splits a reference into, the scheme left out.
///
/// Each part borrows from the text it was split from (the `'a`), so nothing is copied.
struct UriParts<'a> {
    authority: Option<&'a str>,
    path: &'a str,
    query: Option<&'a str>,
    fragment: Option<&'a str>,
}

/// Split off the fragment at the first `#`, then the query at the first `?`, then the scheme,
/// then an authority introduced by `//` and running to the next `/`; what is left is the path.
fn split_uri(text: &str) -> UriParts<'_> {
    let (text, fragment) = match text.split_once('#') {
        Some((before, after)) => (before, Some(after)),
        None => (text, None),
    };
    let (text, query) = match text.split_once('?') {
        Some((before, after)) => (before, Some(after)),
        None => (text, None),
    };
    let text = match uri_scheme(text) {
        Some(scheme) => &text[scheme.len() + 1..],
        None => text,
    };
    let (authority, path) = match text.strip_prefix("//") {
        Some(rest) => {
            let end = rest.find('/').unwrap_or(rest.len());
            (Some(&rest[..end]), &rest[end..])
        }
        None => (None, text),
    };
    UriParts {
        authority,
        path,
        query,
        fragment,
    }
}

/// `URI.getScheme`.
///
/// The text before the first `:` when no `/`, `?` or `#` comes first, and only if it is a
/// letter followed by letters, digits, `+`, `-` or `.` (RFC 3986's scheme grammar).
/// `.then_some(scheme)` is `cond ? Optional.of(scheme) : Optional.empty()`.
pub fn uri_scheme(text: &str) -> Option<&str> {
    let end = text.find([':', '/', '?', '#'])?;
    let scheme = &text[..end];
    (text.as_bytes()[end] == b':'
        && scheme
            .chars()
            .next()
            .is_some_and(|first| first.is_ascii_alphabetic())
        && scheme
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '+' | '-' | '.')))
    .then_some(scheme)
}

/// `new SimpleInterval(String)`, for `--intervals`.
///
/// Split at the LAST colon: `contig` alone is the whole contig (`1` to `Integer.MAX_VALUE`),
/// `contig:start-end` a range, `contig:start+` from `start` to the end, and `contig:start` the
/// one base. The `(String, i32, i32)` is a tuple: contig, start, end.
pub fn interval(text: &str) -> Result<(String, i32, i32), String> {
    let position = |digits: &str| {
        SimpleInterval::parse_position(digits).ok_or_else(|| {
            format!("Problem parsing start/end value in interval string. Value was: {digits}")
        })
    };
    let Some((contig, range)) = text.rsplit_once(':') else {
        return Ok((text.to_string(), 1, i32::MAX));
    };
    let (start, end) = match range.split_once('-') {
        Some((start, end)) => (position(start)?, position(end)?),
        None => match range.strip_suffix('+') {
            Some(start) => (position(start)?, i32::MAX),
            None => {
                let start = position(range)?;
                (start, start)
            }
        },
    };
    Ok((contig.to_string(), start, end))
}

/// Whether `new SimpleInterval(text)` succeeds: not empty, and its positions in order.
pub fn java_interval_parses(text: &str) -> bool {
    !text.is_empty()
        && interval(text).is_ok_and(|(_, start, end)| SimpleInterval::is_valid(start, end))
}

/// Whether `new URI(text)` parses: only the characters RFC 2396 allows, and each escape whole.
///
/// Also refused: an empty authority with an empty path (`scheme://`), a colon that ends a
/// would-be scheme which is not one, a scheme with nothing after its colon, and a second `#`.
/// Non-ASCII characters that are neither control nor space characters pass.
pub fn java_uri_parses(text: &str) -> bool {
    if split_uri(text).authority == Some("") && split_uri(text).path.is_empty() {
        return false;
    }
    if let Some(colon) = text.find(':') {
        let before = &text[..colon];
        if !before.contains(['/', '?', '#']) && uri_scheme(text).is_none() {
            return false;
        }
        if colon + 1 == text.len() && !before.contains(['/', '?', '#']) {
            return false;
        }
    }
    let mut fragments = 0;
    let mut characters = text.chars().peekable();
    while let Some(character) = characters.next() {
        match character {
            '%' => {
                for _ in 0..2 {
                    if !characters.next().is_some_and(|c| c.is_ascii_hexdigit()) {
                        return false;
                    }
                }
            }
            '#' => {
                fragments += 1;
                if fragments > 1 {
                    return false;
                }
            }
            c if c.is_ascii_alphanumeric() => {}
            '-' | '_' | '.' | '!' | '~' | '*' | '\'' | '(' | ')' | ';' | '/' | '?' | ':' | '@'
            | '&' | '=' | '+' | '$' | ',' | '[' | ']' => {}
            c if !c.is_ascii() && !c.is_control() && !c.is_whitespace() => {}
            _ => return false,
        }
    }
    true
}

/// Which part of the URI is being encoded: a path segment allows more characters unescaped than a
/// query parameter does, and only a query parameter turns a space into `+`.
#[derive(Clone, Copy, PartialEq)]
enum Component {
    Path,
    QueryParameter,
}

/// `UriComponent.contextualEncode`: a whole escape stays as it is, and the rest as UTF-8 escapes.
fn encode(text: &str, component: Component) -> String {
    // The braces stay because they are a template's, which is the builder's next complaint.
    let allowed = |c: char| {
        c.is_ascii_alphanumeric()
            || matches!(c, '-' | '_' | '.' | '~' | '{' | '}')
            || (component == Component::Path
                && matches!(
                    c,
                    '!' | '$'
                        | '\''
                        | '('
                        | ')'
                        | '*'
                        | ','
                        | ';'
                        | ':'
                        | '@'
                        | '/'
                        | '&'
                        | '='
                        | '+'
                ))
    };
    // `char_indices` yields each character with its byte offset, so a `%` can look at the two
    // bytes after it. A character outside the allowed set is written as the `%XX` of each of its
    // UTF-8 bytes, upper-case hex.
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    for (at, character) in text.char_indices() {
        if character == '%'
            && bytes.get(at + 1).is_some_and(u8::is_ascii_hexdigit)
            && bytes.get(at + 2).is_some_and(u8::is_ascii_hexdigit)
        {
            out.push('%');
        } else if character == ' ' && component == Component::QueryParameter {
            out.push('+');
        } else if allowed(character) {
            out.push(character);
        } else {
            let mut buffer = [0u8; 4];
            for byte in character.encode_utf8(&mut buffer).bytes() {
                out.push_str(&format!("%{byte:02X}"));
            }
        }
    }
    out
}

/// One of `urls`: `HtsgetResponse.Block`.
///
/// `url` is `None` where the JSON held null. `headers` keeps the JSON object's order, and a null
/// header value becomes the empty string. The block's `class` is checked (see `block`) but not
/// kept: nothing reads it.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Block {
    pub url: Option<String>,
    pub headers: Option<Vec<(String, String)>>,
}

/// A body that is the JSON literal `null`, which Jackson reads as no object at all, whatever
/// follows it.
///
/// `null` must be a whole token: `nullx` is not the literal.
pub fn is_null_literal(body: &str) -> bool {
    body.trim_start()
        .strip_prefix("null")
        .is_some_and(|rest| !rest.starts_with(|c: char| c.is_alphanumeric() || c == '_'))
}

/// `HtsgetResponse`. `None` for `blocks` is a response that named no `urls`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Response {
    pub blocks: Option<Vec<Block>>,
    pub md5: Option<String>,
}

/// `HtsgetErrorResponse`: its `@JsonProperty("error")` and `@JsonProperty("message")`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ErrorResponse {
    pub error: Option<String>,
    pub message: Option<String>,
}

/// `UNWRAP_ROOT_VALUE`: an object whose first and only key is `htsget`, then nothing it reads.
///
/// `@JsonRootName(value = "htsget")` on both response classes names the key. Only the first
/// JSON value of the body is parsed (`into_iter().next()`), so trailing text after it is not
/// read, and every failure is the one [`Unreadable`].
fn root(body: &str) -> Result<serde_json::Map<String, Value>, Unreadable> {
    let value = serde_json::Deserializer::from_str(body)
        .into_iter::<Value>()
        .next()
        .ok_or(Unreadable)?
        .map_err(|_| Unreadable)?;
    let Value::Object(wrapper) = value else {
        return Err(Unreadable);
    };
    if wrapper.len() != 1 {
        return Err(Unreadable);
    }
    match wrapper.into_iter().next() {
        Some((name, Value::Object(inner))) if name == "htsget" => Ok(inner),
        _ => Err(Unreadable),
    }
}

/// `StringDeserializer`: a scalar is its text, and a structure is refused.
///
/// `true` reads as `"true"` and a number as serde_json renders it (`"md5":1` reads as `"1"`,
/// see the tests).
fn text_of(value: &Value) -> Result<Option<String>, Unreadable> {
    match value {
        Value::Null => Ok(None),
        Value::String(text) => Ok(Some(text.clone())),
        Value::Bool(flag) => Ok(Some(flag.to_string())),
        Value::Number(number) => Ok(Some(number.to_string())),
        _ => Err(Unreadable),
    }
}

/// An enum constant, by its exact spelling or null.
///
/// Case-insensitivity covers property NAMES only: `"format":"bam"` is refused (see the tests).
fn constant(value: &Value, allowed: &[&str]) -> Result<(), Unreadable> {
    match value {
        Value::Null => Ok(()),
        Value::String(text) if allowed.contains(&text.as_str()) => Ok(()),
        _ => Err(Unreadable),
    }
}

/// One `HtsgetResponse.Block`: `url` parsed as a `java.net.URI` (a text `new URI` refuses fails
/// the whole body), `headers` as a string map, `class` as an `HtsgetClass`, and any other key
/// refused, since the mapper does not ignore unknown properties.
fn block(value: &Value) -> Result<Block, Unreadable> {
    let Value::Object(fields) = value else {
        return Err(Unreadable);
    };
    let mut block = Block::default();
    for (name, value) in fields {
        match name.to_lowercase().as_str() {
            "url" => {
                block.url = match text_of(value)? {
                    Some(text) if !java_uri_parses(&text) => return Err(Unreadable),
                    text => text,
                }
            }
            "headers" => {
                block.headers = match value {
                    Value::Null => None,
                    Value::Object(headers) => Some(
                        headers
                            .iter()
                            .map(|(name, value)| {
                                Ok((name.clone(), text_of(value)?.unwrap_or_default()))
                            })
                            .collect::<Result<_, Unreadable>>()?,
                    ),
                    _ => return Err(Unreadable),
                }
            }
            "class" => constant(value, &["body", "header"])?,
            _ => return Err(Unreadable),
        }
    }
    Ok(block)
}

/// Jackson's refusal of a body, whose text is never read.
#[derive(Debug)]
pub struct Unreadable;

/// `mapper.readValue(jsonBody, HtsgetResponse.class)`. The refusal's text is never read.
///
/// The three properties are `@JsonProperty("format")` (an `HtsgetFormat`), `"urls"` (a list of
/// blocks) and `"md5"`; any other name is refused. A null `urls` leaves `blocks` as `None`, which
/// is the Java's null list.
pub fn read_response(body: &str) -> Result<Response, Unreadable> {
    let mut response = Response::default();
    for (name, value) in root(body)? {
        match name.to_lowercase().as_str() {
            "format" => constant(&value, &["BAM", "CRAM", "VCF", "BCF"])?,
            "urls" => {
                response.blocks = match &value {
                    Value::Null => None,
                    Value::Array(blocks) => Some(
                        blocks
                            .iter()
                            .map(block)
                            .collect::<Result<_, Unreadable>>()?,
                    ),
                    _ => return Err(Unreadable),
                }
            }
            "md5" => response.md5 = text_of(&value)?,
            _ => return Err(Unreadable),
        }
    }
    Ok(response)
}

/// `mapper.readValue(jsonBody, HtsgetErrorResponse.class)`: `error` and `message` as text, any
/// other property refused.
pub fn read_error_response(body: &str) -> Result<ErrorResponse, Unreadable> {
    let mut response = ErrorResponse::default();
    for (name, value) in root(body)? {
        match name.to_lowercase().as_str() {
            "error" => response.error = text_of(&value)?,
            "message" => response.message = text_of(&value)?,
            _ => return Err(Unreadable),
        }
    }
    Ok(response)
}

/// The `data:` branch of `getData`, with what `IllegalArgumentException` says when it throws.
pub enum DataFailure {
    /// `"data URI must be base64 encoded: " + dataUri`, the `UserException`.
    NotBase64,
    /// The decoder's `IllegalArgumentException`, with its message.
    Decode(String),
}

/// `getData` on a `data:` URI.
///
/// ```java
/// if (!dataUri.matches("^data:.*;base64,.*")) {
/// ```
///
/// then `dataUri.replaceFirst("^data:.*;base64,", "")`. The `.*` is greedy, so what is decoded
/// is everything after the LAST `;base64,`, hence `rfind`. `.filter(..)` on the `Option` keeps
/// the match only when the URI also starts with `data:`.
pub fn decode_data_uri(uri: &str) -> Result<Vec<u8>, DataFailure> {
    let marker = ";base64,";
    let Some(at) = uri.rfind(marker).filter(|_| uri.starts_with("data:")) else {
        return Err(DataFailure::NotBase64);
    };
    java_base64_decode(&uri.as_bytes()[at + marker.len()..]).map_err(DataFailure::Decode)
}

/// `Base64.getDecoder().decode`, refusing what it refuses and saying what it says.
///
/// **How**: four characters make three bytes, accumulated six bits at a time into `bits` from
/// bit 18 down; `=` ends the input, valid only as the third character of a unit followed by a
/// second `=`, or as the fourth. A partial last unit of two or three characters yields one or
/// two bytes; one character alone is refused. Anything after the padding is refused with its
/// offset. The messages are the ones the tests assert. `byte as i8` reads the offending byte as
/// Java's signed `byte`, which is how the message prints a byte above 127 (as `-` and hex).
pub fn java_base64_decode(source: &[u8]) -> Result<Vec<u8>, String> {
    if source.is_empty() {
        return Ok(Vec::new());
    }
    if source.len() < 2 {
        return Err("Input byte[] should at least have 2 bytes for base64 bytes".to_string());
    }
    let value = |byte: u8| -> i32 {
        match byte {
            b'A'..=b'Z' => i32::from(byte - b'A'),
            b'a'..=b'z' => i32::from(byte - b'a') + 26,
            b'0'..=b'9' => i32::from(byte - b'0') + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' => -2,
            _ => -1,
        }
    };
    let mut out = Vec::with_capacity(source.len() / 4 * 3);
    let (mut bits, mut shift, mut at) = (0i32, 18i32, 0usize);
    while at < source.len() {
        let byte = source[at];
        at += 1;
        let digit = value(byte);
        if digit < 0 {
            if digit == -2 {
                let second = shift == 6
                    && (at == source.len() || {
                        at += 1;
                        source[at - 1] != b'='
                    });
                if second || shift == 18 {
                    return Err("Input byte array has wrong 4-byte ending unit".to_string());
                }
                break;
            }
            let signed = i32::from(byte as i8);
            let text = if signed < 0 {
                format!("-{:x}", -signed)
            } else {
                format!("{signed:x}")
            };
            return Err(format!("Illegal base64 character {text}"));
        }
        bits |= digit << shift;
        shift -= 6;
        if shift < 0 {
            out.extend_from_slice(&[(bits >> 16) as u8, (bits >> 8) as u8, bits as u8]);
            shift = 18;
            bits = 0;
        }
    }
    match shift {
        6 => out.push((bits >> 16) as u8),
        0 => out.extend_from_slice(&[(bits >> 16) as u8, (bits >> 8) as u8]),
        12 => return Err("Last unit does not have enough valid bits".to_string()),
        _ => {}
    }
    if at < source.len() {
        return Err(format!(
            "Input byte array has incorrect ending byte at {at}"
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bam_id_refuses_every_format() {
        let request = Request {
            id: "a.bam".to_string(),
            format: Some("BAM".to_string()),
            ..Request::default()
        };
        assert_eq!(
            request.validate().unwrap_err(),
            "Specified format: BAM is incompatible with id's file extension"
        );
    }

    #[test]
    fn base64_says_what_java_says() {
        assert_eq!(java_base64_decode(b"QUJD").unwrap(), b"ABC");
        assert_eq!(java_base64_decode(b"QUI").unwrap(), b"AB");
        assert_eq!(java_base64_decode(b"QQ==").unwrap(), b"A");
        assert_eq!(
            java_base64_decode(b"Q").unwrap_err(),
            "Input byte[] should at least have 2 bytes for base64 bytes"
        );
        assert_eq!(
            java_base64_decode(b"QQ=A").unwrap_err(),
            "Input byte array has wrong 4-byte ending unit"
        );
        assert_eq!(
            java_base64_decode(b"Q Q=").unwrap_err(),
            "Illegal base64 character 20"
        );
        assert_eq!(
            java_base64_decode(b"QUJDQ").unwrap_err(),
            "Last unit does not have enough valid bits"
        );
    }

    #[test]
    fn jackson_refuses_what_it_does_not_know() {
        assert!(read_response(r#"{"htsget":{"urls":[],"unknown":1}}"#).is_err());
        assert!(read_response(r#"{"other":{"urls":[]}}"#).is_err());
        assert!(read_response(r#"{"htsget":{"format":"bam"}}"#).is_err());
        let read = read_response(r#"{"htsget":{"URLS":[{"URL":"data:;base64,QQ=="}],"md5":1}}"#);
        assert_eq!(read.unwrap().md5.as_deref(), Some("1"));
    }
}
