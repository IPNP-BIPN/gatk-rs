//! Kryo 4's `Output`, as far as the PathSeq serializers use it.
//!
//! # Why this exists
//!
//! `PathSeqBuildKmers` writes its whole answer through `PSKmerUtils.writeKryoObject`, which is
//! `new Kryo()` and `kryo.writeObject(output, obj)`:
//!
//! ```java
//! final Output output = new Output(BucketUtils.createFile(uri));
//! final Kryo kryo = new Kryo();
//! kryo.writeObject(output, obj);
//! ```
//!
//! `PathSeqBuildReferenceTaxonomy` writes through `writeTaxonomyDatabase` in
//! `PSBuildReferenceTaxonomyUtils` instead, which is the same with references switched off first:
//!
//! ```java
//! final Kryo kryo = new Kryo();
//! kryo.setReferences(false);
//! Output output = new Output(new FileOutputStream(filePath));
//! kryo.writeObject(output, taxonomyDatabase);
//! ```
//!
//! A runner that writes anything else answers a different file, and the covering array compares a
//! binary output by its digest, so there is no comparison to be had without these bytes
//! (IPNP-BIPN/gatk-rs#1181).
//!
//! # What is measured and what is implemented
//!
//! Kryo's stream is a library's, not a specification, so every rule here comes from the
//! `kryo-stream` golden rather than from reading Kryo's source: the same rule
//! `docs/an-unspecified-order-that-reaches-the-output.md` states for an iteration order. GATK pins
//! Kryo `strictly [4,5)`, so this is Kryo 4 and the golden says so by construction.
//!
//! The four encodings the serializers reach for:
//!
//! - **a fixed `int` is four bytes, most significant first.** `writeInt(-1)` is `ffffffff` and
//!   `writeInt(1)` is `00000001`, which is what tells this apart from the varint below;
//! - **a fixed `long` is eight bytes**, the same way round;
//! - **`writeInt(value, true)` is a varint**, seven bits to a byte, least significant group first,
//!   with the top bit set on every byte but the last. It is written for a value the caller says is
//!   usually positive, and a negative one costs the full five bytes: `-300` is `d4fdffff0f`;
//! - **a string is one of two shapes.** More than one character and every one of them below
//!   `0x80`: the bytes themselves, with the top bit set on the LAST one and no length at all
//!   (`131567` is `3133313536b7`). Otherwise a length of the CHARACTER count plus one, its low six
//!   bits in a byte with `0x80` set, then UTF-8:
//!   a null string is `80`, an empty one `81`, `"A"` is `8241`, and `"Crème"` is `86` and six
//!   bytes for five characters.
//!
//! Nothing here writes a class or a schema, because `writeObject` is told the class by the reader.
//! What it does write is a **reference marker**: Kryo 4 has references ON by default, so the first
//! object in a stream is preceded by `01`. That byte is the object writer's (`write_object`), not
//! the primitives'.
//!
//! # Reading this file without reading Rust
//!
//! The Kryo API is `com.esotericsoftware.kryo.io.Output` (`writeInt(int)`, `writeLong(long)`,
//! `writeInt(int, boolean)`, `writeString(String)`) and `Kryo.writeObject(Output, Object)`. Kryo
//! has no local source, so every byte layout below is the Rust's, checked against the golden.
//!
//! - `&mut self` on a method means it may change the `Output` it is called on; `&self` only reads;
//! - `x as u32`, `x as u8`, `x as u64` are casts that never fail: `i32 as u32` keeps the same 32
//!   bits (so `-1` becomes `0xffffffff`, Java's `>>>` view of an `int`), and `as u8` keeps only
//!   the low eight bits, which is Java's `(byte)` cast;
//! - `Option<&str>` is a string that may be absent, which is how Java's nullable `String` arrives.

/// The bytes one Kryo `Output` has been handed, in order.
///
/// `Output` buffers and flushes; nothing here needs that distinction, because every stream these
/// tools write is built whole and then handed to a file.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Output {
    bytes: Vec<u8>,
}

impl Output {
    pub fn new() -> Self {
        Self::default()
    }

    /// `Output.writeInt(int)`: four bytes, most significant first.
    ///
    /// Big-endian, the opposite of the little-endian order BAM and BCF use: golden row `int-one` is
    /// `00000001`. `to_be_bytes` is the standard library's "this integer as big-endian bytes".
    pub fn write_int(&mut self, value: i32) {
        self.bytes.extend_from_slice(&value.to_be_bytes());
    }

    /// `Output.writeLong(long)`: eight bytes, most significant first.
    pub fn write_long(&mut self, value: i64) {
        self.bytes.extend_from_slice(&value.to_be_bytes());
    }

    /// `Output.writeInt(int, true)`: a varint over the value's THIRTY-TWO bits.
    ///
    /// The `true` is Kryo's "optimize positive", and it optimises nothing else: a negative value is
    /// its two's complement read as unsigned, so every one of them takes five bytes.
    ///
    /// Golden rows: `varint-300` is `ac02` (`0x2c | 0x80`, then `300 >> 7 = 2`) and
    /// `varint-minus-300` is `d4fdffff0f`.
    pub fn write_varint(&mut self, value: i32) {
        // The bits reinterpreted as unsigned, so the shift below is Java's `>>>`, not `>>`: a
        // sign-extending shift would never reach zero for a negative value.
        let mut remaining = value as u32;
        // `loop` is Java's `while (true)`; the only exit is the `return`.
        loop {
            if remaining < 0x80 {
                self.bytes.push(remaining as u8);
                return;
            }
            self.bytes.push((remaining as u8 & 0x7f) | 0x80);
            remaining >>= 7;
        }
    }

    /// `Output.writeString(String)`, both of its shapes.
    ///
    /// The ASCII shape needs at least two characters: a one-character ASCII string takes the length
    /// form, which is why `"A"` is `8241` and `"131567"` carries no length at all.
    ///
    /// What the golden covers: ASCII strings of one to forty characters (`string-long-ascii` is
    /// the longest) and one non-ASCII string whose characters are each one or two UTF-8 bytes
    /// (`string-utf8`). An ASCII string of sixty-four characters or more, and a character outside
    /// the Basic Multilingual Plane (one Rust `char`, two Java `char`s), are not measured.
    pub fn write_string(&mut self, value: Option<&str>) {
        // `let Some(text) = value else { ... }` unwraps the string or takes the `else` exit,
        // which here is the Java `null` case.
        let Some(text) = value else {
            // A null string is the length varint of zero, marked as the last byte.
            self.bytes.push(0x80);
            return;
        };
        // Unicode characters, not bytes: a Rust `&str` is UTF-8, so `text.len()` would count
        // `è` as two.
        let characters = text.chars().count();
        if characters == 0 {
            self.bytes.push(0x81);
            return;
        }
        if characters > 1 && text.is_ascii() {
            // Every byte but the last as it is, then the last with its top bit set: the reader
            // knows the string has ended when it meets a byte at or above `0x80`.
            // `&raw[..n]` is a view of the first `n` bytes, not a copy.
            let raw = text.as_bytes();
            self.bytes.extend_from_slice(&raw[..raw.len() - 1]);
            self.bytes.push(raw[raw.len() - 1] | 0x80);
            return;
        }
        // The count is of CHARACTERS and the bytes that follow are UTF-8, so a string whose
        // characters are not one byte each is longer than its own count says.
        self.write_utf8_length(characters);
        self.bytes.extend_from_slice(text.as_bytes());
    }

    /// The length form's prefix, which is `count + 1` and NOT an ordinary varint.
    ///
    /// The first byte carries the low SIX bits with `0x80` set, so a null string is `80`, an empty
    /// one `81` and a one-character one `82`. A value that does not fit in six bits sets `0x40` as
    /// well and continues in seven-bit groups.
    ///
    /// Only the one-byte form is measured: every string these serializers write is either ASCII,
    /// which takes the other shape, or a name far shorter than sixty-three characters. A longer
    /// non-ASCII string would follow the continuation below, and it has to be measured before it is
    /// believed.
    fn write_utf8_length(&mut self, count: usize) {
        // The `+ 1` is what frees `0` for null: null is `80`, empty `81`.
        let value = count as u64 + 1;
        if value >> 6 == 0 {
            self.bytes.push((value as u8) | 0x80);
            return;
        }
        self.bytes.push((value as u8 & 0x3f) | 0x40 | 0x80);
        let mut remaining = value >> 6;
        loop {
            if remaining < 0x80 {
                self.bytes.push(remaining as u8);
                return;
            }
            self.bytes.push((remaining as u8 & 0x7f) | 0x80);
            remaining >>= 7;
        }
    }

    /// `kryo.writeObject(output, obj)` with Kryo 4's default settings, which have references ON.
    ///
    /// The marker is the varint `1`, written BEFORE the serializer's own bytes, and it is the same
    /// `01` for every object written for the first time rather than a running count: the
    /// `kmer-set` golden row carries `01` before the `PSKmerSet` and `01` again before the nested
    /// set (`010000001f00000000000fffff01...`). What a second write of the SAME object would
    /// produce is not measured, and no serializer here writes one object twice. A serializer that
    /// turns references off for its own nested writes, as `PSTaxonomyDatabase` does, suppresses
    /// them for what it writes and not for itself.
    ///
    /// `referenced` is whether the `Kryo` has references on at this call, which the caller knows
    /// from the Java it ports. `write` is the serializer body, passed as a closure:
    /// `F: FnOnce(&mut Output)` reads "any function that takes the `Output` and is called once",
    /// Java's `Consumer<Output>`.
    pub fn write_object<F>(&mut self, referenced: bool, write: F)
    where
        F: FnOnce(&mut Output),
    {
        if referenced {
            self.write_varint(1);
        }
        write(self);
    }

    /// The stream so far.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// The stream as lower-case hex, which is what the golden carries.
    ///
    /// `.iter().map(..).collect()` is Java's `stream().map(..).collect(joining())`: each byte to
    /// two hex digits, then concatenated into one `String`.
    pub fn hex(&self) -> String {
        self.bytes
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }
}

// `#[cfg(test)]` compiles this module only for `cargo test`. The cases are the golden's own rows,
// so a change to the encoding fails here before it fails the conformance suite.
#[cfg(test)]
mod tests {
    use super::*;

    fn hex_of(write: impl FnOnce(&mut Output)) -> String {
        let mut output = Output::new();
        write(&mut output);
        output.hex()
    }

    #[test]
    fn a_fixed_int_is_four_bytes_most_significant_first() {
        assert_eq!(hex_of(|out| out.write_int(0)), "00000000");
        assert_eq!(hex_of(|out| out.write_int(1)), "00000001");
        assert_eq!(hex_of(|out| out.write_int(-1)), "ffffffff");
        assert_eq!(hex_of(|out| out.write_int(i32::MAX)), "7fffffff");
        assert_eq!(hex_of(|out| out.write_int(i32::MIN)), "80000000");
    }

    #[test]
    fn a_varint_costs_five_bytes_for_a_negative_value() {
        assert_eq!(hex_of(|out| out.write_varint(300)), "ac02");
        assert_eq!(hex_of(|out| out.write_varint(-300)), "d4fdffff0f");
    }

    #[test]
    fn a_string_takes_the_ascii_shape_only_above_one_character() {
        assert_eq!(hex_of(|out| out.write_string(None)), "80");
        assert_eq!(hex_of(|out| out.write_string(Some(""))), "81");
        assert_eq!(hex_of(|out| out.write_string(Some("A"))), "8241");
        assert_eq!(
            hex_of(|out| out.write_string(Some("131567"))),
            "3133313536b7"
        );
    }
}
