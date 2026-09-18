//! Kryo 4's `Output`, as far as the PathSeq serializers use it.
//!
//! # Why this exists
//!
//! `PathSeqBuildKmers` and `PathSeqBuildReferenceTaxonomy` write their whole answer through
//! `PSKmerUtils.writeKryoObject`, which is `new Kryo()` and `kryo.writeObject(output, obj)`. A
//! runner that writes anything else answers a different file, and the covering array compares a
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
//! object in a stream is preceded by `01`. That byte is the object writer's ([`write_object`]), not
//! the primitives'.

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
    pub fn write_varint(&mut self, value: i32) {
        let mut remaining = value as u32;
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
    pub fn write_string(&mut self, value: Option<&str>) {
        let Some(text) = value else {
            // A null string is the length varint of zero, marked as the last byte.
            self.bytes.push(0x80);
            return;
        };
        let characters = text.chars().count();
        if characters == 0 {
            self.bytes.push(0x81);
            return;
        }
        if characters > 1 && text.is_ascii() {
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
    /// The marker is a varint of the object's position in the stream's reference list, starting at
    /// one, and it is written BEFORE the serializer's own bytes. A serializer that turns references
    /// off for its own nested writes, as `PSTaxonomyDatabase` does, suppresses them for what it
    /// writes and not for itself.
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
    pub fn hex(&self) -> String {
        self.bytes
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }
}

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
