//! `DumpTabixIndex`, ported from `org.broadinstitute.hellbender.tools.DumpTabixIndex`
//! (GATK 4.6.2.0).
//!
//! The smallest tool of the reporting-walker archetype: it reads a `.tbi`, which is a gzipped
//! little-endian structure, and prints it as text. The whole observable is that text.
//!
//! # The bin number is turned back into a range by a ladder
//!
//! ```java
//! } else if ( binNo <= 584 ) {
//!     final int binStart = binNo - 73;
//!     writer.print(binNo + "\t" + tigName + ":" + binStart + "M-" + (binStart + 1) + "M\t");
//! } else if ( binNo <= 4680 ) {
//! ```
//!
//! Six cases, whose boundaries are 1, 8, 72, 584 and 4680, and whose unit changes from `M` to `K`
//! at 585. The arithmetic is the reference's and is not derived from the BAI specification here.
//!
//! # The pseudobin summary prints the wrong field
//!
//! ```java
//! "\tend=" + Long.toHexString(chunkEnd >>> 16) + ":" + Long.toHexString(chunkStart & 0xffff)
//! ```
//!
//! The low half of `end=` comes from **chunkStart**, where every other place uses the matching
//! value. It reaches the output, so this port reproduces it: a suite comparing text would fail on
//! the corrected version, and the tool's users read this text.

use std::fmt::Write as _;

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TabixError {
    /// The first four bytes were not `TBI\1`.
    WrongMagic,
    /// A pseudobin whose chunk count is not two.
    PseudobinChunks(i32),
    /// The contig names block did not have its declared length.
    ContigNameLength,
    /// Data after what should have been the end.
    TrailingData,
    /// The stream ended in the middle of a field.
    PastEof,
}

impl TabixError {
    /// The message `UserException` carries, or the `IOException`'s for the last one.
    pub fn message(&self) -> String {
        match self {
            TabixError::WrongMagic => "Incorrect magic number for tabix index".to_string(),
            TabixError::PseudobinChunks(count) => format!("pseudobin has {count} chunks"),
            TabixError::ContigNameLength => {
                "Contig names didn't have the correct length.".to_string()
            }
            TabixError::TrailingData => "Unexpected data follows index.".to_string(),
            TabixError::PastEof => "Tried to read past EOF".to_string(),
        }
    }
}

/// A cursor over the decompressed index, little-endian and a byte at a time.
struct Cursor<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Cursor<'a> {
    fn byte(&mut self) -> Result<i32, TabixError> {
        let value = *self.bytes.get(self.position).ok_or(TabixError::PastEof)?;
        self.position += 1;
        Ok(value as i32)
    }

    fn short(&mut self) -> Result<i32, TabixError> {
        Ok(self.byte()? | (self.byte()? << 8))
    }

    fn int(&mut self) -> Result<i32, TabixError> {
        Ok(self.short()? | (self.short()? << 16))
    }

    /// `readLong`: the low int masked to unsigned and the high one **not**, so a value above 2^63
    /// comes out negative, exactly as Java's `long` does.
    fn long(&mut self) -> Result<i64, TabixError> {
        let low = self.int()? as i64 & 0xffff_ffff;
        let high = self.int()? as i64;
        Ok(low | (high << 32))
    }

    fn at_end(&self) -> bool {
        self.position >= self.bytes.len()
    }
}

/// `Long.toHexString`, which prints an unsigned hexadecimal with no leading zeroes.
fn hex(value: i64) -> String {
    format!("{:x}", value as u64)
}

/// `dumpTabixIndex(is, writer)`: the whole text, from the **decompressed** index.
///
/// The gzip layer is the caller's, because the tool's two failures live on different sides of it:
/// a file that is not gzipped fails inside `java.util.zip` before the magic is ever checked.
pub fn dump_tabix_index(bytes: &[u8]) -> Result<String, TabixError> {
    let mut cursor = Cursor { bytes, position: 0 };
    let mut out = String::new();

    // The fourth byte is compared against the NUMBER 1, not the character.
    if cursor.byte()? != b'T' as i32
        || cursor.byte()? != b'B' as i32
        || cursor.byte()? != b'I' as i32
        || cursor.byte()? != 1
    {
        return Err(TabixError::WrongMagic);
    }

    let n_tigs = cursor.int()?;
    let format = cursor.int()?;
    let seq_col = cursor.int()?;
    let beg_col = cursor.int()?;
    let end_col = cursor.int()?;
    let meta = cursor.int()? as u8 as char;
    let skip = cursor.int()?;
    let names_len = cursor.int()?;

    out.push_str("#tigs\tformat\tseqCol\tbegCol\tendCol\tmetaChr\tskip\n");
    let _ = writeln!(
        out,
        "{n_tigs}\t{format}\t{seq_col}\t{beg_col}\t{end_col}\t{meta}\t{skip}\n"
    );

    let names = read_contig_names(&mut cursor, n_tigs, names_len)?;
    for name in &names {
        let _ = writeln!(out, "{name} binned index:");
        let mut bins = cursor.int()?;
        while bins > 0 {
            bins -= 1;
            let bin_no = cursor.int()?;
            let mut chunks = cursor.int()?;

            if bin_no > 37448 {
                if chunks != 2 {
                    return Err(TabixError::PseudobinChunks(chunks));
                }
                let chunk_start = cursor.long()?;
                let chunk_end = cursor.long()?;
                let mapped = cursor.long()?;
                let unmapped = cursor.long()?;
                // `end=` takes its low half from chunk_START, which is the reference's own slip.
                let _ = writeln!(
                    out,
                    "{name} summary: mapped={mapped}\tplaced={unmapped}\tstart={}:{}\tend={}:{}",
                    hex(((chunk_start as u64) >> 16) as i64),
                    hex(chunk_start & 0xffff),
                    hex(((chunk_end as u64) >> 16) as i64),
                    hex(chunk_start & 0xffff),
                );
                continue;
            }

            out.push_str(&bin_range(bin_no, name));
            while chunks > 0 {
                chunks -= 1;
                let chunk_start = cursor.long()?;
                let chunk_end = cursor.long()?;
                let _ = write!(
                    out,
                    "\t{}:{}->{}:{}",
                    hex(((chunk_start as u64) >> 16) as i64),
                    hex(chunk_start & 0xffff),
                    hex(((chunk_end as u64) >> 16) as i64),
                    hex(chunk_end & 0xffff),
                );
            }
            out.push('\n');
        }

        let mut intervals = cursor.int()?;
        let mut kilobases = 0;
        out.push('\n');
        let _ = writeln!(out, "{name} linear index:");
        while intervals > 0 {
            intervals -= 1;
            let offset = cursor.long()?;
            let _ = writeln!(
                out,
                "{kilobases}K\t{}:{}",
                hex(((offset as u64) >> 16) as i64),
                hex(offset & 0xffff)
            );
            kilobases += 16;
        }
    }

    // Whatever follows the last contig is one long, and anything after that is refused.
    if !cursor.at_end() {
        let unplaced = cursor.long()?;
        let _ = writeln!(out, "{unplaced} unplaced reads.");
        if !cursor.at_end() {
            return Err(TabixError::TrailingData);
        }
    }
    Ok(out)
}

/// The ladder that turns a bin number back into a printed range, with its trailing tab.
fn bin_range(bin_no: i32, name: &str) -> String {
    if bin_no == 0 {
        format!("{bin_no}\t{name}:0M-512M\t")
    } else if bin_no <= 8 {
        let start = (bin_no - 1) * 64;
        format!("{bin_no}\t{name}:{start}M-{}M\t", start + 64)
    } else if bin_no <= 72 {
        let start = (bin_no - 9) * 8;
        format!("{bin_no}\t{name}:{start}M-{}M\t", start + 8)
    } else if bin_no <= 584 {
        let start = bin_no - 73;
        format!("{bin_no}\t{name}:{start}M-{}M\t", start + 1)
    } else if bin_no <= 4680 {
        let start = (bin_no - 585) * 128;
        format!("{bin_no}\t{name}:{start}K-{}K\t", start + 128)
    } else {
        let start = (bin_no - 4681) * 16;
        format!("{bin_no}\t{name}:{start}K-{}K\t", start + 16)
    }
}

/// `readContigNames`: NUL-terminated names whose **total** length is checked.
fn read_contig_names(
    cursor: &mut Cursor<'_>,
    count: i32,
    names_len: i32,
) -> Result<Vec<String>, TabixError> {
    let mut names = Vec::with_capacity(count.max(0) as usize);
    let mut read = 0;
    for _ in 0..count {
        let mut name = String::new();
        loop {
            let next = cursor.byte()?;
            if next == 0 {
                break;
            }
            name.push(next as u8 as char);
            read += 1;
        }
        read += 1;
        names.push(name);
    }
    if read != names_len {
        return Err(TabixError::ContigNameLength);
    }
    Ok(names)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_ladder_changes_unit_at_five_hundred_and_eighty_five() {
        assert_eq!(bin_range(0, "chr1"), "0\tchr1:0M-512M\t");
        assert_eq!(bin_range(1, "chr1"), "1\tchr1:0M-64M\t");
        assert_eq!(bin_range(9, "chr1"), "9\tchr1:0M-8M\t");
        assert_eq!(bin_range(73, "chr1"), "73\tchr1:0M-1M\t");
        assert_eq!(bin_range(584, "chr1"), "584\tchr1:511M-512M\t");
        assert_eq!(bin_range(585, "chr1"), "585\tchr1:0K-128K\t");
        assert_eq!(bin_range(4681, "chr1"), "4681\tchr1:0K-16K\t");
        assert_eq!(bin_range(4685, "chr1"), "4685\tchr1:64K-80K\t");
    }

    #[test]
    fn the_fourth_magic_byte_is_a_number_and_not_a_character() {
        let letter = [b'T', b'B', b'I', b'1', 0, 0, 0, 0];
        assert_eq!(
            dump_tabix_index(&letter).unwrap_err(),
            TabixError::WrongMagic
        );
    }

    #[test]
    fn a_long_is_two_ints_of_which_only_the_low_one_is_unsigned() {
        // 0xFFFFFFFF as the low half and 0 as the high half is 4294967295, not -1.
        let bytes = [0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0];
        let mut cursor = Cursor {
            bytes: &bytes,
            position: 0,
        };
        assert_eq!(cursor.long().expect("eight bytes"), 4_294_967_295);
    }
}
