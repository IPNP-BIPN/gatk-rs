#!/usr/bin/env python3
"""Write a BAM of about N GiB the way `crates/gatk-cli/tests/bounded_memory.rs` does.

    python3 tools/speed/scale_bam.py OUT.bam 2.5

One header block, then one compressed unit of 4096 records repeated: BGZF blocks stand alone, so
the file is written at the speed of the disk rather than of a deflater. The records are the test's:
150 M at chr1:1000, read group rg1, xorshift bases and qualities that hardly compress. Prints the
read count.
"""

import struct
import sys
import zlib

READS_PER_UNIT = 4096
MASK = (1 << 64) - 1
HEADER_TEXT = b"@SQ\tSN:chr1\tLN:1000000\n@RG\tID:rg1\tSM:sample1\n"
TERMINATOR = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
SEQ_CODE = {ord("A"): 1, ord("C"): 2, ord("G"): 4, ord("T"): 8}


def bgzf(data):
    """BGZF blocks of at most 64 KiB - 256 B input each, as htsjdk writes them."""
    out = bytearray()
    for start in range(0, len(data), 0xFF00):
        chunk = data[start : start + 0xFF00]
        packer = zlib.compressobj(5, zlib.DEFLATED, -15)
        deflated = packer.compress(chunk) + packer.flush()
        out += struct.pack("<4BI2BH2BHH", 0x1F, 0x8B, 8, 4, 0, 0, 0xFF, 6, 66, 67, 2, len(deflated) + 25)
        out += deflated + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return bytes(out)


def record(index, state):
    def nxt():
        state[0] ^= (state[0] << 13) & MASK
        state[0] ^= state[0] >> 7
        state[0] ^= (state[0] << 17) & MASK
        return state[0]

    bases = bytes(b"ACGT"[nxt() % 4] for _ in range(150))
    quals = bytes(nxt() % 41 for _ in range(150))
    name = f"r{index}".encode() + b"\0"
    seq = bytes(SEQ_CODE[bases[i]] << 4 | SEQ_CODE[bases[i + 1]] for i in range(0, 150, 2))
    # refID, pos (0-based), l_read_name, mapq, bin (reg2bin(999, 1149)), n_cigar, flag, l_seq,
    # mate refID, mate pos, tlen
    body = struct.pack("<iiBBHHHIiii", 0, 999, len(name), 60, 4681, 1, 0, 150, -1, -1, 0)
    body += name + struct.pack("<I", 150 << 4) + seq + quals + b"RGZrg1\0"
    return struct.pack("<i", len(body)) + body


def header():
    ref = b"chr1\0"
    return (
        b"BAM\1" + struct.pack("<i", len(HEADER_TEXT)) + HEADER_TEXT
        + struct.pack("<ii", 1, len(ref)) + ref + struct.pack("<i", 1_000_000)
    )


def main(path, gib):
    state = [0x9E3779B97F4A7C15]
    unit = bgzf(b"".join(record(i, state) for i in range(READS_PER_UNIT)))
    repeats = int(gib * (1 << 30)) // len(unit) + 1
    with open(path, "wb") as out:
        out.write(bgzf(header()))
        for _ in range(repeats):
            out.write(unit)
        out.write(TERMINATOR)
    print(repeats * READS_PER_UNIT)


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]))
