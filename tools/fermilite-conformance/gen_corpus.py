#!/usr/bin/env python3
"""Seeded corpus of small assemblies for comparing gatk-fermilite with the bundled libfml.

    gen_corpus.py main OUT [--seed N] [--scale N]   defined inputs, native C and JNI both apply
    gen_corpus.py edge OUT [--seed N]               high bytes, NUL, quals too long: JNI only

Format, big endian: "FMLC1\n", then per assembly
    u32 len, name; i32 flag; u32 nreads; per read u32 len, bases; u32 len, quals.

The generator has its own xorshift so that the corpus does not depend on the Python version.
"""
import argparse
import struct
import sys

DEFAULT = 0xC0
AGGRESSIVE = 0x60
FLAGS = (DEFAULT, AGGRESSIVE)


class Rng:
    def __init__(self, seed):
        self.s = (seed * 0x9E3779B97F4A7C15 + 0x1234567) & (2**64 - 1) or 1
        for _ in range(16):
            self.u64()

    def u64(self):
        x = self.s
        x ^= (x << 13) & (2**64 - 1)
        x ^= x >> 7
        x ^= (x << 17) & (2**64 - 1)
        self.s = x
        return x

    def below(self, n):
        return (self.u64() >> 11) % n

    def between(self, a, b):
        return a + self.below(b - a + 1)

    def real(self):
        return (self.u64() >> 11) / float(1 << 53)

    def choice(self, xs):
        return xs[self.below(len(xs))]


COMP = bytes.maketrans(b"ACGTacgtN", b"TGCAtgcaN")


def rc(s):
    return s.translate(COMP)[::-1]


def genome(r, n, gc=0.5):
    out = bytearray(n)
    for i in range(n):
        if r.real() < gc:
            out[i] = 67 if r.below(2) else 71
        else:
            out[i] = 65 if r.below(2) else 84
    return bytes(out)


def mutate_haplotype(r, g, rate, indels=True):
    out = bytearray()
    i = 0
    while i < len(g):
        if r.real() < rate:
            kind = r.below(3) if indels else 0
            if kind == 0:
                out.append(r.choice([c for c in b"ACGT" if c != g[i]]))
                i += 1
            elif kind == 1:
                i += r.between(1, 12)
            else:
                out += genome(r, r.between(1, 12))
        else:
            out.append(g[i])
            i += 1
    return bytes(out)


def qual_model(r, name, n, errs):
    """A qual string for a read of n bases, `errs` the positions with an error."""
    if name == "phred33":
        q = bytearray(r.between(30, 40) + 33 for _ in range(n))
        for i in errs:
            q[i] = r.between(2, 15) + 33
    elif name == "illumina_tail":
        tail = r.below(n // 3 + 1)
        q = bytearray((r.between(32, 41) if i < n - tail else r.between(2, 20)) + 33 for i in range(n))
        for i in errs:
            q[i] = r.between(2, 12) + 33
    elif name == "raw_phred":
        # what a BAM's unconverted qualities look like: small numbers, never 0 here since 0 is the NUL
        q = bytearray(r.between(1, 45) for _ in range(n))
        for i in errs:
            q[i] = r.between(1, 5)
    elif name == "two":
        q = bytearray(2 for _ in range(n))
    elif name == "bang":
        q = bytearray(33 for _ in range(n))
    elif name == "all_high":
        q = bytearray(r.between(60, 126) for _ in range(n))
    elif name == "top":
        q = bytearray(r.between(1, 127) for _ in range(n))
    elif name == "const40":
        q = bytearray(40 + 33 for _ in range(n))
    elif name == "lowonly":
        q = bytearray(r.between(0, 19) + 33 for _ in range(n))
    else:
        raise ValueError(name)
    return bytes(q)


def sample_read(r, g, lens, sub, ins, dele, qual, circular=False, both=True, start=None):
    n = lens if isinstance(lens, int) else r.choice(lens)
    if circular:
        if start is None:
            start = r.below(len(g))
        s = (g + g + g)[start : start + n]
    else:
        if n > len(g):
            n = len(g)
        if start is None:
            start = r.below(len(g) - n + 1)
        s = g[start : start + n]
    if both and r.below(2):
        s = rc(s)
    out = bytearray()
    errs = []
    i = 0
    while i < len(s):
        x = r.real()
        if x < sub:
            errs.append(len(out))
            out.append(r.choice([c for c in b"ACGT" if c != s[i]]))
        elif x < sub + ins:
            errs.append(len(out))
            out.append(s[i])
            out.append(r.choice(b"ACGT"))
        elif x < sub + ins + dele:
            pass
        else:
            out.append(s[i])
        i += 1
    errs = [e for e in errs if e < len(out)]
    return bytes(out), qual_model(r, qual, len(out), errs) if out else b""


def sample_depth(r, g, depth, lens, sub=0.0, ins=0.0, dele=0.0, qual="phred33", circular=False):
    mean = lens if isinstance(lens, int) else sum(lens) // len(lens)
    n = max(1, int(depth * len(g) / max(mean, 1)))
    return [sample_read(r, g, lens, sub, ins, dele, qual, circular) for _ in range(n)]


QUALS = ["phred33", "illumina_tail", "raw_phred", "two", "bang", "all_high", "top", "const40", "lowonly"]
LENS = [30, 36, 50, 75, 100, 101, 125, 150, 151, 250, [75, 100, 151], [30, 151, 250], [36, 50, 76]]


def rand_lens(r):
    return r.choice(LENS)


def rand_err(r):
    return r.choice([(0, 0, 0), (0.001, 0, 0), (0.005, 0.0005, 0.0005), (0.01, 0.001, 0.001), (0.03, 0.002, 0.002)])


def rand_depth(r):
    return r.choice([2, 4, 8, 12, 20, 30, 50, 80, 120, 200])


def g_random(r):
    n = r.choice([100, 200, 400, 800, 1500])
    g = genome(r, n, r.choice([0.3, 0.5, 0.65]))
    d = rand_depth(r)
    if d * n > 60000:
        d = 60000 // n
    e = rand_err(r)
    return sample_depth(r, g, d, rand_lens(r), *e, qual=r.choice(QUALS))


def g_bubble(r):
    g = genome(r, r.choice([300, 600, 1000]))
    h = mutate_haplotype(r, g, r.choice([0.002, 0.005, 0.01, 0.03]), indels=bool(r.below(2)))
    d = r.choice([10, 20, 40, 60])
    e = rand_err(r)
    q = r.choice(["phred33", "illumina_tail", "const40"])
    return sample_depth(r, g, d, rand_lens(r), *e, qual=q) + sample_depth(r, h, d, rand_lens(r), *e, qual=q)


def g_tandem(r):
    unit = genome(r, r.choice([1, 2, 3, 5, 8, 13, 21, 40, 80]))
    left = genome(r, 200)
    right = genome(r, 200)
    g = left + unit * r.between(4, 2000 // len(unit) + 4) + right
    return sample_depth(r, g, r.choice([10, 30, 60]), rand_lens(r), *rand_err(r), qual=r.choice(QUALS[:3]))


def g_inverted(r):
    arm = genome(r, r.choice([40, 80, 150, 400]))
    spacer = genome(r, r.choice([0, 1, 2, 10, 30, 100]))
    g = genome(r, 150) + arm + spacer + rc(arm) + genome(r, 150)
    return sample_depth(r, g, r.choice([10, 30, 60]), rand_lens(r), *rand_err(r), qual=r.choice(QUALS[:3]))


def g_palindrome(r):
    """Reads that are their own reverse complement reach the is_rev_same branch of fml_fmi_gen."""
    arm = genome(r, r.choice([20, 40, 60, 75, 125]))
    pal = arm + rc(arm)
    reads = []
    q = r.choice(["const40", "phred33"])
    for _ in range(r.choice([5, 10, 20, 40])):
        reads.append((pal, qual_model(r, q, len(pal), [])))
    g = genome(r, 100) + pal + genome(r, 100)
    return reads + sample_depth(r, g, 30, r.choice([50, 75, 100]), qual=q)


def g_chimera(r):
    a = genome(r, r.choice([400, 800]))
    b = genome(r, r.choice([400, 800]))
    pa = r.between(100, len(a) - 1)
    pb = r.between(0, len(b) - 100)
    bb = rc(b) if r.below(2) else b
    micro = genome(r, r.choice([0, 0, 1, 3, 10, 30]))
    junction = a[:pa] + micro + bb[pb:]
    d = r.choice([10, 20, 40])
    e = rand_err(r)
    q = r.choice(["phred33", "illumina_tail"])
    reads = sample_depth(r, a, d, rand_lens(r), *e, qual=q)
    reads += sample_depth(r, b, d, rand_lens(r), *e, qual=q)
    reads += sample_depth(r, junction, d, rand_lens(r), *e, qual=q)
    return reads


def g_circle(r):
    g = genome(r, r.choice([100, 150, 300, 600, 1000]))
    return sample_depth(r, g, r.choice([10, 20, 40, 80]), r.choice([50, 75, 100, 151]), *rand_err(r),
                        qual=r.choice(QUALS[:3]), circular=True)


def g_low_complexity(r):
    n = r.choice([150, 300, 800])
    kind = r.below(5)
    if kind == 0:
        g = b"A" * n
    elif kind == 1:
        g = genome(r, 100) + b"A" * n + genome(r, 100)
    elif kind == 2:
        g = (b"AT" * n)[:n]
    elif kind == 3:
        g = bytes(r.choice(b"AAAAAAAC") for _ in range(n))
    else:
        g = genome(r, 100) + (b"GGC" * n)[:n] + b"T" * 40 + genome(r, 100)
    return sample_depth(r, g, r.choice([5, 20, 60, 120]), rand_lens(r), *rand_err(r), qual=r.choice(QUALS))


def dirty(r, reads, p_n, alphabet):
    out = []
    for b, q in reads:
        b = bytearray(b)
        for i in range(len(b)):
            if r.real() < p_n:
                b[i] = r.choice(alphabet)
        out.append((bytes(b), q))
    return out


def g_dirty(r):
    g = genome(r, r.choice([200, 500, 1000]))
    reads = sample_depth(r, g, r.choice([10, 30, 60]), rand_lens(r), *rand_err(r), qual=r.choice(QUALS[:3]))
    kind = r.below(4)
    if kind == 0:
        return dirty(r, reads, r.choice([0.001, 0.01, 0.1]), b"N")
    if kind == 1:
        return [(b.lower() if r.below(2) else b, q) for b, q in reads]
    if kind == 2:
        return dirty(r, reads, r.choice([0.001, 0.01]), b"RYSWKMBDHVNnryswkm-.*")
    return dirty(r, reads, 0.002, bytes(c for c in range(1, 128)))


def g_trivial(r):
    kind = r.below(7)
    if kind == 0:
        return []
    if kind == 1:
        g = genome(r, r.choice([5, 30, 100, 151]))
        return [(g, qual_model(r, "phred33", len(g), []))]
    if kind == 2:
        return [(b"", b"")] * r.between(1, 4)
    if kind == 3:
        return [(genome(r, r.between(1, 12)), b"I" * 12) for _ in range(r.between(1, 30))]
    if kind == 4:
        g = genome(r, r.choice([33, 34, 50, 100, 151]))
        q = qual_model(r, "phred33", len(g), [])
        return [(g, q)] * r.choice([2, 3, 4, 5, 8, 16, 100, 300])
    if kind == 5:
        g = genome(r, 200)
        reads = sample_depth(r, g, 10, r.choice([30, 36, 40]))
        return reads + [(b"", b"")] * 3
    return [(b"N" * 100, b"I" * 100)] * r.between(1, 20)


def g_pairs(r):
    """Interleaved mates, as FindBreakpointEvidence's FASTQs have them."""
    g = genome(r, r.choice([800, 1500, 3000]))
    n = int(r.choice([10, 20, 40]) * len(g) / 300)
    out = []
    e = rand_err(r)
    q = r.choice(["phred33", "illumina_tail", "raw_phred"])
    for _ in range(n):
        ins = r.between(250, 500)
        p = r.below(max(1, len(g) - ins))
        out.append(sample_read(r, g, 151, *e, q, both=False, start=p))
        m, mq = sample_read(r, g, 151, *e, q, both=False, start=p + ins - 151)
        out.append((rc(m), mq[::-1]))
    return out


def g_duplicates(r):
    g = genome(r, r.choice([200, 500]))
    reads = sample_depth(r, g, r.choice([4, 10, 30]), rand_lens(r), *rand_err(r), qual="phred33")
    return reads * r.between(2, 6)


def g_typical(r):
    """Hundreds of 151 bp reads around an SV: the shape GATK actually assembles."""
    a = genome(r, r.choice([600, 1000]))
    b = genome(r, r.choice([600, 1000]))
    j = a[: r.between(200, len(a) - 1)] + rc(b)[r.between(0, len(b) - 200) :]
    h = mutate_haplotype(r, j, 0.002)
    reads = []
    e = (0.005, 0.0005, 0.0005)
    for src, d in ((a, 15), (b, 15), (j, 20), (h, 10)):
        reads += sample_depth(r, src, d, 151, *e, qual="illumina_tail")
    n = len(reads)
    order = list(range(n))
    for i in range(n - 1, 0, -1):
        k = r.below(i + 1)
        order[i], order[k] = order[k], order[i]
    return [reads[i] for i in order]


def g_big(r):
    n_reads = r.choice([5000, 6000, 8000, 12000])
    g = genome(r, r.choice([5000, 20000, 50000]))
    h = mutate_haplotype(r, g, 0.002)
    reads = []
    for i in range(n_reads):
        reads.append(sample_read(r, g if i & 1 else h, 151, 0.005, 0.0005, 0.0005, "illumina_tail"))
    return reads


CATEGORIES = [
    ("random", g_random, 900),
    ("bubble", g_bubble, 500),
    ("tandem", g_tandem, 300),
    ("inverted", g_inverted, 250),
    ("palindrome", g_palindrome, 200),
    ("chimera", g_chimera, 400),
    ("circle", g_circle, 250),
    ("lowcomplex", g_low_complexity, 300),
    ("dirty", g_dirty, 400),
    ("trivial", g_trivial, 150),
    ("pairs", g_pairs, 200),
    ("dups", g_duplicates, 150),
    ("typical", g_typical, 150),
]


def write_record(out, name, flag, reads):
    nb = name.encode()
    out.write(struct.pack(">I", len(nb)) + nb + struct.pack(">iI", flag, len(reads)))
    for b, q in reads:
        out.write(struct.pack(">I", len(b)) + b + struct.pack(">I", len(q)) + q)


def main_corpus(out, seed, scale):
    n = 0
    for ci, (name, fn, count) in enumerate(CATEGORIES):
        for i in range(count * scale):
            r = Rng(seed * 1000003 + ci * 7919 + i)
            reads = fn(r)
            flag = FLAGS[i & 1] if i % 5 else None
            for f in FLAGS if flag is None else (flag,):
                write_record(out, f"{name}-{i}", f, reads)
                n += 1
    for i in range(6):
        r = Rng(seed * 1000003 + 1000 + i)
        reads = g_big(r)
        write_record(out, f"big-{i}", FLAGS[i & 1], reads)
        n += 1
    return n


def edge_corpus(out, seed):
    """What only the JNI can answer, since the bytes before seq_nt6_table are the .so's."""
    n = 0

    def emit(name, reads):
        nonlocal n
        for f in FLAGS:
            write_record(out, name, f, reads)
            n += 1

    for i in range(40):
        r = Rng(seed * 31 + i)
        g = genome(r, 400)
        reads = sample_depth(r, g, 30, 75, qual="const40")
        # bytes whose table entry is a zero of .data, whichever way it is loaded
        for byte in (0xD0, 0xD5, 0xDF, 0xE0, 0xE8, 0xF0, 0xFF):
            emit(f"highbyte-{byte:02x}-{i}", dirty(r, reads, r.choice([0.001, 0.01]), bytes([byte])))
        # a NUL in the bases: the JNI's two strings go out of step
        b0, q0 = reads[0]
        mid = len(b0) // 2
        emit(f"nulbase-{i}", [(b0[:mid] + b"\0" + b0[mid + 1 :], q0)] + reads[1:])
        # quals shorter by a NUL at the very end are refused by us: see nulqual below
        # quals longer than the bases, with room left by nothing: BufferOverflowException for the last
        emit(f"longqual-{i}", reads[:-1] + [(reads[-1][0], reads[-1][1] + b"II")])
        # quals longer on a read followed by one with shorter quals: the buffer still has room
        b1, q1 = reads[1]
        emit(f"quals-shift-{i}", [(reads[0][0], reads[0][1] + b"IIII"), (b1, q1[:-4])] + reads[2:])
    # a Phred 0 in FastqRead's raw qualities: truly undefined, so this says what it does in practice.
    # Two identical runs of each, to see whether the heap decides.
    for i in range(24):
        r = Rng(seed * 37 + i)
        g = genome(r, 500)
        reads = sample_depth(r, g, 40, 100, qual="raw_phred")
        share = (0.002, 0.05, 1.0)[i % 3]
        out_reads = []
        for b, q in reads:
            if r.real() < share:
                q = bytearray(q)
                q[r.below(len(q))] = 0
                q = bytes(q)
            out_reads.append((b, q))
        for run in (1, 2):
            emit(f"nulqual-{i}-run{run}", out_reads)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["main", "edge"])
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--scale", type=int, default=1)
    a = ap.parse_args()
    with open(a.out, "wb") as out:
        out.write(b"FMLC1\n")
        n = main_corpus(out, a.seed, a.scale) if a.kind == "main" else edge_corpus(out, a.seed)
    print(f"{n} assemblies written to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
