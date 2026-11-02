#!/usr/bin/env python3
"""Compare the port's TrainVariantAnnotationsModel and ScoreVariantAnnotations with the reference's.

Usage: compare_port.py <reference directory> <port directory>

The h5py files and the VCF must match byte for byte. The scorer must match with the digits of the
reference's temporary script name zeroed; those are random in every run of the reference, and one
name in ten is a digit shorter, which only moves the frames. score.scores.hdf5 goes through Java's
HDF5 library in the reference, which holds the time and lays the file out its own way: its doubles
are compared by the isolation_forest_measured test.
"""
import os
import re
import sys

EXACT = [
    "score.vcf",
    "train.snp.trainingScores.hdf5",
    "train.snp.calibrationScores.hdf5",
    "train.indel.trainingScores.hdf5",
    "train.indel.calibrationScores.hdf5",
]
PICKLES = ["train.snp.scorer.pkl", "train.indel.scorer.pkl"]
SCRIPT_NAME = re.compile(rb"/tmp/isolation-forest\.(\d+)\.py")
# The digits of `File.createTempFile`'s name, as many as the reference happened to draw.
PORT_DIGITS = 20


def read(directory, name):
    with open(os.path.join(directory, name), "rb") as stream:
        return stream.read()


def first_difference(a, b):
    for index, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return index
    return min(len(a), len(b))


def masked(data):
    return SCRIPT_NAME.sub(lambda m: b"/tmp/isolation-forest." + b"0" * len(m.group(1)) + b".py", data)


def unframed(data):
    """The opcodes without the FRAME headers, which only a frame's own length can place."""
    out, at = bytearray(data[:2]), 2
    while at < len(data):
        if data[at] == 0x95:
            length = int.from_bytes(data[at + 1:at + 9], "little")
            out += data[at + 9:at + 9 + length]
            at += 9 + length
        else:
            out += data[at:]
            break
    return bytes(out)


def normal(data):
    """The script's name, and the length before it, as one string."""
    return re.sub(rb".\/tmp/isolation-forest\.\d+\.py", b"?NAME", data, count=1, flags=re.S)


def main():
    reference, port = sys.argv[1], sys.argv[2]
    failed = 0
    for name in EXACT:
        a, b = read(reference, name), read(port, name)
        same = a == b
        print(f"{name}: {'identical' if same else f'differs at {first_difference(a, b)}'} ({len(a)}, {len(b)})")
        failed += not same
    for name in PICKLES:
        a, b = read(reference, name), read(port, name)
        digits = len(SCRIPT_NAME.search(a).group(1))
        if digits == PORT_DIGITS:
            a = masked(a)
        else:
            # A shorter name moves the frames along with it: compare what is in them.
            print(f"{name}: the reference drew {digits} digits; comparing without the frames")
            a, b = normal(unframed(a)), normal(unframed(b))
        same = a == b
        print(f"{name}: {'identical but the name' if same else f'differs at {first_difference(a, b)}'} ({len(a)}, {len(b)})")
        failed += not same
    print("FAILED" if failed else "OK")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
