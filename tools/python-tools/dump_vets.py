#!/usr/bin/env python
"""What measure_vets.sh left, as the lines crates/gatk-cli/tests/python_tools_end_to_end.rs reads.

Usage: dump_vets.py <fixture directory> (<pass> <directory>)...

    input  <file>         <hex>        a fixture file
    file   <pass> <file>  <hex>        a file that is the same in every run of the reference
    scorer <pass> <file>  <sha256>     the pickle, without its frames and its temporary script name
    scores <pass> <file>  <hex>        the doubles of a file Java wrote, which holds the time

The name in a scorer is random in every run, and one name in ten is a digit shorter, which moves the
frames along with it; the rest of the pickle is a function of the fit. CPU features go on a comment
line: scores can depend on them, and a golden should say what it was made on.
"""
import hashlib
import os
import re
import sys

import h5py
import numpy as np


def read(*path):
    with open(os.path.join(*path), "rb") as stream:
        return stream.read()


def normal_scorer(data):
    out, at = bytearray(data[:2]), 2
    while at < len(data):
        if data[at] == 0x95:
            length = int.from_bytes(data[at + 1 : at + 9], "little")
            out += data[at + 9 : at + 9 + length]
            at += 9 + length
        else:
            out += data[at:]
            break
    return re.sub(rb"./tmp/isolation-forest\.\d+\.py", b"?NAME", bytes(out), count=1, flags=re.S)


def main():
    features = np.core._multiarray_umath.__cpu_features__
    print("# cpu " + " ".join(name for name, on in features.items() if on))
    fixture = sys.argv[1]
    for name in ("calls.vcf", "resource.vcf"):
        print(f"input\t{name}\t{read(fixture, name).hex()}")
    passes = sys.argv[2:]
    for tag, directory in zip(passes[::2], passes[1::2]):
        for name in ("score.vcf",) + tuple(
            f"train.{kind}.{what}Scores.hdf5"
            for kind in ("snp", "indel")
            for what in ("training", "calibration")
        ):
            print(f"file\t{tag}\t{name}\t{read(directory, name).hex()}")
        print(f"input\t{tag}.extract.annot.hdf5\t{read(directory, 'extract.annot.hdf5').hex()}")
        for kind in ("snp", "indel"):
            digest = hashlib.sha256(normal_scorer(read(directory, f"train.{kind}.scorer.pkl"))).hexdigest()
            print(f"scorer\t{tag}\ttrain.{kind}.scorer.pkl\t{digest}")
        with h5py.File(os.path.join(directory, "score.scores.hdf5"), "r") as stream:
            values = stream["/data/scores"][()].astype("<f8").tobytes().hex()
        print(f"scores\t{tag}\tscore.scores.hdf5\t{values}")


if __name__ == "__main__":
    main()
