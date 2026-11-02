#!/usr/bin/env python3
"""What the reference's Python makes of the awkward corners, in the reference's image.

Prints one line for each: an annotation written some way in a header and a record, as what
`encode_annotations` gives or the exception it ends on; `random_generator.Random` over a run of
bounds, as java.util.Random is not what NumPy's integers make of it where they wrap; and the two
numbers `base_quality_to_p_hot_array` takes from `10.0 ** exponent`, as float bits, for all 256
qualities. crates/gatk-tools/src/nv_score_variants/encode.rs pins these.

    python3 probe_edges.py <scorevariants' parent directory> <work directory>
"""
import hashlib
import struct
import sys

sys.path.insert(0, sys.argv[1])
import numpy as np
from pysam import VariantFile
from scorevariants.encoders import encode_annotations
from scorevariants.random_generator import Random

HEAD = "##fileformat=VCFv4.2\n##contig=<ID=c1,length=2600>\n"
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

CASES = [
    ("Number=1,Type=Float", "1.5", "A"),
    ("Number=1,Type=Float", ".", "A"),
    ("Number=1,Type=Float", "", "A"),
    ("Number=1,Type=Float", "1,2", "A"),
    ("Number=1,Type=Float", "nan", "A"),
    ("Number=1,Type=Float", "inf", "A"),
    ("Number=1,Type=Float", "-inf", "A"),
    ("Number=1,Type=Float", "x", "A"),
    ("Number=1,Type=Integer", "7", "A"),
    ("Number=1,Type=Integer", ".", "A"),
    ("Number=1,Type=Integer", "1,2", "A"),
    ("Number=1,Type=Integer", "1.5", "A"),
    ("Number=2,Type=Integer", "1,2", "A"),
    ("Number=2,Type=Integer", "1", "A"),
    ("Number=A,Type=Float", "1.5", "A"),
    ("Number=A,Type=Float", ".", "A"),
    ("Number=A,Type=Float", "1.5,2.5", "A,C"),
    ("Number=R,Type=Float", "1.5", "A"),
    ("Number=G,Type=Float", "1.5", "A"),
    ("Number=.,Type=Float", "1.5", "A"),
    ("Number=.,Type=Float", ".", "A"),
    ("Number=.,Type=Float", "1.5,2.5", "A"),
    ("Number=.,Type=Integer", "3", "A"),
    ("Number=0,Type=Flag", None, "A"),
    ("Number=1,Type=String", "1.5", "A"),
    ("Number=1,Type=String", "7", "A"),
    ("Number=1,Type=String", "x", "A"),
    ("Number=1,Type=String", ".", "A"),
    ("Number=1,Type=String", "", "A"),
    ("Number=1,Type=String", "nan", "A"),
    ("Number=1,Type=String", " 2.5", "A"),
    ("Number=1,Type=String", "1_0", "A"),
    ("Number=.,Type=String", "1.5", "A"),
    ("Number=.,Type=String", "1.5,2", "A"),
    (None, "3.5", "A"),
    (None, None, "A"),
    (None, "", "A"),
    (None, "1,2", "A"),
]


def annotation_cases(work):
    for at, (definition, value, alts) in enumerate(CASES):
        path = f"{work}/edge{at}.vcf"
        with open(path, "w") as out:
            out.write(HEAD)
            if definition is not None:
                out.write(f'##INFO=<ID=MQ,{definition},Description="x">\n')
            out.write(COLUMNS)
            info = "MQ" if value is None else f"MQ={value}"
            out.write(f"c1\t500\t.\tG\t{alts}\t.\t.\t{info}\n")
        variant = next(iter(VariantFile(path)))
        try:
            encoded = encode_annotations(variant, ["MQ"])[0]
            shown = struct.pack("<f", np.float32(encoded)).hex()
        except Exception as error:
            shown = type(error).__name__
        print(f"annotation\t{definition}\t{value}\t{alts}\t{shown}")


def random_cases():
    Random.multiplier = np.int64(0x5DEECE66D)
    Random.addend = np.int64(0xB)
    Random.mask = np.int64((1 << 48) - 1)
    random = Random(np.int64(47382911))
    answers = [int(random.nextInt(n)) for n in range(129, 20129)]
    digest = hashlib.sha256(",".join(map(str, answers)).encode()).hexdigest()
    print("random\tfirst\t" + ",".join(map(str, answers[:8])))
    print("random\tat256\t" + str(answers[256 - 129]))
    print("random\tminimum\t" + str(min(answers)))
    print("random\tsha256\t" + digest)


def phred_cases():
    table = []
    for quality in range(256):
        exponent = float(-quality) / 10.0
        p = 1.0 - (10.0 ** exponent)
        not_p = (1.0 - p) / 3.0
        table.append(struct.pack("<ff", np.float32(p), np.float32(not_p)))
    print("phred\tsha256\t" + hashlib.sha256(b"".join(table)).hexdigest())


if __name__ == "__main__":
    annotation_cases(sys.argv[2])
    random_cases()
    phred_cases()
