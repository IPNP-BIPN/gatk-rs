#!/usr/bin/env python3
"""One shot: the two shipped models' tensors, out of the pickled `nn.Module`s, into plain files.

`NVScoreVariants` ships `1d_cnn_mix_train_full_bn.pt` and `small_2d.pt` inside the GATK jar. Each is
a pickled module that only opens with the `scorevariants` classes on the path, so no Rust code can
read it. This writes what the port needs and nothing else, little-endian:

    b"NVSW1\\n"
    u32 number of tensors
    for each, in `state_dict()` order, dropping `num_batches_tracked`:
        u32 name length, name, u32 dimensions, u64 each dimension, 32-bit floats in C order

The two files come to about 8 MB, so they are not committed: the port checks each against the
SHA-256 in `crates/gatk-tools/src/nv_score_variants/weights.rs`, and this script prints those hashes
along with the hash of the model it read, so a file that moved is refused rather than scored with.

    docker run --rm --platform linux/amd64 -v "$PWD/tools/nvscorevariants:/tool:ro" -v "$OUT:/out" \\
        broadinstitute/gatk:4.6.2.0 python3 /tool/extract_weights.py /out
"""

import hashlib
import io
import struct
import sys
import zipfile

import torch

JAR = "/gatk/gatk-package-4.6.2.0-local.jar"
MODELS = {
    "1d_cnn_mix_train_full_bn.pt": "nvscore-1d.weights",
    "small_2d.pt": "nvscore-2d.weights",
}


def main():
    out = sys.argv[1]
    with zipfile.ZipFile(JAR) as jar:
        for model, target in MODELS.items():
            shipped = jar.read(f"large/nvscorevariants/{model}")
            print(f"{hashlib.sha256(shipped).hexdigest()}  {model}")
            tensors = torch.load(io.BytesIO(shipped)).state_dict()
            named = [(k, v) for k, v in tensors.items() if not k.endswith("num_batches_tracked")]
            payload = bytearray(b"NVSW1\n")
            payload += struct.pack("<I", len(named))
            for name, tensor in named:
                if tensor.dtype != torch.float32:
                    sys.exit(f"{model}: {name} is {tensor.dtype}, not 32-bit floats")
                encoded = name.encode()
                payload += struct.pack("<I", len(encoded)) + encoded
                payload += struct.pack("<I", tensor.dim())
                for size in tensor.shape:
                    payload += struct.pack("<Q", size)
                payload += tensor.numpy().astype("<f4").tobytes(order="C")
            with open(f"{out}/{target}", "wb") as file:
                file.write(payload)
            print(f"{hashlib.sha256(payload).hexdigest()}  {target}  {len(payload)} bytes")


if __name__ == "__main__":
    main()
