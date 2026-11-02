#!/usr/bin/env python3
"""Run the reference's own forward pass and write every tensor it passes through, bit for bit.

The same steps as `nvscorevariants.py`, with the same model file, the same `ReferenceDataset` and
the same `DataLoader`, without Lightning's progress bars: a forward hook on every module records its
output, so the port is compared against each layer rather than only against the printed score.

Each record of the dump is, little-endian:

    u32 name length, name, u32 kind (0 f32, 1 i64), u32 dimensions, u64 each dimension, values

and the names are `<batch>/<layer>`. The script also writes the temporary score file the reference
would, so a run of this harness can be checked against a run of the tool.

    dump_layers.py <reference|read_tensor> <vcf> <fasta> <bam or -> <batch size> <dump> <scores>
"""

import io
import struct
import sys
import zipfile

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

from scorevariants.dataset import ReferenceDataset
from scorevariants.models.common import predictions_to_score
from scorevariants.readers import ReferenceTensorReader, TensorReader

JAR = "/gatk/gatk-package-4.6.2.0-local.jar"


def write_tensor(out, name, tensor):
    tensor = tensor.detach()
    encoded = name.encode()
    out.write(struct.pack("<I", len(encoded)) + encoded)
    if tensor.dtype == torch.float32:
        kind, data = 0, tensor.numpy().astype("<f4").tobytes(order="C")
    elif tensor.dtype == torch.int64:
        kind, data = 1, tensor.numpy().astype("<i8").tobytes(order="C")
    else:
        sys.exit(f"{name}: unexpected {tensor.dtype}")
    out.write(struct.pack("<II", kind, tensor.dim()))
    for size in tensor.shape:
        out.write(struct.pack("<Q", size))
    out.write(data)


def main():
    kind, vcf, fasta, bam, batch_size, dump, scores = sys.argv[1:8]
    torch.manual_seed(724)
    if kind == "reference":
        model_file = "1d_cnn_mix_train_full_bn.pt"
        reader = ReferenceTensorReader.from_files(vcf, fasta)
    else:
        model_file = "small_2d.pt"
        reader = TensorReader.from_files(vcf, fasta, bam)
    with zipfile.ZipFile(JAR) as jar:
        model = torch.load(io.BytesIO(jar.read(f"large/nvscorevariants/{model_file}")))
    model.eval()

    out = open(dump, "wb")
    current = {"batch": 0}

    def hook(name):
        def record(_module, _inputs, output):
            write_tensor(out, f"{current['batch']}/{name}", output)

        return record

    for name, module in model.named_modules():
        if name and not list(module.children()) and "loss_fn" not in name:
            module.register_forward_hook(hook(name))

    loader = DataLoader(ReferenceDataset(reader), batch_size=int(batch_size))
    with open(scores, "w") as tmp, torch.inference_mode():
        for number, batch in enumerate(loader):
            current["batch"] = number
            for key in ("type", "best_practices", "reference", "read_tensor"):
                if key in batch:
                    write_tensor(out, f"{number}/input.{key}", batch[key])
            predictions = model(batch)
            write_tensor(out, f"{number}/logits", predictions)
            probabilities = F.softmax(predictions, dim=1)
            write_tensor(out, f"{number}/softmax", probabilities)
            final = predictions_to_score(probabilities, batch["type"])
            write_tensor(out, f"{number}/score", final)
            for i in range(len(final)):
                tmp.write(
                    "%s\t%d\t%s\t[%s]\t%.3f\n"
                    % (batch["chrom"][i], batch["pos"][i], batch["ref"][i], batch["alt"][i], final[i])
                )
    out.close()


if __name__ == "__main__":
    main()
