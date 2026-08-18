# The ML surface cannot be bit-identical, and it is one tool

Measured against GATK 4.6.2.0's own shipped model files, before any port was written. Recorded here
because the conclusion changes what a milestone can ask for, and a conclusion that lives only in an
issue comment is a conclusion the next person re-derives.

## CNNScoreVariants is gone, and TensorFlow with it

`DeprecatedToolsRegistry` retires it as of 4.6.1.0:

```java
deprecatedTools.put("CNNScoreVariants", Pair.of("4.6.1.0",
    "Please use the replacement tool NVScoreVariants instead, which produces virtually identical results"));
```

`CNNVariantTrain` is retired with it, and `scripts/gatkcondaenv.yml.template` contains **no TensorFlow
and no Keras**. The ML surface of the version this programme reproduces is one tool,
`NVScoreVariants`, on PyTorch 2.1.0 with an MKL BLAS. Any plan naming two tools is naming one that
does not ship.

## The reference does not reproduce itself

The same shipped model (`small_2d.pt`), the same deterministic input, `eval()` so dropout is off:

| comparison | floats differing | worst gap |
|---|---|---|
| CPU torch 2.1.0 against CPU torch 2.10.0 | 6 of 8 | 1 ULP |
| CPU against MPS, torch 2.1.0 | 7 of 8 | 4 ULP |
| CPU against MPS, torch 2.10.0 | 7 of 8 | 4 ULP |
| CPU rerun, same version | 0 of 8 | 0 |

Each run is deterministic: the last row is what rules out noise. What moves the bits is the kernel,
so ten minor versions of PyTorch move six outputs of eight, and a different accelerator moves four
bits of seven.

**A bit-identity claim against this tool is therefore unavailable at any effort**, including the
effort of staying inside PyTorch on the pinned version, because the claim would then be against one
build of one runtime on one accelerator rather than against the reference.

These numbers were taken on this laptop rather than in the pinned container, which for once is the
right place: the finding *is* that the answer depends on the machine, and a container would have
hidden the comparison that shows it.

## The shipped artefacts are pickled Python modules

`1d_cnn_mix_train_full_bn.pt` (6.4 MB) and `small_2d.pt` (2.0 MB) are zip archives holding
`archive/data.pkl` and raw tensors, with **no `code/` directory**: they are neither TorchScript nor
state dicts. `torch.load` returns a live `GATK_CNN_2D` object, which unpickles only because the class
is importable. The artefact references Python code, not only numbers.

| route | what the pickle forces |
|---|---|
| `tch` and rebuild in `tch::nn` | extract the tensors from the pickle, then transcribe the architecture |
| `tch` and `CModule::load` | impossible as shipped: there is no TorchScript in the file |
| `ort` or `tract` | an offline ONNX conversion, which is a Python build step |
| a pure-Rust framework | the same extraction, and kernels further from the reference |

No route reads the file as it stands.

## What follows

1. **Set a tolerance and justify it from the table above, or quarantine the tool** and report it as
   bio-identical, which the milestone text already allows for. What is not available is a byte claim.
2. **Choose the runtime on other grounds** -- dependency weight, licence, whether an ONNX conversion
   step is acceptable -- because bit-identity cannot discriminate between them.
3. **Scope it as one tool.**

The probe that produced the table is not committed; it becomes a suite if the tolerance route is
taken, and there is nothing to compare against until then.
