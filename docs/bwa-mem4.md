# BWA through bwa-mem4, against the JNI GATK runs

Issue #1390.

**Verdict: adopted, and already landed in #1397.** bwa-mem4 at rev 41b4612 (MIT) replaced
gatk-bwamem-jni for every tool that needs BWA-MEM:

- `BwaSpark` and `BwaAndMarkDuplicatesPipelineSpark`;
- `PathSeqBwaSpark`, `PathSeqPipelineSpark` and `PathSeqFilterSpark --filter-bwa-image`;
- `ReadsPipelineSpark --align`.

Measured against the reference's JNI, every field `applyAlignment` builds agrees, except `MC`, which
bwa-mem4 adds and the port drops. Still open: FilterAlignmentArtifacts' realignment, and the
contig alignment of the SV pipeline. Both are separate work that would call the same aligner.

## Measurements

| What | Result | Where |
| --- | --- | --- |
| JNI vs bwa-mem4 on the fields `BwaMemAlignmentUtils.applyAlignment` sets: flag, contig, start, end, MAPQ, CIGAR, NM, AS, XS, MD, XA, mate, TLEN | Identical on 6,000 paired and 6,205 single-end records | `crates/gatk-tools/src/bwa_aligner.rs` module doc |
| `BwaSpark` against the reference on three seeds, paired and `--single-end-alignment` | 66,000 pairs, no difference on any SAM field or tag other than `MC` | #1397 |
| `BwaSpark` covering array, t=2 | 15/15 tuples (100%), 2,200 probes, 0/8 diverge | `docs/VALIDATION.md` |
| `BwaAndMarkDuplicatesPipelineSpark` | 19/19 (100%), 496 probes, 0/4 diverge | same |
| `PathSeqBwaSpark` | 22/22 (100%), 1,512 probes, 0/6 diverge | same |
| `BwaMemIndexImageCreator` | 15/15 (100%), 1,038 probes, 0/8 diverge, on every byte that depends on the FASTA | same |
| Thread independence (bwa-mem4's pool is the only parallel path in the Spark tools) | 1 and 4 threads give the same exit, stdout and file digests on every row | `docs/spark-parallel.md`, `spark-threads` CI job |

## The issue's questions, answered

- **The batch.** Pair statistics are computed per batch, and the JNI's batch is a partition. `-K` is
  set beyond anything a partition can hold, so bwa-mem4's batch is the partition too.
- **The `.img` index.** The port reads the reference's image directly. `fasta_from_image` rebuilds
  the FASTA from the image's `.pac` and annotations. bwa-mem4 then indexes it in a temporary
  directory, and its index files never leave the run.
- **Which API.** bwa-mem4's `mem` command is called as a library, through its own clap argument
  types, and gives SAM. The alignment's own fields are read from that SAM. Everything
  `applyAlignment` takes from the input (the name, the bases as stored, the qualities, the read
  group) comes from the input record, so nothing is lost in the round trip through SAM text.
- **Provenance.** bwa-mem4 is IPNP-BIPN's MIT reimplementation of bwa-mem2 (MIT), which is itself
  derived from lh3/bwa. Transcribing lh3/bwa (GPL-3.0) stays ruled out.
- **Pin.** A git dependency at `41b46129606d0f996e72a91557b1c2f0f0853121` in the workspace
  `Cargo.toml`, with default features off: the binary's allocator, htslib output and FASTQ
  decompressors are not needed by a library call.

## Known limits

- `MC` is dropped. Only bwa-mem4's CLI emits it, and the JNI has no field for it.
- An image that is not an image kills the reference's JVM, so no row measures that refusal.
- `PathSeqPipelineSpark` adds scores one partition per input. Past
  `--pipeline-reads-per-partition`, a score could differ in its last bit, and the `psp-batches` rows
  measure that.
- No speed comparison with the JNI was made. The criterion was byte identity, and the JNI runs only
  inside the oracle's container.
