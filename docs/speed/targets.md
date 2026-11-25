# Targets, chosen from the baseline

Milestone S, S.5 (#112). A path is on this list because the baseline ([baseline.md](baseline.md))
says it is slow, not because it looks slow. Each entry carries its measured ratio, where the time
goes according to a profile, the optimisation and the argument for why it cannot move a byte
(written before the change, not after), and the suite that proves it did not.

The baseline names four tools whose steady ratio reaches 1. Two are targets; two are not yet.

## 1. `jmath`'s correctly rounded logarithm (`ModelSegments`, `VariantRecalibrator`)

**Ratio.** `ModelSegments` 18.4 steady, 4.83 cold (the only tool slower than the reference even
with the JVM's start-up charged to it); `VariantRecalibrator` 4.63 steady.

**Where the time goes.** A native profile of `ModelSegments`' baseline row (macOS `sample`): 94% of
samples in `jmath::log::ln_dd`, the 22-term double-double series behind `jmath::math::log`;
`VariantRecalibrator`'s row is dominated by the same function. Both call `Math.log` where the
reference has a hardware intrinsic, a dozen times per het per sampler evaluation in the first
case, per Gaussian per variant in the second.

**The optimisation.** In htsjdk-rs (IPNP-BIPN/htsjdk-rs#243, decision 0045 there): a fast phase,
a table, an exact reduction and a short polynomial with an error bound, answering only when
twice the bound rounds to one `f64`, and falling back to the series otherwise.

**Why it cannot move a byte.** `Math.log` is correctly rounded (decision 0006), so it has one
answer per input; the fast phase returns an answer only where it is provably that one, and the
series, unchanged, decides everything else. Checked on 200 million inputs (99.94% decided fast,
none differing) and on two million more on every CI run.

**The suites that prove it.** The `jmath` suite and its hard-to-round corpus in htsjdk-rs; here,
every conformance suite and covering array on the bump's PR, `model-segments` and
`VariantRecalibrator`'s array among them, all green; locally, the four tools' baseline rows give
identical outputs before and after.

**Result, natively (M-series, median of three, identical outputs).**

| Tool | before | after | |
|---|---:|---:|---:|
| `ModelSegments` | 5.19 s | 0.46 s | 11.3x |
| `VariantRecalibrator` | 0.358 s | 0.114 s | 3.1x |
| `HaplotypeBasedVariantRecaller` | 68 ms | 59 ms | 1.15x |
| `FlowPairHMMAlignReadsToHaplotypes` | 87 ms | 87 ms | 1.0x |

**Result on real x86-64** (`Speed` run 37923376066 on this change, `tools/speed/current.json`,
now STATUS.md's cost column; the baseline stays in `baseline.json`):

| Tool | steady, before | steady, after | cold, before | cold, after | port, before | port, after |
|---|---:|---:|---:|---:|---:|---:|
| `ModelSegments` | 18.39 | **1.50** | 4.83 | **0.38** | 9.67 s | 1.26 s |
| `VariantRecalibrator` | 4.63 | **0.75** | 0.37 | 0.06 | 1.00 s | 0.15 s |
| `HaplotypeBasedVariantRecaller` | 1.30 | 1.05 | 0.07 | 0.06 | 158 ms | 125 ms |
| `FlowPairHMMAlignReadsToHaplotypes` | 1.05 | 0.98 | 0.09 | 0.09 | 120 ms | 183 ms |

`ModelSegments` is now cheaper than the reference cold and within a factor of 1.5 of it warm;
`VariantRecalibrator` is cheaper than a warm JVM. The last two rows moved inside the noise, in
both directions, which is what an unchanged path does. The other 131 tools moved a median 7% cold
and 10% steady, at most 45%, the same spread two runs of one commit showed in the baseline.

## Scale N. `VcfToIntervalList` held the VCF

**Ratio.** At scale (gnomAD chr20 sites, 61 MB `.vcf.gz`, linux/arm64, 6 GiB container): GATK
7.4 s CPU, 0.51 GiB peak RSS; the port killed at the limit (9.4x).

**Where the memory goes.** `vcf_text` inflated the whole file into one string, `read_vcf` decoded
every record of it, and the intervals, the merged intervals and the output text were each a
whole list on top.

**The optimisation.** The header is read alone and the records after it one at a time
(`vcf_stream`, #1500). Each is named, merged and written as it comes, to a file beside `OUTPUT`
that is renamed onto it at the end. The merger is `IntervalMergerIterator` a step at a time
(`Merger`), which the whole-list `merge_intervals` now drains too. It keeps the group's least
start and names rather than the group, and writes each line without a string of its own.

**Why it cannot move a byte.** The text is the one `vcf_text` reads. A file that starts with the
gzip magic is inflated by `BgzfReader`, which is what `decompress_all` runs serially and refuses
the same files with. Both are decoded lossily, line by line here, which is the same text because
a newline is never part of a multi-byte sequence. The records are `read_vcf`'s, in its order. On
any refusal (a record, the inflater, no contig lines, an output that cannot be created) nothing is
written, the scratch file is removed, and the whole-text path runs and refuses as before.

**The suites that prove it.** `vcf-to-interval-list` and its golden, plus a new test for the one
case the golden leaves out: a group overlapping backwards, where `CONCAT_ALL` takes the least
start and `USE_FIRST` the first member's. Natively, old against new, the same list, exit and
message were checked:
- on gnomAD chr20 by default, with `INCLUDE_FILTERED` and with `USE_FIRST`;
- on small files: an unsorted file with touching groups and repeated names, CRLF lines, a
  malformed record in the middle, no contig lines, no `#CHROM` line, no records, a Latin-1 byte
  in plain text and in BGZF;
- on BGZF with trailing garbage, truncated BGZF, and plain gzip;
- on an output in a missing directory.

**Result.** macOS native, gnomAD chr20 (4.5 million intervals):

| | main | this | GATK |
|---|---:|---:|---:|
| RSS | 6.22 GB | 16 MB | 0.51 GiB |
| CPU | 4.85 s | 3.87 s | 7.4 s |

The RSS is 0.03x GATK's and the CPU 0.52x. Of what is left, 80% is `htsjdk-vcf`'s
`decode_line` (`parse_alleles` growing its `Vec`, the INFO split) and the single-threaded
inflater, both of which belong to the decoder (htsjdk-io).

## Scale O. The Picard VCF writers held the VCF (`MakeSitesOnlyVcf`, `MergeVcfs`, `SplitVcfs`, `UpdateVcfSequenceDictionary`)

**Ratio.** At scale (linux/arm64, 6 GiB container), each killed at the limit:

| Tool | input | GATK CPU | GATK RSS |
|---|---|---:|---:|
| `MakeSitesOnlyVcf` | dbsnp138 chr20 | 8.2 s | 0.50 GiB |
| `MergeVcfs` | gnomAD chr20 + dbsnp138 chr20 | 17.3 s | 0.51 GiB |
| `SplitVcfs` (`STRICT false`) | gnomAD chr20 | 13.0 s | 0.55 GiB |
| `UpdateVcfSequenceDictionary` | gnomAD chr20, `chr20.dict` | 18.3 s | 0.63 GiB |

**Where the memory goes.** Each input was inflated whole into one string, decoded whole, then
written into a second whole string. `write_vcf` (the Picard writer) compressed that string in
memory before the file was written.

**The optimisation.** A Picard writer a record at a time (`picard::VcfStream`). It is #1500's
`VariantStreamWriter` given `CREATE_INDEX` and the `.idx` dictionary by its caller. As with
`SelectVariants`, a `.vcf.gz` is tabix indexed as it is written. The inputs are read with
`vcf_lines` (Scale N).
- `MakeSitesOnlyVcf` subsets each record's genotypes and writes it.
- `MergeVcfs`' library is split into the steps the whole run takes:
  - `HeaderChecks` checks each input's header in input order.
  - `finish` builds the merged header and the comparator.
  - `merge_sorted` is the same Java `PriorityQueue` of iterators, over the peeked heads of
    record streams.
  The whole-list `merge` now runs these same steps.
- `SplitVcfs` writes each record to its file as it is read. Some cases are left to the
  whole-text run: a `.vcf.gz` SNP file indexed beside a text INDEL file (they share one tabix
  creator), and two outputs of one name. The dictionary is read the way the whole run reads it,
  and `vcf_text` no longer runs first.
- `UpdateVcfSequenceDictionary` writes each record under the updated header (`updated_header`,
  which the whole run's `update` now builds too). It runs once the dictionary is read; the
  refusals of a missing dictionary come first, as before.

**Why it cannot move a byte.**
- The header and records are encoded as `write_vcf` encodes them, through the writer
  `write_variant_file` uses.
- The refusals the Picard writer makes before writing a byte are made before anything is written:
  a Tribble index with no dictionary, and a tabix the records' order refuses.
- `merge_sorted` checks each stream's order against the record it last gave, which is the
  whole-list check over each pair before the merge. Its queue makes the same comparisons in the
  same order.
- An input's header is its text up to `#CHROM`, which holds every line the contig checks, the
  dictionary and `with_contigs` read.
- On any refusal (a record, an input's header, the order, the writer) nothing is written and the
  whole-text path runs, refusing in its own words and order.

**The suites that prove it.** `make-sites-only-vcf`, `merge-vcfs` and their covering arrays.
Natively, old against new, the outputs and their `.tbi` are byte-identical, and every `.idx` too
except its timestamp. On the files above, and on small files with samples:
- `MakeSitesOnlyVcf`:
  - to `.vcf` and `.vcf.gz`, with and without `CREATE_INDEX`, with `SAMPLE` (repeated, and one
    the file lacks);
  - an unsorted input (tabix refusal), a malformed record, no contig lines, an undeclared FORMAT
    key, and truncated BGZF and plain gzip input.
- `MergeVcfs`:
  - two and three inputs in both orders, a `.list`, `COMMENT`, `CREATE_INDEX false`;
  - `SEQUENCE_DICTIONARY` with and without contig lines in the inputs;
  - unsorted, malformed, different samples, incompatible contigs, no dictionary, a missing input.
- `SplitVcfs`, over every pair of `.vcf` and `.vcf.gz` outputs including the shared index:
  - with `STRICT` and without, `CREATE_INDEX false`, and `SEQUENCE_DICTIONARY`;
  - a file of mixed, MNP and symbolic records, a clean one, an unsorted one, a malformed record,
    no contig lines, two outputs of one name, and an output name with no VCF extension.
- `UpdateVcfSequenceDictionary`, to `.vcf` and `.vcf.gz`, with and without `CREATE_INDEX`, with
  the dictionary from a `.dict`, a VCF, a VCF with no contig lines, a FASTA and a missing file,
  over a good input, a BGZF one, a malformed record, no contig lines and truncated BGZF.

**Result.** macOS native:

| Tool | main RSS | RSS | main CPU | CPU | RSS vs GATK | CPU vs GATK |
|---|---:|---:|---:|---:|---:|---:|
| `MakeSitesOnlyVcf` | 5.66 GB | 21 MB | 10.0 s | 6.3 s | 0.04x | 0.76x |
| `MergeVcfs` | 13.5 GB | 28 MB | 74.0 s | 52.4 s | 0.05x | 3.0x |
| `SplitVcfs` | 7.82 GB | 22 MB | 21.0 s | 7.7 s | 0.04x | 0.59x |
| `UpdateVcfSequenceDictionary` | 5.25 GB | 20 MB | 6.7 s | 6.8 s | 0.03x | 0.37x |

`MergeVcfs`' CPU is 70% `VcfEncoder::encode_into` checking each INFO, FORMAT and FILTER key
against the header. htsjdk-vcf does that with a linear scan of every header line, and gnomAD's
header has hundreds. Sent to htsjdk-io as a per-encoder set of ids, which answers the same
booleans; it is every VCF writer's cost, not this tool's.

A text `.vcf` with `CREATE_INDEX` read its whole output back for the `.idx` (650 MB for
`MakeSitesOnlyVcf` on gnomAD to `.vcf`); Scale P streams it.

## Scale P. The last Picard VCF writers (`VcfFormatConverter`, `FixVcfHeader`) and a text `.idx`

**Ratio.** At scale (gnomAD chr20 sites, linux/arm64, 6 GiB container), each killed at the limit:

| Tool | output | GATK CPU | GATK RSS |
|---|---|---:|---:|
| `VcfFormatConverter` | `.vcf` (`.idx`) | 24.2 s | 1.27 GiB |
| `FixVcfHeader` | `.vcf.gz` (`.tbi`) | 17.8 s | 0.98 GiB |

**Where the memory goes.** Each tool read, decoded and wrote the whole file, as in Scale O. A text
output with an index had one more copy: Scale O's writer read the whole output back for the
`.idx`.

**The optimisation.**
- `on_the_fly_index` is now `OnTheFlyIndex` fed the written text line by line. The whole-text
  path feeds it the text once; `VariantStreamWriter` feeds it each line as it writes. A text
  output's `.idx` no longer reads the file back.
- `VcfFormatConverter` writes each record as it is read, once the index check passes. A missing
  index is left to the whole run.
- `FixVcfHeader` reads the file twice, as `doWork` does. The first pass finds the lines the records
  use and the header lacks (`Invented`, which `fix` now runs too); the second writes the records
  under the fixed header. Every record is decoded in the first pass, whatever
  `CHECK_FIRST_N_RECORDS` says, as the whole read decoded them. `HEADER` is left to the whole run.

**Why it cannot move a byte.**
- `OnTheFlyIndex` sees the same lines at the same offsets and is finished at the same length, with
  the same dictionary, path, size and timestamp sources.
- The records and the header are Scale O's: `write_vcf`'s encoding through `write_variant_file`'s
  writer.
- On any refusal nothing is written and the whole run refuses as before.

**The suites that prove it.** `vcf-format-converter`, `fix-vcf-header` and every Scale O suite.
Natively, old against new, 117 runs give the same output, index and message (an `.idx` differs
only in its timestamp):
- `FixVcfHeader` with and without `CHECK_FIRST_N_RECORDS` and `HEADER`;
- `VcfFormatConverter` with `REQUIRE_INDEX` and `CREATE_INDEX` on and off;
- to `.vcf`, `.vcf.gz` and a name with no VCF extension;
- over an indexed input, an unindexed one, a BGZF one, undeclared INFO, FORMAT and FILTER keys,
  a malformed record, no contig lines, an unsorted file and truncated BGZF.
The 162 Scale O runs were repeated with this build and still match.

**Result.** macOS native:

| Tool | main RSS | RSS | main CPU | CPU | RSS vs GATK | CPU vs GATK |
|---|---:|---:|---:|---:|---:|---:|
| `VcfFormatConverter` | 5.25 GB | 34 MB | 7.0 s | 7.5 s | 0.03x | 0.31x |
| `FixVcfHeader` | 9.49 GB | 23 MB | 22.7 s | 12.0 s | 0.02x | 0.67x |
| `MakeSitesOnlyVcf` to `.vcf` (gnomAD) | 8.17 GB | 33 MB | 8.0 s | 7.6 s | | |

`FixVcfHeader`'s CPU is two decodes of the file, as GATK's is; what is left is the decoder's
(htsjdk-io).

## Scale Q. `LeftAlignAndTrimVariants` held the VCF

**Ratio.** At scale (gnomAD chr20 sites to `.vcf.gz`, linux/arm64, 6 GiB container), GATK takes
15.2 s CPU and 0.86 GiB; the port was killed at the limit (4.5 GiB, as far as it got).

**Where the memory goes.** The run read the whole file into one string, decoded every record, kept
every written record in a `Vec`, then encoded the whole output as one more string: four copies of
the file. Natively the peak was 9.30 GB.

**Where the time goes.** A `sample` profile of the streamed run: most of it in the engine bridge,
building per record a JEXL context that no expression reads (a hash map of every INFO key), and
rendering every attribute twice to bring it back. Under that, the decoder and the encoder.

**The optimisation.**
- A run with no `-L` streams: each record is read, aligned and written before the next is read,
  through `VariantStreamWriter`. A run with `-L` (whose traversal is the index's), GenomicsDB, or
  any refusal is left to the whole run, which refuses in its own words and order.
- The header and the per-record body (`LeftAligner::align`) are shared by both paths.
- The bridge is asked for no JEXL context (`to_engine_for(.., false)`).
- A record with no genotypes that the alignment leaves where it is (not split, and too long, not
  an indel, or with an empty window) is written as read, with no round trip through the engine.

**Why it cannot move a byte.**
- The streamed records are the whole run's, in its order (no `-L` means the file's order), aligned
  by the same body with the same bound from the record before, filtered by the same
  `variant_output_filter` and written by the same encoder and writer.
- The JEXL context is only read by `SelectVariants`' expressions.
- For a record with no genotypes and an unchanged place, `from_engine(to_engine(record))` is the
  record: the alleles come back by identity, every INFO value finds its own original (the decoder
  keeps one value per key), and the id, QUAL and filters are copied. The decision is the slow
  path's own: `largest_indel_length`, `is_indel` and the window on the same alleles.

**The suites that prove it.** `left-align-and-trim-variants`, `left-align-and-trim` and
`select-variants-output`. Natively, old against new, 243 runs give the same output, index and
message: nine inputs (plain, BGZF with an index, a malformed record, a contig the reference lacks,
an undeclared INFO key, unsorted, no records, missing, two samples with genotypes) by nine option
sets (none, `--split-multi-allelics`, with `--keep-original-ac`, `--sites-only-vcf-output`,
`--dont-trim-alleles`, `--max-indel-length 5`, `-L`, `-XL`, `--max-leading-bases 3`) by three
outputs (`.vcf`, `.vcf.gz`, no VCF extension). gnomAD chr20 gives the same `.vcf.gz` and `.tbi`.

**Result.** macOS native, gnomAD chr20:

| | RSS | CPU | RSS vs GATK | CPU vs GATK |
|---|---:|---:|---:|---:|
| main | 9.30 GB | 17.6 s | | |
| streamed | 152 MB | 17.3 s | | |
| no JEXL, unchanged records as read | 152 MB | 10.0 s | 0.16x | 0.66x |

The 152 MB is mostly chr20's bases, read once for the alignment, as GATK's reference cache holds
them. What is left of the CPU is the decoder and the encoder (htsjdk-io), and the BGZF codecs.

## Scale R. `ValidateVariants` held the VCF

**Ratio.** At scale (gnomAD chr20 sites against chr20, linux/arm64, 6 GiB container), GATK takes
9.1 s CPU and 1.19 GiB; the port took 5.36 GiB.

**Where the memory goes.** The run read the whole file into one string and decoded every record
before checking the first. Natively the peak was 4.94 GB.

**Where the time goes.** Once streamed, a `sample` profile put the system time (1.5 s of 7.4 s) in
`pread`: one read of the FASTA per record, for the REF check.

**The optimisation.**
- A run with no `-L` streams: each record is read and checked before the next is read. Nothing
  is written, so any refusal, the reader's included, goes to the whole run, which refuses in its
  own words and order. A run with `-L` or GenomicsDB is left to it as well.
- Both paths share the checks (`validate_variant_records`), handed the records where the whole run
  traversed them, after the arguments and the reference are checked.
- `ReferenceFileSource::query_windowed` reads a mebibase ahead on a miss and answers the queries
  inside that window from it.

**Why it cannot move a byte.**
- The streamed records are the whole run's, in its order, through the same checks; the run writes
  no file, and its only output is the refusal, which the streamed path never makes.
- The reference's flattening (upper case, IUPAC to `N`) is base by base, so a slice of a longer
  query is the shorter query; a query `query` refuses, or a read ahead that fails, goes to `query`
  itself. `reference_query`'s golden test now asks every query again through the window, forwards
  and then backwards.

**The suites that prove it.** `validate-variants`, and `reference_query`. Natively, old against
new, 144 runs give the same message and status: twelve inputs (plain, BGZF with an index, a
malformed record, a contig the reference lacks, an undeclared INFO key, unsorted, no records,
missing, two samples with genotypes, a GVCF, an overlapping GVCF, a wrong REF) by twelve option
sets (no reference, `-R`, `--warn-on-errors`, `--validation-type-to-exclude ALL` and `REF`,
`--validate-GVCF` alone, with `-L` and with `--fail-gvcf-on-overlap`, `-L`, `-XL`,
`--do-not-validate-filtered-records`).

**Result.** macOS native, gnomAD chr20:

| | RSS | CPU | RSS vs GATK | CPU vs GATK |
|---|---:|---:|---:|---:|
| main | 4.94 GB | 7.8 s | | |
| streamed | 17 MB | 7.4 s | | |
| reference read ahead | 22 MB | 5.4 s | 0.02x | 0.59x |

What is left of the CPU is the decoder (htsjdk-io).

## Scale S. `SortVcf` held every input

**Ratio.** At scale (gnomAD chr20 sites to `.vcf.gz`, linux/arm64, 6 GiB container), GATK takes
33.8 s CPU and 3.57 GiB (its `SortingCollection` spills to disk); the port was killed at the
limit.

**Where the memory goes.** Each input was read into one string and decoded whole, every record
held for one in-memory sort, and the output encoded as one more string. Natively the peak was
6.64 GB for one gnomAD chromosome and 11.5 GB for two inputs.

**The optimisation.** An external sort over the encoded lines.
- The headers are read first, all of them, and checked as before (`Checks`, which `sort` now runs
  too); the merged header is written at once.
- Each record is read, encoded under that header and kept with its key (contig index, start).
  Every 128 MiB of text the run is sorted stably and written beside the output; at the end the
  runs are merged, the lowest run winning a tie, and each line written as it comes out. A run
  that never fills one is sorted and written from memory.
- Any refusal (a reader, a dictionary, the samples, a contig the dictionary lacks, the encoder,
  the writer) leaves nothing written and the whole run refuses in its own words and order.

**Why it cannot move a byte.**
- `write_vcf` is the header followed by each record's encoding, so encoding first and ordering
  the lines after is the same text.
- The order is the whole run's stable sort: a run holds only records read after every earlier
  run's, each run is sorted stably, and the merge breaks a tie by run, so equal keys keep their
  input order across runs. The key is the comparator's (its contig lookup, `contig_indexes`).
- The output goes through Scale O's writer and index.

**The suites that prove it.** `sort-vcf`, and a unit test merging runs of one line, of a few and
of all of them over repeated keys. Natively, old against new, 204 runs give the same output,
index (but its timestamp) and message: one to three inputs (plain, unsorted, two samples,
unsorted samples, no records, a malformed record, a contig the dictionary lacks, an undeclared
INFO key, missing, no contig lines, truncated BGZF, two inputs either way round, different
samples, one with and one without a dictionary) by `SEQUENCE_DICTIONARY` (none, the reference's,
another) and `CREATE_INDEX false`, by three outputs. gnomAD chr20 gives the same `.vcf.gz` and
`.tbi`; a shuffled copy of it sorted with the original (five runs, every key tied across them)
gives the same `.vcf`.

**Result.** macOS native:

| input | main RSS | RSS | main CPU | CPU | RSS vs GATK | CPU vs GATK |
|---|---:|---:|---:|---:|---:|---:|
| gnomAD chr20 | 6.64 GB | 244 MB | 18.6 s | 8.5 s | 0.06x | 0.25x |
| it and a shuffled copy | 11.5 GB | 315 MB | 16.3 s | 16.1 s | | |

## Not targets yet

`HaplotypeBasedVariantRecaller` (1.30 steady) and `FlowPairHMMAlignReadsToHaplotypes` (1.05) are
at or near one, and the baseline's own noise says a single tool can move by up to half between two
runs on different runners. Neither is "slow" by a margin that survives that. Their rows take 60 to
90 ms natively, too short for a sampling profiler to attribute, so there is also no profile to
choose an optimisation from. They re-enter this list when a larger input, still covered by a
golden, gives them a run long enough to profile and a ratio that clears the noise.

## Out of scope, as S.2 says

Anything that would reorder a reduction, contract an FMA or swap a `jmath` call for the host's,
however tempting a profile makes it look. Target 1 is admissible precisely because it changes how
the correctly rounded value is found and never which value it is.
