# A second reader: noodles reads what htsjdk-rs reads

Issue #1393.

**Verdict: keep it in CI as a gate.** noodles 0.119.0 is a dev-dependency of `gatk-corpus`. The test
`noodles_reads_what_htsjdk_rs_reads` reads every BAM, BAI, BGZF stream, tabix index and VCF that
the 368 goldens carry. It reads each with htsjdk-rs and with noodles, and fails on any disagreement
nobody has judged. It also fails if a judged disagreement stops happening, so a fix cannot leave a
stale entry in the list. It found one real bug, in htsjdk-rs's genotype model, and that bug leads
to a SelectVariants divergence from the oracle. The other three disagreements are noodles being
stricter than the reference on purpose. The gate costs about 8 s of CPU per `cargo test` run. Nothing
reaches the byte path, which follows `docs/when-a-dependency-is-cheaper-than-a-port.md`.

## What it reads

The goldens are the port's output: every port output is byte-equal to its golden, so reading the
files embedded in the goldens (base64 BGZF and BAI, escaped text VCF) means reading what the port
writes.

| Read by both | Count |
| --- | --- |
| BGZF streams, inflated by both and compared byte for byte | 582 (4 refused by both, 13 neither BAM, VCF nor index) |
| BAM records: name, flag, reference, position, MAPQ, CIGAR, mate, TLEN, bases, qualities, every tag | 88,224 |
| BAI: bins and chunks, linear index, pseudo-bin metadata, unplaced count | 228 (1 whose BAM the golden does not carry) |
| BAI queries through noodles, each record looked for over its own span | 50,599 measured once; 4,598 in the gate |
| Tabix indices, each record of the VCF they index queried through noodles | 13, 112 queries |
| VCF files (text and BGZF) / records | 527 (110 refused by both) / 2,599 |

Each VCF record is compared field by field: CHROM, POS, ID, alleles, QUAL as `f32` bits, FILTER,
INFO keys and values, sample names, GT allele indices and phasing, GQ, DP, AD, PL, FT and the other
FORMAT fields. A value is compared after flattening it on commas, and floats are compared as `f32`
bits. That flattening removes a representation difference: htsjdk keeps a FORMAT extended value as
one string, and noodles splits it by `Number`.

The BAI pairing is part of the test. A dump puts its `.bai` before its BAM as often as after it,
and several BAMs in a golden can share a dictionary and record counts. So the BAMs a BAI could
describe are chosen by the per-reference aligned and unaligned counts and the unplaced count. The
one the index answers for is then picked by querying each reference's first and last mapped record.

## The disagreements, judged

| Where | What | Whose | Count |
| --- | --- | --- | --- |
| `check_terminator_block.txt.gz` | noodles refuses a block: `block data checksum mismatch` | **noodles is stricter than the reference.** The dump corrupts a payload on purpose. htsjdk's `BlockGunzipper` checks CRCs only when asked, and htsjdk-bgzf does the same on purpose | 1 |
| `apply_vqsr_site_filtering.txt.gz` | `InvalidValue("VQSLOD", InvalidFloat)` | **An input malformed on purpose** (`VQSLOD=nonsense`). htsjdk decodes INFO lazily, and the golden holds the reference's answer | 1 |
| `mutect2_gvcf.txt.gz` | `InvalidValue("AS_MQ", EmptyValue)` | **A gap in the specification, which GATK fills itself.** It writes `AS_MQ=57.52,`, an empty value for `<NON_REF>`. The specification has no empty element, so noodles refuses it; the port writes the reference's bytes | 1 |
| 15 VCF records | FT: htsjdk-rs `[""]`, noodles `PASS` | **htsjdk-rs's bug.** A genotype `FT=PASS` reads as `Some("")`, where htsjdk's `GenotypeBuilder.filters(emptyList)` leaves `null`. `is_filtered()` agrees with htsjdk; the field does not | 15 |

The test does not judge these with a tolerance. Each one is an entry in `KNOWN` (where it is, what it
says, the verdict above). An entry that matches nothing fails the test, so fixing the FT bug means
deleting its entry.

### The FT bug reaches a tool's output

The oracle settled whether the field matters. On a three-record VCF (`FT=PASS`, `FT=lowQ`, no FT),
the pinned GATK 4.7.0.0 gives these SelectVariants `--select-genotype` results:

| Expression | GATK 4.7.0.0 keeps | The port (x567 build) |
| --- | --- | --- |
| `FT == "PASS"` | `pass`, `noft` | `Invalid JEXL expression detected for select-genotype-0` |
| `isPassFT == "1"` | `pass`, `noft` | the same error |

On a VCF where every genotype has FT, GATK keeps `pass` for both expressions. The port gives an
empty output for `FT == "PASS"` and the error for `isPassFT`. The genotype JEXL context the port
builds has neither `FT` nor `isPassFT`, and its FT would be `""` where the reference's is `PASS`.
That fix belongs to SelectVariants and htsjdk-rs, not to this test. It is left as a follow-up: new
rows in the SelectVariants filter dump, a golden frozen in CI, then deleting the `KNOWN` entry.

## What it costs

`/usr/bin/time -l`, release build, test binary already built, on an arm64 Mac:

| Version | Real | User |
| --- | --- | --- |
| First: every mapped record queried through its BAI, candidates found by full queries | — | 151 s |
| Candidates picked by the first and last record of each reference, then every record queried | 101 s | 78.8 s |
| **Kept:** each reference's ends plus every `n/64`-th record, so about 64 queries per index | 16.0–17.1 s | 7.2–8.5 s |

Timing each phase showed the BAI queries took 81 of the 82 s: 1.6 ms a query, because noodles
decodes each chunk's records to answer one. Inflating, the BAM and VCF comparisons, and tabix take
under 1 s together. The full query of all 50,599 mapped records returned every one of them. The
sample keeps every structural comparison (bins, linear index, metadata) and queries each reference's
ends, where an off-by-one in a bin or linear-index entry shows first. A comment on `QUERIES_PER_INDEX`
says when to query every record again: when the index writer or either reader changes.

## Pinning

`noodles = { version = "0.119.0", features = [...] }`, plus the 23 packages it adds to `Cargo.lock`.
Its rust-version is 1.91, below the pinned 1.97.1. Its README calls the API experimental, so an
upgrade can break the test's code. That is a test-only change, and the `KNOWN` list says what the
upgraded reader must still refuse.
