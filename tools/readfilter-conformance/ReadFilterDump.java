/*
 * The read filters' decisions, taken by the reference itself.
 *
 * Each filter is instantiated from ReadFilterLibrary and applied to a fixed corpus through
 * SAMRecordToGATKReadAdapter, which is the adapter every GATK tool uses. That last part is the
 * point: a filter's predicate is written against GATKRead, and three of GATKRead's accessors are
 * not the flag test they look like (isUnmapped is three criteria, isFirstOfPair requires pairing).
 * Comparing against the flags would compare against a reimplementation of the adapter rather than
 * against the reference.
 *
 * Output, one line per filter:
 *
 *     filter\t<Name>\t<one character per record: 1 kept, 0 filtered out>
 *
 * plus the corpus itself, so the port judges exactly the records the reference judged:
 *
 *     record\t<index>\tname|flags|refIdx|start|mapq|cigar|mateRef|mateStart|isize|bases|quals
 *     tag\t<index>\t<NAME>\t<value>
 *
 * Fields, not a SAM line. The corpus deliberately contains a record whose flags say mapped while
 * its reference is absent, because that is one of the three criteria of GATKRead.isUnmapped, and
 * htsjdk's *reader* rejects exactly that record ("RNAME is not specified but flags indicate
 * mapped"; htsjdk-rs decision 0015 is about the writer emitting what the reader refuses). Routing
 * the corpus through SAM text would therefore drop the one case the filter most needs.
 *
 * Usage: ReadFilterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class ReadFilterDump {

    static final int CHR1 = 2_000;
    static final int CHR2 = 1_000;

    /**
     * The filters ported so far, by the name GATK exposes on the command line.
     *
     * Kept in step by hand with gatk-readfilter's PORTED list; the Rust test fails if the golden
     * carries a filter it cannot evaluate, so a name added on one side and not the other is loud.
     */
    static Map<String, ReadFilter> filters() {
        final Map<String, ReadFilter> map = new LinkedHashMap<>();
        map.put("AllowAllReadsReadFilter", ReadFilterLibrary.ALLOW_ALL_READS);
        map.put("FirstOfPairReadFilter", ReadFilterLibrary.FIRST_OF_PAIR);
        map.put("HasReadGroupReadFilter", ReadFilterLibrary.HAS_READ_GROUP);
        map.put("MappedReadFilter", ReadFilterLibrary.MAPPED);
        map.put("MappingQualityAvailableReadFilter", ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        map.put("MappingQualityNotZeroReadFilter", ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        map.put("MatchingBasesAndQualsReadFilter", ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        map.put("MateDifferentStrandReadFilter", ReadFilterLibrary.MATE_DIFFERENT_STRAND);
        map.put("MateOnSameContigOrNoMappedMateReadFilter",
                ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        map.put("NonZeroFragmentLengthReadFilter", ReadFilterLibrary.NONZERO_FRAGMENT_LENGTH_READ_FILTER);
        map.put("NotDuplicateReadFilter", ReadFilterLibrary.NOT_DUPLICATE);
        map.put("NotProperlyPairedReadFilter", ReadFilterLibrary.NOT_PROPERLY_PAIRED);
        map.put("NotSecondaryAlignmentReadFilter", ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        map.put("NotSupplementaryAlignmentReadFilter", ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);
        map.put("PairedReadFilter", ReadFilterLibrary.PAIRED);
        map.put("PassesVendorQualityCheckReadFilter", ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        map.put("SecondOfPairReadFilter", ReadFilterLibrary.SECOND_OF_PAIR);
        map.put("ProperlyPairedReadFilter", ReadFilterLibrary.PROPERLY_PAIRED);
        map.put("CigarContainsNoNOperator", new ReadFilterLibrary.CigarContainsNoNOperator());
        map.put("GoodCigarReadFilter", ReadFilterLibrary.GOOD_CIGAR);
        map.put("NonZeroReferenceLengthAlignmentReadFilter",
                ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        map.put("PrimaryLineReadFilter", ReadFilterLibrary.PRIMARY_LINE);
        map.put("ReadLengthEqualsCigarLengthReadFilter",
                ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH);
        map.put("SeqIsStoredReadFilter", ReadFilterLibrary.SEQ_IS_STORED);
        map.put("ValidAlignmentStartReadFilter", ReadFilterLibrary.VALID_ALIGNMENT_START);
        map.put("ValidAlignmentEndReadFilter", ReadFilterLibrary.VALID_ALIGNMENT_END);
        map.put("MateUnmappedAndUnmappedReadFilter",
                ReadFilterLibrary.MATE_UNMAPPED_AND_UNMAPPED_READ_FILTER);
        map.put("NonChimericOriginalAlignmentReadFilter",
                ReadFilterLibrary.NON_CHIMERIC_ORIGINAL_ALIGNMENT_READ_FILTER);
        return map;
    }

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = header();
        final List<SAMRecord> corpus = corpus(header);

        System.out.println("# ReadFilterDump: the reference's own decision per filter, per record");

        for (int i = 0; i < corpus.size(); i++) {
            final SAMRecord record = corpus.get(i);
            System.out.printf("record\t%d\t%s%n", i, fields(record));
            // Tags get their own rows rather than a delimited column: an OA value ends with a
            // semicolon of its own, so any in-line separator collides with the data it carries.
            for (final SAMRecord.SAMTagAndValue tag : record.getAttributes()) {
                System.out.printf("tag\t%d\t%s\t%s%n", i, tag.tag, String.valueOf(tag.value));
            }
        }

        for (final Map.Entry<String, ReadFilter> entry : filters().entrySet()) {
            final ReadFilter filter = entry.getValue();
            filter.setHeader(header);
            final StringBuilder decisions = new StringBuilder();
            for (final SAMRecord record : corpus) {
                final GATKRead read = new SAMRecordToGATKReadAdapter(record);
                boolean kept;
                try {
                    kept = filter.test(read);
                } catch (final Exception e) {
                    // A filter that throws on a record is a third outcome, not a silent false:
                    // mateIsUnmapped asserts pairing, for instance. The port has to match that.
                    decisions.append('E');
                    continue;
                }
                decisions.append(kept ? '1' : '0');
            }
            System.out.printf("filter\t%s\t%s%n", entry.getKey(), decisions);
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader h = new SAMFileHeader();
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", CHR1));
        dict.addSequence(new SAMSequenceRecord("chr2", CHR2));
        h.setSequenceDictionary(dict);
        h.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("sample1");
        rg.setLibrary("lib1");
        rg.setPlatform("ILLUMINA");
        h.addReadGroup(rg);
        return h;
    }

    /**
     * A corpus built to make the filters disagree with each other.
     *
     * Every record here exists to separate two predicates that a careless port would collapse: the
     * unmapped-by-flag from the unmapped-by-missing-reference, the 0x40-without-0x1 from the real
     * first-of-pair, the absent qualities from the mismatched ones.
     */
    static List<SAMRecord> corpus(final SAMFileHeader header) {
        final List<SAMRecord> out = new ArrayList<>();

        out.add(read(header, "plain_mapped", 0, 0, 100, 60, "10M", 0, 200, 100, true));

        SAMRecord r = read(header, "flag_unmapped", 0x4, 0, 100, 0, "*", 0, 0, 0, true);
        out.add(r);

        // The three criteria of GATKRead.isUnmapped, one at a time: this one has the flag clear.
        r = read(header, "no_reference", 0, -1, 100, 30, "10M", -1, 0, 0, true);
        out.add(r);

        r = read(header, "zero_start", 0, 0, 0, 30, "10M", 0, 0, 0, true);
        out.add(r);

        // 0x40 without 0x1: first-of-pair by flag, not by GATKRead.
        r = read(header, "second_flag_unpaired", 0x40, 0, 300, 60, "10M", 0, 0, 0, true);
        out.add(r);

        out.add(read(header, "first_of_pair", 0x1 | 0x40 | 0x2, 0, 400, 60, "10M", 0, 500, 150, true));
        out.add(read(header, "second_of_pair", 0x1 | 0x80 | 0x2, 0, 500, 60, "10M", 0, 400, -150, true));

        // Mate on another contig, and mate unmapped, which the same filter treats differently.
        out.add(read(header, "mate_other_contig", 0x1 | 0x40, 0, 600, 60, "10M", 1, 100, 0, true));
        out.add(read(header, "mate_unmapped", 0x1 | 0x40 | 0x8, 0, 700, 60, "10M", 0, 700, 0, true));

        // Strands: same and different, both ends mapped.
        out.add(read(header, "mate_same_strand", 0x1 | 0x40, 0, 800, 60, "10M", 0, 900, 100, true));
        out.add(read(header, "mate_diff_strand", 0x1 | 0x40 | 0x20, 0, 900, 60, "10M", 0, 1000, 100, true));

        out.add(read(header, "duplicate", 0x400, 0, 1000, 60, "10M", 0, 0, 0, true));
        out.add(read(header, "secondary", 0x100, 0, 1100, 60, "10M", 0, 0, 0, true));
        out.add(read(header, "supplementary", 0x800, 0, 1200, 60, "10M", 0, 0, 0, true));
        out.add(read(header, "vendor_fail", 0x200, 0, 1300, 60, "10M", 0, 0, 0, true));
        out.add(read(header, "mapq_zero", 0, 0, 1400, 0, "10M", 0, 0, 0, true));
        out.add(read(header, "mapq_unavailable", 0, 0, 1500, 255, "10M", 0, 0, 0, true));

        // No read group at all.
        out.add(read(header, "no_read_group", 0, 0, 1600, 60, "10M", 0, 0, 0, false));

        // Cigars, one per rule the good-cigar test applies. Each is legal SAM text, so they reach
        // the filter rather than being rejected by the record's own setter.
        out.add(read(header, "cigar_with_n", 0, 0, 1750, 60, "4M2N4M", 0, 0, 0, true));
        out.add(read(header, "cigar_consecutive_indels", 0, 0, 1760, 60, "3M1I1D5M", 0, 0, 0, true));
        out.add(read(header, "cigar_leading_deletion", 0, 0, 1770, 60, "2S1D8M", 0, 0, 0, true));
        out.add(read(header, "cigar_trailing_deletion", 0, 0, 1780, 60, "8M1D2S", 0, 0, 0, true));
        out.add(read(header, "cigar_all_insertion", 0, 0, 1790, 60, "10I", 0, 0, 0, true));
        out.add(read(header, "cigar_soft_inside_hard", 0, 0, 1800, 60, "2H2S6M2S2H", 0, 0, 0, true));

        // Read length disagreeing with the cigar's read length, on a *mapped* read: the filter
        // lets every unmapped read through, so an unmapped one would prove nothing.
        r = read(header, "cigar_length_mismatch", 0, 0, 1810, 60, "5M", 0, 0, 0, true);
        out.add(r);

        // OA and XM: same contig, then different, then only one of the two present.
        //
        // XM, not MC. The first corpus used MC, which is not the tag AddOriginalAlignmentTags
        // writes, so the filter saw no tags and passed everything: the golden agreed with a port
        // that had made the same wrong assumption, and neither was tested.
        r = read(header, "oa_same_contig", 0, 0, 1820, 60, "10M", 0, 0, 0, true);
        r.setAttribute("OA", "chr1,100,+,10M,60,0;");
        r.setAttribute("XM", "chr1");
        out.add(r);

        r = read(header, "oa_other_contig", 0, 0, 1830, 60, "10M", 0, 0, 0, true);
        r.setAttribute("OA", "chr2,100,+,10M,60,0;");
        r.setAttribute("XM", "chr1");
        out.add(r);

        r = read(header, "oa_without_mate_contig", 0, 0, 1840, 60, "10M", 0, 0, 0, true);
        r.setAttribute("OA", "chr2,100,+,10M,60,0;");
        out.add(r);

        // Qualities absent: getBaseQualityCount() is 0 while the read has ten bases.
        r = read(header, "no_qualities", 0, 0, 1700, 60, "10M", 0, 0, 0, true);
        r.setBaseQualities(SAMRecord.NULL_QUALS);
        out.add(r);

        return out;
    }

    static SAMRecord read(
            final SAMFileHeader header,
            final String name,
            final int flags,
            final int refIndex,
            final int start,
            final int mapq,
            final String cigar,
            final int mateRef,
            final int mateStart,
            final int insertSize,
            final boolean readGroup) {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName(name);
        r.setFlags(flags);
        r.setReferenceIndex(refIndex);
        r.setAlignmentStart(start);
        r.setMappingQuality(mapq);
        r.setCigarString(cigar);
        r.setMateReferenceIndex(mateRef);
        r.setMateAlignmentStart(mateStart);
        r.setInferredInsertSize(insertSize);
        r.setReadBases("ACGTACGTAC".getBytes());
        r.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        if (readGroup) r.setAttribute("RG", "rg1");
        return r;
    }

    /** The record as fields, in the order the port reads them. */
    static String fields(final SAMRecord r) {
        final byte[] quals = r.getBaseQualities();
        final StringBuilder qualText = new StringBuilder();
        for (final byte q : quals) qualText.append(qualText.length() == 0 ? "" : ",").append(q);
        return String.join("|",
                r.getReadName(),
                String.valueOf(r.getFlags()),
                String.valueOf(r.getReferenceIndex()),
                String.valueOf(r.getAlignmentStart()),
                String.valueOf(r.getMappingQuality()),
                r.getCigarString(),
                String.valueOf(r.getMateReferenceIndex()),
                String.valueOf(r.getMateAlignmentStart()),
                String.valueOf(r.getInferredInsertSize()),
                new String(r.getReadBases()),
                qualText.toString());
    }
}
