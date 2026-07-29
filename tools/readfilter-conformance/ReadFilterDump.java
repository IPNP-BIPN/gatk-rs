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
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.AmbiguousBaseReadFilter;
import org.broadinstitute.hellbender.engine.filters.ExcessiveEndClippedReadFilter;
import org.broadinstitute.hellbender.engine.filters.LibraryReadFilter;
import org.broadinstitute.hellbender.engine.filters.MetricsReadFilter;
import org.broadinstitute.hellbender.engine.filters.NotOpticalDuplicateReadFilter;
import org.broadinstitute.hellbender.engine.filters.OverclippedReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadGroupBlackListReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadGroupReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadTagValueFilter;
import org.broadinstitute.hellbender.engine.filters.SoftClippedReadFilter;
import org.broadinstitute.hellbender.engine.filters.PlatformReadFilter;
import org.broadinstitute.hellbender.engine.filters.PlatformUnitReadFilter;
import org.broadinstitute.hellbender.engine.filters.SampleReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.filters.FragmentLengthReadFilter;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.MateDistantReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadNameReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadStrandFilter;
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

        // The parameterised filters, instantiated by the reference itself. The label carries the
        // parameters, and the port rebuilds the instance from that label, so there is one list of
        // instantiations rather than two that could drift.
        //
        // Values are chosen to sit *inside* the corpus: a threshold no record can reach makes a
        // filter that always answers the same thing, which is a row of identical characters and no
        // evidence at all.
        map.put("MappingQualityReadFilter(min=10,max=null)", new MappingQualityReadFilter(10));
        map.put("MappingQualityReadFilter(min=30,max=60)", new MappingQualityReadFilter(30, 60));
        map.put("ReadLengthReadFilter(min=1,max=9)", new ReadLengthReadFilter(1, 9));
        map.put("ReadLengthReadFilter(min=10,max=100)", new ReadLengthReadFilter(10, 100));
        map.put("FragmentLengthReadFilter(min=0,max=120)", fragmentLength(0, 120));
        map.put("FragmentLengthReadFilter(min=150,max=1000000)", fragmentLength(150, 1000000));
        map.put("MateDistantReadFilter(threshold=50)", new MateDistantReadFilter(50));
        map.put("MateDistantReadFilter(threshold=1000)", new MateDistantReadFilter(1000));
        map.put("ReadNameReadFilter(names=plain_mapped+duplicate)",
                readName("plain_mapped", "duplicate"));
        map.put("ReadStrandFilter(keepReverse=false)", readStrand(false));
        map.put("ReadStrandFilter(keepReverse=true)", readStrand(true));
        map.put("AmbiguousBaseReadFilter(maxBases=0,maxFraction=0.05)",
                new AmbiguousBaseReadFilter(0));
        map.put("AmbiguousBaseReadFilter(maxBases=null,maxFraction=0.05)",
                new AmbiguousBaseReadFilter());

        // The header-dependent family. Each is given the header by setHeader() in main, exactly as
        // the engine gives it to them.
        map.put("HasReadGroupWithHeader()", ReadFilterLibrary.HAS_READ_GROUP);
        map.put("AlignmentAgreesWithHeaderReadFilter()", new AlignmentAgreesWithHeaderReadFilter());
        map.put("WellformedReadFilter()", new WellformedReadFilter());
        map.put("LibraryReadFilter(keep=lib1)", library("lib1"));
        map.put("SampleReadFilter(keep=sample2)", sample("sample2"));
        map.put("PlatformReadFilter(names=ILLUM)", platform("ILLUM"));
        map.put("PlatformUnitReadFilter(blacklist=unit-rg2)", platformUnit("unit-rg2"));
        map.put("PlatformUnitReadFilter(blacklist=)", platformUnit());

        // The clipping family. Thresholds are chosen to sit between corpus records rather than
        // outside them: OverclippedReadFilter with both instances at minAlignedBases=8 separates
        // only because one of them requires soft clips at both ends, and 5S5M has one.
        map.put("SoftClippedReadFilter(ratio=0.4,leadingTrailingRatio=null)",
                configure(new SoftClippedReadFilter(), "maximumSoftClippedRatio", 0.4d));
        map.put("SoftClippedReadFilter(ratio=0.1,leadingTrailingRatio=null)",
                configure(new SoftClippedReadFilter(), "maximumSoftClippedRatio", 0.1d));
        map.put("SoftClippedReadFilter(ratio=null,leadingTrailingRatio=0.4)",
                configure(new SoftClippedReadFilter(),
                        "maximumLeadingTrailingSoftClippedRatio", 0.4d));
        map.put("SoftClippedReadFilter(ratio=null,leadingTrailingRatio=0.2)",
                configure(new SoftClippedReadFilter(),
                        "maximumLeadingTrailingSoftClippedRatio", 0.2d));
        map.put("OverclippedReadFilter(minAlignedBases=8,dontRequireBothEnds=false)",
                configure(new OverclippedReadFilter(), "minAlignedBases", 8));
        map.put("OverclippedReadFilter(minAlignedBases=8,dontRequireBothEnds=true)",
                configure(new OverclippedReadFilter(), "minAlignedBases", 8,
                        "doNotRequireSoftClipsOnBothEnds", true));
        map.put("ExcessiveEndClippedReadFilter(maxClippedBases=4)",
                configure(new ExcessiveEndClippedReadFilter(), "maxClippedBases", 4));
        map.put("ExcessiveEndClippedReadFilter(maxClippedBases=2)",
                configure(new ExcessiveEndClippedReadFilter(), "maxClippedBases", 2));

        // The tag and read-group family. Three of these can throw rather than decide, which the
        // dump records as 'E': ReadGroupReadFilter dereferences a read group that may be absent,
        // and both tag filters parse a value that may not be a number.
        map.put("NotOpticalDuplicateReadFilter()", new NotOpticalDuplicateReadFilter());
        map.put("ReadGroupReadFilter(keep=rg1)", readGroup("rg1"));
        map.put("ReadGroupReadFilter(keep=rg2)", readGroup("rg2"));
        map.put("MetricsReadFilter(pfReadOnly=true,alignedReadsOnly=true)",
                new MetricsReadFilter(true, true));
        map.put("MetricsReadFilter(pfReadOnly=false,alignedReadsOnly=false)",
                new MetricsReadFilter(false, false));
        map.put("ReadTagValueFilter(tag=TV,op=EQUAL,value=0.0)",
                new ReadTagValueFilter("TV", 0.0f, ReadTagValueFilter.Operator.EQUAL));
        map.put("ReadTagValueFilter(tag=TV,op=NOT_EQUAL,value=0.0)",
                new ReadTagValueFilter("TV", 0.0f, ReadTagValueFilter.Operator.NOT_EQUAL));
        map.put("ReadTagValueFilter(tag=TV,op=LESS,value=0.0)",
                new ReadTagValueFilter("TV", 0.0f, ReadTagValueFilter.Operator.LESS));
        map.put("ReadTagValueFilter(tag=TV,op=GREATER_OR_EQUAL,value=5.0)",
                new ReadTagValueFilter("TV", 5.0f, ReadTagValueFilter.Operator.GREATER_OR_EQUAL));

        // Exact values, not substrings, despite the argument documentation saying <TAG>:<SUBSTRING>.
        map.put("ReadGroupBlackListReadFilter(blacklist=PU:unit-rg1)", blackList("PU:unit-rg1"));
        map.put("ReadGroupBlackListReadFilter(blacklist=PU:unit)", blackList("PU:unit"));
        map.put("ReadGroupBlackListReadFilter(blacklist=ID:rg2+PL:ILLUMINA)",
                blackList("ID:rg2", "PL:ILLUMINA"));
        return map;
    }

    /**
     * Set argument fields the reference declares private, by their declared names.
     *
     * The clipping filters expose their parameters only through @Argument fields and
     * package-private test constructors, neither of which this class can reach. Naming the field
     * is deliberate: if a rename upstream makes the name wrong, this throws rather than leaving the
     * filter on its default, which would produce a full golden describing a filter nobody asked
     * for.
     */
    static <T extends ReadFilter> T configure(final T filter, final Object... nameThenValue) {
        for (int i = 0; i < nameThenValue.length; i += 2) {
            final String name = (String) nameThenValue[i];
            try {
                final java.lang.reflect.Field field = filter.getClass().getDeclaredField(name);
                field.setAccessible(true);
                field.set(filter, nameThenValue[i + 1]);
            } catch (final ReflectiveOperationException e) {
                throw new IllegalStateException(
                        filter.getClass().getSimpleName() + " has no argument field '" + name + "'",
                        e);
            }
        }
        return filter;
    }

    static ReadGroupReadFilter readGroup(final String keep) {
        final ReadGroupReadFilter filter = new ReadGroupReadFilter();
        filter.readGroup = keep;
        return filter;
    }

    static ReadGroupBlackListReadFilter blackList(final String... entries) {
        final ReadGroupBlackListReadFilter filter = new ReadGroupBlackListReadFilter();
        filter.blackList = new ArrayList<>(java.util.Arrays.asList(entries));
        return filter;
    }

    static LibraryReadFilter library(final String... libraries) {
        final LibraryReadFilter filter = new LibraryReadFilter();
        filter.libraryToKeep = new java.util.LinkedHashSet<>(java.util.Arrays.asList(libraries));
        return filter;
    }

    static SampleReadFilter sample(final String... samples) {
        final SampleReadFilter filter = new SampleReadFilter();
        filter.samplesToKeep = new java.util.LinkedHashSet<>(java.util.Arrays.asList(samples));
        return filter;
    }

    static PlatformReadFilter platform(final String... names) {
        final PlatformReadFilter filter = new PlatformReadFilter();
        filter.PLFilterNames = new java.util.LinkedHashSet<>(java.util.Arrays.asList(names));
        return filter;
    }

    static PlatformUnitReadFilter platformUnit(final String... lanes) {
        final PlatformUnitReadFilter filter = new PlatformUnitReadFilter();
        filter.blackListedLanes = new java.util.LinkedHashSet<>(java.util.Arrays.asList(lanes));
        return filter;
    }

    /** The two filters whose parameters have no constructor: set the public fields directly. */
    static FragmentLengthReadFilter fragmentLength(final int min, final int max) {
        final FragmentLengthReadFilter filter = new FragmentLengthReadFilter();
        filter.minFragmentLength = min;
        filter.maxFragmentLength = max;
        return filter;
    }

    static ReadNameReadFilter readName(final String... names) {
        final ReadNameReadFilter filter = new ReadNameReadFilter();
        filter.readNames = new java.util.LinkedHashSet<>(java.util.Arrays.asList(names));
        return filter;
    }

    static ReadStrandFilter readStrand(final boolean keepOnlyReverse) {
        final ReadStrandFilter filter = new ReadStrandFilter();
        filter.keepOnlyReverse = keepOnlyReverse;
        return filter;
    }

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = header();
        final List<SAMRecord> corpus = corpus(header);

        System.out.println("# ReadFilterDump: the reference's own decision per filter, per record");

        // The header travels too: the resolved filters read the library, sample, platform and
        // contig lengths out of it, so a port given a different header would be answering a
        // different question.
        for (final SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
            System.out.printf("sq\t%d\t%s\t%d%n",
                    seq.getSequenceIndex(), seq.getSequenceName(), seq.getSequenceLength());
        }
        for (final SAMReadGroupRecord group : header.getReadGroups()) {
            System.out.printf("rg\t%s\tLB=%s\tSM=%s\tPL=%s\tPU=%s%n",
                    group.getId(), group.getLibrary(), group.getSample(),
                    group.getPlatform(), group.getPlatformUnit());
        }

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
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        rg1.setSample("sample1");
        rg1.setLibrary("lib1");
        rg1.setPlatform("ILLUMINA");
        rg1.setPlatformUnit("unit-rg1");
        h.addReadGroup(rg1);

        // A second group differing in every resolved field, so the library, sample, platform and
        // platform-unit filters each have something to separate.
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        rg2.setSample("sample2");
        rg2.setLibrary("lib2");
        rg2.setPlatform("PACBIO");
        rg2.setPlatformUnit("unit-rg2");
        h.addReadGroup(rg2);
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

        // A short read and a reverse-strand read. Without them ReadLengthReadFilter and
        // ReadStrandFilter answer the same thing for every record, which is a row of identical
        // characters and no evidence: a filter that never changes its mind tests nothing.
        r = read(header, "short_read", 0, 0, 1850, 60, "5M", 0, 0, 0, true);
        r.setReadBases("ACGTA".getBytes());
        r.setBaseQualities(new byte[] {30, 30, 30, 30, 30});
        out.add(r);

        out.add(read(header, "reverse_strand", 0x10, 0, 1860, 60, "10M", 0, 0, 0, true));

        // Ambiguous bases, and the one that looks ambiguous and is not: BaseUtils maps '*' to A's
        // index ("the wildcard character counts as an A"), so a read of '*' has none.
        r = read(header, "two_n_bases", 0, 0, 1870, 60, "10M", 0, 0, 0, true);
        r.setReadBases("ACGTNNGTAC".getBytes());
        out.add(r);

        r = read(header, "wildcard_bases", 0, 0, 1880, 60, "10M", 0, 0, 0, true);
        r.setReadBases("**********".getBytes());
        out.add(r);

        // The header path: a read in the second group, a read naming a group the header does not
        // declare (which resolves to null, so "has a read group" is false even though the tag is
        // there), and a read aligned past the end of its contig.
        r = read(header, "second_read_group", 0, 0, 1890, 60, "10M", 0, 0, 0, false);
        r.setAttribute("RG", "rg2");
        out.add(r);

        r = read(header, "undeclared_read_group", 0, 0, 1900, 60, "10M", 0, 0, 0, false);
        r.setAttribute("RG", "rg_absent");
        out.add(r);

        r = read(header, "past_contig_end", 0, 1, CHR2 + 50, 60, "10M", 0, 0, 0, true);
        out.add(r);

        // The clipping family. Every one of these keeps ten read bases, so they separate the
        // clipping filters without also moving ReadLengthEqualsCigarLengthReadFilter.
        //
        // Note what the denominator of SoftClippedReadFilter actually is: the sum of *every* cigar
        // element's length, hard clips, deletions and skips included, not the read's length. The
        // documentation says "total bases in read". A port written from the documentation would
        // divide by the read length and agree on every record here except the hard-clipped ones.
        out.add(read(header, "soft_clip_both_ends", 0, 0, 1910, 60, "3S4M3S", 0, 0, 0, true));
        out.add(read(header, "soft_clip_one_end", 0, 0, 1920, 60, "5S5M", 0, 0, 0, true));

        // Consecutive soft clips: two elements, one block. OverclippedReadFilter counts blocks,
        // not elements, and this is the record that tells the two apart.
        out.add(read(header, "soft_clip_consecutive", 0, 0, 1930, 60, "1S2S5M2S", 0, 0, 0, true));

        // Hard clips count towards ExcessiveEndClippedReadFilter and towards SoftClippedReadFilter's
        // denominator, but never towards its numerator.
        out.add(read(header, "hard_clip_front", 0, 0, 1940, 60, "9H1S9M", 0, 0, 0, true));
        out.add(read(header, "hard_clip_back", 0, 0, 1950, 60, "9M1S9H", 0, 0, 0, true));

        // Tag-reading filters. OD is an integer, TV a float, and one of each is text that does not
        // parse, which is how both filters reach their exception rather than a decision.
        r = read(header, "optical_duplicate", 0, 0, 1960, 60, "10M", 0, 0, 0, true);
        r.setAttribute("OD", 3);
        out.add(r);

        r = read(header, "optical_duplicate_zero", 0, 0, 1970, 60, "10M", 0, 0, 0, true);
        r.setAttribute("OD", 0);
        out.add(r);

        r = read(header, "optical_duplicate_not_a_number", 0, 0, 1980, 60, "10M", 0, 0, 0, true);
        r.setAttribute("OD", "many");
        out.add(r);

        // Positive and negative zero, separately: EQUAL and NOT_EQUAL go through Float.equals,
        // which compares floatToIntBits and so calls -0.0 different from 0.0. Every other operator
        // unboxes to a primitive comparison, where the two are equal. One record cannot show that;
        // two can.
        r = read(header, "tag_zero", 0, 0, 1990, 60, "10M", 0, 0, 0, true);
        r.setAttribute("TV", 0.0f);
        out.add(r);

        r = read(header, "tag_negative_zero", 0, 0, 2000, 60, "10M", 0, 0, 0, true);
        r.setAttribute("TV", -0.0f);
        out.add(r);

        r = read(header, "tag_nan", 0, 0, 1690, 60, "10M", 0, 0, 0, true);
        r.setAttribute("TV", Float.NaN);
        out.add(r);

        r = read(header, "tag_high", 0, 0, 1680, 60, "10M", 0, 0, 0, true);
        r.setAttribute("TV", 7.5f);
        out.add(r);

        r = read(header, "tag_not_a_number", 0, 0, 1670, 60, "10M", 0, 0, 0, true);
        r.setAttribute("TV", "high");
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
