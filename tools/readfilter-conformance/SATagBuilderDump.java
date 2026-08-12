/*
 * SATagBuilder, taken from the reference.
 *
 * The SA tag is how a split read says where its other pieces went, and SplitNCigarReads is what
 * needs it: `repairSupplementaryTags` marks every piece after the first as supplementary and calls
 * `setReadsAsSupplemental`, which is this class. Everything below is measured because the class is
 * a builder with an order and a set of refusals, and neither is visible from the field list.
 *
 * Eight behaviours this is built to catch.
 *
 *   - A UNIT IS SIX FIELDS SPLIT WITH `split(",", -1)`, so a trailing empty NM counts as a field and
 *     `chr1,10,+,10M,60,` parses while `chr1,10,+,10M,60` (five fields) is a GATKException. The -1
 *     limit is the whole difference: without it the trailing empty string would be dropped;
 *   - THE THREE VALIDATIONS ARE POS, CIGAR AND MAPQ, and each one accepts `*` as a special case.
 *     A negative POS and a negative MAPQ are refused; the cigar must match `\*|([0-9]+[MIDNSHPX=])+`,
 *     which accepts a concatenation with no length checking and refuses an empty string;
 *   - NM IS NOT VALIDATED AT ALL. Any text survives the round trip, including text that is not a
 *     number;
 *   - A UNIT BUILT FROM A READ TAKES `*` FOR AN ABSENT NM. `toString` also carries defaults for a
 *     null contig, position, cigar and mapping quality, but only the contig's is reachable from a
 *     read: an unmapped read comes out `*,0,+,*,0,*;`, where the position and the mapping quality
 *     are the read's own zeroes and the cigar is what an empty Cigar prints as. The `0` and the
 *     `255` in the source are never the ones that reach the output;
 *   - THE STRAND IS NORMALISED ON THE WAY OUT: anything that is not exactly `-` prints as `+`, so a
 *     unit parsed with a strand of `x` comes back as `+`;
 *   - addTag PUTS A NON-SUPPLEMENTARY READ AT THE FRONT AND A SUPPLEMENTARY ONE AT THE BACK, which
 *     is what makes the primary alignment the first unit of everyone's tag;
 *   - setReadsAsSupplemental MARKS EVERY READ BUT THE FIRST AS SUPPLEMENTARY BEFORE BUILDING ANY
 *     TAG, so the order the units come out in follows the marking and not the argument order. It is
 *     an all-pairs loop: each read gets a unit for every other read, and never for itself;
 *   - AND EXISTING SA TAGS ARE PRESERVED, because the builder parses the read's own tag at
 *     construction and the new units are added to that list. A read that already claimed a piece
 *     keeps the claim.
 *
 * `setSATag` on a builder with no units writes nothing at all, so a read with no SA tag and no
 * additions keeps having no SA tag rather than gaining an empty one.
 *
 * Output:
 *
 *     parse\t<label>\t<the tag as it comes back out>
 *     error\t<label>\t<exception class>:<message>
 *     unit\t<label>\t<the unit built from a read>
 *     group\t<label>\t<read name>\t<flags>\t<the SA tag, or absent>
 *
 * Usage: SATagBuilderDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.SATagBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SATagBuilderDump {

    public static void main(final String[] args) {
        System.out.println("# SATagBuilderDump: the SA tag SplitNCigarReads writes");

        final SAMFileHeader header = header();

        // What survives a parse and comes back out unchanged, and what does not parse at all.
        final String[][] tags = {
                {"one", "chr1,100,+,10M,60,2;"},
                {"two", "chr1,100,+,10M,60,2;chr2,200,-,5M5S,30,0;"},
                {"empty-nm", "chr1,100,+,10M,60,;"},
                {"star-pos", "chr1,*,+,10M,60,2;"},
                {"star-mapq", "chr1,100,+,10M,*,2;"},
                {"star-cigar", "chr1,100,+,*,60,2;"},
                // The strand is normalised on the way out, so this one does not round trip.
                {"odd-strand", "chr1,100,x,10M,60,2;"},
                // NM is never validated, so this one does round trip.
                {"text-nm", "chr1,100,+,10M,60,not-a-number;"},
                // A cigar the aligner would never write, but the regex accepts.
                {"odd-cigar", "chr1,100,+,1M1M1M,60,2;"},
                {"five-fields", "chr1,100,+,10M,60;"},
                {"seven-fields", "chr1,100,+,10M,60,2,extra;"},
                {"negative-pos", "chr1,-1,+,10M,60,2;"},
                {"negative-mapq", "chr1,100,+,10M,-1,2;"},
                {"bad-cigar", "chr1,100,+,10Z,60,2;"},
                {"empty-cigar", "chr1,100,+,,60,2;"},
        };
        for (final String[] pair : tags) {
            roundTrip(header, pair[0], pair[1]);
        }

        // The unit a read itself becomes, which is what every group tag is made of.
        unit(header, "plain", read(header, "plain", "chr1", 100, "10M", 0, 60, null));
        unit(header, "with-nm", read(header, "with-nm", "chr1", 100, "10M", 0, 60, 3));
        unit(header, "reverse", read(header, "reverse", "chr1", 100, "10M", 0x10, 60, 3));
        unit(header, "zero-mapq", read(header, "zero-mapq", "chr1", 100, "10M", 0, 0, null));
        // No contig at all, which is where toString's own defaults become reachable.
        unit(header, "unmapped", unmapped(header, "unmapped"));

        // Two reads, then three, which is the shape SplitNCigarReads produces for a read with one
        // N and with two.
        group(header, "pair", Arrays.asList(
                read(header, "piece", "chr1", 100, "5M5S", 0, 60, 1),
                read(header, "piece", "chr1", 200, "5S5M", 0, 60, 2)));
        group(header, "triple", Arrays.asList(
                read(header, "piece", "chr1", 100, "3M7S", 0, 60, 1),
                read(header, "piece", "chr1", 200, "3S3M4S", 0, 60, 2),
                read(header, "piece", "chr1", 300, "6S4M", 0, 60, 3)));
        // A read that already carries an SA tag: the existing unit is kept and the new ones follow.
        final GATKRead primaryWithTag = read(header, "piece", "chr1", 100, "5M5S", 0, 60, 1);
        primaryWithTag.setAttribute("SA", "chr9,999,-,4M,20,0;");
        group(header, "existing-tag", Arrays.asList(
                primaryWithTag,
                read(header, "piece", "chr1", 200, "5S5M", 0, 60, 2)));
        // One where the primary is already marked supplementary, which the method does not undo.
        group(header, "primary-already-supplementary", Arrays.asList(
                read(header, "piece", "chr1", 100, "5M5S", 0x800, 60, 1),
                read(header, "piece", "chr1", 200, "5S5M", 0, 60, 2)));
        // And a single read, where the all-pairs loop has no pairs and nothing is written.
        group(header, "single", Arrays.asList(read(header, "alone", "chr1", 100, "10M", 0, 60, 1)));
    }

    /** A tag parsed and written straight back out, or the refusal it raised. */
    static void roundTrip(final SAMFileHeader header, final String label, final String tag) {
        final GATKRead read = read(header, "carrier", "chr1", 1, "10M", 0, 60, null);
        read.setAttribute("SA", tag);
        try {
            final SATagBuilder builder = new SATagBuilder(read);
            builder.setSATag();
            System.out.printf("parse\t%s\t%s%n", label, read.getAttributeAsString("SA"));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }

    /** The single unit a read becomes, read off a second read's tag. */
    static void unit(final SAMFileHeader header, final String label, final GATKRead read) {
        final GATKRead carrier = read(header, "carrier", "chr1", 1, "10M", 0, 60, null);
        final SATagBuilder builder = new SATagBuilder(carrier);
        builder.addTag(read);
        builder.setSATag();
        System.out.printf("unit\t%s\t%s%n", label, carrier.getAttributeAsString("SA"));
    }

    /** One family of reads set as supplemental to each other, with the tag each one ends up with. */
    static void group(final SAMFileHeader header, final String label, final List<GATKRead> family) {
        final List<GATKRead> reads = new ArrayList<>(family);
        final GATKRead primary = reads.remove(0);
        SATagBuilder.setReadsAsSupplemental(primary, reads);

        final List<GATKRead> all = new ArrayList<>();
        all.add(primary);
        all.addAll(reads);
        for (final GATKRead read : all) {
            final String tag = read.getAttributeAsString("SA");
            System.out.printf("group\t%s\t%s:%d\t%d\t%s%n", label, read.getName(), read.getStart(),
                    read.convertToSAMRecord(header).getFlags(), tag == null ? "absent" : tag);
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 5000),
                new SAMSequenceRecord("chr2", 5000),
                new SAMSequenceRecord("chr9", 5000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        return header;
    }

    static GATKRead read(final SAMFileHeader header, final String name, final String contig,
                         final int start, final String cigar, final int flags, final int mapq,
                         final Integer nm) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases("ACGTACGTAC".getBytes(StandardCharsets.UTF_8));
        record.setBaseQualities(new byte[] {35, 35, 35, 35, 35, 35, 35, 35, 35, 35});
        record.setMappingQuality(mapq);
        record.setAttribute("RG", "rg1");
        if (nm != null) {
            record.setAttribute("NM", nm);
        }
        return new SAMRecordToGATKReadAdapter(record);
    }

    static GATKRead unmapped(final SAMFileHeader header, final String name) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(0x4);
        record.setReadUnmappedFlag(true);
        record.setReadBases("ACGTACGTAC".getBytes(StandardCharsets.UTF_8));
        record.setBaseQualities(new byte[] {35, 35, 35, 35, 35, 35, 35, 35, 35, 35});
        record.setAttribute("RG", "rg1");
        return new SAMRecordToGATKReadAdapter(record);
    }
}
