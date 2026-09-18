//! What `PathSeqBuildKmers` leaves on disk: a `PSKmerSet` through Kryo.
//!
//! The arithmetic of the k-mer set is [`crate::pathseq_kmers`]; this is the file it becomes. Every
//! rule here is the `kryo-stream` golden's, measured in the pinned container, and the encoding
//! itself is [`gatk_engine::kryo`].
//!
//! # The shape, outside in
//!
//! `PSKmerUtils.writeKryoObject` is `new Kryo()` and one `writeObject`, and Kryo 4 has references
//! ON by default, so the file opens with the reference marker `01`. Then:
//!
//! - **`PSKmerSet`** writes `kmerSize` as a fixed int, the mask as a fixed long, and the set;
//! - **`LargeLongHopscotchSet`** turns references OFF around its own body, writes the number of
//!   partitions as a fixed int, and then each partition. So the partitions carry no marker of
//!   their own -- but the set ITSELF does, because `PSKmerSet` never turned them off, which is why
//!   the file carries `01` twice;
//! - **`LongHopscotchSet`** writes its capacity and its size, then the chain heads in bucket order
//!   and the squatters after them.
//!
//! The capacity is the one fact a port cannot shortcut: it is the legal size above the request, so
//! a set asked for eight and one asked for sixty-four both write 251 and the same seventy-three
//! bytes.

use gatk_engine::hopscotch::{LargeLongHopscotchSet, LongHopscotchSet};
use gatk_engine::kryo::Output;

/// `LongHopscotchSet.Serializer.write`, without a reference marker of its own.
pub fn write_hopscotch_set(output: &mut Output, set: &LongHopscotchSet) {
    output.write_int(set.capacity());
    output.write_int(set.size());
    for entry in set.serialized_entries() {
        output.write_long(entry);
    }
}

/// `LargeLongHopscotchSet.serialize`, which writes its partitions with references off.
pub fn write_large_hopscotch_set(output: &mut Output, set: &LargeLongHopscotchSet) {
    output.write_int(set.sets().len() as i32);
    for partition in set.sets() {
        write_hopscotch_set(output, partition);
    }
}

/// `PSKmerSet.serialize`: the k-mer size, the mask as a long, and the set.
///
/// The set is written through `kryo.writeObject`, and `PSKmerSet` does NOT turn references off, so
/// the nested object carries a marker of its own. `LargeLongHopscotchSet` turns them off for what
/// IT writes, which is why its partitions carry none. Measured: `kmer-set` has `01` twice and
/// `large-hopscotch` has it once.
pub fn write_kmer_set(
    output: &mut Output,
    kmer_size: i32,
    kmer_mask: i64,
    set: &LargeLongHopscotchSet,
) {
    output.write_int(kmer_size);
    output.write_long(kmer_mask);
    output.write_object(true, |inner| write_large_hopscotch_set(inner, set));
}

/// The whole file `PathSeqBuildKmers` writes, reference marker included.
pub fn kmer_set_file(kmer_size: i32, kmer_mask: i64, set: &LargeLongHopscotchSet) -> Vec<u8> {
    let mut output = Output::new();
    output.write_object(true, |inner| {
        write_kmer_set(inner, kmer_size, kmer_mask, set)
    });
    output.bytes().to_vec()
}
