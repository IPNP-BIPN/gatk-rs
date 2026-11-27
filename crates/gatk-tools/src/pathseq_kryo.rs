//! What the two PathSeq builders leave on disk: a `PSKmerSet` and a `PSTaxonomyDatabase`, through
//! Kryo.
//!
//! The arithmetic is [`crate::pathseq_kmers`] and [`crate::pathseq_taxonomy`]; this is the file
//! each becomes. Every rule here is the `kryo-stream` golden's, measured in the pinned container,
//! and the encoding itself is [`gatk_engine::kryo`].
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

//!
//! # The taxonomy, outside in
//!
//! `PSTaxonomyDatabase` is NOT written the same way: `writeTaxonomyDatabase` builds its own `Kryo`
//! and turns references off before the first write, so the file carries no marker at all and opens
//! with the root's string. Measured: every run in `pathseq-taxonomy-kryo` opens `8231`, where the
//! hand-built case in `kryo-stream`, written through a default `Kryo`, opens `01`. Then:
//!
//! - **`PSTree`** writes its root as a STRING, the node count as a fixed int, and each node as its
//!   id (a fixed int) followed by the node;
//! - **`PSTreeNode`** writes its name and its rank as strings, its parent as a fixed int, its
//!   length as a fixed long, the number of children as a fixed int, and each child as a STRING;
//! - **the database** then writes the map's size as a fixed int and each entry as two strings, the
//!   contig name and the taxon id in decimal.
//!
//! Three of those are written in a `HashMap`'s order -- the nodes, each node's children and the map
//! -- so the bytes are only as right as [`crate::pathseq_taxonomy`]'s tables, and a table whose
//! order was never measured is refused here rather than written.
//!
//! # Reading this file without reading Rust
//!
//! - `output: &mut Output` is the Kryo `Output` borrowed for writing; the caller keeps it;
//! - `|inner| write_kmer_set(inner, ...)` is a closure, Java's `inner -> writeKmerSet(inner, ...)`;
//! - `n as i32` converts a `usize` count to Java's `int` without a range check (a count above
//!   `i32::MAX` would wrap rather than fail);
//! - `?` after a call returns its error to the caller at once, Java's `throws` at one call site;
//! - `node.name.as_deref()` turns an owned optional `String` into the borrowed `Option<&str>`
//!   [`Output::write_string`] takes, so a Java `null` name stays null.

use crate::pathseq_taxonomy::{PsTree, TreeNode};
use gatk_engine::hopscotch::{LargeLongHopscotchSet, LongHopscotchSet};
use gatk_engine::java_hash::{HashOrderError, JavaHashMap};
use gatk_engine::kryo::Output;

/// `LongHopscotchSet.Serializer.write`, without a reference marker of its own.
///
/// `LongHopscotchSet.serialize`, two fixed ints and then every entry as a fixed long:
///
/// ```java
/// stream.writeInt(capacity);
/// stream.writeInt(size);
///
/// // write the chain heads, and then the squatters
/// ```
///
/// The order of the longs is the table's, which [`LongHopscotchSet::serialized_entries`]
/// reproduces (chain heads in bucket order, then the other occupied buckets in bucket order), so
/// two sets holding the same values in different tables write different bytes. The reference's
/// closing `count != size` check (an `IllegalStateException`) is an assertion inside
/// `serialized_entries`.
pub fn write_hopscotch_set(output: &mut Output, set: &LongHopscotchSet) {
    output.write_int(set.capacity());
    output.write_int(set.size());
    for entry in set.serialized_entries() {
        output.write_long(entry);
    }
}

/// `LargeLongHopscotchSet.serialize`, which writes its partitions with references off.
///
/// ```java
/// kryo.setReferences(false);
///
/// stream.writeInt(sets.size());
/// for (final LongHopscotchSet set : sets) {
///     kryo.writeObject(stream, set);
/// }
/// ```
///
/// So each partition is a `writeObject` with no marker, which is a plain call here.
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
///
/// ```java
/// output.writeInt(kmerSize);
/// output.writeLong(kmerMask.getLong());
/// kryo.writeObject(output, kmerSet);
/// ```
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
///
/// `PSKmerUtils.writeKryoObject` uses a default `new Kryo()`, so references are on and the file
/// opens with `01`. This is the `else` branch of `PathSeqBuildKmers` (`bloomFpp` zero); the
/// `PSKmerBloomFilter` it writes when `bloomFpp > 0` is not written by this file.
pub fn kmer_set_file(kmer_size: i32, kmer_mask: i64, set: &LargeLongHopscotchSet) -> Vec<u8> {
    let mut output = Output::new();
    output.write_object(true, |inner| {
        write_kmer_set(inner, kmer_size, kmer_mask, set)
    });
    output.bytes().to_vec()
}

/// `PSTreeNode.serialize`.
///
/// ```java
/// output.writeString(name);
/// output.writeString(rank);
/// output.writeInt(parent);
/// output.writeLong(length);
/// output.writeInt(children.size());
/// for (final int child : children) {
///     output.writeString(String.valueOf(child));
/// }
/// ```
///
/// The child ids are written as decimal STRINGS, not ints, so a child `2` is the two bytes `8232`
/// and a child `131567` is six. `children` is a `HashSet<Integer>`, which the port models as a
/// `JavaHashMap` with empty values, so `.keys()` walks it in the table's bucket order.
pub fn write_tree_node(output: &mut Output, node: &TreeNode) {
    output.write_string(node.name.as_deref());
    output.write_string(node.rank.as_deref());
    output.write_int(node.parent);
    output.write_long(node.length);
    output.write_int(node.children.len() as i32);
    for child in node.children.keys() {
        output.write_string(Some(&child.to_string()));
    }
}

/// `PSTree.serialize`, whose nested nodes carry no marker because it turned references off.
///
/// ```java
/// output.writeString(String.valueOf(root));
/// output.writeInt(tree.size());
/// for (final int key : tree.keySet()) {
///     output.writeInt(key);
///     kryo.writeObject(output, tree.get(key));
/// }
/// ```
///
/// The root is a string and every node id an int, the same number in two encodings: golden row
/// `tree-three-nodes` opens `01 8231 00000003 00000001`, the root `"1"` and then node id `1`.
pub fn write_tree(output: &mut Output, tree: &PsTree) {
    output.write_string(Some(&tree.root().to_string()));
    output.write_int(tree.node_ids().len() as i32);
    for (id, node) in tree.nodes() {
        output.write_int(*id);
        write_tree_node(output, node);
    }
}

/// `PSTaxonomyDatabase.serialize`, over the map's entries in the order they are handed in.
///
/// The entries are an iterator rather than the map so that the `LinkedHashMap` case the golden
/// measures, which pins the encoding apart from any order, is the same function.
///
/// ```java
/// kryo.writeObject(output, tree);
/// output.writeInt(accessionToTaxId.size());
/// for (final String key : accessionToTaxId.keySet()) {
///     output.writeString(key);
///     output.writeString(String.valueOf(accessionToTaxId.get(key)));
/// }
/// ```
///
/// The tree is written with references already off, so it has no marker. The taxon id, an
/// `Integer`, is written as its decimal STRING, like the child ids.
///
/// `impl ExactSizeIterator<Item = ...>` accepts any iterator that knows its length up front,
/// which is what lets the size be written before the entries. `'a` is a lifetime: it says the
/// borrowed names and ids live at least as long as the call, and changes nothing at run time.
pub fn write_taxonomy_database<'a>(
    output: &mut Output,
    tree: &PsTree,
    entries: impl ExactSizeIterator<Item = (&'a String, &'a i32)>,
) {
    write_tree(output, tree);
    output.write_int(entries.len() as i32);
    for (name, tax_id) in entries {
        output.write_string(Some(name));
        output.write_string(Some(&tax_id.to_string()));
    }
}

/// The whole file `PathSeqBuildReferenceTaxonomy` writes, which has no reference marker.
///
/// Refused if any table on the way crowded a bucket past what the probes measured.
///
/// `writeTaxonomyDatabase` calls `kryo.setReferences(false)` before its one `writeObject`, so the
/// file opens with the root's string and not with `01`: every run in `pathseq-taxonomy-kryo`
/// opens `8231`. The `Result` is this port's refusal, not a reference exception: the Java always
/// writes, and the port declines when it cannot vouch for the `HashMap` order.
pub fn taxonomy_database_file(
    tree: &PsTree,
    map: &JavaHashMap<String, i32>,
) -> Result<Vec<u8>, HashOrderError> {
    // Both orders are checked before any byte is built, so a refusal returns no partial stream.
    tree.check_order()?;
    map.check()?;
    // Collected into a `Vec` so the iterator handed on knows its length.
    let entries: Vec<(&String, &i32)> = map.iter().collect();
    let mut output = Output::new();
    output.write_object(false, |inner| {
        write_taxonomy_database(inner, tree, entries.into_iter())
    });
    Ok(output.bytes().to_vec())
}
