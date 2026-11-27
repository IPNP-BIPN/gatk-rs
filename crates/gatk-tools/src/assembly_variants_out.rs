//! `--debug-assembly-variants-out`, ported from
//! `ReadThreadingAssembler.addAssembledVariantsToEventMapOutput` and the writer HaplotypeCaller's
//! and Mutect2's engines open for it (GATK 4.7.0.0; `ReadThreadingAssembler.java` and
//! `Mutect2Engine.java` are unchanged since 4.6.2.0, and `HaplotypeCallerEngine.java`'s changes
//! since then do not touch this writer).
//!
//! Every assembled event of every region, before trimming, as a sites-only record, goes into a
//! `PriorityQueue` ordered by contig and start. Once the queue holds 300, each new event first
//! writes the queue's first 200 **in its array order** and removes none of them, so they are
//! written again later; at shutdown the whole queue is written, again in array order. A heap's
//! array is sorted only while events arrive in order, which the regions mostly do.
//!
//! ```java
//! if (queue.size() >= 300) {
//!     queue.stream().limit(200).forEachOrdered(vc ->
//!             assembledEventMapVcfOutputWriter.get().add(vc));
//! }
//! queue.add(event.convertToVariantContext("assembly"));}));
//! ```
//!
//! (`ReadThreadingAssembler.java:884-887`; the `forEachOrdered` line is broken in two here to fit
//! the line limit.) `queue.stream()` walks the `PriorityQueue`'s backing array, not its priority
//! order, which is why the heap is kept here as a plain `Vec` with the sift written out: the
//! array, not only its minimum, is what is written.
//!
//! # Reading this file without reading Rust
//!
//! - `Vec<T>` is Java's `ArrayList<T>`; `&self` / `&mut self` is `this`, read-only or writable.
//! - `.iter().take(200).cloned().collect()` is `stream().limit(200)` with each element copied
//!   into a new list.
//! - `|vc| ..` is a closure, Java's lambda; `Option<usize>` is "an index or nothing".

use htsjdk_vcf::variant::{Value, VariantContext};

/// The queue and what has been written ahead of it.
#[derive(Debug, Clone, Default)]
pub struct AssemblyVariantsOut {
    /// The contigs in dictionary order, for `VariantContextComparator`.
    contigs: Vec<String>,
    /// The `PriorityQueue`'s backing array, a binary min-heap under `compare`.
    heap: Vec<VariantContext>,
    /// Every record written ahead of shutdown, in the order the writer received them.
    written: Vec<VariantContext>,
}

impl AssemblyVariantsOut {
    /// The engines' queue, with no record yet written:
    ///
    /// ```java
    /// new PriorityQueue<>(200,
    ///         new VariantContextComparator(readsHeader.getSequenceDictionary()))
    /// ```
    ///
    /// (`HaplotypeCallerEngine.java:335`, broken in two here.) The 200 is the initial capacity
    /// and bounds nothing.
    pub fn new(contigs: Vec<String>) -> AssemblyVariantsOut {
        AssemblyVariantsOut {
            contigs,
            heap: Vec::new(),
            written: Vec::new(),
        }
    }

    /// `VariantContextComparator` (htsjdk 5.0.0): the contig's index, then the start.
    ///
    /// ```java
    /// return contigCompare != 0 ? contigCompare
    ///         : firstVariantContext.getStart() - secondVariantContext.getStart();
    /// ```
    ///
    /// A contig missing from the dictionary is a `NullPointerException` in the reference and
    /// sorts first here (`None` before any `Some`); assembled events are always on a dictionary
    /// contig, so neither is reached.
    fn compare(&self, a: &VariantContext, b: &VariantContext) -> std::cmp::Ordering {
        let index = |vc: &VariantContext| self.contigs.iter().position(|c| *c == vc.contig);
        index(a).cmp(&index(b)).then(a.start.cmp(&b.start))
    }

    /// `queue.add(event.convertToVariantContext("assembly"))`, after the early write.
    pub fn add_event(&mut self, event: &gatk_engine::event_map::Event) {
        // The early write: the first 200 of the backing array, none removed from the queue.
        if self.heap.len() >= 300 {
            let first: Vec<VariantContext> = self.heap.iter().take(200).cloned().collect();
            self.written.extend(first);
        }
        let mut vc = VariantContext::new(
            event.contig(),
            i64::from(event.start()),
            vec![event.ref_allele().clone(), event.alt_allele().clone()],
        );
        vc.stop = i64::from(event.end());
        // `convertToVariantContext` copies the event's attributes, and the only one it can carry
        // is `EventMap`'s `EXT_COLLAPSED_TAG`, set to "1".
        if event.is_collapsed() {
            vc.attributes
                .push(("XC".to_string(), Value::Str("1".to_string())));
        }
        // `PriorityQueue.offer`: appended, then sifted up past every greater parent. A parent
        // that compares equal stops the sift, so an event tied on contig and start stays below the
        // one that arrived first.
        self.heap.push(vc);
        let mut k = self.heap.len() - 1;
        while k > 0 {
            let parent = (k - 1) / 2;
            if self.compare(&self.heap[k], &self.heap[parent]) != std::cmp::Ordering::Less {
                break;
            }
            self.heap.swap(k, parent);
            k = parent;
        }
    }

    /// What the file holds at shutdown: the early writes, then the queue in its array order.
    pub fn records(&self) -> Vec<VariantContext> {
        let mut out = self.written.clone();
        out.extend(self.heap.iter().cloned());
        out
    }
}
