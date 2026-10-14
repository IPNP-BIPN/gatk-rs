//! `--debug-assembly-variants-out`, ported from `ReadThreadingAssembler.addAssembledVariantsToEventMapOutput`
//! and the writer HaplotypeCaller's and Mutect2's engines open for it (GATK 4.6.2.0).
//!
//! Every assembled event of every region, before trimming, as a sites-only record, goes into a
//! `PriorityQueue` ordered by contig and start. Once the queue holds 300, each new event first
//! writes the queue's first 200 **in its array order** and removes none of them, so they are
//! written again later; at shutdown the whole queue is written, again in array order. A heap's
//! array is sorted only while events arrive in order, which the regions mostly do.

use htsjdk_vcf::variant::{Value, VariantContext};

/// The queue and what has been written ahead of it.
#[derive(Debug, Clone, Default)]
pub struct AssemblyVariantsOut {
    /// The contigs in dictionary order, for `VariantContextComparator`.
    contigs: Vec<String>,
    heap: Vec<VariantContext>,
    written: Vec<VariantContext>,
}

impl AssemblyVariantsOut {
    pub fn new(contigs: Vec<String>) -> AssemblyVariantsOut {
        AssemblyVariantsOut {
            contigs,
            heap: Vec::new(),
            written: Vec::new(),
        }
    }

    /// `VariantContextComparator`: the contig's index, then the start.
    fn compare(&self, a: &VariantContext, b: &VariantContext) -> std::cmp::Ordering {
        let index = |vc: &VariantContext| self.contigs.iter().position(|c| *c == vc.contig);
        index(a).cmp(&index(b)).then(a.start.cmp(&b.start))
    }

    /// `queue.add(event.convertToVariantContext("assembly"))`, after the early write.
    pub fn add_event(&mut self, event: &gatk_engine::event_map::Event) {
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
        if event.is_collapsed() {
            vc.attributes
                .push(("XC".to_string(), Value::Str("1".to_string())));
        }
        // `PriorityQueue.offer`: appended, then sifted up past every greater parent.
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
