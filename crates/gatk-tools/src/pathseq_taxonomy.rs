//! `PathSeqBuildReferenceTaxonomy`, ported from the tool, `PSBuildReferenceTaxonomyUtils`,
//! `PSTree` and `PSPathogenReferenceTaxonProperties` (GATK 4.6.2.0).
//!
//! The NCBI taxonomy dump and one or two accession catalogs turned into the tree PathSeq scores
//! against, trimmed to the taxa the reference actually holds. The tool writes a Kryo serialisation,
//! which is not what is ported: what is here is the tree it holds and the map beside it.
//!
//! # The map is keyed by the contig name, not by the accession
//!
//! ```java
//! addReferenceAccessionToTaxon(taxId, nameAndLength._1, nameAndLength._2, taxIdToProperties);
//! ```
//!
//! `nameAndLength._1` is the record's NAME. The accession is only ever the key the catalog is
//! searched with; what lands in the taxon's list, and therefore in the map the tool writes, is the
//! whole contig name. An entry looks up as `ref|NC_BACT.1|`, not as `NC_BACT.1`, and anything
//! reading that map has to know the reference's naming rather than NCBI's.
//!
//! # A reference name is read for `taxid|` first
//!
//! ```java
//! for (int i = 0; i < tokens.length - 1 && recordTaxId == PSTree.NULL_NODE; i++) {
//!     if (tokens[i].equals("ref")) { recordAccession = tokens[i + 1]; }
//!     else if (tokens[i].equals("taxid")) { recordTaxId = parseTaxonId(tokens[i + 1]); }
//! }
//! ```
//!
//! The loop stops at the first taxon id, so a name carrying both `ref|` and `taxid|` is placed by
//! its taxon and its accession is never looked up. A name carrying neither falls back to the first
//! WORD of the first bar-delimited token.
//!
//! # A blank line ends a catalog
//!
//! `while ((line = reader.readLine()) != null && !line.isEmpty())` stops on the first empty line, so
//! everything after it is invisible and its contigs are merely reported as not found. A line with
//! too few columns, on the other hand, is a refusal, and its message says GenBank whatever format
//! was being read.
//!
//! # The length filter is on the map, not the tree
//!
//! `--min-non-virus-contig-length` decides which contigs reach the map; the tree's node lengths
//! were summed before it ran and still count the contigs it drops. So a node can carry a length no
//! reachable contig accounts for. A contig whose path holds the virus node, 10239, is never
//! dropped whatever its length.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// `PSTree.NULL_NODE`.
pub const NULL_NODE: i32 = 0;

/// `PSTaxonomyConstants.ROOT_ID`.
pub const ROOT_ID: i32 = 1;

/// `PSTaxonomyConstants.VIRUS_ID`, the node the length filter is exempted by.
pub const VIRUS_ID: i32 = 10239;

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaxonomyError {
    /// `parseTaxonId` on anything that is not an integer.
    NotAnInteger { value: String },
    /// A catalog line with fewer columns than the format needs. The message names GenBank whatever
    /// the format was.
    TooFewColumns {
        expected: usize,
        found: usize,
        line: u64,
    },
    /// Neither catalog given.
    NoCatalog,
    /// Nothing in the reference landed on a taxon.
    NoRelevantTaxa,
    /// A names or nodes line with too few columns.
    MalformedDump { file: &'static str, found: usize },
}

impl TaxonomyError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    }

    /// The message as `UserException.BadInput` prints it, its "Bad input: " prefix included.
    pub fn message(&self) -> String {
        let body = match self {
            TaxonomyError::NotAnInteger { value } => {
                format!("Expected taxonomy ID to be an integer but found \"{value}\"")
            }
            TaxonomyError::TooFewColumns {
                expected,
                found,
                line,
            } => format!(
                "Expected at least {expected} tab-delimited columns in GenBank catalog file, \
                 but only found {found} on line {line}"
            ),
            TaxonomyError::NoCatalog => {
                "At least one of --refseq-catalog or --genbank-catalog must be specified"
                    .to_string()
            }
            TaxonomyError::NoRelevantTaxa => {
                "Did not find any taxa corresponding to reference sequence names.\n\n\
                 Check that reference names follow one of the required formats:\n\n\
                 \t...|ref|<accession.version>|...\n\
                 \t...|taxid|<taxonomy_id>|...\n\
                 \t<accession.version><mask>..."
                    .to_string()
            }
            TaxonomyError::MalformedDump { file, found } => format!(
                "Expected at least {} columns in tax dump {file} file but found {found}",
                if *file == "names" { 4 } else { 3 }
            ),
        };
        format!("Bad input: {body}")
    }
}

/// `parseTaxonId`.
fn parse_taxon_id(value: &str) -> Result<i32, TaxonomyError> {
    value
        .parse::<i32>()
        .map_err(|_| TaxonomyError::NotAnInteger {
            value: value.to_string(),
        })
}

/// `String.split("\\s*\\|\\s*")`, whose trailing empty fields Java drops.
fn split_bars(text: &str) -> Vec<String> {
    let mut fields: Vec<String> = Vec::new();
    let mut current = String::new();
    let bytes: Vec<char> = text.chars().collect();
    let mut index = 0;
    while index < bytes.len() {
        // A separator is any run of whitespace around exactly one bar.
        let mut probe = index;
        while probe < bytes.len() && bytes[probe].is_whitespace() {
            probe += 1;
        }
        if probe < bytes.len() && bytes[probe] == '|' {
            probe += 1;
            while probe < bytes.len() && bytes[probe].is_whitespace() {
                probe += 1;
            }
            fields.push(std::mem::take(&mut current));
            index = probe;
            continue;
        }
        current.push(bytes[index]);
        index += 1;
    }
    fields.push(current);
    while fields.last().is_some_and(|field| field.is_empty()) {
        fields.pop();
    }
    fields
}

/// `PSPathogenReferenceTaxonProperties`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TaxonProperties {
    pub name: Option<String>,
    pub rank: Option<String>,
    pub parent: i32,
    /// The running total, which `addAccession` adds to whether or not the name was already there.
    pub length: i64,
    pub accessions: BTreeMap<String, i64>,
}

impl TaxonProperties {
    fn named(name: &str) -> Self {
        TaxonProperties {
            name: Some(name.to_string()),
            ..Default::default()
        }
    }

    fn add_accession(&mut self, name: &str, length: i64) {
        self.accessions.insert(name.to_string(), length);
        self.length += length;
    }
}

fn add_reference_accession(
    properties: &mut BTreeMap<i32, TaxonProperties>,
    tax_id: i32,
    name: &str,
    length: i64,
) {
    properties
        .entry(tax_id)
        .or_default()
        .add_accession(name, length);
}

/// One node of `PSTree`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TreeNode {
    pub name: Option<String>,
    pub rank: Option<String>,
    pub parent: i32,
    pub length: i64,
    pub children: BTreeSet<i32>,
}

/// `PSTree`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsTree {
    root: i32,
    nodes: BTreeMap<i32, TreeNode>,
}

impl PsTree {
    pub fn new(root: i32) -> Self {
        let mut nodes = BTreeMap::new();
        nodes.insert(
            root,
            TreeNode {
                name: Some("root".to_string()),
                rank: Some("root".to_string()),
                parent: NULL_NODE,
                length: 0,
                children: BTreeSet::new(),
            },
        );
        PsTree { root, nodes }
    }

    /// `addNode`, which creates a placeholder for a parent it has not seen.
    pub fn add_node(&mut self, id: i32, name: &str, parent: i32, length: i64, rank: &str) {
        let node = self.nodes.entry(id).or_default();
        node.name = Some(name.to_string());
        node.parent = parent;
        node.length = length;
        node.rank = Some(rank.to_string());
        self.nodes.entry(parent).or_default().children.insert(id);
    }

    pub fn node_ids(&self) -> Vec<i32> {
        self.nodes.keys().copied().collect()
    }

    pub fn has_node(&self, id: i32) -> bool {
        self.nodes.contains_key(&id)
    }

    pub fn name_of(&self, id: i32) -> Option<&str> {
        self.nodes.get(&id).and_then(|node| node.name.as_deref())
    }

    pub fn rank_of(&self, id: i32) -> Option<&str> {
        self.nodes.get(&id).and_then(|node| node.rank.as_deref())
    }

    pub fn parent_of(&self, id: i32) -> i32 {
        self.nodes.get(&id).map_or(NULL_NODE, |node| node.parent)
    }

    pub fn length_of(&self, id: i32) -> i64 {
        self.nodes.get(&id).map_or(0, |node| node.length)
    }

    /// `getPathOf`, from the node up to the root.
    pub fn path_of(&self, id: i32) -> Vec<i32> {
        let mut path = Vec::new();
        let mut seen = BTreeSet::new();
        let mut current = id;
        while current != NULL_NODE {
            if !seen.insert(current) {
                break;
            }
            match self.nodes.get(&current) {
                Some(node) => {
                    path.push(current);
                    current = node.parent;
                }
                None => break,
            }
        }
        path
    }

    /// `traverse`, breadth first from the root.
    fn reachable(&self) -> BTreeSet<i32> {
        let mut visited = BTreeSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(self.root);
        while let Some(id) = queue.pop_front() {
            if visited.contains(&id) {
                continue;
            }
            if let Some(node) = self.nodes.get(&id) {
                queue.extend(node.children.iter().copied());
            }
            visited.insert(id);
        }
        visited
    }

    /// `removeUnreachableNodes`.
    pub fn remove_unreachable_nodes(&mut self) -> BTreeSet<i32> {
        let reachable = self.reachable();
        let unreachable: BTreeSet<i32> = self
            .nodes
            .keys()
            .copied()
            .filter(|id| !reachable.contains(id))
            .collect();
        self.retain_nodes(&reachable);
        unreachable
    }

    /// `retainNodes`, which also cuts the pointers that would dangle.
    pub fn retain_nodes(&mut self, keep: &BTreeSet<i32>) {
        let mut kept: BTreeMap<i32, TreeNode> = BTreeMap::new();
        for (id, node) in &self.nodes {
            if !keep.contains(id) {
                continue;
            }
            let mut copy = node.clone();
            copy.children.retain(|child| keep.contains(child));
            if !keep.contains(&copy.parent) {
                copy.parent = NULL_NODE;
            }
            kept.insert(*id, copy);
        }
        self.nodes = kept;
    }
}

/// `parseReferenceRecords`: the contigs that name a taxon go straight to the properties, the rest
/// go to a map from accession to the contig's name and length.
#[allow(clippy::type_complexity)]
pub fn parse_reference_records(
    records: &[(String, i64)],
    properties: &mut BTreeMap<i32, TaxonProperties>,
) -> Result<BTreeMap<String, (String, i64)>, TaxonomyError> {
    let mut by_accession = BTreeMap::new();
    for (name, length) in records {
        let tokens = split_bars(name);
        let mut accession: Option<String> = None;
        let mut tax_id = NULL_NODE;
        let mut index = 0;
        while index + 1 < tokens.len() && tax_id == NULL_NODE {
            if tokens[index] == "ref" {
                accession = Some(tokens[index + 1].clone());
            } else if tokens[index] == "taxid" {
                tax_id = parse_taxon_id(&tokens[index + 1])?;
            }
            index += 1;
        }
        if tax_id == NULL_NODE {
            let accession = accession.unwrap_or_else(|| {
                // The default accession is the first word of the first token.
                tokens
                    .first()
                    .map(|token| token.split(' ').next().unwrap_or("").to_string())
                    .unwrap_or_default()
            });
            by_accession.insert(accession, (name.clone(), *length));
        } else {
            add_reference_accession(properties, tax_id, name, *length);
        }
    }
    Ok(by_accession)
}

/// Which columns a catalog holds its taxon and its accession in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatalogFormat {
    RefSeq,
    GenBank,
}

impl CatalogFormat {
    fn tax_id_column(&self) -> usize {
        match self {
            CatalogFormat::RefSeq => 0,
            CatalogFormat::GenBank => 6,
        }
    }

    fn accession_column(&self) -> usize {
        match self {
            CatalogFormat::RefSeq => 2,
            CatalogFormat::GenBank => 1,
        }
    }
}

/// `parseCatalog`, which returns the accessions still not found.
pub fn parse_catalog(
    text: &str,
    format: CatalogFormat,
    by_accession: &BTreeMap<String, (String, i64)>,
    properties: &mut BTreeMap<i32, TaxonProperties>,
    not_found_in: Option<&BTreeSet<String>>,
) -> Result<BTreeSet<String>, TaxonomyError> {
    let mut not_found: BTreeSet<String> = match not_found_in {
        Some(previous) => previous.clone(),
        None => by_accession.keys().cloned().collect(),
    };
    let tax_column = format.tax_id_column();
    let accession_column = format.accession_column();
    let min_columns = tax_column.max(accession_column) + 1;
    let mut line_number: u64 = 1;
    #[allow(clippy::explicit_counter_loop)]
    for line in text.split('\n') {
        // The read loop stops at the first empty line, so everything after it is invisible.
        if line.is_empty() {
            break;
        }
        let tokens = split_with_limit(line.trim(), '\t', min_columns + 1);
        if tokens.len() < min_columns {
            return Err(TaxonomyError::TooFewColumns {
                expected: min_columns,
                found: tokens.len(),
                line: line_number,
            });
        }
        let tax_id = parse_taxon_id(&tokens[tax_column])?;
        let accession = &tokens[accession_column];
        if let Some((name, length)) = by_accession.get(accession) {
            add_reference_accession(properties, tax_id, name, *length);
            not_found.remove(accession);
        }
        line_number += 1;
    }
    Ok(not_found)
}

/// `String.split(regex, limit)` for a plain separator: the last field keeps everything left.
fn split_with_limit(text: &str, separator: char, limit: usize) -> Vec<String> {
    let mut fields: Vec<String> = Vec::new();
    let mut rest = text;
    while fields.len() + 1 < limit {
        match rest.find(separator) {
            Some(position) => {
                fields.push(rest[..position].to_string());
                rest = &rest[position + separator.len_utf8()..];
            }
            None => break,
        }
    }
    fields.push(rest.to_string());
    // Java drops trailing empty fields only when the limit is zero, which it is not here.
    fields
}

/// `parseNamesFile`, which keeps only the scientific names.
pub fn parse_names(
    text: &str,
    properties: &mut BTreeMap<i32, TaxonProperties>,
) -> Result<(), TaxonomyError> {
    for line in text.lines() {
        let tokens = split_bars(line);
        if tokens.len() < 4 {
            return Err(TaxonomyError::MalformedDump {
                file: "names",
                found: tokens.len(),
            });
        }
        if tokens[3] != "scientific name" {
            continue;
        }
        let tax_id = parse_taxon_id(&tokens[0])?;
        let name = tokens[1].clone();
        match properties.get_mut(&tax_id) {
            Some(existing) => existing.name = Some(name),
            None => {
                properties.insert(tax_id, TaxonProperties::named(&name));
            }
        }
    }
    Ok(())
}

/// `parseNodesFile`, which returns the taxa it had never seen.
pub fn parse_nodes(
    text: &str,
    properties: &mut BTreeMap<i32, TaxonProperties>,
) -> Result<Vec<i32>, TaxonomyError> {
    let mut not_found = Vec::new();
    for line in text.lines() {
        let tokens = split_bars(line);
        if tokens.len() < 3 {
            return Err(TaxonomyError::MalformedDump {
                file: "nodes",
                found: tokens.len(),
            });
        }
        let tax_id = parse_taxon_id(&tokens[0])?;
        let parent = parse_taxon_id(&tokens[1])?;
        let rank = tokens[2].clone();
        let node = properties.entry(tax_id).or_insert_with(|| {
            // A node the reference and the names file never mentioned is named after its id.
            not_found.push(tax_id);
            TaxonProperties::named(&format!("tax_{tax_id}"))
        });
        node.rank = Some(rank);
        if tax_id != ROOT_ID {
            // The root's parent stays unset.
            node.parent = parent;
        }
    }
    Ok(not_found)
}

/// `buildTaxonomicTree`.
pub fn build_taxonomic_tree(
    properties: &BTreeMap<i32, TaxonProperties>,
) -> Result<PsTree, TaxonomyError> {
    let mut tree = PsTree::new(ROOT_ID);
    for (tax_id, taxon) in properties {
        if *tax_id == ROOT_ID {
            continue;
        }
        match (&taxon.name, taxon.parent, &taxon.rank) {
            (Some(name), parent, Some(rank)) if parent != NULL_NODE => {
                tree.add_node(*tax_id, name, parent, taxon.length, rank);
            }
            // A node missing a name, a parent or a rank is simply left out.
            _ => {}
        }
    }
    tree.remove_unreachable_nodes();

    let mut relevant: BTreeSet<i32> = BTreeSet::new();
    for (tax_id, taxon) in properties {
        if !taxon.accessions.is_empty() && tree.has_node(*tax_id) {
            relevant.extend(tree.path_of(*tax_id));
        }
    }
    if relevant.is_empty() {
        return Err(TaxonomyError::NoRelevantTaxa);
    }
    tree.retain_nodes(&relevant);
    Ok(tree)
}

/// `removeUnusedTaxIds`.
pub fn remove_unused_tax_ids(properties: &mut BTreeMap<i32, TaxonProperties>, tree: &PsTree) {
    let kept: BTreeSet<i32> = tree.node_ids().into_iter().collect();
    properties.retain(|tax_id, _| kept.contains(tax_id));
}

/// `buildAccessionToTaxIdMap`, which is the map keyed by contig name.
pub fn build_accession_to_tax_id(
    properties: &BTreeMap<i32, TaxonProperties>,
    tree: &PsTree,
    min_non_virus_contig_length: i64,
) -> BTreeMap<String, i32> {
    let mut map = BTreeMap::new();
    for (tax_id, taxon) in properties {
        let is_virus = tree.path_of(*tax_id).contains(&VIRUS_ID);
        for (name, length) in &taxon.accessions {
            if is_virus || *length >= min_non_virus_contig_length {
                map.insert(name.clone(), *tax_id);
            }
        }
    }
    map
}

/// The whole run: the reference's contigs, the catalogs and the two dump files, in.
///
/// The third element is the accessions no catalog accounted for, which the reference only logs.
#[allow(clippy::type_complexity)]
pub fn build(
    contigs: &[(String, i64)],
    refseq_catalog: Option<&str>,
    genbank_catalog: Option<&str>,
    names: &str,
    nodes: &str,
    min_non_virus_contig_length: i64,
) -> Result<(PsTree, BTreeMap<String, i32>, BTreeSet<String>), TaxonomyError> {
    if refseq_catalog.is_none() && genbank_catalog.is_none() {
        return Err(TaxonomyError::NoCatalog);
    }
    let mut properties: BTreeMap<i32, TaxonProperties> = BTreeMap::new();
    let by_accession = parse_reference_records(contigs, &mut properties)?;
    let mut not_found: Option<BTreeSet<String>> = None;
    if let Some(text) = refseq_catalog {
        not_found = Some(parse_catalog(
            text,
            CatalogFormat::RefSeq,
            &by_accession,
            &mut properties,
            None,
        )?);
    }
    if let Some(text) = genbank_catalog {
        not_found = Some(parse_catalog(
            text,
            CatalogFormat::GenBank,
            &by_accession,
            &mut properties,
            not_found.as_ref(),
        )?);
    }
    parse_names(names, &mut properties)?;
    parse_nodes(nodes, &mut properties)?;
    let tree = build_taxonomic_tree(&properties)?;
    remove_unused_tax_ids(&mut properties, &tree);
    let map = build_accession_to_tax_id(&properties, &tree, min_non_virus_contig_length);
    Ok((tree, map, not_found.unwrap_or_default()))
}
