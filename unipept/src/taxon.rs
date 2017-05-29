//! Defines operations and data structures over taxons.

use std;
use std::fmt;
use std::io;
use std::io::BufRead;
use std::fs::File;
use std::collections::HashMap;
use std::collections::HashSet;

use std::str::FromStr;

/// Represents the rank of a [Taxon](struct.Taxon.html)
#[allow(missing_docs)]
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum Rank {
    NoRank, Superkingdom, Kingdom, Subkingdom, Superphylum, Phylum, Subphylum, Superclass, Class,
    Subclass, Infraclass, Superorder, Order, Suborder, Infraorder, Parvorder, Superfamily, Family,
    Subfamily, Tribe, Subtribe, Genus, Subgenus, SpeciesGroup, SpeciesSubgroup, Species,
    Subspecies, Varietas, Forma
}

impl Rank {
    /// Converts a rank to a usize.
    pub fn index(&self) -> usize {
        *self as usize
    }
}

impl FromStr for Rank {
    type Err = Error;
    fn from_str(rank: &str) -> Result<Self> {
        match rank {
            "no rank"          => Ok(Rank::NoRank),
            "superkingdom"     => Ok(Rank::Superkingdom),
            "kingdom"          => Ok(Rank::Kingdom),
            "subkingdom"       => Ok(Rank::Subkingdom),
            "superphylum"      => Ok(Rank::Superphylum),
            "phylum"           => Ok(Rank::Phylum),
            "subphylum"        => Ok(Rank::Subphylum),
            "superclass"       => Ok(Rank::Superclass),
            "class"            => Ok(Rank::Class),
            "subclass"         => Ok(Rank::Subclass),
            "infraclass"       => Ok(Rank::Infraclass),
            "superorder"       => Ok(Rank::Superorder),
            "order"            => Ok(Rank::Order),
            "suborder"         => Ok(Rank::Suborder),
            "infraorder"       => Ok(Rank::Infraorder),
            "parvorder"        => Ok(Rank::Parvorder),
            "superfamily"      => Ok(Rank::Superfamily),
            "family"           => Ok(Rank::Family),
            "subfamily"        => Ok(Rank::Subfamily),
            "tribe"            => Ok(Rank::Tribe),
            "subtribe"         => Ok(Rank::Subtribe),
            "genus"            => Ok(Rank::Genus),
            "subgenus"         => Ok(Rank::Subgenus),
            "species group"    => Ok(Rank::SpeciesGroup),
            "species subgroup" => Ok(Rank::SpeciesSubgroup),
            "species"          => Ok(Rank::Species),
            "subspecies"       => Ok(Rank::Subspecies),
            "varietas"         => Ok(Rank::Varietas),
            "forma"            => Ok(Rank::Forma),
            _                  => Err(ErrorKind::UnknownRank(rank.to_string()).into())
        }
    }
}

impl fmt::Display for Rank {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let stringified = match *self {
            Rank::NoRank          => "no rank",
            Rank::Superkingdom    => "superkingdom",
            Rank::Kingdom         => "kingdom",
            Rank::Subkingdom      => "subkingdom",
            Rank::Superphylum     => "superphylum",
            Rank::Phylum          => "phylum",
            Rank::Subphylum       => "subphylum",
            Rank::Superclass      => "superclass",
            Rank::Class           => "class",
            Rank::Subclass        => "subclass",
            Rank::Infraclass      => "infraclass",
            Rank::Superorder      => "superorder",
            Rank::Order           => "order",
            Rank::Suborder        => "suborder",
            Rank::Infraorder      => "infraorder",
            Rank::Parvorder       => "parvorder",
            Rank::Superfamily     => "superfamily",
            Rank::Family          => "family",
            Rank::Subfamily       => "subfamily",
            Rank::Tribe           => "tribe",
            Rank::Subtribe        => "subtribe",
            Rank::Genus           => "genus",
            Rank::Subgenus        => "subgenus",
            Rank::SpeciesGroup    => "species group",
            Rank::SpeciesSubgroup => "species subgroup",
            Rank::Species         => "species",
            Rank::Subspecies      => "subspecies",
            Rank::Varietas        => "varietas",
            Rank::Forma           => "forma",
        };
        write!(f, "{}", stringified)
    }
}

/// A unique identifier for a [Taxon](struct.Taxon.html).
pub type TaxonId = usize;

/// Represents a group of organisms with similar qualities.
#[derive(Clone, PartialEq, Debug)]
pub struct Taxon {
    /// The taxon's unique id
    pub id: TaxonId,
    /// The taxon's name
    pub name: String,
    /// The rank of the taxon
    pub rank: Rank,
    /// The taxon's parent
    pub parent: TaxonId,
    /// Whether the taxon is valid. `false` taxons are discarded in some calculations
    pub valid: bool
}

impl Taxon {
    /// Creates a taxon from the given arguments.
    pub fn new(id: TaxonId, name: String, rank: Rank, parent: TaxonId, valid: bool) -> Taxon {
        Taxon { id: id, name: name, rank: rank, parent: parent, valid: valid }
    }

    /// Creates a taxon from the given arguments (using a &str instead of a String).
    pub fn from_static(id: TaxonId, name: &str, rank: Rank, parent: TaxonId, valid: bool) -> Taxon {
        Taxon::new(id, name.to_string(), rank, parent, valid)
    }
}

impl FromStr for Taxon {
    type Err = Error;

    /// Parses a taxon from the given string.
    ///
    /// # Fields
    /// A line is defined by 5 columns, separated with a tab.
    /// Note that all fields are required, in the following order:
    /// `id`,`name`,`rank`,`parent`,`valid`.
    ///
    /// The `valid`-field will be parsed as true for `"\x01"` and false for `"\x00"`.
    ///
    /// # Examples
    /// ```
    /// use unipept::taxon::Taxon;
    /// let taxon = "1\tFelis catus\tspecies\t4\t\x01".parse::<Taxon>();
    /// // Will return: Taxon {
    /// //                  id: 1,
    /// //                  name: "Felis catus",
    /// //                  rank: Rank::Species,
    /// //                  parent: 4,
    /// //                  valid: true
    /// //              }
    /// ```
    fn from_str(s: &str) -> Result<Self> {
        let split: Vec<&str> = s.trim_right().split('\t').collect();

        if split.len() != 5 { bail!("Taxon requires five fields"); }
        match (
            split[0].parse::<TaxonId>(),
            split[1].to_string(),
            split[2].parse::<Rank>(),
            split[3].parse::<TaxonId>(),
            split[4]
        ) {
            (Ok(id), name, Ok(rank), Ok(parent), "\x01") => Ok(Taxon::new(id, name, rank, parent, true)),
            (Ok(id), name, Ok(rank), Ok(parent), "\x00") => Ok(Taxon::new(id, name, rank, parent, false)),
            (Err(e), _,    _,        _,          _)      => Err(e.into()),
            (_,      _,    Err(e),   _,          _)      => Err(e.into()),
            (_,      _,    _,        Err(e),     _)      => Err(e.into()),
            _                                            => bail!("Couldn't parse the valid byte")
        }
    }
}

/// Reads in a file where each line can be parsed as a taxon.
///
/// See [Taxon::from_str()](struct.Taxon.html#method.from_str) for more details on the line format.
pub fn read_taxa_file(filename: &str) -> Result<Vec<Taxon>> {
    let file   = try!(File::open(filename).chain_err(|| "Failed opening taxon file."));
    let reader = io::BufReader::new(file);
    let mut taxa = Vec::new();
    for mline in reader.lines() {
        let line = try!(mline.map_err(|_| "Failed to read all lines."));
        taxa.push(try!(line.parse::<Taxon>()));
    }
    Ok(taxa)
}


/// A newtype definition for a (pretty dense) list of taxons by ID.
pub struct TaxonList(Vec<Option<Taxon>>);

impl TaxonList {
    /// Groups a list of taxons by their TaxonId.
    pub fn new(taxons: Vec<Taxon>) -> Self {
        // new vec, with at least the length of the current one
        let max_id = taxons.iter().map(|t: &Taxon| t.id).max().unwrap_or(0);
        let mut by_id = Vec::with_capacity(max_id + 1);
        by_id.resize(max_id + 1, None);
        for taxon in taxons {
            let id = taxon.id;
            by_id[id] = Some(taxon);
        }
        TaxonList(by_id)
    }


    /// Constructs a vector mapping a TaxonId to the id of its parent, if it has one.
    pub fn ancestry(&self) -> Vec<Option<TaxonId>> {
        self.0.iter()
            .map(|opttaxon| opttaxon.as_ref().map(|taxon| taxon.parent))
            .collect()
    }

    /// Retrieve a taxon from the taxon list by id.
    pub fn get(&self, index: TaxonId) -> Option<&Taxon> {
        if index >= self.0.len() { None }
        else { self.0[index].as_ref() }
    }
}

/// Represents a taxonomy tree. Each node is a [Taxon](struct.Taxon.html) represented by its
/// TaxonId.
pub struct TaxonTree {
    /// The root taxon
    pub root: TaxonId,
    /// A map that maps each taxon to its children
    pub children: HashMap<TaxonId, Vec<TaxonId>>,
    max: TaxonId
}

impl TaxonTree {
    /// Creates a taxon tree from the given taxons.
    pub fn new(taxons: &Vec<Taxon>) -> TaxonTree {
        let mut map = HashMap::with_capacity(taxons.len());
        let mut max = taxons[0].id;
        let mut roots: HashSet<TaxonId> = taxons.into_iter().map(|t| t.id).collect();
        for taxon in taxons {
            if taxon.id > max { max = taxon.id }
            if taxon.id == taxon.parent { continue; }
            let siblings = map.entry(taxon.parent).or_insert(Vec::new());
            siblings.push(taxon.id);
            roots.remove(&taxon.id);
        }
        if roots.len() > 1 {
            panic!("More than one root!");
        }
        TaxonTree {
            root: roots.into_iter().next().expect("There's no root!"),
            max: max,
            children: map
        }
    }

    // Takes a (mutable) vector of taxons indexed by their id, and adds the current taxon if it
    // passes the filter.
    fn with_filtered<F>(&self, mut ancestors: &mut Vec<Option<TaxonId>>, current: TaxonId, ancestor: Option<TaxonId>, filter: &F)
    where F: Fn(TaxonId) -> bool {
        let ancestor = if filter(current) { Some(current) } else { ancestor };
        ancestors[current] = ancestor;
        if let Some(children) = self.children.get(&current) {
            for child in children {
                self.with_filtered(&mut ancestors, *child, ancestor, filter);
            }
        }
    }

    /// Returns a filtered list of taxons (or more specifically, their identifiers)
    pub fn filter_ancestors<F>(&self, filter: F) -> Vec<Option<TaxonId>>
    where F: Fn(TaxonId) -> bool {
        let mut valid_ancestors = (0..self.max + 1).map(|_| None).collect();
        self.with_filtered(&mut valid_ancestors, self.root, Some(self.root), &filter);
        valid_ancestors
    }

    /// Returns the amount of children a given taxon has in this tree.
    pub fn child_count(&self, whose: TaxonId) -> usize {
        self.children.get(&whose).map(|v| v.len()).unwrap_or(0)
    }

    /// Converts a list of taxons into their respective taxon id's for this tree. Replaces each
    /// invalid (or unranked) taxon with it's first valid (and ranked) ancestor.
    ///
    /// # Arguments:
    /// * `taxons`: a vector of taxons, indexed by their TaxonId.
    /// * `ranked_only`: whether to include only taxons with a rank.
    pub fn snapping(&self, taxons: &TaxonList, ranked_only: bool) -> Vec<Option<TaxonId>> {
        self.filter_ancestors(|i: TaxonId| {
            taxons.get(i)
                .map(|t| t.valid && (!ranked_only || t.rank != Rank::NoRank))
                .unwrap_or(false)
        })
    }
}

/// The depth in the tree
pub type Depth = usize;

/// An iterator that takes a [Euler tour](https://en.wikipedia.org/wiki/Euler_tour_technique)
/// through a [TaxonTree](struct.TaxonTree.html).
pub struct EulerIterator {
    tree: TaxonTree,
    path: Vec<(TaxonId, usize, usize)>,
    current: TaxonId,
    currentn: usize,
    children: usize
}

impl EulerIterator {
    fn new(tree: TaxonTree) -> EulerIterator {
        let child_count = tree.child_count(tree.root);
        let TaxonTree { root, children, max } = tree;
        EulerIterator {
            tree: TaxonTree { root: root, children: children, max: max },
            path: Vec::new(),
            current: root,
            currentn: 0,
            children: child_count
        }
    }
}

impl Iterator for EulerIterator {
    type Item = (TaxonId, Depth);
    fn next(&mut self) -> Option<(TaxonId, Depth)> {
        if self.currentn > self.children {
            match self.path.pop() {
                None => None,
                Some((parent, currentn, children)) => {
                    self.current = parent;
                    self.currentn = currentn;
                    self.children = children;
                    self.next()
                }
            }
        } else if self.currentn == self.children {
            let current   = self.current;
            match self.path.pop() {
                None => {
                    self.current = 0;
                    self.currentn = 1;
                    self.children = 0;
                    Some((self.tree.root, 0))
                },
                Some((parent, currentn, children)) => {
                    self.current = parent;
                    self.currentn = currentn;
                    self.children = children;
                    Some((current, self.path.len() + 1))
                }
            }
        } else {
            let current = self.current;
            // there must be unvisited children, as currentn < children
            let children = self.tree.children.get(&current).unwrap();
            self.path.push((self.current, self.currentn + 1, self.children));
            self.current = children[self.currentn];
            self.currentn = 0;
            self.children = self.tree.child_count(self.current);
            Some((current, self.path.len() - 1))
        }
    }
}

impl IntoIterator for TaxonTree {
    type Item = (TaxonId, Depth);
    type IntoIter = EulerIterator;
    fn into_iter(self) -> Self::IntoIter {
        EulerIterator::new(self)
    }
}

error_chain! {
    foreign_links {
        InvalidID(std::num::ParseIntError) #[doc = "Indicates failure to parse a Taxon ID"];
    }
    errors {
        /// Unrecognized taxon rank
        UnknownRank(rank: String) {
            description("Unrecognized taxon rank")
            display("Unknown rank: {}", rank)
        }
        /// Encountered an unknown taxon ID
        UnknownTaxon(tid: TaxonId) {
            description("Encountered an unknown taxon ID")
            display("Unknown Taxon ID: {}", tid)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fixtures;

    #[test]
    fn test_taxon_parsing() {
        assert_eq!(Taxon::from_static(1, "root", Rank::NoRank, 1,  true),  "1	root	no rank	1	\x01".parse().unwrap());
        assert_eq!(Taxon::from_static(1, "root", Rank::Family, 1,  true),  "1	root	family	1	\x01".parse().unwrap());
        assert_eq!(Taxon::from_static(1, "root", Rank::NoRank, 22, true),  "1	root	no rank	22	\x01".parse().unwrap());
        assert_eq!(Taxon::from_static(1, "root", Rank::NoRank, 1,  false), "1	root	no rank	1	\x00".parse().unwrap());

        assert_matches!(*"hello world".parse::<Taxon>().unwrap_err().kind(), ErrorKind::Msg(_));
        assert_matches!(*"a	root	no_rank	1	\x00".parse::<Taxon>().unwrap_err().kind(),  ErrorKind::InvalidID(_));
        assert_matches!(*"1	root	no_rank	1	\x00".parse::<Taxon>().unwrap_err().kind(),  ErrorKind::UnknownRank(_));
        assert_matches!(*"1	root	no rank	#	\x00".parse::<Taxon>().unwrap_err().kind(),  ErrorKind::InvalidID(_));
        assert_matches!(*"1	root	no rank		\x00".parse::<Taxon>().unwrap_err().kind(),  ErrorKind::InvalidID(_));
        assert_matches!(*"1	root	no rank	7	hello".parse::<Taxon>().unwrap_err().kind(), ErrorKind::Msg(_));
    }

    #[test]
    fn test_euler_tour() {
        let euler: Vec<(TaxonId, Depth)> = fixtures::tree().into_iter().collect();
        assert_eq!(
            vec![
                (1, 0), (2, 1),
                (1, 0), (10239, 1),
                (1, 0), (12884, 1), (185751, 2),
                        (12884, 1), (185752, 2),
                        (12884, 1),
                (1, 0)
            ],
            euler
        );
    }

    #[test]
    fn test_taxon_list() {
        let list  = fixtures::taxon_list();
        let by_id = fixtures::by_id();
        assert_eq!(Some(&list[0]), by_id.get(1));
        assert_eq!(Some(&list[1]), by_id.get(2));
        assert_eq!(None,           by_id.get(3));

        let ancestry = by_id.ancestry();
        assert_eq!(Some(1),     ancestry[1]);
        assert_eq!(Some(1),     ancestry[2]);
        assert_eq!(Some(1),     ancestry[10239]);
        assert_eq!(Some(1),     ancestry[12884]);
        assert_eq!(Some(12884), ancestry[185751]);
        assert_eq!(Some(12884), ancestry[185752]);
        assert_eq!(None,        ancestry[3]);
    }
}
