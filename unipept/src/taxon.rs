use std::fmt;
use std::io;
use std::io::BufRead;
use std::fs::File;
use std::collections::HashMap;
use std::collections::HashSet;

use std::str::FromStr;

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum Rank {
    NoRank, Superkingdom, Kingdom, Subkingdom, Superphylum, Phylum, Subphylum, Superclass, Class,
    Subclass, Infraclass, Superorder, Order, Suborder, Infraorder, Parvorder, Superfamily, Family,
    Subfamily, Tribe, Subtribe, Genus, Subgenus, SpeciesGroup, SpeciesSubgroup, Species,
    Subspecies, Varietas, Forma
}

impl Rank {
    pub fn index(&self) -> usize {
        *self as usize
    }
}

impl FromStr for Rank {
    type Err = &'static str;
    fn from_str(rank: &str) -> Result<Self, Self::Err> {
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
            _                  => Err("Unknown taxon rank.")
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

pub type TaxonId = usize;

#[derive(Clone, PartialEq, Debug)]
pub struct Taxon {
    pub id: TaxonId,
    pub name: String,
    pub rank: Rank,
    pub parent: TaxonId,
    pub valid: bool
}

impl Taxon {
    pub fn new(id: TaxonId, name: String, rank: Rank, parent: TaxonId, valid: bool) -> Taxon {
        Taxon { id: id, name: name, rank: rank, parent: parent, valid: valid }
    }

    pub fn from_static(id: TaxonId, name: &str, rank: Rank, parent: TaxonId, valid: bool) -> Taxon {
        Taxon::new(id, name.to_string(), rank, parent, valid)
    }
}

impl FromStr for Taxon {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let split: Vec<&str> = s.trim_right().split('\t').collect();

        if split.len() != 5 { return Err("Taxon requires five fields"); }
        match (
            split[0].parse::<TaxonId>(),
            split[1].to_string(),
            split[2].parse::<Rank>(),
            split[3].parse::<TaxonId>(),
            split[4]
        ) {
            (Ok(id), name, Ok(rank), Ok(parent), "\x01") => Ok(Taxon::new(id, name, rank, parent, true)),
            (Ok(id), name, Ok(rank), Ok(parent), "\x00") => Ok(Taxon::new(id, name, rank, parent, false)),
            (Err(_), _,    _,        _,          _)      => Err("Couldn't parse the ID"),
            (_,      _,    Err(_),   _,          _)      => Err("Couldn't parse the Rank"),
            (_,      _,    _,        Err(_),     _)      => Err("Couldn't parse the parent ID"),
            _                                            => Err("Couldn't parse the valid byte")
        }
    }
}

pub fn read_taxa_file(filename: &str) -> Result<Vec<Taxon>, &'static str> {
    let file   = try!(File::open(filename).map_err(|_| "Failed opening taxon file."));
    let reader = io::BufReader::new(file);
    let mut taxa = Vec::new();
    for mline in reader.lines() {
        let line = try!(mline.map_err(|_| "Failed to read all lines."));
        taxa.push(try!(line.parse::<Taxon>()));
    }
    Ok(taxa)
}

pub fn taxa_vector_by_id(taxons: Vec<Taxon>) -> Vec<Option<Taxon>> {
    // new vec, with at least the length of the current one
    let max_id = taxons.iter().map(|t: &Taxon| t.id).max().unwrap_or(0);
    let mut by_id = Vec::with_capacity(max_id + 1);
    by_id.resize(max_id + 1, None);
    for taxon in taxons {
        let id = taxon.id;
        by_id[id] = Some(taxon);
    }
    by_id
}

pub struct TaxonTree {
    pub root: TaxonId,
    pub children: HashMap<TaxonId, Vec<TaxonId>>,
    max: TaxonId
}

impl TaxonTree {
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

    pub fn filter_ancestors<F>(&self, filter: F) -> Vec<Option<TaxonId>>
    where F: Fn(TaxonId) -> bool {
        let mut valid_ancestors = (0..self.max + 1).map(|_| None).collect();
        self.with_filtered(&mut valid_ancestors, self.root, Some(self.root), &filter);
        valid_ancestors
    }

    pub fn child_count(&self, whose: TaxonId) -> usize {
        self.children.get(&whose).map(|v| v.len()).unwrap_or(0)
    }

    pub fn snapping(&self, taxons: &Vec<Option<Taxon>>, ranked_only: bool) -> Vec<Option<TaxonId>> {
        self.filter_ancestors(|i: TaxonId| {
            taxons[i]
                .as_ref()
                .map(|t| t.valid && (!ranked_only || t.rank != Rank::NoRank))
                .unwrap_or(false)
        })
    }
}

pub type Depth = usize;

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

