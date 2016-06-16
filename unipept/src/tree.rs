
use std::ops::Add;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use taxon;
use taxon::{Taxon, TaxonId};
use agg::Aggregator;

pub struct TreeBasedAggregator {
    root: TaxonId,
    taxons: Vec<Option<Taxon>>,
    snapping: Vec<Option<TaxonId>>,
}

struct SubTree<T: Default + Copy> {
    root: TaxonId,
    value: T,
    children: Vec<SubTree<T>>,
}

impl<T: Default + Copy> SubTree<T> {
    fn new(aggregator: &TreeBasedAggregator, taxons: HashMap<TaxonId, T>) -> Result<Self, String> {
        let mut tree: HashMap<TaxonId, HashSet<TaxonId>> = HashMap::with_capacity(taxons.len());
        let mut queue: VecDeque<TaxonId> = taxons.keys().map(|t| *t).collect();
        while let Some(id) = queue.pop_front() {
            let parent = try!(aggregator.taxons[id].as_ref().ok_or("Error: unknown taxon in input.")).parent;
            if id == parent { continue; }
            if !tree.contains_key(&parent) { queue.push_back(parent); }
            let siblings = tree.entry(parent).or_insert(HashSet::new());
            siblings.insert(id);
        }
        Ok(SubTree::create(aggregator.root, &tree, &taxons))
    }

    fn create(root: TaxonId, children: &HashMap<TaxonId, HashSet<TaxonId>>, taxons: &HashMap<TaxonId, T>) -> Self {
        SubTree {
            root: root,
            value: *taxons.get(&root).unwrap_or(&T::default()),
            children: children.get(&root).map(|set| set.iter().map(|&tid| SubTree::create(tid, children, taxons)).collect()).unwrap_or(Vec::new())
        }
    }

    fn collapse<F>(&self, combine: &F) -> Self 
    where F: Fn(T, T) -> T {
        let mut value = self.value;
        let mut new   = self;
        while new.children.len() == 1 {
            new   = &new.children[0];
            value = combine(value, new.value);
        }
        SubTree {
            root: new.root,
            value: value,
            children: new.children.iter().map(|c| c.collapse(combine)).collect()
        }
    }
}

impl<T: Default + Copy + ToString> SubTree<T> {
    fn _print(&self, depth: usize) {
        let mut string = "".to_string();
        for _ in 0..depth { string = string + " " }
        println!("{}-> {} ({})", string, self.root, self.value.to_string());
        for child in self.children.iter() {
            child._print(depth + 1);
        }
    }
}

impl TreeBasedAggregator {
    pub fn new(taxons: Vec<Taxon>, ranked_only: bool) -> Self {
        let tree     = taxon::TaxonTree::new(&taxons);
        let by_id    = taxon::taxa_vector_by_id(taxons);
        let snapping = tree.snapping(&by_id, ranked_only);
        TreeBasedAggregator {
            root:     tree.root,
            taxons:   by_id,
            snapping: snapping,
        }
    }
}

impl Aggregator for TreeBasedAggregator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        // Count the occurences
        let mut counts = HashMap::new();
        for taxon in taxons {
            *counts.entry(*taxon).or_insert(0) += 1;
        }

        // Build a subtree of the taxonomy
        let subtree = SubTree::new(self, counts).unwrap().collapse(&Add::add);

        // TODO (for now, just an LCA implementation
        let lca = self.snapping[subtree.root].expect("LCA does not exist");
        self.taxons[lca].as_ref().expect("No LCA* found.")
    }
}
