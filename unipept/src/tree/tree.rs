
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use taxon::{Taxon, TaxonId};

pub struct SubTree<T: Default + Copy> {
    pub root: TaxonId,
    pub value: T,
    pub children: Vec<SubTree<T>>,
}

impl<T: Default + Copy> SubTree<T> {
    pub fn new(root: TaxonId, taxonomy: &Vec<Option<Taxon>>, taxons: HashMap<TaxonId, T>) -> Result<Self, String> {
        let mut tree: HashMap<TaxonId, HashSet<TaxonId>> = HashMap::with_capacity(taxons.len());
        let mut queue: VecDeque<TaxonId> = taxons.keys().map(|t| *t).collect();
        while let Some(id) = queue.pop_front() {
            let parent = try!(taxonomy[id].as_ref().ok_or("Error: unknown taxon in input.")).parent;
            if id == parent { continue; }
            if !tree.contains_key(&parent) { queue.push_back(parent); }
            let siblings = tree.entry(parent).or_insert(HashSet::new());
            siblings.insert(id);
        }
        Ok(SubTree::create(root, &tree, &taxons))
    }

    fn create(root: TaxonId, children: &HashMap<TaxonId, HashSet<TaxonId>>, taxons: &HashMap<TaxonId, T>) -> Self {
        SubTree {
            root: root,
            value: *taxons.get(&root).unwrap_or(&T::default()),
            children: children.get(&root).map(|set| set.iter().map(|&tid| SubTree::create(tid, children, taxons)).collect()).unwrap_or(Vec::new())
        }
    }

    pub fn collapse<F>(&self, combine: &F) -> Self 
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

    pub fn aggregate<F>(&self, combine: &F) -> Self
    where F: Fn(T, T) -> T {
        let children: Vec<SubTree<T>> = self.children.iter().map(|c| c.aggregate(combine)).collect();
        let value                     = children.iter().map(|c| c.value).fold(self.value, combine);
        SubTree {
            root: self.root,
            value: value,
            children: children
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

