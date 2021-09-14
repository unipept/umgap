//! Defines aggregation operations using only the tree structure.

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use crate::taxon;
use crate::taxon::TaxonId;

/// A recursive tree of TaxonId's and a label.
pub struct Tree<T: Default + Copy> {
    /// The root of the (sub)tree.
    pub root: TaxonId,
    /// The label of this node.
    pub value: T,
    /// The children of this node.
    pub children: Vec<Tree<T>>,
}

impl<T: Default + Copy> Tree<T> {
    /// Create a new Tree.
    ///
    /// `root` will be the root taxon of the tree. In `parents`, on index i should be the parent of
    /// taxon with id i, if any. `taxons` links taxa to their labels. Only taxon ids included in
    /// `taxons` and their ancestors will be in the tree.
    pub fn new(
        root: TaxonId,
        parents: &[Option<TaxonId>],
        taxons: &HashMap<TaxonId, T>,
    ) -> taxon::Result<Self> {
        let mut tree: HashMap<TaxonId, HashSet<TaxonId>> = HashMap::with_capacity(taxons.len());
        let mut queue: VecDeque<TaxonId> = taxons.keys().copied().collect();
        while let Some(id) = queue.pop_front() {
            let parent = parents[id].ok_or(taxon::ErrorKind::UnknownTaxon(id))?;
            if id == parent {
                continue;
            }
            if !tree.contains_key(&parent) {
                queue.push_back(parent);
            }
            let siblings = tree.entry(parent).or_insert_with(HashSet::new);
            siblings.insert(id);
        }
        Ok(Tree::create(root, &tree, taxons))
    }

    fn create(
        root: TaxonId,
        children: &HashMap<TaxonId, HashSet<TaxonId>>,
        taxons: &HashMap<TaxonId, T>,
    ) -> Self {
        Tree {
            root,
            value: *taxons.get(&root).unwrap_or(&T::default()),
            children: children
                .get(&root)
                .map(|set| {
                    set.iter()
                        .map(|&tid| Tree::create(tid, children, taxons))
                        .collect()
                })
                .unwrap_or_else(Vec::new),
        }
    }

    /// Collapses the label of a parent with a single child into that child using the `combine`
    /// function, as long as there are parents with single children.
    pub fn collapse<F>(&self, combine: &F) -> Self
    where
        F: Fn(T, T) -> T,
    {
        let mut value = self.value;
        let mut new = self;
        while new.children.len() == 1 {
            new = &new.children[0];
            value = combine(value, new.value);
        }
        Tree {
            root: new.root,
            value,
            children: new.children.iter().map(|c| c.collapse(combine)).collect(),
        }
    }

    /// Replaces every label in the tree with the aggregate label of all their descendants, with
    /// given `combine` function. The aggregation of children happens in pre-order.
    pub fn aggregate<F>(&self, combine: &F) -> Self
    where
        F: Fn(T, T) -> T,
    {
        let children: Vec<Tree<T>> = self.children.iter().map(|c| c.aggregate(combine)).collect();
        let value = children.iter().map(|c| c.value).fold(self.value, combine);
        Tree {
            root: self.root,
            value,
            children,
        }
    }
}

impl<T: Default + Copy + ToString> Tree<T> {
    fn _print(&self, depth: usize) {
        let mut string = "".to_string();
        for _ in 0..depth {
            string += " "
        }
        println!("{}-> {} ({})", string, self.root, self.value.to_string());
        for child in self.children.iter() {
            child._print(depth + 1);
        }
    }
}
