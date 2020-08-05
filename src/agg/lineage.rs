//! A lineage-based method of Taxon aggregation

use std::iter::{Iterator, Peekable};

use crate::rank::RANK_COUNT;
use crate::taxon::{Lineage, Taxon, TaxonId, TaxonList};

/// An Iterator which yields the aggegrated pairs.
pub struct LineageAggregator<T: Iterator<Item = (String, TaxonId)>>(Peekable<LineageIterator<T>>);

impl<T: Iterator<Item = (String, TaxonId)>> LineageAggregator<T> {
    /// Aggregates the given records using a taxonomy derived from taxa.
    pub fn new(records: T, taxa: Vec<Taxon>) -> Self {
        LineageAggregator(LineageIterator::new(records.into_iter(), taxa).peekable())
    }
}

struct LineageIterator<T: Iterator<Item = (String, TaxonId)>> {
    records: T,
    taxon_list: TaxonList,
}

impl<T: Iterator<Item = (String, TaxonId)>> LineageIterator<T> {
    pub fn new(records: T, taxa: Vec<Taxon>) -> Self {
        let taxon_list = TaxonList::new_with_unknown(taxa, true);
        LineageIterator {
            records,
            taxon_list,
        }
    }
}

impl<T: Iterator<Item = (String, TaxonId)>> Iterator for LineageIterator<T> {
    type Item = (String, Lineage);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((sequence, taxon_id)) = self.records.next() {
            if let Ok(lineage) = self.taxon_list.lineage(taxon_id) {
                Some((sequence, lineage))
            } else {
                None
            }
        } else {
            None
        }
    }
}

impl<T: Iterator<Item = (String, TaxonId)>> Iterator for LineageAggregator<T> {
    type Item = (String, TaxonId);

    fn next(&mut self) -> Option<(String, TaxonId)> {
        if let Some((_sequence, _aggregate)) = self.0.next() {
            let _join_rank = RANK_COUNT;
            None
        } else {
            None
        }
    }
}
