//! A ranked-based method of Taxon aggregation

use std::cmp::min;
use std::iter::{Iterator, Peekable};

use crate::rank::Rank;
use crate::taxon::{Taxon, TaxonId, TaxonList, TaxonTree};

/// An Iterator which yields the aggregated pairs.
pub struct RankAggregator<T: Iterator<Item = (String, TaxonId)>> {
	records: Peekable<T>,
	ancestors: Vec<Option<TaxonId>>,
	ranks: Vec<Option<Rank>>,
}

impl<T: Iterator<Item = (String, TaxonId)>> RankAggregator<T> {
	/// Uses a taxonomy constructed from taxa to aggregate. Assumes the
	/// taxon with TaxonId 0 is absent or represents an "unknown" taxon.
	pub fn new(records: T, taxa: Vec<Taxon>) -> Self {
		let taxon_tree = TaxonTree::new(&taxa);
		let taxon_list = TaxonList::new_with_unknown(taxa, true);
		RankAggregator { records: records.peekable(),
		                 ancestors: taxon_tree.snapping(&taxon_list, true),
		                 ranks: taxon_list.0
		                                  .iter()
		                                  .map(|mt| mt.as_ref().map(|t| t.rank))
		                                  .collect() }
	}

	fn raise_to_rank(&self, taxon: TaxonId, target: Rank) -> Option<TaxonId> {
		let mut ancestor = Some(taxon);
		while ancestor.map_or(false, |a| self.ranks[a].map_or(true, |r| target < r)) {
			ancestor = ancestor.and_then(|a| self.ancestors[a]);
		}
		ancestor
	}

	fn with_rank(&self, taxon: TaxonId) -> Option<(TaxonId, Rank)> {
		if let Some(rank) = self.ranks[taxon] {
			Some((taxon, rank))
		} else if let Some(ancestor) = self.ancestors[taxon] {
			self.ranks[ancestor].map(|r| (ancestor, r))
		} else {
			None
		}
	}
}

impl<T: Iterator<Item = (String, TaxonId)>> Iterator for RankAggregator<T> {
	type Item = (String, TaxonId);

	fn next(&mut self) -> Option<(String, TaxonId)> {
		if let Some((sequence, initial_tid)) = self.records.next() {
			let mut join_rank: Option<Rank> = None;
			let (mut aggregate, mut aggregate_rank) =
				self.with_rank(initial_tid).expect("reeeeeeee");

			while self.records.peek().map_or(false, |(s, _)| *s == sequence) {
				let (_, next) = self.records.next().expect("we just peeked at it");
				let (next_taxon, next_rank) = self.with_rank(next).expect("reeeeeee");

				let compare_rank = min(next_rank, join_rank.unwrap_or(aggregate_rank));
				let mut raised_aggregate = self.raise_to_rank(aggregate, compare_rank)
				                               .expect("none if not?");
				let mut raised_next = self.raise_to_rank(next_taxon, compare_rank)
				                          .expect("none if not?");
				if raised_aggregate != raised_next {
					// samen stijgen tot gelijk, dit is nieuwe aggregate en join rank
					while raised_aggregate != raised_next {
						raised_aggregate = self.ancestors[raised_aggregate].expect("eeeeeeeeee");
						raised_next = self.ancestors[raised_next].expect("eeeeeeeeee");
					}
					aggregate = raised_aggregate;
					aggregate_rank = self.ranks[aggregate].expect("reeeeeeeeeeeeeeeee");
					join_rank = Some(aggregate_rank);
				} else if join_rank.is_none() && compare_rank != next_rank {
					aggregate = next_taxon; // && aggregate rank updaten
					aggregate_rank = next_rank;
				}
			}
			Some((sequence, aggregate))
		} else {
			None
		}
	}
}
