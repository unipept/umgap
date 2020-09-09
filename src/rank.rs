//! Represents the rank of a [Taxon](struct.Taxon.html)
#![allow(missing_docs)]

use std::cmp::Ordering;
use strum::IntoEnumIterator;

#[cfg_attr(rustfmt, rustfmt_skip)]
#[derive(PartialEq, Eq, Debug, Clone, Copy, Display, EnumString, EnumIter)]
pub enum Rank {
    #[strum(serialize="no rank")]          NoRank,
    #[strum(serialize="superkingdom")]     Superkingdom,
    #[strum(serialize="kingdom")]          Kingdom,
    #[strum(serialize="subkingdom")]       Subkingdom,
    #[strum(serialize="superphylum")]      Superphylum,
    #[strum(serialize="phylum")]           Phylum,
    #[strum(serialize="subphylum")]        Subphylum,
    #[strum(serialize="superclass")]       Superclass,
    #[strum(serialize="class")]            Class,
    #[strum(serialize="subclass")]         Subclass,
    #[strum(serialize="infraclass")]       Infraclass,
    #[strum(serialize="superorder")]       Superorder,
    #[strum(serialize="order")]            Order,
    #[strum(serialize="suborder")]         Suborder,
    #[strum(serialize="infraorder")]       Infraorder,
    #[strum(serialize="parvorder")]        Parvorder,
    #[strum(serialize="superfamily")]      Superfamily,
    #[strum(serialize="family")]           Family,
    #[strum(serialize="subfamily")]        Subfamily,
    #[strum(serialize="tribe")]            Tribe,
    #[strum(serialize="subtribe")]         Subtribe,
    #[strum(serialize="genus")]            Genus,
    #[strum(serialize="subgenus")]         Subgenus,
    #[strum(serialize="species group")]    SpeciesGroup,
    #[strum(serialize="species subgroup")] SpeciesSubgroup,
    #[strum(serialize="species")]          Species,
    #[strum(serialize="subspecies")]       Subspecies,
    #[strum(serialize="varietas")]         Varietas,
    #[strum(serialize="forma")]            Forma,
}

pub const RANK_COUNT: usize = 29;

static RANKS: &[&str] = &[
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species group",
    "species subgroup",
    "species",
    "subspecies",
    "varietas",
    "forma",
];

impl Rank {
    pub fn index(&self) -> usize {
        *self as usize
    }

    #[cfg_attr(rustfmt, rustfmt_skip)]
    pub fn score(&self) -> Option<usize> {
        if self < &Rank::Species { Some(10) }
        else if self < &Rank::SpeciesGroup { Some(9) }
        else if self < &Rank::Genus {        Some(8) }
        else if self < &Rank::Tribe {        Some(7) }
        else if self < &Rank::Superfamily {  Some(6) }
        else if self < &Rank::Superorder {   Some(5) }
        else if self < &Rank::Superclass {   Some(4) }
        else if self < &Rank::Superphylum {  Some(3) }
        else if self < &Rank::Superkingdom { Some(2) }
        else { None }
    }

    /// Iterator over all the real ranks (NoRank is skipped)
    pub fn ranks() -> impl Iterator<Item = Rank> {
        Self::iter().skip(1)
    }

    pub fn variants() -> &'static [&'static str] {
        RANKS
    }
}

impl PartialOrd for Rank {
    fn partial_cmp(&self, other: &Rank) -> Option<Ordering> {
        if self == &Rank::NoRank || other == &Rank::NoRank {
            None
        } else {
            Some(self.index().cmp(&other.index()))
        }
    }
}

impl Ord for Rank {
    fn cmp(&self, other: &Rank) -> Ordering {
        self.index().cmp(&other.index())
    }
}
