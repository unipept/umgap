//! Represents the rank of a [Taxon](struct.Taxon.html)
#![allow(missing_docs)]

use plain_enum::*;

use std::fmt;
use std::str::FromStr;

plain_enum_mod!{taxon, Rank {
    NoRank, Superkingdom, Kingdom, Subkingdom, Superphylum, Phylum, Subphylum, Superclass, Class,
    Subclass, Infraclass, Superorder, Order, Suborder, Infraorder, Parvorder, Superfamily, Family,
    Subfamily, Tribe, Subtribe, Genus, Subgenus, SpeciesGroup, SpeciesSubgroup, Species,
    Subspecies, Varietas, Forma,
}}

pub const RANK_SIZE : usize = Rank::SIZE;

impl Rank {
	/// Converts a rank to a usize.
	pub fn index(&self) -> usize {
		*self as usize
	}

	/// Converts a rank to its score, if it has one
    #[cfg_attr(rustfmt, rustfmt_skip)]
    pub fn score(&self) -> Option<usize> {
        if      self < &Rank::Species       { Some(10) }
        else if self < &Rank::SpeciesGroup  { Some(9) }
        else if self < &Rank::Genus         { Some(8) }
        else if self < &Rank::Tribe         { Some(7) }
        else if self < &Rank::Superfamily   { Some(6) }
        else if self < &Rank::Superorder    { Some(5) }
        else if self < &Rank::Superclass    { Some(4) }
        else if self < &Rank::Superphylum   { Some(3) }
        else if self < &Rank::Superkingdom  { Some(2) }
        else { None }
    }
}

impl FromStr for Rank {
	type Err = Error;

	#[cfg_attr(rustfmt, rustfmt_skip)]
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
	#[cfg_attr(rustfmt, rustfmt_skip)]
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

error_chain!{
	errors {
		/// Unrecognized taxon rank
		UnknownRank(rank: String) {
			description("Unrecognized taxon rank")
			display("Unknown rank: {}", rank)
		}
	}
}

//impl PartialOrd for Rank {
//	fn partial_cmp(&self, other: &Rank) -> Option<Ordering> {
//		if self == &Rank::NoRank || other == &Rank::NoRank {
//			None
//		} else {
//			Some(self.index().cmp(&other.index()))
//		}
//	}
//}
