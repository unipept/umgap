//! Defines the possible errors in Unipept.

use std::io;
use std::num;

use csv;

use fst;

use regex;

use crate::agg;
use crate::dna::translation;
use crate::taxon;

error_chain! {
	links {
		Taxon(taxon::Error, taxon::ErrorKind) #[doc = "Error propagated from Taxon"];
		Translation(translation::Error, translation::ErrorKind) #[doc = "Error propagated from Translation"];
		Aggregation(agg::Error, agg::ErrorKind) #[doc = "Error propagated from Aggregation"];
	}
	foreign_links {
		Csv(csv::Error) #[doc = "CSV"];
		Fst(fst::Error) #[doc = "Fst"];
		Io(io::Error) #[doc = "IO"];
		ParseI(num::ParseIntError) #[doc = "Parse Integer"];
		ParseF(num::ParseFloatError) #[doc = "Parse Float"];
		Regex(regex::Error) #[doc = "Regex"];
	}
	errors {
		/// Invalid invocation
		InvalidInvocation(message: String) {
			description("Invalid invocation")
			display("Invalid invocation: {}", message)
		}
	}
}
