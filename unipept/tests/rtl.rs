
extern crate unipept;

use unipept::taxon::*;
use unipept::agg::*;
use unipept::rtl::*;

fn taxon_list() -> Vec<Taxon> {
    vec![
        Taxon::from_static(1,      "root",          Rank::NoRank,       1,     true),
        Taxon::from_static(2,      "Bacteria",      Rank::Superkingdom, 1,     true),
        Taxon::from_static(10239,  "Viruses",       Rank::Superkingdom, 1,     true),
        Taxon::from_static(12884,  "Viroids",       Rank::Superkingdom, 1,     true),
        Taxon::from_static(185751, "Pospiviroidae", Rank::Family,       12884, true),
        Taxon::from_static(185752, "Avsunviroidae", Rank::Family,       12884, true),
    ]
}

#[test]
fn test_all_on_same_path() {
    let taxon_list = taxon_list();
    let calculator = RTLCalculator::new(taxon_list, false);
    assert_eq!(1, calculator.aggregate(&vec![1]).id);
    assert_eq!(12884, calculator.aggregate(&vec![1, 12884]).id);
    assert_eq!(185751, calculator.aggregate(&vec![1, 12884, 185751]).id);
}

#[test]
fn favouring_root() {
    let taxon_list = taxon_list();
    let calculator = RTLCalculator::new(taxon_list, false);
    assert_eq!(185751, calculator.aggregate(&vec![1, 1, 1, 185751, 1, 1]).id);
}

#[test]
fn leaning_close() {
    let taxon_list = taxon_list();
    let calculator = RTLCalculator::new(taxon_list, false);
    assert_eq!(185751, calculator.aggregate(&vec![1, 1, 185752, 185751, 185751, 1]).id);
}

#[test]
fn non_deterministic() {
    let taxon_list = taxon_list();
    let calculator = RTLCalculator::new(taxon_list, false);
    assert!(vec![185751, 185752].contains(&calculator.aggregate(&vec![1, 1, 185752, 185751, 1]).id))
}


