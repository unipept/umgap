
extern crate unipept;

use unipept::agg::*;
use unipept::rtl::*;

mod fixtures;

use fixtures::taxon_list;

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


