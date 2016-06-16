
extern crate unipept;

use unipept::agg::Aggregator;
use unipept::tree::lca::LCACalculator;

mod fixtures;

#[test]
fn test_two_on_same_path() {
    let taxon_list = fixtures::taxon_list();
    let calculator = LCACalculator::new(taxon_list, false);
    assert_eq!(185752, calculator.aggregate(&vec![12884, 185752]).id);
    assert_eq!(185752, calculator.aggregate(&vec![185752, 12884]).id);
    assert_eq!(2, calculator.aggregate(&vec![1, 2]).id);
    assert_eq!(2, calculator.aggregate(&vec![2, 1]).id);
}

#[test]
fn test_two_on_fork() {
    let calculator = LCACalculator::new(fixtures::taxon_list(), false);
    assert_eq!(1, calculator.aggregate(&vec![2, 10239]).id);
    assert_eq!(1, calculator.aggregate(&vec![10239, 2]).id);

    let calculator = LCACalculator::new(fixtures::taxon_list(), true);
    assert_eq!(12884, calculator.aggregate(&vec![185751, 185752]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185752, 185751]).id);
}

#[test]
fn test_three_on_triangle() {
    let taxon_list = fixtures::taxon_list();
    let calculator = LCACalculator::new(taxon_list, false);
    assert_eq!(12884, calculator.aggregate(&vec![12884, 185751, 185752]).id);
    assert_eq!(12884, calculator.aggregate(&vec![12884, 185752, 185751]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185751, 12884, 185752]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185752, 12884, 185751]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185751, 185752, 12884]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185752, 185751, 12884]).id);
}
