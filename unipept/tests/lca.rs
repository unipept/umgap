
extern crate unipept;

use unipept::taxon::*;
use unipept::agg::*;
use unipept::rmq::lca::*;

mod fixtures;

use fixtures::taxon_list;

#[test]
fn test_two_on_same_path() {
    let taxon_list = taxon_list();
    let calculator = LCACalculator::new(taxon_list, false);
    assert_eq!(185752, calculator.aggregate(&vec![12884, 185752]).id);
    assert_eq!(185752, calculator.aggregate(&vec![185752, 12884]).id);
    assert_eq!(2, calculator.aggregate(&vec![1, 2]).id);
    assert_eq!(2, calculator.aggregate(&vec![2, 1]).id);
}

#[test]
fn test_two_on_fork() {
    let calculator = LCACalculator::new(taxon_list(), false);
    assert_eq!(1, calculator.aggregate(&vec![2, 10239]).id);
    assert_eq!(1, calculator.aggregate(&vec![10239, 2]).id);

    let calculator = LCACalculator::new(taxon_list(), true);
    assert_eq!(12884, calculator.aggregate(&vec![185751, 185752]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185752, 185751]).id);
}

#[test]
fn test_three_on_triangle() {
    let taxon_list = taxon_list();
    let calculator = LCACalculator::new(taxon_list, false);
    assert_eq!(12884, calculator.aggregate(&vec![12884, 185751, 185752]).id);
    assert_eq!(12884, calculator.aggregate(&vec![12884, 185752, 185751]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185751, 12884, 185752]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185752, 12884, 185751]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185751, 185752, 12884]).id);
    assert_eq!(12884, calculator.aggregate(&vec![185752, 185751, 12884]).id);
}

fn taxon(id: TaxonId, parent: TaxonId) -> Taxon {
    Taxon::from_static(id, "", Rank::NoRank, parent, true)
}

fn large_taxon_list() -> Vec<Taxon> {
    vec![
        taxon(1, 1),
          taxon(2, 1),
            taxon(5, 2),
            taxon(6, 2),
          taxon(3, 1),
            taxon(7, 3),
              taxon(10, 7),
                taxon(13, 10),
                  taxon(14, 13),
              taxon(15, 3),
            taxon(8, 3),
              taxon(11, 8),
              taxon(12, 8),
            taxon(9, 3),
          taxon(4, 1)
    ]
}

#[test]
fn test_with_deeper_interns() {
    let taxon_list = large_taxon_list();
    let calculator = LCACalculator::new(taxon_list, false);
    assert_eq!(3, calculator.aggregate(&vec![9, 7]).id);
    assert_eq!(3, calculator.aggregate(&vec![9, 10]).id);
    assert_eq!(3, calculator.aggregate(&vec![7, 9]).id);
    assert_eq!(3, calculator.aggregate(&vec![14, 8]).id);
    assert_eq!(3, calculator.aggregate(&vec![14, 8]).id);
}

