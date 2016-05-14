
extern crate unipept;

use unipept::taxon::*;
use unipept::lca::*;

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
fn test_two_on_same_path() {
    let taxon_list = taxon_list();
    let calculator = LCACalculator::new(taxon_list);
    assert_eq!(185752, calculator.calc_lca(&vec![12884, 185752], false));
    assert_eq!(185752, calculator.calc_lca(&vec![185752, 12884], false));
    assert_eq!(2, calculator.calc_lca(&vec![1, 2], false));
    assert_eq!(2, calculator.calc_lca(&vec![2, 1], false));
}

#[test]
fn test_two_on_fork() {
    let taxon_list = taxon_list();
    let calculator = LCACalculator::new(taxon_list);
    assert_eq!(1, calculator.calc_lca(&vec![2, 10239], false));
    assert_eq!(1, calculator.calc_lca(&vec![10239, 2], false));
    assert_eq!(12884, calculator.calc_lca(&vec![185751, 185752], true));
    assert_eq!(12884, calculator.calc_lca(&vec![185752, 185751], true));
}

#[test]
fn test_three_on_triangle() {
    let taxon_list = taxon_list();
    let calculator = LCACalculator::new(taxon_list);
    assert_eq!(12884, calculator.calc_lca(&vec![12884, 185751, 185752], false));
    assert_eq!(12884, calculator.calc_lca(&vec![12884, 185752, 185751], false));
    assert_eq!(12884, calculator.calc_lca(&vec![185751, 12884, 185752], false));
    assert_eq!(12884, calculator.calc_lca(&vec![185752, 12884, 185751], false));
    assert_eq!(12884, calculator.calc_lca(&vec![185751, 185752, 12884], false));
    assert_eq!(12884, calculator.calc_lca(&vec![185752, 185751, 12884], false));
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
    let calculator = LCACalculator::new(taxon_list);
    assert_eq!(3, calculator.calc_lca(&vec![9, 7], false));
    assert_eq!(3, calculator.calc_lca(&vec![9, 10], false));
    assert_eq!(3, calculator.calc_lca(&vec![7, 9], false));
    assert_eq!(3, calculator.calc_lca(&vec![14, 8], false));
    assert_eq!(3, calculator.calc_lca(&vec![14, 8], false));
}

