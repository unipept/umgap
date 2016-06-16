
extern crate num_rational;
use num_rational::Ratio;

extern crate unipept;
use unipept::taxon::*;
use unipept::agg::*;
use unipept::rmq::mix::*;

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
fn test_full_rtl() {
    let taxon_list = taxon_list();
    let calculator = MixCalculator::new(taxon_list, false, Ratio::new(0, 1));
    assert_eq!(185751, calculator.aggregate(&vec![12884, 185751]).id);
    assert_eq!(185752, calculator.aggregate(&vec![12884, 185751, 185752, 185752]).id);
    assert_eq!(10239, calculator.aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]).id);
}

#[test]
fn test_full_lca() {
    let taxon_list = taxon_list();
    let calculator = MixCalculator::new(taxon_list, false, Ratio::new(1, 1));
    assert_eq!(12884, calculator.aggregate(&vec![12884, 185751]).id);
    assert_eq!(12884, calculator.aggregate(&vec![12884, 185751, 185752, 185752]).id);
    assert_eq!(1, calculator.aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]).id);
}

/* TODO: third example might fail because 12884 and 185751 have the same score. */
#[test]
#[ignore]
fn test_one_half() {
    let taxon_list = taxon_list();
    let calculator = MixCalculator::new(taxon_list, false, Ratio::new(1, 2));
    assert_eq!(185751, calculator.aggregate(&vec![12884, 185751]).id);
    assert_eq!(185751, calculator.aggregate(&vec![12884, 185751]).id);
    assert_eq!(185751, calculator.aggregate(&vec![1, 12884, 12284, 185751]).id);
}
