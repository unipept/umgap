
extern crate unipept;

use unipept::taxon::*;

mod fixtures;

use fixtures::taxon_list;

#[test]
fn test_taxon_parsing() {
    assert_eq!(Ok(Taxon::from_static(1, "root", Rank::NoRank, 1,  true)),  "1	root	no rank	1	\x01".parse());
    assert_eq!(Ok(Taxon::from_static(1, "root", Rank::Family, 1,  true)),  "1	root	family	1	\x01".parse());
    assert_eq!(Ok(Taxon::from_static(1, "root", Rank::NoRank, 22, true)),  "1	root	no rank	22	\x01".parse());
    assert_eq!(Ok(Taxon::from_static(1, "root", Rank::NoRank, 1,  false)), "1	root	no rank	1	\x00".parse());

    assert_eq!(Err("Taxon requires five fields"),    "hello world".parse::<Taxon>());
    assert_eq!(Err("Couldn't parse the ID"),         "a	root	no_rank	1	\x00".parse::<Taxon>());
    assert_eq!(Err("Couldn't parse the Rank"),       "1	root	no_rank	1	\x00".parse::<Taxon>());
    assert_eq!(Err("Couldn't parse the parent ID"),  "1	root	no rank	#	\x00".parse::<Taxon>());
    assert_eq!(Err("Couldn't parse the parent ID"),  "1	root	no rank		\x00".parse::<Taxon>());
    assert_eq!(Err("Couldn't parse the valid byte"), "1	root	no rank	7	hello".parse::<Taxon>());
}

#[test]
fn test_euler_tour() {
    let taxon_list = taxon_list();
    let tree = TaxonTree::new(&taxon_list);
    let euler: Vec<(TaxonId, Depth)> = tree.into_iter().collect();
    assert_eq!(
        vec![
            (1, 0), (2, 1),
            (1, 0), (10239, 1),
            (1, 0), (12884, 1), (185751, 2),
                    (12884, 1), (185752, 2),
                    (12884, 1),
            (1, 0)
        ],
        euler
    );
}
