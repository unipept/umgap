extern crate rmq;

use rmq::*;

#[test]
fn test_block_minima() {
    assert_eq!(
        vec![3, 1],
        RMQInfo::block_min(&vec![12, 17, 23, 2, 20, 4, 8, 27, 26, 19, 31, 22, 28, 16, 24, 14, 5, 29, 32, 11, 7, 9, 25, 30, 21, 13, 6, 18, 15, 33, 10, 3, /**/ 33, 1])
    );
}

#[test]
fn test_sparse() {
    assert_eq!(
        vec![vec![3, 5, 15, 16, 16, 26, 31, 33],
             vec![3, 5, 16, 16, 31, 33],
             vec![3, 33]],
        RMQInfo::sparse(
            &vec![12, 17, 23, 2,
                  20, 4,  8,  27,
                  26, 19, 31, 22,
                  28, 16, 24, 14,
                  5,  29, 32, 11,
                  7,  9,  25, 30,
                  21, 13, 6,  18,
                  15, 33, 10, 3,
                  33, 1],
            &vec![3, 5, 9, 15, 16, 20, 26, 31, 33]
        )
    )
}

fn array() -> Vec<usize> {
    vec![
             /* 0   1   2   3   4   5   6   7   8   9 */
        /* 0 */ 39, 60, 15, 94, 25, 3,  88, 94, 71, 68,
        /* 1 */ 17, 15, 73, 32, 59, 89, 25, 36, 12, 85,
        /* 2 */ 80, 94, 56, 30, 62, 3,  10, 58, 69, 56,
        /* 3 */ 10, 8,  48, 25, 34, 5,  61, 22, 99, 64,
        /* 4 */ 22, 49, 80, 28, 13, 71, 17, 38, 40, 61,
        /* 5 */ 55, 20, 55, 43, 82, 49, 78, 24, 8,  47,
        /* 6 */ 12, 50, 87, 61, 8,  21, 66, 69, 76, 66,
        /* 7 */ 65, 98, 47, 77, 58, 60, 81, 76, 98, 21,
        /* 8 */ 69, 85, 73, 25, 29, 88, 74, 7,  12, 14,
        /* 9 */ 87, 25, 97, 74, 86, 5,  28, 84, 6,  4,
             /* 0   1   2   3   4   5   6   7   8   9 */
    ]
}

#[test]
fn test_rmq_single_block() {
    let array = array();
    let info = RMQInfo::new(&array);
    assert_eq!(5, info.query(0, 9));
    assert_eq!(18, info.query(10, 19));
}

#[test]
fn test_rmq_two_blocks() {
    let array = array();
    let info = RMQInfo::new(&array);
    assert_eq!(5, info.query(0, 39));
}

#[test]
fn test_rmq_three_blocks() {
    let array = array();
    let info = RMQInfo::new(&array);
    assert_eq!(5, info.query(0, 69));
    assert_eq!(99, info.query(40, 99));
}

#[test]
fn test_rmq_more_blocks() {
    let array = array();
    let info = RMQInfo::new(&array);
    assert_eq!(5, info.query(0, 99));
    assert_eq!(25, info.query(10, 99));
    assert_eq!(99, info.query(30, 99));
}

