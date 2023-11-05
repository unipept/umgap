//! Implements aggregation operations using
//! [Range Minimum Query](https://en.wikipedia.org/wiki/Range_minimum_query) (RMQ)

pub mod lca;
pub mod mix;
pub mod rtl;

use std::fmt::Display;
use std::mem::size_of;

/// Represents a Range Minimum Query (RMQ), which can efficiently return the minimal value in a
/// given range of an array.
pub struct RMQ<T: Ord + Display> {
    /// The full array.
    pub array: Vec<T>,
    /// The absolute position (i.e. the index in array) of the minimum for each block.
    pub block_min: Vec<usize>,
    /// `sparse[i][j]` is the position of the minimum in `block[i]` to `block[i + 2^(j+1) - 1]`.
    pub sparse: Vec<Vec<usize>>,
    /// The j'th bit of `labels[i]` is 1 iff j is the first position (in the block) left of i where
    /// `array[j] < array[i]`.
    pub labels: Vec<usize>,
}

/* clear the least significant x - 1 bits */
fn clearbits(n: usize, x: usize) -> usize {
    (n >> x) << x
}

fn size() -> usize {
    size_of::<usize>() * 8
}

fn intlog2(n: usize) -> usize {
    (n.leading_zeros() as usize) ^ (size() - 1)
}

impl<T: Ord + Display> RMQ<T> {
    /// Constructs an RMQ for the given array.
    pub fn new(array: Vec<T>) -> RMQ<T> {
        let block_min = RMQ::<T>::block_min(&array);
        let sparse = RMQ::<T>::sparse(&array, &block_min);
        let labels = RMQ::<T>::labels(&array);
        RMQ {
            array,
            block_min,
            sparse,
            labels,
        }
    }

    /// Calculates the position of each block's minimum.
    pub fn block_min(array: &[T]) -> Vec<usize> {
        array
            .chunks(size())
            .enumerate()
            .map(|(i, c)| {
                c.iter()
                    .enumerate()
                    .min_by_key(|&(_, val)| val)
                    .expect("So, it has come to this.")
                    .0
                    + i * size()
            })
            .collect()
    }

    fn aggregate_minima(array: &[T], shift: usize, minima: &[usize]) -> Vec<usize> {
        minima
            .iter()
            .zip(minima.iter().skip(shift))
            .map(|(&l, &r)| if array[l] < array[r] { l } else { r })
            .collect()
    }

    /// Calculate the values of the sparse field.
    pub fn sparse(array: &[T], block_min: &[usize]) -> Vec<Vec<usize>> {
        let length = intlog2(block_min.len());
        let mut sparse = Vec::with_capacity(length);
        sparse.push(RMQ::<T>::aggregate_minima(array, 1, block_min));
        for i in 1..length {
            let minima = RMQ::<T>::aggregate_minima(array, 1 << i, &sparse[i - 1]);
            sparse.push(minima);
        }
        sparse
    }

    /// Calculate the values of the label field.
    pub fn labels(array: &[T]) -> Vec<usize> {
        let mut gstack = Vec::with_capacity(size());
        let mut labels = Vec::with_capacity(array.len());
        for i in 0..array.len() {
            if i % size() == 0 {
                gstack.clear();
            }
            labels.push(0);
            while !gstack.is_empty() && array[i] < array[gstack[gstack.len() - 1]] {
                gstack.pop();
            }
            if !gstack.is_empty() {
                let g = gstack[gstack.len() - 1];
                labels[i] = labels[g] | ((1_usize) << (g % size()));
            }
            gstack.push(i);
        }
        labels
    }

    /// Returns the position of the minimal value in a given block
    fn min_in_block(labels: &[usize], left: usize, right: usize) -> usize {
        let v = clearbits(labels[right], left % size());
        if v == 0 {
            right
        } else {
            clearbits(left, intlog2(size())) + (v.trailing_zeros() as usize)
        }
    }

    /// Returns the position of the minimal value in a given sublist.
    #[rustfmt::skip]
    pub fn query(&self, start: usize, end: usize) -> usize {
        if start == end { return start; }
        let (left, right)  = if start < end { (start, end) } else { (end, start) };
        let (log2, size)   = (intlog2(size()), size());
        let block_diff     = (right >> log2) - (left >> log2);
        match block_diff {
            0 => {
                /* one inblock query, in left_block from (l % 32) to (r % 32) */
                RMQ::<T>::min_in_block(&self.labels, left, right)
            },
            1 => {
                /* two inblock queries:
                 *   - in left_block from (l % 32) to 31
                 *   - in right_block from 0 to (r % 32)
                 * minimum is the minimum of these two
                 */
                let l = RMQ::<T>::min_in_block(&self.labels, left, clearbits(left, log2) + size - 1);
                let r = RMQ::<T>::min_in_block(&self.labels, clearbits(right, log2), right);
                if self.array[l] <= self.array[r] { l } else { r }
            },
            _ => {
                let l = RMQ::<T>::min_in_block(&self.labels, left, clearbits(left, log2) + size - 1);
                let r = RMQ::<T>::min_in_block(&self.labels, clearbits(right, log2), right);
                let m = if block_diff == 2 {
                    self.block_min[(left >> log2) + 1]
                } else {
                    let k = intlog2(block_diff - 1) - 1;
                    let t1 = self.sparse[k][(left >> log2) + 1];
                    let t2 = self.sparse[k][(right >> log2) - (1 << (k + 1))];
                    if self.array[t1] <= self.array[t2] { t1 } else { t2 }
                };
                let ex = if self.array[l] <= self.array[m] { l } else { m };
                if self.array[ex] <= self.array[r] { ex } else { r }
            }
        }
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn test_block_minima() {
        assert_eq!(
            if size() == 32 { vec![3, 33] } else { vec![33] },
            RMQ::block_min(&[12, 17, 23, 2, 20, 4, 8, 27, 26, 19, 31, 22, 28, 16, 24, 14, 5, 29, 32, 11, 7, 9, 25, 30, 21, 13, 6, 18, 15, 33, 10, 3, /**/ 33, 1])
        );
    }

    #[test]
    fn test_sparse() {
        assert_eq!(
            vec![vec![3, 5, 15, 16, 16, 26, 31, 33],
                 vec![3, 5, 16, 16, 31, 33],
                 vec![3, 33]],
            RMQ::sparse(
                &[12, 17, 23, 2,
                  20, 4,  8,  27,
                  26, 19, 31, 22,
                  28, 16, 24, 14,
                  5,  29, 32, 11,
                  7,  9,  25, 30,
                  21, 13, 6,  18,
                  15, 33, 10, 3,
                  33, 1],
                &[3, 5, 9, 15, 16, 20, 26, 31, 33]
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
        let info = RMQ::new(array);
        assert_eq!(5, info.query(0, 9));
        assert_eq!(18, info.query(10, 19));
    }

    #[test]
    fn test_rmq_two_blocks() {
        let array = array();
        let info = RMQ::new(array);
        assert_eq!(5, info.query(0, 39));
    }

    #[test]
    fn test_rmq_three_blocks() {
        let array = array();
        let info = RMQ::new(array);
        assert_eq!(5, info.query(0, 69));
        assert_eq!(99, info.query(40, 99));
    }

    #[test]
    fn test_rmq_more_blocks() {
        let array = array();
        let info = RMQ::new(array);
        assert_eq!(5, info.query(0, 99));
        assert_eq!(25, info.query(10, 99));
        assert_eq!(99, info.query(30, 99));
        assert_eq!(105, info.query(30, 140));
    }

    #[test]
    fn test_wave_of_33() {
        let array = vec![
            1, 2, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 2, 1
        ];
        let info = RMQ::new(array);
        assert_eq!(2, info.query(2, 64));
    }

    #[test]
    fn test_wave_of_65() {
        let array = vec![
            1, 2, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
            3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 2, 1
        ];
        let info = RMQ::new(array);
        assert_eq!(2, info.query(2, 64));
    }
}
