
use std::fmt::Display;

pub struct RMQ<T: Ord + Display> {
    pub array: Vec<T>, // The original array
    pub block_min: Vec<usize>, // The position (index in array) of the minimum for each block
    pub sparse: Vec<Vec<usize>>, // RMQ of sequences of blocks. sparse[i][j] is the position of the minimum in block[i] to block[i + 2**{j+1} - 1]
    pub labels: Vec<usize> // The j'th bit of labels[i] is 1 iff j is the first position (in the block) left of i where array[j] < array[i]
}

/* clear the least significant x - 1 bits */
fn clearbits(n: usize, x: usize) -> usize { (n >> x) << x }

fn intlog2(n: usize) -> usize {
    ((n as u32).leading_zeros() ^ 31) as usize
}

impl<T: Ord + Display> RMQ<T> {
    pub fn block_min(array: &Vec<T>) -> Vec<usize> {
        array.chunks(32)
             .enumerate()
             .map(|(i, c)| c.iter().enumerate()
                            .min_by_key(|&(_, val)| val)
                            .expect("So, it has come to this.")
                            .0 + i * 32)
             .collect()
    }

    fn aggregate_minima(array: &Vec<T>, shift: usize, minima: &Vec<usize>) -> Vec<usize> {
        minima.iter().zip(minima.iter().skip(shift)).map(|(&l, &r)|
            if array[l] < array[r] { l } else { r }
        ).collect()
    }

    pub fn sparse(array: &Vec<T>, block_min: &Vec<usize>) -> Vec<Vec<usize>> {
        let length = intlog2(block_min.len());
        let mut sparse = Vec::with_capacity(length);
        sparse.push(RMQ::<T>::aggregate_minima(array, 1, block_min));
        for i in 1..length {
            let minima = RMQ::<T>::aggregate_minima(array, 1 << i, &sparse[i - 1]);
            sparse.push(minima);
        }
        sparse
    }

    pub fn labels(array: &Vec<T>) -> Vec<usize> {
        let mut gstack = Vec::with_capacity(32);
        let mut labels = Vec::with_capacity(array.len());
        for i in 0..array.len() {
            if i % 32 == 0 {
                gstack.clear();
            }
            labels.push(0);
            while !gstack.is_empty() && array[i] < array[gstack[gstack.len() - 1]] {
                gstack.pop();
            }
            if !gstack.is_empty() {
                let g = gstack[gstack.len() - 1];
                labels[i] = labels[g] | ((1 as usize) << (g % 32));
            }
            gstack.push(i);
        }
        labels
    }

    pub fn new(array: Vec<T>) -> RMQ<T> {
        let block_min = RMQ::<T>::block_min(&array);
        let sparse    = RMQ::<T>::sparse(&array, &block_min);
        let labels    = RMQ::<T>::labels(&array);
        RMQ {
            array:     array,
            block_min: block_min,
            sparse:    sparse,
            labels:    labels
        }
    }


    fn min_in_block(labels: &Vec<usize>, left: usize, right: usize) -> usize {
        let v = clearbits(labels[right], left % 32);
        if v == 0 {
            right
        } else {
            clearbits(left, 5) + (v.trailing_zeros() as usize)
        }
    }

    pub fn query(&self, start: usize, end: usize) -> usize {
        if start == end { return start; }
        let (left, right)  = if start < end { (start, end) } else { (end, start) };
        let block_diff     = (right >> 5) - (left >> 5);
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
                let l = RMQ::<T>::min_in_block(&self.labels, left, clearbits(left, 5) + 31);
                let r = RMQ::<T>::min_in_block(&self.labels, clearbits(right, 5), right);
                if self.array[l] <= self.array[r] { l } else { r }
            },
            _ => {
                let l = RMQ::<T>::min_in_block(&self.labels, left, clearbits(left, 5) + 31);
                let r = RMQ::<T>::min_in_block(&self.labels, clearbits(right, 5), right);
                let m = if block_diff == 2 {
                    self.block_min[(left >> 5) + 1]
                } else {
                    let k = intlog2(block_diff - 1) - 1;
                    let t1 = self.sparse[k][(left >> 5) + 1];
                    let t2 = self.sparse[k][(right >> 5) - (1 << (k + 1))];
                    if self.array[t1] <= self.array[t2] { t1 } else { t2 }
                };
                let ex = if self.array[l] <= self.array[m] { l } else { m };
                if self.array[ex] <= self.array[r] { ex } else { r }
            }
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_block_minima() {
        assert_eq!(
            vec![3, 33],
            RMQ::block_min(&vec![12, 17, 23, 2, 20, 4, 8, 27, 26, 19, 31, 22, 28, 16, 24, 14, 5, 29, 32, 11, 7, 9, 25, 30, 21, 13, 6, 18, 15, 33, 10, 3, /**/ 33, 1])
        );
    }

    #[test]
    fn test_sparse() {
        assert_eq!(
            vec![vec![3, 5, 15, 16, 16, 26, 31, 33],
                 vec![3, 5, 16, 16, 31, 33],
                 vec![3, 33]],
            RMQ::sparse(
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
}
