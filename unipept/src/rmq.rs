
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
                //println!("case >2 ({})", block_diff);
                let l = RMQ::<T>::min_in_block(&self.labels, left, clearbits(left, 5) + 31);
                let r = RMQ::<T>::min_in_block(&self.labels, clearbits(right, 5), right);
                let m = if block_diff == 2 {
                    //println!("length: {}", self.array.len());
                    //println!("block: {}", (left >> 5) + 1);
                    //println!("[{}, {}, {}, ...]", self.array[33], self.array[34], self.array[35]);
                    //for min in self.block_min.iter() { println!("{}", min); }
                    self.block_min[(left >> 5) + 1]
                } else {
                    let k = intlog2(block_diff - 1) - 1;
                    let t1 = self.sparse[k][(left >> 5) + 1];
                    let t2 = self.sparse[k][(right >> 5) - (1 << (k + 1))];
                    if self.array[t1] <= self.array[t2] { t1 } else { t2 }
                };
                //println!("l={} ({}), m={} ({}), r={} ({})", l, self.array[l], m, self.array[m], r, self.array[r]);
                let ex = if self.array[l] <= self.array[m] { l } else { m };
                if self.array[ex] <= self.array[r] { ex } else { r }
            }
        }
    }

}
