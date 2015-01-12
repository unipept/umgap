/*  
 * Implementation of a simple range minimum query algorithm described in
 * S. Alstrup, C. Gavoille, H. Kaplan, T. Rauhe.
 * Nearest common ancestors: a survey and a new distributed algorithm,
 * In Proc. 14th annual ACM symposium on Parallel algorithms and architectures, 258-264, 2002.
 *
 * Copyright (C) 2005 Hideo Bannai (http://tlas.i.kyushu-u.ac.jp/~bannai/)
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef RMQ_H_
#define RMQ_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Type for array index.
 * Assumes at leat 32 bits.
 * Hasn't been tested with 64 bits */
typedef unsigned int INT;

/* Type for values in the array. */
typedef unsigned int VAL;	

/* compare the value inside the array */
#define VAL_LT(x,y) x < y

/* a struct to hold preprocessed information */

struct rmqinfo {
  INT alen;     // length of original array
  VAL * array;  // pointer to original array
  INT ** sparse;
  INT * block_min;
  INT * labels;
};

/*
 * Return the position in array which gives the minimum value
 * in the subarray a[x..y] using a naive algorithm
 * When there are multiple positions with the same minimum value,
 * in the range, the smallest index is returned.
 */
INT rm_query_naive(VAL * a, INT x, INT y);

/*
 * Preprocess array in linear time so that range minimum
 * queries can be conducted in constant time.
 */
struct rmqinfo * rm_query_preprocess(VAL * a, INT alen);


/*
 * Return the position in array which gives the minimum value 
 * in the subarray rmqinfo.array[x..y] using preprocessed information.
 * When there are multiple positions with the same minimum value,
 * in the range, the smallest index is returned.
 */
INT rm_query(const struct rmqinfo * info, INT x, INT y);

/*
 * Free all memory associated with rmqinfo EXCEPT the original array
 */
void rm_free(struct rmqinfo * info);

#ifdef __cplusplus
}
#endif

#endif /*RMQ_H_*/
