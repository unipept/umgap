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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rmq.h"

#define N 1000

int main(int argc, char * argv[]){
  INT i, j;
  VAL a[N];
  struct rmqinfo * ri;
  srandom(time());
  for(i = 0; i < N; i++)
    a[i] = random() % 10;
  ri = rm_query_preprocess(a, N);
  for(i = 0; i < N; i++){
    for(j = i; j < N; j++){
      INT x,y;
      x = rm_query_naive(a,i,j);
      y = rm_query(ri, i, j);
      // assert(a[x] == a[y]);
      assert(x == y);
    }
  }
  rm_free(ri);
  return(0);
}
