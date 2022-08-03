#include "MPIX_interface_proposal.h"

#include <math.h>
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

/** Proposal for improvements of the MPI API related to application aware
 * hardware topology mapping.
 *
 * @author Christoph Niethammer and Rolf Rabenseifner, 
 * High Performance Computing Center Stuttgart (HLRS), 
 * University of Stuttgart, Germany, 2018.
 *
 * Copyright (c) 2019, HLRS, University of Stuttgart, Germany.
 * All rights reserved.
 *
 * This software is made available under the BSD 3-clause license ("BSD License 2.0"):
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the copyright holder nor the names of its contributors
 *       may be used to endorse or promote products derived from this software 
 *       without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/** Return list of divisors (internal routine only)
 *
 * @param[in] p number for which divisors shall be computed
 * @param[out] divisors list of divisors in descending order
 */
static int MPIX_calc_divisors(int n, int *divisors) {
  int f, j, i = 0;
  divisors[i++] = 1;
  if (n % 2 == 0) {
    divisors[i++] = 2;
    for (f = 3; f <= floor(sqrt(n)); f++) {
      if (n % f == 0) {
        divisors[i++] = f;
      }
    }

  } else {
    for (f = 3; f <= floor(sqrt(n)); f += 2) {
      if (n % f == 0) {
        divisors[i++] = f;
      }
    }
  }

  for (j = i - 1; j >= 0; j--) {
    divisors[i] = n / divisors[j];
    int tmp = divisors[j];
    divisors[j] = divisors[i];
    divisors[i] = tmp;
    i++;
  }
  return MPI_SUCCESS;
}

/** recursive factorization optimization looping over (remaining) divisors for
 * each factor (internal routine only)
 * @param[in] p number to factorize
 * @param[in] divisors list of divisors to check, entry 1 mark list end
 * @param[in] ndims number of factors
 * @param[in] i index of current dimension
 * @param[in] dim_weights weight factors for dimensions, must be sorted in
 * increasing order
 * @param[in,out] dims[ndims] current array for dimensions
 * @param[in,out] min_sum current minimal sum
 * @param[out] min_dims[ndims] computed optimal values for dimensions
 * @param[in,out] min_diff current minimal diff
 */
void MPIX_optdims(const int p, int *divisors, const int ndims, int i,
             double *dim_weights, int *dims, double *min_sum, int *min_dims,
             int *min_diff) {
  if (p == 1) {
    int k;
    for (k = i; k < ndims; k++) {
      dims[k] = 1;
    }
    i = ndims - 1; /* no need for recursion, all remaining factors computed, so
                      set to end index */
  } else {
    if (i < ndims - 1) {
      int *divisor;
      for (divisor = divisors; *divisor > 1; divisor++) {
        int k = *divisor;
        if (p % k != 0) {
          continue;
        }
        if (k < round(pow(p / k, 1. / (ndims - i - 1)))) {
          break;
        }
        dims[i] = k;
        MPIX_optdims(p / k, divisor, ndims, i + 1, dim_weights, dims, min_sum,
                min_dims, min_diff);
      }
    } else {
      dims[i] = p;
    }
  }

  double sum = 0.0;
  int min = ndims * p;
  int max = 0;
  if (i == ndims - 1) {
    int k;
    for (k = 0; k < ndims; k++) {
      sum += dim_weights[k] * dims[k];
      if (dims[k] < min) min = dims[k];
      if (dims[k] > max) max = dims[k];
#ifndef NDEBUG
      printf("%d ", dims[k]);
#endif
    }
    int diff = max - min;
#ifndef NDEBUG
    printf(" => sum: %lf, delta: %d\n", sum, diff);
#endif
    if ((sum < *min_sum) || ((sum == *min_sum) && (diff < *min_diff)) ||
        ((sum == *min_sum) && (diff == *min_diff) && (p < dims[0]))) {
      int k;
      for (k = 0; k < ndims; k++) {
        min_dims[k] = dims[k];
      }
      *min_sum = sum;
      *min_diff = diff;
    }
  }
}

/** Compute dimensions based on weights
 * @param[in] nnodes number of nodes in a grid
 * @param[in] ndims number of Cartesian dimensions
 * @param[in] dim_weights[ndims] weight factors for dimensions
 * or MPI_WEIGHTS_EQUAL
 * @param[in,out] dims[ndims] nodes in each dimension
 *
 * @todo implement fixed dimensions for dims[i] > 0
 */

int MPIX_Dims_weighted_create(int nnodes, int ndims, double *dim_weights,
                             int *dims) {

  if (ndims < 0) {
    return MPI_ERR_DIMS;
  }
  if (ndims == 0) {
    return MPI_SUCCESS;
  }
  if(nnodes == 1) {
    int i;
    for(i = 0; i < ndims; i++) {dims[i] = 1;}
    return MPI_SUCCESS;
  }


  int divisors[1344]; /* maximum number of divisors for int32 according to http://oeis.org/A066150 is 1344 */

  if(nnodes == 1) {
          int i;
	  for(i = 0; i < ndims; i++) {dims[i] = 1;}
	  return MPI_SUCCESS;
  }

  double tmp_dim_weights[ndims];
  if(dim_weights == MPIX_WEIGHTS_EQUAL) {
        int i;
	for(i = 0; i < ndims; i++) { tmp_dim_weights[i] = 1.0; }
  } else {
    int i;
    for(i = 0; i < ndims; i++) { tmp_dim_weights[i] = dim_weights[i]; }
  }

  int i;
  int permutation[ndims];
  for (i = 0; i < ndims; i++) {
    permutation[i] = i;
  }
  for (i = 0; i < ndims; i++) {
    int j;
    for (j = i + 1; j < ndims; j++) {
      if (tmp_dim_weights[i] > tmp_dim_weights[j]) {
        double tmp = tmp_dim_weights[i];
        tmp_dim_weights[i] = tmp_dim_weights[j];
        tmp_dim_weights[j] = tmp;
        int ti = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = ti;
      }
    }
  }

  double min_sum = ((double)nnodes) * ndims * tmp_dim_weights[ndims - 1];
  int min_diff = nnodes - 1;
  int min_dims[ndims];
  for(i = 0; i < ndims; i++) { min_dims[i] = 1; }
  int tmp_dims[ndims];
  for(i = 0; i < ndims; i++) { tmp_dims[i] = 1; }

  MPIX_calc_divisors(nnodes, divisors);
  MPIX_optdims(nnodes, divisors, ndims, 0, tmp_dim_weights, tmp_dims, &min_sum, min_dims, &min_diff);

  for (i = 0; i < ndims; i++) {
    dims[permutation[i]] = min_dims[i];
  }

  return MPI_SUCCESS;
}

