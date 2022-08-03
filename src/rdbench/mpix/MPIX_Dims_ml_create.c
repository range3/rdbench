#include <mpi.h>

#include "MPIX_interface_proposal.h"

/** Proposal for improvements of the MPI API related to application aware
 * hardware topology mapping.
 *
 * @author Christoph Niethammer and Rolf Rabenseifner, 
 * High Performance Computing Center Stuttgart (HLRS), 
 * University of Stuttgart, Germany, 2018.
 *
 * Copyright (c) 2019, HLRS, University of Stuttgart, Germany.
 * Copyright 2022 range3 ( https://github.com/range3 )
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


/** Compute multi-level dimensions based on weights
 * @param[in] nnodes number of nodes
 * @param[in] ndims number of Cartesian dimensions
 * @param[in] dim_weights[ndims] weight factors for dimensions
 * or MPI_WEIGHTS_EQUAL
 * @param[in] nlevels number of hardware topology levels
 * @param[in] level_sizes[nlevels] size of each hardware topology level, from coarsest to
 * finest. The product of all level_sizes must be equal to the number of processes
 * nnodes. Example:  nlevels = 3 and level_sizes[0] == number of nodes,
 * level_sizes[1] == number of sockets per node, level_sizes[2] == number of
 * cores per socket.
 * @param[out] dims_ml[ndims][nlevels] computed optimal values for dimensions in the nlevels
 * hardware topology levels and ndims application dimensions. For each level l,
 * the product of dims_ml[*][l] == level_sizes[l]. For each dimension d the product
 * dims_ml[d][*] == the i-th dimension of a Cartesian communicator when dims_ml
 * is used as input for MPI_Cart_multilevel_create.
 *
 * Advice to users: Optimization expects that coarse grained levels may have
 * slower communication capabilities than fine grained levels.  The input array
 * level_sizes can be obtained from MPI_Comm_multilevel_split_type, in which
 * case the order of the output dims_ml is corresponding and suitable for use
 * with MPI_Cart_multilevel_create [advice to implemetors!?]
 *
 * @todo implement fixed dimensions for dims[i] > 0
 */
int MPIX_Dims_ml_create(int nnodes, int ndims, double *dim_weights,
                               int nlevels, int *level_sizes,
                               int *dims_ml) {
  int i, l, level, d, prod = 1;
  for(l = 0; l < nlevels; l++) { prod *= level_sizes[l]; }
  if(prod != nnodes) {
    return MPI_ERR_DIMS;
  }

  double tmp_dim_weights[ndims];
  if(dim_weights == MPIX_WEIGHTS_EQUAL) {
	for(i = 0; i < ndims; i++) { tmp_dim_weights[i] = 1.0; }
  } else {
    for(i = 0; i < ndims; i++) { tmp_dim_weights[i] = dim_weights[i]; }
  }

  for(level = 0; level < nlevels; level++) {
    int p = level_sizes[level];
    int dims[ndims];
    MPIX_Dims_weighted_create(p, ndims, tmp_dim_weights, dims);
    for(d = 0; d < ndims; d++) {
      dims_ml[d * nlevels + level] = dims[d];
      tmp_dim_weights[d] = tmp_dim_weights[d] * dims[d];
    }
  }
  return MPI_SUCCESS;
}
