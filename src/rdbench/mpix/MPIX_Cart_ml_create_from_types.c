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


/** Create a Cartesian communicator, which takes into account hardware topology
 * levels provided as an array of communicators.
 * @param[in]   comm_old      old communicator, will be duplicated into new_comms[0]
 * @param[in]   ntype_levels  number of levels to split the communicator into
 * @param[in]   type_levels[ntype_levels]   ntype_levels keys as for MPI_Comm_split_type
 * @param[in]   ndims         number of dimensions
 * @param[in]   dim_weights[ndims]   weight factors for dimensions
 * or MPI_WEIGHTS_EQUAL
 * @param[in]   periods[ndims]       logical array of size ndims specifying wether the
 *                            grid is periodic (true) or not (false) in each dimension
 * @param[in]   info          additional hints via MPI info
 * @param[out]  dims[ndims]          integer array of size ndims specifiying the number
 *                            of nodes in each dimension
 * @param[out]  comm_cart     communicator with new_Cartesian topology (handle)
 *
 * @return      MPI_SUCCESS if successful, otherwise state of dims[] and comm_cart
 * is undefined
 *
 * @todo implement fixed dimensions for dims[i] > 0
 */
int MPIX_Cart_ml_create_from_types(const MPI_Comm comm_old, const int ntype_levels,
                                  const int *type_levels,
                                  int ndims, double *dim_weights,
                                  int *periods, const MPI_Info info,
                                  int *dims, MPI_Comm *comm_cart) {
  int nlevels, level;
  nlevels = ntype_levels+1;
  MPI_Comm level_comms[nlevels];
  level_comms[0] = comm_old;
  for(level = 0; level < ntype_levels; ++level) {
    if (type_levels[level] == MPIX_COMM_TYPE_HALFNODE) {
      int size, my_rank, color;
      MPI_Comm_size(level_comms[level], &size);
      MPI_Comm_rank(level_comms[level], &my_rank);
      if ((size % 2) == 0) {
        color = my_rank / (size/2);
      } else {
        color = 0;
      }
      MPI_Comm_split(level_comms[level], color, 0, &level_comms[level+1]);
    } else { 
      MPI_Comm_split_type(level_comms[level], type_levels[level], 0, info, &level_comms[level+1]);
    }
  }

  MPIX_Cart_ml_create_from_comms(nlevels, level_comms, ndims, dim_weights, periods, info, dims, comm_cart);

  for(level = 0; level < ntype_levels; ++level) {
    MPI_Comm_free(&level_comms[level+1]);
  }
  
  return MPI_SUCCESS;
}
