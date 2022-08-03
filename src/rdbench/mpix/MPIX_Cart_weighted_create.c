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


/** Makes a new communicator to which Cartesian topology information has been attached.
 * @param[in]   comm_old      old communicator
 * @param[in]   ndims         number of dimensions
 * @param[in]   dim_weights[ndims]   weight factors for dimensions
 * or MPI_WEIGHTS_EQUAL
 * @param[in]   periods[ndims] logical array of size ndims specifying wether the
 * grid is periodic (true) or not (false) in each dimension
 * @param[in]   info          additional hints via MPI info
 * @param[in,out]  dims[ndims]   integer array of size ndims specifiying the number
 * of nodes in each dimension
 * @param[out]  comm_cart     communicator with new_Cartesian topology (handle)
 *
 * @return      MPI_SUCCESS if successful
 *
 * @todo implement fixed dimensions for dims[i] > 0
 */
int MPIX_Cart_weighted_create(const MPI_Comm comm_old, int ndims,
                               double *dim_weights, int *periods, const MPI_Info info,
                               int *dims, MPI_Comm *comm_cart) {
  int type_levels[] = {
    MPI_COMM_TYPE_SHARED,
#ifdef OPEN_MPI
    OMPI_COMM_TYPE_NUMA,
#else
    MPIX_COMM_TYPE_HALFNODE,
#endif
  };
  int ntype_levels = sizeof(type_levels) / sizeof(type_levels[0]);
  MPIX_Cart_ml_create_from_types(comm_old, ntype_levels, type_levels, ndims, dim_weights, periods, info, dims, comm_cart);
  
  return MPI_SUCCESS;
}
