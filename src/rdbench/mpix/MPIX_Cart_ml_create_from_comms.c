#include <mpi.h>

#include "MPIX_interface_proposal.h"

#ifdef MPIX_debug
#include <stdio.h>
#include <stdlib.h>
#endif

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


/** Create a Cartesian communicator, which takes into account hardware topology
 * levels provided as an array of communicators.
 * @param[in]   nlevels       number of hardware topology levels
 * @param[in]   level_comms[nlevels]   Array of communicators for the nlevels hardware
 * topology levels.  The communicator comms[0] includes all processes for the
 * new Cartesian communicator. Each communicator comms[i+1] must be a subset of
 * comms[i].  The array comms can be created with the
 * MPI_Comm_multilevel_split_type function.
 * @param[in]   ndims         number of dimensions
 * @param[in]   dim_weights[ndims]   weight factors for dimensions
 * or MPI_WEIGHTS_EQUAL
 * @param[in]   periods[ndims]       logical array of size ndims specifying wether the
 * grid is periodic (true) or not (false) in each dimension
 * @param[in]   info          additional hints via MPI info
 * @param[out]  dims[ndims]   integer array of size ndims specifiying the number
 * of nodes in each dimension
 * @param[out]  comm_cart     communicator with new_Cartesian topology (handle)
 *
 * @return      MPI_SUCCESS if successful, otherwise state of dims[] and comm_cart
 * is undefined
 *
 * @todo implement fixed dimensions for dims[i] > 0
 */
int MPIX_Cart_ml_create_from_comms(int nlevels, MPI_Comm *level_comms, 
                                   int ndims, double *dim_weights, int *periods,
                                   const MPI_Info info,
                                   int *dims,  MPI_Comm *comm_cart) {
  int total_sizes[nlevels], sizes[nlevels], my_rank[nlevels], my_rank_on_level[nlevels];
  int level;
  int global_equal = 1;
  for(level = 0; level < nlevels; ++level) {
    MPI_Comm_size(level_comms[level], &total_sizes[level]);
    MPI_Comm_rank(level_comms[level], &my_rank[level]);
    int sendbuf[2] = {-total_sizes[level], total_sizes[level]};
    int recvbuf[2] = {-1, -1};
    MPI_Allreduce(sendbuf, recvbuf, 2, MPI_INT, MPI_MAX, level_comms[0]);
    if( (-1*recvbuf[0]) != recvbuf[1] ) global_equal = 0;
  }
  if (global_equal) {
    int d, dd, total_dim, new_rank;
    int dims_ml[ndims*nlevels+nlevels], my_coord_per_level[ndims][nlevels], my_coords[ndims]; 
    MPI_Comm new_comm;
    for (level = 0; level < nlevels-1; ++level)
      sizes[level] = total_sizes[level] / total_sizes[level+1];
    sizes[nlevels-1] = total_sizes[nlevels-1];  
#ifdef MPIX_debug
    if (my_rank[0]==0) {
      printf("MPIX_Cart_ml_create_from_comms: global_equal==1, nlevels=%d, ndims=%d \n",nlevels,ndims);
      printf("MPIX_Cart_ml_create_from_comms: total_sizes=");
         for (level = 0; level < nlevels; ++level) printf(" %d",total_sizes[level]); printf("\n");
      printf("MPIX_Cart_ml_create_from_comms:       sizes=");
         for (level = 0; level < nlevels; ++level) printf(" %d",sizes[level]); printf("\n");
      if (dim_weights == MPIX_WEIGHTS_EQUAL) { printf("MPIX_Cart_ml_create_from_comms: dim_weights == MPIX_WEIGHTS_EQUAL\n");
      } else {
        printf("MPIX_Cart_ml_create_from_comms: dim_weights="); 
         for (d = 0; d<ndims; ++d) printf(" %lf",dim_weights[d]); printf("\n");
      }
    }
#endif
    MPIX_Dims_ml_create(total_sizes[0], ndims, dim_weights, nlevels, sizes, dims_ml);
#ifdef MPIX_debug
    if (my_rank[0]==0) {
      for (d = 0; d<ndims; ++d) {
        printf("MPIX_Cart_ml_create_from_comms: dims_ml[d=%d]=",d);
          for (level = 0; level < nlevels; ++level) printf(" %3d",dims_ml[d*nlevels+level]); printf("\n");
      }
    }    
#endif
    /* Computation of the ranks for each level:
       E.g., with 3 levels, all processes within the same ccNUMA node 
       (i.e. together with MPI_Split_type(MPI_COMM_TYPE_SHARED))
       will have the same my_rank_on_level[0], and
       all processes within same NUMA domain will have same my_rank_on_level[1].
       On most inner level, my_rank_on_level[level] := my_rank[level] */ 
    for(level = 0; level < nlevels-1; ++level) {
      MPI_Comm heads_comm;
      MPI_Comm_split(level_comms[level], (my_rank[level+1]==0 ? 0 : MPI_UNDEFINED), 0, &heads_comm);
      if (heads_comm != MPI_COMM_NULL) { 
        MPI_Comm_rank(heads_comm, &my_rank_on_level[level]);
        MPI_Comm_free(&heads_comm);
      }
      MPI_Bcast(&my_rank_on_level[level], 1, MPI_INT, 0, level_comms[level+1]);
    }
    my_rank_on_level[nlevels-1] = my_rank[nlevels-1];

    /* now, for each level, based on dims_ml[*][level], the coords for each level must be computed: */
    for(level = 0; level < nlevels; ++level) {
      for (d = 0; d<ndims; ++d) my_coord_per_level[d][level] = my_rank_on_level[level];
      for (d = 0; d<ndims-1; ++d) 
        for (dd = d+1; dd<ndims; ++dd) my_coord_per_level[d][level] = my_coord_per_level[d][level] / dims_ml[dd*nlevels+level];
      for (d = 0; d<ndims; ++d)  my_coord_per_level[d][level] =  my_coord_per_level[d][level] % dims_ml[d*nlevels+level]; 
    }

    /* now, the total coordinates are calculated for each dimension: */
    for (d = 0; d<ndims; ++d) { 
      total_dim = dims_ml[d*nlevels+nlevels-1];
      my_coords[d] =  my_coord_per_level[d][nlevels-1]; 
      for(level = nlevels-2; level >= 0; --level) {
        my_coords[d] = my_coords[d] + my_coord_per_level[d][level]*total_dim;
        total_dim = total_dim * dims_ml[d*nlevels+level]; 
      }
      dims[d] = total_dim;
    }

    /* and this is the base for the new rank: */
    new_rank = my_coords[ndims-1];
    total_dim = dims[ndims-1];
    for (d = ndims-2; d >= 0; --d) { 
      new_rank = new_rank + my_coords[d] * total_dim;
      total_dim = total_dim * dims[d];
    }
#ifdef MPIX_debug
    if (total_dim != total_sizes[0]) printf("MPIX_Cart_ml_create_from_comms: total_dim=%d != total_sizes[0]=%d\n",total_dim,total_sizes[0]);
    if (my_rank[0]==0) {
      printf("MPIX_Cart_ml_create_from_comms:        dims=");
         for (d = 0; d<ndims; ++d) printf(" %d",dims[d]); printf("\n");
    }
# if MPIX_debug > 1
    printf("MPIX_Cart_ml_create_from_comms: my_rank=");
        for (level = 0; level < nlevels; ++level) printf(" %3d",my_rank[level]);
    printf("  my_rank_on_level="); 
        for (level = 0; level < nlevels; ++level) {
          printf(" %3d (",my_rank_on_level[level]); 
            for (d = 0; d<ndims; ++d) printf(" %2d",my_coord_per_level[d][level]);
          printf(")");
        } 
    printf("  coords="); 
        for (d = 0; d<ndims; ++d) printf(" %3d",my_coords[d]);
    printf(" new= %3d\n",new_rank);
# endif
#endif
  
    MPI_Comm_split(level_comms[0], /*color*/ 0, new_rank, &new_comm);
    MPI_Cart_create(new_comm, ndims, dims, periods, 0, comm_cart);
    MPI_Comm_free(&new_comm);

  } else {
    /* hard to optimize - future work */
    MPIX_Dims_weighted_create( total_sizes[0], ndims, dim_weights, dims);
    MPI_Cart_create(level_comms[0], ndims, dims, periods, 1, comm_cart);
  }
  return MPI_SUCCESS;
}
