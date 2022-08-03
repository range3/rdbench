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
 * Copyright 2022, range3 ( https://github.com/range3 )
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
 * levels provided in form of level sizes.
 * @param[in]   comm_old      old communicator. The communicator must be
 * sequentially ordered through the hardware topologies, i.e., rank numbers 0,
 * 1, ..., (level_sizes[nlevels-1]-1) are contiguous in the finest hardware and
 * progress with the next hardware level (nlevels-2), and so on. For
 * level_sizes[level] see definition of dims below.
 * @param[in]   ndims         number of dimensions
 * @param[in]   periods[ndims] logical array of size ndims specifying wether the
 * grid is periodic (true) or not (false) in each dimension
 * @param[in]   nlevels       number of hardware topology levels
 * @param[in]   dims_ml[ndims][nlevels]       values for dimensions in the nlevels hardware
 * topology levels and ndims application dimensions. For each level l, the
 * product of dims[l][*] reflects the size of the hardware topology level, i.e.,
 * level_sizes[l]. For each dimension d the product dims[*][d] == the i-th dimension
 * of the returned Cartesian communicator comm_cart.
 * @param[in]   info          additional hints via MPI info
 * @param[out]  dims[ndims]   integer array of size ndims specifiying the number
 * of nodes in each dimension
 * @param[out]  comm_cart     communicator with new_Cartesian topology (handle)
 *
 * @return      MPI_SUCCESS if successful, otherwise state of dims[] and comm_cart
 * is undefined
 *
 * @todo Return status in case of internal errors from MPI_* calls.
 */
int MPIX_Cart_ml_create(const MPI_Comm comm_old, int ndims, int *periods,
                        // int nlevels, int dims_ml[ndims][nlevels],
                        int nlevels, int *dims_ml,
                        const MPI_Info info,
                        int *dims,  MPI_Comm *comm_cart) {

  /* This routine is only useful if comm_old is sequentially ranked !!! */ 

/* define MPIX__Cart_ml_create_TESTCODE */
#ifdef MPIX__Cart_ml_create_TESTCODE
  {
    /* MPIX__Cart_ml_create() TEST source code, only for nlevels <= 3 and ndims <= 3: */
    int d, level, fixed_dims_ml[3][3];
    int my_old_rank;
    int inner_d0, mid_d0, outer_d0, dim0; 
    int inner_d1, mid_d1, outer_d1, dim1; 
    int inner_d2, mid_d2, outer_d2, dim2; 
    int idim, mdim, odim, whole_size, size, *ranks, old_rank, new_rank;
    int oc0,mc0,ic0, oc1,mc1,ic1, oc2,mc2,ic2, c0,c1,c2;
    MPI_Group world_group, new_group;
    MPI_Comm new_comm;

    MPI_Comm_size(comm_old, &size);
    MPI_Comm_rank(comm_old, &my_old_rank);

    for (d = 0; d< 3   ; ++d) for (level = 0; level <  3     ; ++level) fixed_dims_ml[d][level]=1;
    for (d = 0; d<ndims; ++d) for (level = 0; level < nlevels; ++level) fixed_dims_ml[d][level]=dims_ml[d*nlevels + level];

    outer_d0=fixed_dims_ml[0][0];  mid_d0=fixed_dims_ml[0][1];  inner_d0=fixed_dims_ml[0][2];
    outer_d1=fixed_dims_ml[1][0];  mid_d1=fixed_dims_ml[1][1];  inner_d1=fixed_dims_ml[1][2];
    outer_d2=fixed_dims_ml[2][0];  mid_d2=fixed_dims_ml[2][1];  inner_d2=fixed_dims_ml[2][2];
  
    dim0=inner_d0*mid_d0*outer_d0;
    dim1=inner_d1*mid_d1*outer_d1;
    dim2=inner_d2*mid_d2*outer_d2;
    idim=inner_d0*inner_d1*inner_d2;
    mdim=mid_d0*mid_d1*mid_d2;
    odim=outer_d0*outer_d1*outer_d2;
    whole_size=dim0*dim1*dim2   /* or  =idim*mdim*odim */;  

    dims[0] = dim0; dims[1] = dim1; dims[2] = dim2;

   if(whole_size != size) 
   { if(my_old_rank == 0)  printf("whole_size=%d != size=%d\n",whole_size,size);
     MPI_Abort(comm_old, MPI_ERR_DIMS);
   }

#ifdef MPIX_debug
    if (my_old_rank==0) {
      printf("MPIX_Cart_ml_create TEST: nlevels=%d, ndims=%d \n",nlevels,ndims);
      for (d = 0; d<ndims; ++d) {
        printf("MPIX_Cart_ml_create TEST: dims_ml[d=%d]=",d);
          for (level = 0; level < nlevels; ++level) printf(" %3d",dims_ml[d*nlevels + level]); printf("\n");
      }
      printf("MPIX_Cart_ml_create TEST:        dims=");
         for (d = 0; d<ndims; ++d) printf(" %d",dims[d]); printf("\n");
    }    
#endif

    ranks= malloc(whole_size*sizeof(int));
    for (oc0=0; oc0<outer_d0; oc0++) /*any sequence of the nested loops works*/
     for (mc0=0; mc0<mid_d0;   mc0++)
      for (ic0=0; ic0<inner_d0; ic0++)
       for (oc1=0; oc1<outer_d1; oc1++)
        for (mc1=0; mc1<mid_d1;   mc1++)
         for (ic1=0; ic1<inner_d1; ic1++)
          for (oc2=0; oc2<outer_d2; oc2++)
           for (mc2=0; mc2<mid_d2;   mc2++)
            for (ic2=0; ic2<inner_d2; ic2++)
            { old_rank =   (ic2 + inner_d2*(ic1 + inner_d1*ic0))
                         + (mc2 +   mid_d2*(mc1 +   mid_d1*mc0))*idim
                         + (oc2 + outer_d2*(oc1 + outer_d1*oc0))*idim*mdim;
              c0 = ic0 + inner_d0*mc0 + inner_d0*mid_d0*oc0;
              c1 = ic1 + inner_d1*mc1 + inner_d1*mid_d1*oc1;
              c2 = ic2 + inner_d2*mc2 + inner_d2*mid_d2*oc2;
              new_rank = c2 + dim2*(c1 + dim1*c0);
              ranks[new_rank] = old_rank;
#ifdef MPIX_debug
# if MPIX_debug > 1
              if (my_old_rank==0) {
                printf("MPIX_Cart_ml_create TEST: old_rank= %3d", old_rank);
                printf("  my_rank_on_level= %2d %2d %2d", (oc2 + outer_d2*(oc1 + outer_d1*oc0)),
                   (mc2 +   mid_d2*(mc1 +   mid_d1*mc0)), (ic2 + inner_d2*(ic1 + inner_d1*ic0)) ); 
                printf("  coords= %3d %3d %3d", c0, c1, c2); 
                printf(" new= %3d\n",new_rank);
              }    
# endif
#endif
            }
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, whole_size, ranks, &new_group);  free(ranks); 
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    MPI_Cart_create(new_comm, ndims, dims, periods, 0, comm_cart);
    MPI_Comm_free(&new_comm);
  }
#else
  { /* MPIX__Cart_ml_create() production source code: */
    int d, dd, total_dim, new_rank;
    int level, my_old_rank;
    int my_coord_per_level[ndims][nlevels], my_coords[ndims]; 
    MPI_Comm new_comm;
    int total_sizes[nlevels], sizes[nlevels], my_rank_on_level[nlevels];

    MPI_Comm_rank(comm_old, &my_old_rank);

    for (level = nlevels-1; level >= 0; --level) {
      sizes[level] = 1;
      for (d = 0; d<ndims; ++d) sizes[level] = sizes[level] * dims_ml[d*nlevels+level];
      total_sizes[level] = (level==(nlevels-1) ? sizes[level] : sizes[level] * total_sizes[level+1]);
    }
#ifdef MPIX_debug
    if (my_old_rank==0) {
      printf("MPIX_Cart_ml_create: nlevels=%d, ndims=%d \n",nlevels,ndims);
      for (d = 0; d<ndims; ++d) {
        printf("MPIX_Cart_ml_create: dims_ml[d=%d]=",d);
          for (level = 0; level < nlevels; ++level) printf(" %3d",dims_ml[d*nlevels+level]); printf("\n");
      }
      printf("MPIX_Cart_ml_create:       sizes=");
         for (level = 0; level < nlevels; ++level) printf(" %d",sizes[level]); printf("\n");
      printf("MPIX_Cart_ml_create: total_sizes=");
         for (level = 0; level < nlevels; ++level) printf(" %d",total_sizes[level]); printf("\n");
    }    
#endif
    /* Computation of the ranks for each level:
       E.g., with 3 levels, all processes within the same ccNUMA node 
       (i.e. together with MPI_Split_type(MPI_COMM_TYPE_SHARED))
       will have the same my_rank_on_level[0], and
       all processes within same NUMA domain will have same my_rank_on_level[1].
       On most inner level, my_rank_on_level[level] := my_old_rank % sizes[nlevels-1] */ 
    my_rank_on_level[nlevels-1] = my_old_rank % sizes[nlevels-1];
    for (level = 0; level < nlevels-1; ++level) my_rank_on_level[level] = my_old_rank / total_sizes[level+1] % sizes[level];

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
    if (total_dim != total_sizes[0]) printf("MPIX_Cart_ml_create: total_dim=%d != total_sizes[0]=%d\n",total_dim,total_sizes[0]);
    if (my_old_rank==0) {
      printf("MPIX_Cart_ml_create:        dims=");
         for (d = 0; d<ndims; ++d) printf(" %d",dims[d]); printf("\n");
    }
# if MPIX_debug > 1
    printf("MPIX_Cart_ml_create: my_old_rank= %3d", my_old_rank);
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
  
    MPI_Comm_split(comm_old, /*color*/ 0, new_rank, &new_comm);
    MPI_Cart_create(new_comm, ndims, dims, periods, 0, comm_cart);
    MPI_Comm_free(&new_comm);
  } /* end of production code */
#endif
  return MPI_SUCCESS;
}
