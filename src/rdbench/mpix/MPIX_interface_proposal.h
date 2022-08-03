#ifndef MPIX_INTERFACE_PROPOSAL_H_
#define MPIX_INTERFACE_PROPOSAL_H_

#include <mpi.h>

#if __cplusplus
extern "C" {
#endif

// #define MPIX_debug 1

#define MPIX_WEIGHTS_EQUAL 0
#define MPIX_COMM_TYPE_HALFNODE -1234

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


/** Compute dimensions
 * @param[in]     nnodes      number of nodes in a grid
 * @param[in]     ndims       number of Cartesian dimensions
 * @param[in,out] dims[ndims] nodes in each dimension
 */
/*
int MPIX_Dims_create(int nnodes, int ndims,
                    int *dims);
*/

                    /** @todo add optimization hints for implementors */




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
                             int *dims);

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
                               int *dims_ml);


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
int MPIX_Cart_ml_create(const MPI_Comm comm_old, int ndims,
                               int *periods, int nlevels,
                              //  int (*dims_ml)[], const MPI_Info info,
                               int *dims_ml, const MPI_Info info,
                               int *dims, MPI_Comm *comm_cart);


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
 */
int MPIX_Cart_weighted_create(const MPI_Comm comm_old, int ndims,
                               double *dim_weights, int *periods, const MPI_Info info,
                               int *dims, MPI_Comm *comm_cart);


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
 */
int MPIX_Cart_ml_create_from_comms(int nlevels, MPI_Comm *level_comms,
                                   int ndims, double *dim_weights, int *periods,
                                   const MPI_Info info,
                                   int *dims, MPI_Comm *comm_cart);

/** Splits a given communicator into a set of sub communicators according to a
 * set of given keys.
 * @param[in]   comm_old      old communicator
 * @param[in]   ntype_levels       number of levels to split the communicator into
 * @param[in]   type_levels[ntype_levels] keys as for MPI_Comm_split_type
 * @param[in]   info          additional hints via MPI info
 * @param[out]  new_comms[ntype_levels]   new communicators for level_comms
 * (note that level_comms combines the comm_old and these new_comms in one array) 
 * @return      MPI_SUCCESS if successful, otherwise state of new_comms is undefined
 *
 *
 * @todo Return status in case of internal errors from MPI_* calls.
 */
/*  THIS ROUTINE MAY BE NOT NEEDED
int MPIX_Comm_ml_split_type(const MPI_Comm comm_old, const int ntype_levels,
                                   const int *type_levels,
                                   const MPI_Info info,
                                   MPI_Comm *new_comms);
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
 */
int MPIX_Cart_ml_create_from_types(const MPI_Comm comm_old, const int ntype_levels,
                                  const int *type_levels,
                                  int ndims, double *dim_weights,
                                  int *periods,
                                  const MPI_Info info,
                                  int *dims, MPI_Comm *comm_cart);





/* Internal call sequences:
   ========================

   MPIX_Cart_weighted_create(...)
    -->  chooses available and useful types for splitting, 
           e.g., {MPI_COMM_TYPE_SHARED,OMPI_COMM_TYPE_NUMA}
         MPIX_Cart_ml_create_from_types(...)

   MPIX_Cart_ml_create_from_types(...)
    -->  loop over MPIX_Comm_split_type  [insted of MPIX_Comm_ml_split_type(...)]
         MPIX_Cart_ml_create_from_comms(...)

   MPIX_Cart_ml_create_from_comms(...)
    --> must calculate level_sizes[nlevels] and wether all levels are equally sized
        if (equally-sized) then
         { MPIX_Dims_ml_create(...)
           Appropriate renumbering based on dims_ml and the level_comms
           Calculation of dims[] and creation of comm_cart
         } else
         { algorithm of Thorsten Hoefler }

   MPIX_Cart_ml_create(...)
    --> MPIX_Dims_ml_create(...)
        Appropriate renumbering based on dims_ml
                 and sequenially ranked comm_old
        Calculation of dims[] and creation of comm_cart

   MPIX_Dims_ml_create(...)
    --> MPIX_Dims_weighted_create(...) on each level

   MPIX_Dims_weighted_create(...)
    --> is the base routine 
*/

#if __cplusplus
}
#endif

#endif /* MPIX_INTERFACE_PROPOSAL_H_ */
