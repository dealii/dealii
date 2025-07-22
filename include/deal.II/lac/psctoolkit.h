// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_psctoolkit_h
#define dealii_psctoolkit_h

#include <deal.II/base/config.h>
#include <deal.II/lac/full_matrix.h>

#ifdef DEAL_II_WITH_PSBLAS

#include "psb_base_cbind.h"
#include "psb_config.h"

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkit
{
    /**
     * Namespace for PSBLAS communicator functions.
     */
    namespace Communicator
    {
        /**
         * Initialize the PSBLAS communicator and return a pointer to the context.
         * This function creates the context over the MPI_COMM_WORLD communicator.
         *
         * @return Pointer to the initialized PSBLAS context.
         */
        psb_c_ctxt *Init();
        /**
         * Initialize the PSBLAS communicator from an existing MPI communicator.
         *
         * @param comm The MPI communicator to initialize the PSBLAS context from.
         * @return Pointer to the initialized PSBLAS context.
         */
        psb_c_ctxt *InitFromMPI(MPI_Comm mpi_comm);
        /**
         * Exit the PSBLAS communicator and free the context.
         *
         * @param cctxt Pointer to the PSBLAS context to be exited and freed.
         */
        void Exit(psb_c_ctxt *cctxt);
        /**
         * Get information about the current process in the PSBLAS communicator.
         *
         * @param cctxt Pointer to the PSBLAS context.
         * @param iam Pointer to store the rank of the current process.
         * @param nproc Pointer to store the total number of processes.
         */
        void Info(psb_c_ctxt *cctxt, int *iam, int *nproc);
        /**
         * Synchronize all processes in the PSBLAS communicator.
         *
         * @param cctxt Pointer to the PSBLAS context.
         */
        void Barrier(psb_c_ctxt *cctxt);
        /**
         * Create a new PSBLAS descriptor.
         *
         * @param index_set The IndexSet to create the descriptor from.
         * @param cctxt PSBLAS context.
         * @return Pointer to the created PSBLAS descriptor.
         */
        psb_c_descriptor *CreateDescriptor(IndexSet &index_set, psb_c_ctxt cctxt);
        /**
         * Assemble the PSBLAS descriptor.
         *
         * @param cd Pointer to the PSBLAS descriptor to be assembled.
         * @return 0 on success, non-zero on failure.
         */
        int DescriptorAssembly(psb_c_descriptor *cd);
        /**
         * Free the PSBLAS descriptor.
         *
         * @param cd Pointer to the PSBLAS descriptor to be freed.
         * @return 0 on success, non-zero on failure.
         */
        int DescriptorFree(psb_c_descriptor *cd);
    }   // namespace Communicator

    /**
     * Namespace for PSBLAS matrix functions.
     */
    namespace Matrix
    {
        /**
         * Create a new PSBLAS sparse matrix and initialize it with the given descriptor.
         * @param cd Pointer to the PSBLAS descriptor to initialize the matrix.
         * @return Pointer to the created PSBLAS sparse matrix.
         */
        psb_c_dspmat *CreateSparseMatrix(psb_c_descriptor *cd);
        /**
         * Insert a value into the PSBLAS sparse matrix.
         *
         * @param nz Number of non-zero entries to insert.
         * @param irw Row indices of the non-zero entries.
         * @param icl Column indices of the non-zero entries.
         * @param val Values of the non-zero entries.
         * @param mh Pointer to the PSBLAS sparse matrix.
         * @param cdh Pointer to the PSBLAS descriptor.
         * @return 0 on success, non-zero on failure.
         */
        int InsertValue(psb_i_t nz, const psb_l_t *irw, const psb_l_t *icl, const psb_d_t *val, psb_c_dspmat *mh, psb_c_descriptor *cdh);
        /**
        * Distribute local indices and values to the PSBLAS sparse matrix.
        *
        * @param local_dof_indices Vector of local dof indices (in global numbering).
        * @param cell_matrix FullMatrix containing the local matrix values.
        * @param mh Pointer to the PSBLAS sparse matrix.
        * @param cdh Pointer to the PSBLAS descriptor.
        * @return 0 on success, non-zero on failure.
        */
        int distribute_local_to_global(const std::vector<types::global_dof_index>& local_dof_indices, const FullMatrix<double>& cell_matrix,psb_c_dspmat *mh, psb_c_descriptor *cdh);
        /**
         * Assemble the PSBLAS sparse matrix.
         * @param mh Pointer to the PSBLAS sparse matrix to be assembled.
         * @param cdh Pointer to the PSBLAS descriptor.
         * @return 0 on success, non-zero on failure.
         **/
        int AssembleSparseMatrix(psb_c_dspmat *mh, psb_c_descriptor *cdh);
        /**
         * Free the PSBLAS sparse matrix.
         * @param mh Pointer to the PSBLAS sparse matrix to be freed.
         * @param cdh Pointer to the PSBLAS descriptor.
         * @return
         * 0 on success, non-zero on failure.
         */
        int FreeSparseMatrix(psb_c_dspmat *mh, psb_c_descriptor *cdh);
    } // namespace Matrix

    namespace Vector {

    } // namespace Vector   

} // namespace PSCToolkit


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif // dealii_psctoolkit_h