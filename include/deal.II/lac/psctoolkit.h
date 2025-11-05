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
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_PSBLAS

#include "psb_base_cbind.h"
#include "psb_krylov_cbind.h"
#include "psb_prec_cbind.h"
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
         * Distribute local indices and values to the PSBLAS sparse matrix.
         *
         * @param local_dof_indices Vector of local dof indices (in global numbering).
         * @param cell_matrix FullMatrix containing the local matrix values.
         * @param cell_rhs Vector containing the local right-hand side values.
         * @param mh Pointer to the PSBLAS sparse matrix.
         * @param vec Pointer to the PSBLAS vector.
         * @param cdh Pointer to the PSBLAS descriptor.
         * @return 0 on success, non-zero on failure.
         */
        int distribute_local_to_global(const std::vector<types::global_dof_index>& local_dof_indices, const FullMatrix<double>& cell_matrix, const Vector<double>& cell_rhs, psb_c_dspmat *mh, psb_c_dvector *vec, psb_c_descriptor *cdh);
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

    namespace PSBVector {
        /**
         * Create a new PSBLAS vector and initialize it with the given descriptor.
         *
         * @param cd Pointer to the PSBLAS descriptor to initialize the vector.
         * @return Pointer to the created PSBLAS vector.
         */
        psb_c_dvector *CreateVector(psb_c_descriptor *cd);
        /**
        * Distribute local indices and values to the PSBLAS vector.
        *
        * @param local_dof_indices Vector of local dof indices (in global numbering).
        * @param cell_rhs Vector containing the local vector values.
        * @param vector Pointer to the PSBLAS vector.
        * @param cdh Pointer to the PSBLAS descriptor.
        */
        void distribute_local_to_global(const std::vector<types::global_dof_index>& local_dof_indices,  const Vector<double>& cell_rhs, psb_c_dvector *vec, psb_c_descriptor *cdh);
         /**
         * Assemble the PSBLAS vector.
         *
         * @param vector Pointer to the PSBLAS vector to be assembled.
         * @param cd Pointer to the PSBLAS descriptor.
         * @return 0 on success, non-zero on failure.
         */
        int AssembleVector(psb_c_dvector *vector, psb_c_descriptor *cd);
        /**
         * Free the PSBLAS vector.
         *
         * @param vector Pointer to the PSBLAS vector to be freed.
         * @param cd Pointer to the PSBLAS descriptor.
         * @return 0 on success, non-zero on failure.
         */
        int FreeVector(psb_c_dvector *vector, psb_c_descriptor *cd);

    } // namespace Vector   

    namespace Solvers 
    {
        
        /**
         * Create and init a new PSBLAS preconditioner.
         *
         * @param cctxt PSBLAS context.
         * @param ptype Type of the preconditioner
         * @return Pointer to the created PSBLAS preconditioner.
         */
        psb_c_dprec* CreateBasePreconditioner(psb_c_ctxt cctxt, const char *ptype);
        
        /**
        * Build the PSBLAS preconditioner.
        *
        * @param ph Pointer to the PSBLAS preconditioner to be built.
        * @param mh Pointer to the PSBLAS sparse matrix.
        * @param cdh Pointer to the PSBLAS descriptor.
        * @return 0 on success, non-zero on failure.
        */
        int BasePrecBuild(psb_c_dprec *ph, psb_c_dspmat *mh, psb_c_descriptor *cdh);

        /**
         * Free the PSBLAS preconditioner.
         *
         * @param ph Pointer to the PSBLAS preconditioner to be freed.
         * @return 0 on success, non-zero on failure.
         */
        int FreeBasePreconditioner(psb_c_dprec *ph);
  
        /**
         * Create a new PSBLAS solver options structure.
         *
         * @return Pointer to the created PSBLAS solver options structure.
         */
        psb_c_SolverOptions *CreateSolverOptions();

        /**
         * Set solver options with custom values.
         *
         * @param opt Pointer to the PSBLAS solver options structure.
         * @param itmax Maximum number of iterations.
         * @param itrace Print info every itrace iterations (0 for no output).
         * @param irst Restart depth for GMRES or BiCGSTAB(L).
         * @param istop Stopping criterion (1: backward error, 2: relative residual).
         * @param eps Stopping tolerance.
         */
        void SetSolverOptions(psb_c_SolverOptions *opt, int itmax, int itrace, 
                            int irst, int istop, double eps);

        /**
         * Print the current solver options to the log.
         *
         * @param opt Pointer to the PSBLAS solver options structure.
         * @param cctxt Pointer to the PSBLAS context for logging.
         */
        void PrintSolverOptions(const psb_c_SolverOptions *opt, psb_c_ctxt *cctxt);

        /**
         * Print the current solver options to an output file.
         *
         * @param opt Pointer to the PSBLAS solver options structure.
         * @param cctxt Pointer to the PSBLAS context for logging.
         * @param output_file Output file stream to write the options.
         */
        void PrintSolverOptions(const psb_c_SolverOptions *opt, psb_c_ctxt *cctxt, std::ofstream &output_file);

        /**
         * Free the PSBLAS solver options structure.
         *
         * @param opt Pointer to the PSBLAS solver options structure to be freed.
         */
        void FreeSolverOptions(psb_c_SolverOptions *opt);

        /**
         * Solve a linear system using the PSBLAS Krylov solver.
         *
         * @param method The name of the Krylov method to use (e.g., "cg", "gmres").
         * @param ah Pointer to the PSBLAS sparse matrix.
         * @param ph Pointer to the PSBLAS preconditioner.
         * @param bh Pointer to the right-hand side vector.
         * @param xh Pointer to the solution vector.
         * @param cdh Pointer to the PSBLAS descriptor.
         * @param opt Pointer to the PSBLAS solver options.
         * @return 0 on success, non-zero on failure.
         */
        int BaseKrylov(const char *method, psb_c_dspmat *ah, psb_c_dprec *ph, 
            psb_c_dvector *bh, psb_c_dvector *xh,
            psb_c_descriptor *cdh, psb_c_SolverOptions *opt);

    }

} // namespace PSCToolkit


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif // dealii_psctoolkit_h