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


#include <deal.II/base/logstream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/psctoolkit.h>


#ifdef DEAL_II_WITH_PSBLAS

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkit
{

    namespace Communicator
    {
    
      /**
       * Initialize the PSBLAS communicator and return a pointer to the context.
       * This function creates the context over the MPI_COMM_WORLD communicator
       *
       * @return Pointer to the initialized PSBLAS context.
       */
      psb_c_ctxt *Init()
      {
        psb_c_ctxt *cctxt = psb_c_new_ctxt();
        psb_c_init(cctxt);
        return cctxt;
      }

      /**
       * Initialize the PSBLAS communicator from an existing MPI communicator.
       *
       * @param comm The MPI communicator to initialize the PSBLAS context from.
       * @return Pointer to the initialized PSBLAS context.
       * 
       * Note: This routine has to be called *if and only if* the
       * communicator has to inherit the MPI communicator initialized by Deal.II
       * Utilities::MPI::MPI_InitFinalize class. Observe that if the communicator was initialized
       * with Init(), then the responsability of exiting the MPI communicator lies
       * on the user and on a call on PSCToolkit::Communicator::Exit() routine.
       */
      psb_c_ctxt *InitFromMPI(MPI_Comm comm)
      {
        // Conver the MPI communicator to a Fortran-style communicator
        // and initialize the PSBLAS context from it.
        if (comm == MPI_COMM_NULL)
        {
          deallog << "Error: MPI_COMM_NULL passed to InitFromMPI." << std::endl;
          return nullptr; // Handle the case where the communicator is null
        }
        // Convert MPI_Comm to Fortran-style communicator
        MPI_Fint f_comm = MPI_Comm_c2f(comm);
        psb_c_ctxt *cctxt = psb_c_new_ctxt();
        psb_c_init_from_fint(cctxt, f_comm);
        
        return cctxt;
      }

      /**
       * Exit the PSBLAS communicator and free the context.
       *
       * @param cctxt Pointer to the PSBLAS context to be exited and freed.
       * 
       * Note: This routine has to be called *if and only if* the
       * communicator was initialized with Init(), if the communicator was
       * initialized with InitFromMPI(), then the responsability of exiting
       * the MPI communicator lies on the Deal.II Utilities::MPI::MPI_InitFinalize
       * class. 
       */
      void Exit(psb_c_ctxt *cctxt)
      {
        psb_c_exit(*cctxt);
        free(cctxt);
      }
      /**
       * Get information about the current process in the PSBLAS communicator.
       *
       * @param cctxt Pointer to the PSBLAS context.
       * @param iam Pointer to store the rank of the current process.
       * @param nproc Pointer to store the total number of processes.
       */
      void Info(psb_c_ctxt *cctxt, int *iam, int* nproc)
      {
        psb_c_info(*cctxt, iam, nproc);
      }

      /**
       * Synchronize all processes in the PSBLAS communicator.
       *
       * @param cctxt Pointer to the PSBLAS context.
       */
      void Barrier(psb_c_ctxt *cctxt)
      {
        psb_c_barrier(*cctxt);
      }
      
      /**
       * Create a PSBLAS descriptor from an IndexSet.
       *
       * @param index_set The IndexSet containing the indices to be included in the descriptor.
       * @param cctxt PSBLAS context.
       * @return Pointer to the created PSBLAS descriptor.
       */
      psb_c_descriptor *CreateDescriptor(IndexSet &index_set, psb_c_ctxt cctxt)
      {
        // Create a new PSBLAS descriptor
        psb_c_descriptor *cd = psb_c_new_descriptor();

        // Use get_index_vector() from IndexSet to get the indexes
        const std::vector<types::global_dof_index> &indexes = index_set.get_index_vector();

        psb_i_t number_of_local_indexes = indexes.size(); // Number of local indexes
        // Copy the indexes into a psb_l_t array called vl
        psb_l_t *vl = (psb_l_t *)malloc(number_of_local_indexes * sizeof(psb_l_t));
        for (psb_i_t i = 0; i < number_of_local_indexes; ++i)
        {
          vl[i] = static_cast<psb_l_t>(indexes[i]);
        }

        // Insert the indexes into the descriptor
        psb_c_cdall_vl(number_of_local_indexes,vl,cctxt,cd);

        // Free the vl array
        free(vl);
        return cd;
      }

      /**
      * Assemble the PSBLAS descriptor.
      *
      * @param cd Pointer to the PSBLAS descriptor to be assembled.
      * @return 0 on success, non-zero on failure.
      */
      int DescriptorAssembly(psb_c_descriptor *cd)
      {
        // Assemble the descriptor
        int result = psb_c_cdasb(cd);
        if (result != 0)
        {
          deallog << "Error assembling PSBLAS descriptor: " << result << std::endl;
        }
        return result;
      }

      /**
       * Free the PSBLAS descriptor.
       *
       * @param cd Pointer to the PSBLAS descriptor to be freed.
       */
      int DescriptorFree(psb_c_descriptor *cd)
      {
        if (cd != nullptr)
        {
          int result = psb_c_cdfree(cd);
          if (result != 0)
          {
            deallog << "Error freeing PSBLAS descriptor: " << result << std::endl;
          }
          return result; // Success
        }
        else
        {
          deallog << "Warning: Attempted to free a null PSBLAS descriptor." << std::endl;
          return -1; // Failure
        }
      }



    }

    namespace Matrix
    {
      /**
       * Create a new PSBLAS sparse matrix and initialize it with the given descriptor.
       *
       * @param cd Pointer to the PSBLAS descriptor to initialize the matrix.
       * @return Pointer to the created PSBLAS sparse matrix.
       */
      psb_c_dspmat* CreateSparseMatrix(psb_c_descriptor *cd)
      {
        // Create a new PSBLAS sparse matrix
        psb_c_dspmat *matrix = psb_c_new_dspmat();
        if (matrix == nullptr)
        {
          deallog << "Error creating PSBLAS sparse matrix." << std::endl;
        }
        // Initialize the sparse matrix with the descriptor
        int info = psb_c_dspall_remote(matrix, cd);
        if (info != 0)
        {
          deallog << "Error initializing PSBLAS sparse matrix: " << info << std::endl;
        }
        return matrix;
      }

      /**
       * Insert a value into the PSBLAS sparse matrix.
       *
       * @param nz Number of non-zero entries to insert.
       * @param irw Array of row indices for the non-zero entries.
       * @param icl Array of column indices for the non-zero entries.
       * @param val Array of values for the non-zero entries.
       * @param mh Pointer to the PSBLAS sparse matrix.
       * @param cdh Pointer to the PSBLAS descriptor.
       * @return 0 on success, non-zero on failure.
       */
      int InsertValue(psb_i_t nz, const psb_l_t *irw, const psb_l_t *icl, const psb_d_t *val, psb_c_dspmat *mh, psb_c_descriptor *cdh)
      {
        // Insert a value into the sparse matrix
        int info = psb_c_dspins(nz,irw,icl,val,mh,cdh);
        if (info != 0)
        {
          deallog << "Error inserting value into PSBLAS sparse matrix: " << info << std::endl;
        }
        return info;
      }

      /**
       * Distribute local indices and values to the PSBLAS sparse matrix.
       *
       * @param local_dof_indices Vector of local dof indices (in global numbering).
       * @param cell_matrix FullMatrix containing the local matrix values.
       * @param mh Pointer to the PSBLAS sparse matrix.
       * @param cdh Pointer to the PSBLAS descriptor.
       * @return 0 on success, non-zero on failure.
       */
      int distribute_local_to_global(const std::vector<types::global_dof_index>& local_dof_indices, const FullMatrix<double>& cell_matrix,psb_c_dspmat *mh, psb_c_descriptor *cdh)
      {
        // Get the number of local dofs
        const unsigned int dofs_per_cell = local_dof_indices.size();
        psb_i_t nz = dofs_per_cell * dofs_per_cell; // Number of non-zero entries

        // Allocate memory for row and column indices and values
        psb_l_t *irw = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
        psb_l_t *icl = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
        psb_d_t *val = (psb_d_t *)malloc(nz * sizeof(psb_d_t));

        // Fill the arrays with local indices and values
        for (unsigned int i = 0; i < dofs_per_cell; ++i)  
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            irw[i * dofs_per_cell + j] = local_dof_indices[i];
            icl[i * dofs_per_cell + j] = local_dof_indices[j];
            val[i * dofs_per_cell + j] = cell_matrix(i, j);
          }
        }

        // Insert the values into the sparse matrix
        int info = psb_c_dspins(nz,irw,icl,val,mh,cdh);

        // Free allocated memory
        free(irw);
        free(icl);
        free(val);

        return info;
      }

      /**
       * Assemble the PSBLAS sparse matrix.
       *
       * @param mh Pointer to the PSBLAS sparse matrix to be assembled.
       * @param cdh Pointer to the PSBLAS descriptor.
       * @return 0 on success, non-zero on failure.
       */
      int AssembleSparseMatrix(psb_c_dspmat *mh, psb_c_descriptor *cdh)
      {
        int info = -1;
        // Check if the sparse matrix is not already assembled
        if (! psb_c_dis_matasb(mh, cdh))
        {
        // Assemble the sparse matrix since it is not already assembled
        info = psb_c_dspasb(mh, cdh);
          if (info != 0)
          {
            deallog << "Error assembling PSBLAS sparse matrix: " << info << std::endl;
          }
        }
        return info;
      }

      /**
       * Free the PSBLAS sparse matrix.
       *
       * @param mh Pointer to the PSBLAS sparse matrix to be freed.
       * @param cdh Pointer to the PSBLAS descriptor.
       * @return 0 on success, non-zero on failure.
       */
      int FreeSparseMatrix(psb_c_dspmat *mh, psb_c_descriptor *cdh)
      {
        if (mh != nullptr && cdh != nullptr)
        {
          int result = psb_c_dspfree(mh,cdh);
          return result; // Success
        }
        else
        {
          deallog << "Warning: Attempted to free a null PSBLAS sparse matrix." << std::endl;
          return -1; // Failure
        }
      }
    
    } // namespace Matrix

  namespace Vector {

  } // namespace Vector

} // namespace PSCToolkit

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
