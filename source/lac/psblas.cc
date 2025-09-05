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
        psb_c_set_index_base(0); // Set index base to 0
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
        psb_c_set_index_base(0); // Set index base to 0
        
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
      int distribute_local_to_global(const std::vector<types::global_dof_index>& local_dof_indices, const FullMatrix<double>& cell_matrix, const Vector<double>& cell_rhs, psb_c_dspmat *mh, psb_c_dvector *vec, psb_c_descriptor *cdh)
      {
        // Get the number of local dofs
        const unsigned int dofs_per_cell = local_dof_indices.size();
        psb_i_t nz = dofs_per_cell * dofs_per_cell; // Number of non-zero entries

        // Allocate memory for row and column indices and values
        psb_l_t *irw_vec = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
        psb_l_t *irw = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
        psb_l_t *icl = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
        psb_d_t *val = (psb_d_t *)malloc(nz * sizeof(psb_d_t));
        psb_d_t *vec_val = (psb_d_t *)malloc(dofs_per_cell * sizeof(psb_d_t));

        // Fill the arrays with local indices and values
        for (unsigned int i = 0; i < dofs_per_cell; ++i)  
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            irw[i * dofs_per_cell + j] = local_dof_indices[i];
            icl[i * dofs_per_cell + j] = local_dof_indices[j];
            val[i * dofs_per_cell + j] = cell_matrix(i, j);
          }
          // Fill the vector values
          irw_vec[i] = local_dof_indices[i];
          vec_val[i] = cell_rhs(i);
        }

        // Insert the values into the sparse matrix
        int info = psb_c_dspins(nz,irw,icl,val,mh,cdh);
        if (info != 0)
        {
          deallog << "Error inserting values into PSBLAS sparse matrix: " << info << std::endl;
        }

        // Insert the values into the vector
        info = psb_c_dgeins(dofs_per_cell, irw_vec, vec_val, vec, cdh);
        if (info != 0)
        {
          deallog << "Error inserting values into PSBLAS vector: " << info << std::endl;
        }

        // Free allocated memory
        free(irw);
        free(icl);
        free(val);
        free(irw_vec);
        free(vec_val);

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

  namespace PSBVector {

      /**
       * Create a new PSBLAS vector and initialize it with the given descriptor.
       *
       * @param cd Pointer to the PSBLAS descriptor to initialize the vector.
       * @return Pointer to the created PSBLAS vector.
       */
      psb_c_dvector* CreateVector(psb_c_descriptor *cd)
      {
        // Create a new PSBLAS vector
        psb_c_dvector *vector = psb_c_new_dvector();
        if (vector == nullptr)
        {
          deallog << "Error creating PSBLAS vector." << std::endl;
        }
        // Initialize the vector with the descriptor
        int info = psb_c_dgeall_remote(vector, cd);
        if (info != 0)
        {
          deallog << "Error initializing PSBLAS vector: " << info << std::endl;
        }
        return vector;
      }

      /** 
       * Assemble the PSBLAS vector.
       * @param vector Pointer to the PSBLAS vector to be assembled.
       * @param cd Pointer to the PSBLAS descriptor.
       * @return 0 on success, non-zero on failure.
       */
      int AssembleVector(psb_c_dvector *vector, psb_c_descriptor *cd)
      {
        if (vector == nullptr || cd == nullptr)
        {
          deallog << "Error assembling PSBLAS vector: null pointer." << std::endl;
          return -1;
        }
        int info = psb_c_dgeasb(vector, cd);
        if (info != 0)
        {
          deallog << "Error assembling PSBLAS vector: " << info << std::endl;
        }
        return info;
      }

      /**
       * Distribute local indices and values to the PSBLAS vector.
       *
       * @param local_dof_indices Vector of local dof indices (in global numbering).
       * @param cell_rhs Vector containing the local vector values.
       * @param vec Pointer to the PSBLAS vector.
       * @param cdh Pointer to the PSBLAS descriptor.
       */
      void distribute_local_to_global(const std::vector<types::global_dof_index>& local_dof_indices,  const Vector<double> &cell_rhs, psb_c_dvector *vec, psb_c_descriptor *cdh)
      {
        // Get the number of local dofs
        const unsigned int dofs_per_cell = local_dof_indices.size();
        psb_i_t nz = dofs_per_cell; // Number of non-zero entries

        // Allocate memory for row indices and values
        psb_l_t *irw = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
        psb_d_t *val = (psb_d_t *)malloc(nz * sizeof(psb_d_t));

        // Fill the arrays with local indices and values
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          irw[i] = local_dof_indices[i];
          val[i] = cell_rhs(i);
        }

        // Insert the values into the vector
        int info = psb_c_dgeins(nz,irw,val,vec,cdh);

        // Free allocated memory
        free(irw);
        free(val);

        if (info != 0)
        {
          deallog << "Error inserting values into PSBLAS vector: " << info << std::endl;
        }
                                   
      }

      /**
       * Free the PSBLAS vector.
       *
       * @param vector Pointer to the PSBLAS vector to be freed.
       * @param cd Pointer to the PSBLAS descriptor.
       * @return 0 on success, non-zero on failure.
       */
      int FreeVector(psb_c_dvector *vector, psb_c_descriptor *cd)
      {
        if (vector != nullptr && cd != nullptr)
        {
          int result = psb_c_dgefree(vector, cd);
          if (result != 0)
          {
            deallog << "Error freeing PSBLAS vector: " << result << std::endl;
          }
          return result; // Success
        }
        else
        {
          deallog << "Warning: Attempted to free a null PSBLAS vector." << std::endl;
          return -1; // Failure
        }
      }

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
    psb_c_dprec* CreateBasePreconditioner(psb_c_ctxt cctxt, const char *ptype)
    {
      psb_c_dprec *ph = psb_c_new_dprec();
      if (ph == nullptr)
      {
        deallog << "Error allocating memory for PSBLAS preconditioner." << std::endl;
        return nullptr;
      }
      // Initialize the preconditioner with default values
      int info = psb_c_dprecinit(cctxt,ph,ptype);
      if (info != 0)
      {
        deallog << "Error initializing PSBLAS preconditioner: " << info << std::endl;
        free(ph);
        return nullptr;
      }
      return ph;
    }

    /**
     * Build the PSBLAS preconditioner.
     *
     * @param ph Pointer to the PSBLAS preconditioner to be built.
     * @param mh Pointer to the PSBLAS sparse matrix.
     * @param cdh Pointer to the PSBLAS descriptor.
     * @return 0 on success, non-zero on failure.
     */
    int BasePrecBuild(psb_c_dprec *ph, psb_c_dspmat *mh, psb_c_descriptor *cdh)
    {
      if (ph == nullptr || mh == nullptr || cdh == nullptr)
      {
        deallog << "Error building PSBLAS preconditioner: null pointer." << std::endl;
        return -1; // Failure
      }
      int info = psb_c_dprecbld(mh, cdh, ph);
      if (info != 0)
      {
        deallog << "Error building PSBLAS preconditioner: " << info << std::endl;
      }
      return info; // Success or failure
    }

    /**
     * Free the PSBLAS preconditioner.
     *
     * @param ph Pointer to the PSBLAS preconditioner to be freed.
     * @return 0 on success, non-zero on failure.
     */
    int FreeBasePreconditioner(psb_c_dprec *ph)
    {
      if (ph != nullptr)
      {
        int result = psb_c_dprecfree(ph);
        if (result != 0)
        {
          deallog << "Error freeing PSBLAS preconditioner: " << result << std::endl;
        }
        return result; // Success
      }
      else
      {
        deallog << "Warning: Attempted to free a null PSBLAS preconditioner." << std::endl;
        return -1; // Failure
      }
    }

    /**
     * Create and initialize a PSBLAS solver options structure with default values.
     *
     * @return Pointer to the initialized PSBLAS solver options structure.
     */
    psb_c_SolverOptions* CreateSolverOptions()
    {
      psb_c_SolverOptions *opt = (psb_c_SolverOptions*)malloc(sizeof(psb_c_SolverOptions));
      if (opt == nullptr)
      {
        deallog << "Error allocating memory for PSBLAS solver options." << std::endl;
        return nullptr;
      }
      
      // Initialize with default values
      opt->iter = 0;          // Initial iteration count
      opt->itmax = 1000;      // Maximum iterations
      opt->itrace = 1;        // Trace output enabled by default
      opt->irst = 30;         // Restart depth for GMRES
      opt->istop = 2;         // Use relative residual norm ||r||_2/||b||_2
      opt->eps = 1.0e-6;      // Convergence tolerance
      opt->err = 0.0;         // Initial error
      
      return opt;
    }

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
                          int irst, int istop, double eps)
    {
      if (opt == nullptr)
      {
        deallog << "Error: null pointer passed to SetSolverOptions." << std::endl;
        return;
      }
      
      opt->itmax = itmax;
      opt->itrace = itrace;
      opt->irst = irst;
      opt->istop = istop;
      opt->eps = eps;
    }

    /**
     * Print the current solver options to the log.
     *
     * @param opt Pointer to the PSBLAS solver options structure.
     * @param cctxt Pointer to the PSBLAS context for logging.
     */
    void PrintSolverOptions(const psb_c_SolverOptions *opt, psb_c_ctxt *cctxt)
    {
      if (opt == nullptr)
      {
        deallog << "Error: null pointer passed to PrintSolverOptions." << std::endl;
        return;
      }
      
      int iam, nproc;
      psb_c_info(*cctxt, &iam, &nproc);
      if (iam == 0) // Only print from rank 0
      {
        deallog << "PSBLAS Solver Options:" << std::endl;
        deallog << "  Maximum iterations (itmax): " << opt->itmax << std::endl;
        deallog << "  Trace frequency (itrace): " << opt->itrace << std::endl;
        deallog << "  Restart depth (irst): " << opt->irst << std::endl;
        deallog << "  Stopping criterion (istop): " << opt->istop << std::endl;
        deallog << "  Tolerance (eps): " << opt->eps << std::endl;
        deallog << "  Iterations performed (iter): " << opt->iter << std::endl;
        deallog << "  Final error (err): " << opt->err << std::endl;
      }
    }

    /**
     * Print the current solver options to an output file.
     *
     * @param opt Pointer to the PSBLAS solver options structure.
     * @param cctxt Pointer to the PSBLAS context for logging.
     * @param output_file Output file stream to write the options.
     */
    void PrintSolverOptions(const psb_c_SolverOptions *opt, psb_c_ctxt *cctxt, std::ofstream &output_file)
    {
      if (opt == nullptr)
      {
        deallog << "Error: null pointer passed to PrintSolverOptions." << std::endl;
        return;
      }
      
      int iam, nproc;
      psb_c_info(*cctxt, &iam, &nproc);
      if (iam == 0) // Only print from rank 0
      {
        output_file << "PSBLAS Solver Options:" << std::endl;
        output_file << "  Maximum iterations (itmax): " << opt->itmax << std::endl;
        output_file << "  Trace frequency (itrace): " << opt->itrace << std::endl;
        output_file << "  Restart depth (irst): " << opt->irst << std::endl;
        output_file << "  Stopping criterion (istop): " << opt->istop << std::endl;
        output_file << "  Tolerance (eps): " << opt->eps << std::endl;
        output_file << "  Iterations performed (iter): " << opt->iter << std::endl;
        output_file << "  Final error (err): " << opt->err << std::endl;
      }
    }

    /**
     * Free the PSBLAS solver options structure.
     *
     * @param opt Pointer to the PSBLAS solver options structure to be freed.
     */
    void FreeSolverOptions(psb_c_SolverOptions *opt)
    {
      if (opt != nullptr)
      {
        free(opt);
      }
      else
      {
        deallog << "Warning: Attempted to free a null PSBLAS solver options structure." << std::endl;
      }
    }

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
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt)
    {
    if (method == nullptr || ah == nullptr || bh == nullptr || xh == nullptr || cdh == nullptr || opt == nullptr)
    {
      deallog << "Error: null pointer passed to Krylov solver." << std::endl;
      return -1; // Failure
    }
    // Call the PSBLAS Krylov solver
    int info = psb_c_dkrylov(method, ah, ph, bh, xh, cdh, opt);
    if (info != 0)
    {
      deallog << "Error solving linear system with PSBLAS Krylov solver: " << info << std::endl;
    }
    return info; // Success or failure
    }

  } // namespace Solvers

} // namespace PSCToolkit

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
