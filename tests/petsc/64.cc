// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// This test should be run on multiple processors. note that this test also
// started to fail with the upgrade to petsc 2.2.1 which required a fix in
// PETScWrappers::MatrixBase::operator=


#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

#include "../testmatrix.h"


template <typename MatrixType>
void
test(MatrixType &m)
{
  m.add(0, 0, 1);
  m.compress(VectorOperation::add);
  m = 0;

  Assert(fabs(m.frobenius_norm()) < 1e-15, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        const unsigned int N      = 20;
        const unsigned int n_dofs = N * N;
        // check
        // PETScWrappers::SparseMatrix
        PETScWrappers::SparseMatrix v1(n_dofs, n_dofs, 5);
        test(v1);

        // check
        // PETScWrappers::MPI::SparseMatrix
        MPI_Comm mpi_communicator(MPI_COMM_WORLD);
        int      n_jobs = 1;
        MPI_Comm_size(mpi_communicator, &n_jobs);
        const unsigned int n_mpi_processes = static_cast<unsigned int>(n_jobs);
        Assert(n_dofs % n_mpi_processes == 0, ExcInternalError());
        PETScWrappers::MPI::SparseMatrix v2;
        {
          FDMatrix               fd_matrix(N, N);
          DynamicSparsityPattern dsp(n_dofs, n_dofs);
          fd_matrix.five_point_structure(dsp);
          dsp.add(0, 0); // be sure that we have this one
          SparsityPattern sparsity_pattern;
          sparsity_pattern.copy_from(dsp);
          IndexSet all_dofs(n_dofs);
          all_dofs.add_range(0, n_dofs);
          v2.reinit(all_dofs, all_dofs, sparsity_pattern, PETSC_COMM_WORLD);
        }
        test(v2);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
