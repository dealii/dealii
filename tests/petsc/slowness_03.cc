// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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



// this is part of a whole suite of tests that checks the relative speed of
// using PETSc for sparse matrices as compared to the speed of our own
// library. the tests therefore may not all actually use PETSc, but they are
// meant to compare it
//
// the tests build the 5-point stencil matrix for a uniform grid of size N*N

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  const unsigned int N      = 200;
  const unsigned int n_dofs = N * N;

  // build the sparse matrix
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int my_id = Utilities::MPI::this_mpi_process(mpi_communicator);
  Assert(n_dofs % n_mpi_processes == 0, ExcInternalError());
  const unsigned int n_local_dofs = n_dofs / n_mpi_processes;
  IndexSet           locally_owned_dofs(n_dofs);
  locally_owned_dofs.add_range(my_id * n_dofs / n_mpi_processes,
                               (my_id + 1) * n_dofs / n_mpi_processes);
  IndexSet               locally_relevant_dofs = locally_owned_dofs;
  DynamicSparsityPattern dsp(n_dofs);
  for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
      {
        const unsigned int global = i * N + j;
        dsp.add(global, global);
        if (j > 0)
          {
            dsp.add(global - 1, global);
            dsp.add(global, global - 1);
          }
        if (j < N - 1)
          {
            dsp.add(global + 1, global);
            dsp.add(global, global + 1);
          }
        if (i > 0)
          {
            dsp.add(global - N, global);
            dsp.add(global, global - N);
          }
        if (i < N - 1)
          {
            dsp.add(global + N, global);
            dsp.add(global, global + N);
          }
      }

  PETScWrappers::MPI::SparseMatrix matrix;
  matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
      {
        const unsigned int global = i * N + j;
        matrix.add(global, global, 4);
        if (j > 0)
          {
            matrix.add(global - 1, global, -1);
            matrix.add(global, global - 1, -1);
          }
        if (j < N - 1)
          {
            matrix.add(global + 1, global, -1);
            matrix.add(global, global + 1, -1);
          }
        if (i > 0)
          {
            matrix.add(global - N, global, -1);
            matrix.add(global, global - N, -1);
          }
        if (i < N - 1)
          {
            matrix.add(global + N, global, -1);
            matrix.add(global, global + N, -1);
          }
      }
  matrix.compress(VectorOperation::add);

  // then do a single matrix-vector
  // multiplication with subsequent formation
  // of the matrix norm
  PETScWrappers::MPI::Vector v1(PETSC_COMM_WORLD, n_dofs, n_dofs);
  PETScWrappers::MPI::Vector v2(PETSC_COMM_WORLD, n_dofs, n_dofs);
  for (unsigned int i = 0; i < n_dofs; ++i)
    v1(i) = i;
  matrix.vmult(v2, v1);

  deallog << v1 * v2 << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test();
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
