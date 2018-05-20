// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Check that all the functionality needed for calling
// ConstraintMatrix::distribute_local_to_global on mutiple combinations of
// PETScWrappers objects is correctly instantiated and works.

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main(int argc, char** argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  const unsigned int n = 4;

  PETScWrappers::SparseMatrix      serial_matrix(n, n, n);
  PETScWrappers::MPI::SparseMatrix mpi_matrix(MPI_COMM_WORLD, n, n, n, n, n);
  PETScWrappers::MPI::Vector       mpi_vector(MPI_COMM_WORLD, n, n);

  FullMatrix<PetscScalar> cell_matrix(n, n);
  Vector<PetscScalar>     cell_rhs(n);

  ConstraintMatrix constraints;
  constraints.close();

  std::vector<types::global_dof_index> gdi(n);
  std::iota(gdi.begin(), gdi.end(), 0);
  std::iota(cell_rhs.begin(), cell_rhs.end(), 0);

  for(unsigned int i = 0; i < n; ++i)
    for(unsigned int j = 0; j < n; ++j)
      cell_matrix(i, j) = i * n + j;

  constraints.distribute_local_to_global(
    cell_matrix, cell_rhs, gdi, serial_matrix, mpi_vector);
  mpi_vector.compress(VectorOperation::add);
  serial_matrix.compress(VectorOperation::add);
  mpi_vector.print(deallog.get_file_stream());
  serial_matrix.print(deallog.get_file_stream());

  constraints.distribute_local_to_global(
    cell_matrix, cell_rhs, gdi, mpi_matrix, mpi_vector);
  mpi_vector.compress(VectorOperation::add);
  mpi_matrix.compress(VectorOperation::add);
  mpi_vector.print(deallog.get_file_stream());
  mpi_matrix.print(deallog.get_file_stream());
}
