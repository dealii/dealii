// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that all the functionality needed for calling
// AffineConstraints<double>::distribute_local_to_global on multiple
// combinations of PETScWrappers objects is correctly instantiated and works.

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

#include <numeric>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  const unsigned int n = 4;

  PETScWrappers::SparseMatrix serial_matrix(n, n, n);

  // PETSc requires us to set up the sparsity pattern ahead of construction
  DynamicSparsityPattern dsp(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      dsp.add(i, j);
  IndexSet all_dofs(n);
  all_dofs.add_range(0, n);

  PETScWrappers::MPI::SparseMatrix mpi_matrix;
  mpi_matrix.reinit(all_dofs, all_dofs, dsp, MPI_COMM_WORLD);

  PETScWrappers::MPI::Vector mpi_vector(MPI_COMM_WORLD, n, n);

  FullMatrix<PetscScalar> cell_matrix(n, n);
  Vector<PetscScalar>     cell_rhs(n);

  AffineConstraints<PetscScalar> constraints;
  constraints.close();

  std::vector<types::global_dof_index> gdi(n);
  std::iota(gdi.begin(), gdi.end(), 0);
  std::iota(cell_rhs.begin(), cell_rhs.end(), 0);

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
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
