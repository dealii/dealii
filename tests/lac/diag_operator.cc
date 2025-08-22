// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check diagonal_operator

#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"


void
check_vmult()
{
  const unsigned int size = 4;
  IndexSet           local_indices(size);
  local_indices.add_range(0, size);
  local_indices.compress();
  TrilinosWrappers::MPI::Vector diagonal(local_indices, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector u(local_indices, MPI_COMM_WORLD);

  // Create identity matrix as the exemplar.
  TrilinosWrappers::SparsityPattern sparsity(local_indices,
                                             local_indices,
                                             MPI_COMM_WORLD);
  for (unsigned int i = 0; i < size; ++i)
    sparsity.add(i, i);
  sparsity.compress();

  TrilinosWrappers::SparseMatrix A;
  A.reinit(sparsity);
  for (unsigned int i = 0; i < size; ++i)
    {
      A.set(i, i, 1.0);
    }
  A.compress(VectorOperation::insert);

  // diagonal=[1, 2, 3, 4]
  for (unsigned int i = 0; i < size; ++i)
    {
      diagonal[i] = i + 1;
    }
  diagonal.compress(VectorOperation::insert);

  // u=[2, 2, 2, 2]
  for (unsigned int i = 0; i < size; ++i)
    {
      u[i] = 2;
    }
  u.compress(VectorOperation::insert);

  TrilinosWrappers::MPI::Vector v(local_indices, MPI_COMM_WORLD);

  auto exemplar = linear_operator<
    TrilinosWrappers::MPI::Vector,
    TrilinosWrappers::MPI::Vector,
    TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>(
    A);
  auto diag_operator = diagonal_operator(exemplar, diagonal);

  diag_operator.vmult(v, u);

  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << std::endl;

  diag_operator.vmult_add(v, u);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << std::endl;

  diag_operator.Tvmult(v, u);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << std::endl;

  diag_operator.Tvmult_add(v, u);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  initlog();
  deallog << std::setprecision(0) << std::fixed;

  check_vmult();
  deallog << "OK" << std::endl;
}
