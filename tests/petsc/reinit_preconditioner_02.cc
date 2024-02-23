// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check re-initializing a preconditioner (serial version)

#include <deal.II/base/index_set.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

template <class PRE>
void
test()
{
  DynamicSparsityPattern csp(5, 5);

  for (unsigned int i = 0; i < 5; ++i)
    csp.add(i, i);

  csp.add(0, 1);
  csp.add(1, 0);

  PETScWrappers::SparseMatrix mat;
  mat.reinit(csp);

  for (unsigned int i = 0; i < 5; ++i)
    mat.set(i, i, 1.0 + i * 2.0);
  mat.set(0, 1, 0.1);
  mat.set(1, 0, 0.1);

  mat.compress(VectorOperation::insert);

  {
    IndexSet indices(5);
    indices.add_range(0, 5);
    PETScWrappers::MPI::Vector src, dst;
    src.reinit(indices, MPI_COMM_WORLD);
    dst.reinit(indices, MPI_COMM_WORLD);
    src(0) = 1.0;
    src(1) = 2.0;
    src.compress(VectorOperation::insert);

    PRE pre;
    pre.initialize(mat);
    pre.vmult(dst, src);
    dst.print(deallog.get_file_stream());

    mat.add(0, 0, 1.0);
    mat.compress(VectorOperation::add);

    pre.initialize(mat);
    pre.vmult(dst, src);
    dst.print(deallog.get_file_stream());
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<PETScWrappers::PreconditionJacobi>();
  test<PETScWrappers::PreconditionBlockJacobi>();
  test<PETScWrappers::PreconditionSOR>();
  test<PETScWrappers::PreconditionSSOR>();
  test<PETScWrappers::PreconditionICC>();
  test<PETScWrappers::PreconditionILU>();
  test<PETScWrappers::PreconditionLU>();
  test<PETScWrappers::PreconditionBoomerAMG>();
  test<PETScWrappers::PreconditionParaSails>();
  test<PETScWrappers::PreconditionNone>();
}
