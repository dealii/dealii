// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Verify that a PETSc error raised by the linear optimizer's internally
// generated objective callback is propagated unchanged.

#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.templates.h>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using Optimizer  = PETScWrappers::LinearOptimizer<VectorType>;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  // We here create vectors x, c of different sizes. This will lead to
  // an error when trying to compute the objective function c^T
  // x. Check below that this expectedly leads to an error:
  VectorType solution(MPI_COMM_SELF, 1, 1);
  VectorType objective(MPI_COMM_SELF, 2, 2);
  Optimizer  optimizer;

  bool did_throw = false;
  try
    {
      optimizer.solve(solution, objective);
    }
  catch (const ExcPETScError &)
    {
      did_throw = true;
    }

  AssertThrow(did_throw, ExcInternalError());
  deallog << "OK" << std::endl;
}
