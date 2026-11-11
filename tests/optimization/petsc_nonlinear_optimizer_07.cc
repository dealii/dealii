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

// Verify that an exception thrown by a nonlinear optimizer callback is
// propagated unchanged.

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#include <stdexcept>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType solution(MPI_COMM_SELF, 1, 1);
  MatrixType hessian;
  Mat        petsc_hessian;
  PetscCallAbort(
    MPI_COMM_SELF,
    MatCreateSeqAIJ(MPI_COMM_SELF, 1, 1, 1, nullptr, &petsc_hessian));
  static_cast<PETScWrappers::MatrixBase &>(hessian).reinit(petsc_hessian);
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_hessian));

  Optimizer optimizer;
  optimizer.objective = [](const VectorType &) -> PetscReal {
    throw std::runtime_error("Objective callback failure.");
  };
  optimizer.gradient = [](const VectorType &, VectorType &) {};
  optimizer.hessian  = [](const VectorType &, MatrixType &, MatrixType &) {};

  bool did_throw = false;
  try
    {
      optimizer.solve(solution, hessian);
    }
  catch (const std::runtime_error &exception)
    {
      did_throw =
        std::string(exception.what()) == "Objective callback failure.";
    }

  AssertThrow(did_throw, ExcInternalError());
  deallog << "OK" << std::endl;
}
