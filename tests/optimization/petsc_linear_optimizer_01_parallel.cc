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

// On four MPI processes, solve min sum_{r=0}^3 (x_{2r}+2x_{2r+1}), subject to
// x_{2r}+x_{2r+1}=3 and x_{2r+1}>=0 for r=0,...,3. This is four independent
// linear programs, each with minimizer (3,0).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const MPI_Comm comm  = MPI_COMM_WORLD;
  const PetscInt rank  = Utilities::MPI::this_mpi_process(comm);
  const PetscInt first = 2 * rank;
  VectorType     solution_vector(comm, 8, 2);
  VectorType     objective(comm, 8, 2);
  VectorType     equality_right_hand_side(comm, 4, 1);
  MatrixType     equality_matrix;
  Mat            petsc_matrix;
  PetscCallAbort(
    comm,
    MatCreateAIJ(comm, 1, 2, 4, 8, 2, nullptr, 0, nullptr, &petsc_matrix));
  static_cast<PETScWrappers::MatrixBase &>(equality_matrix)
    .reinit(petsc_matrix);
  PetscCallAbort(comm, MatDestroy(&petsc_matrix));

  solution_vector      = 0.;
  objective(first)     = 1.;
  objective(first + 1) = 2.;
  objective.compress(VectorOperation::insert);
  equality_matrix.set(rank, first, 1.);
  equality_matrix.set(rank, first + 1, 1.);
  equality_matrix.compress(VectorOperation::insert);
  equality_right_hand_side(rank) = 3.;
  equality_right_hand_side.compress(VectorOperation::insert);

  std::optional<Optimizer::Bounds> bounds;
  bounds.emplace();
  bounds->lower.reinit(comm, 8, 2);
  bounds->upper.reinit(comm, 8, 2);
  bounds->lower(first)     = PETSC_NINFINITY;
  bounds->lower(first + 1) = 0.;
  bounds->upper(first)     = PETSC_INFINITY;
  bounds->upper(first + 1) = PETSC_INFINITY;
  bounds->lower.compress(VectorOperation::insert);
  bounds->upper.compress(VectorOperation::insert);

  PETScWrappers::LinearOptimizer<VectorType, MatrixType> solver;
  const auto iterations = solver.solve(solution_vector,
                                       objective,
                                       &equality_matrix,
                                       &equality_right_hand_side,
                                       nullptr,
                                       nullptr,
                                       bounds);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;
  AssertThrow(std::abs(solution_vector(first) - 3.) < 1e-6 &&
                std::abs(solution_vector(first + 1)) < 1e-6,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
