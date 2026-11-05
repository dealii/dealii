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

// Solve min x_0 + 2 x_1, subject to x_0 + x_1 = 3 and x_1 >= 0. The
// minimizer is (3, 0).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>

#include <cmath>
#include <cstring>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType solution_vector(MPI_COMM_SELF, 2, 2);
  VectorType objective(MPI_COMM_SELF, 2, 2);
  VectorType equality_right_hand_side(MPI_COMM_SELF, 1, 1);
  MatrixType equality_matrix;
  Mat        petsc_matrix;
  PetscCallAbort(
    MPI_COMM_SELF,
    MatCreateSeqAIJ(MPI_COMM_SELF, 1, 2, 2, nullptr, &petsc_matrix));
  static_cast<PETScWrappers::MatrixBase &>(equality_matrix)
    .reinit(petsc_matrix);
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_matrix));

  solution_vector = 0.;
  objective(0)    = 1.;
  objective(1)    = 2.;
  objective.compress(VectorOperation::insert);
  equality_matrix.set(0, 0, 1.);
  equality_matrix.set(0, 1, 1.);
  equality_matrix.compress(VectorOperation::insert);
  equality_right_hand_side(0) = 3.;
  equality_right_hand_side.compress(VectorOperation::insert);

  std::optional<Optimizer::Bounds> bounds;
  bounds.emplace();
  bounds->lower.reinit(MPI_COMM_SELF, 2, 2);
  bounds->upper.reinit(MPI_COMM_SELF, 2, 2);
  bounds->lower(0) = PETSC_NINFINITY;
  bounds->lower(1) = 0.;
  bounds->upper(0) = PETSC_INFINITY;
  bounds->upper(1) = PETSC_INFINITY;
  bounds->lower.compress(VectorOperation::insert);
  bounds->upper.compress(VectorOperation::insert);

  PETScWrappers::LinearOptimizer<VectorType, MatrixType> solver;
  Tao                                                    subsolver;
  const char                                            *subsolver_type;
  PetscCallAbort(MPI_COMM_SELF,
                 TaoALMMGetSubsolver(solver.petsc_tao(), &subsolver));
  PetscCallAbort(MPI_COMM_SELF, TaoGetType(subsolver, &subsolver_type));
  AssertThrow(std::strcmp(subsolver_type, TAOBNCG) == 0, ExcInternalError());
  const auto iterations = solver.solve(solution_vector,
                                       objective,
                                       &equality_matrix,
                                       &equality_right_hand_side,
                                       nullptr,
                                       nullptr,
                                       bounds);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  AssertThrow(std::abs(solution_vector(0) - 3.) < 1e-6 &&
                std::abs(solution_vector(1)) < 1e-6,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
