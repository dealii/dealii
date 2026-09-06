// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test if Anasazi can obtain the correct eigenvalues of a diagonal matrix
// for different number types.

#include <deal.II/lac/anasazi_solver.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <typename Number>
void
run_test()
{
  unsigned int                                        matrows = 10;
  LinearAlgebra::TpetraWrappers::SparseMatrix<Number> K(matrows, matrows, 1);

  for (int i = 0; i < matrows; i++)
    {
      K_double.set(i, i, i + 1);
    }

  K.compress(VectorOperation::add);

  SolverControl solver_control(1000, 1e-9);

  LinearAlgebra::TpetraWrappers::AnasaziSolverBlockKrylovSchur<Number>
    anasazi_eigensolver(solver_control);

  // define the remaining inputs
  std::vector<Number>                                        eigenvalues;
  std::vector<LinearAlgebra::TpetraWrappers::Vector<Number>> eigenvectors;

  anasazi_eigensolver.solve(K, eigenvalues, eigenvectors, matrows);

  for (int j = 0; j < matrows; j++)
    {
      deallog << "eigenvalue" << j + 1 << ": " << eigenvalues[j] << std::endl;
    }

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  // run with both double and float
  deallog << "Run with number-type double" << std::endl;
  run_test<double>();
  deallog << "Run with number-type float" << std::endl;
  run_test<float>();
}