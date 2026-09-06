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

// Test if Anasazi can correctly compute the eigenvalues of a FiniteDifference
// five-point matrix. Basically the same as slepc/solve_01.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/anasazi_solver.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "../testmatrix.h"
#include "testmatrix.h"



int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(7);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  SolverControl control(5000, 1e-9, false, false);

  const unsigned int size = 46;
  unsigned int       dim  = (size - 1) * (size - 1);

  const unsigned n_eigenvalues = 4;

  deallog << "Size " << size << " Unknowns " << dim << std::endl << std::endl;

  FDMatrix     testproblem(size, size);
  FDDiagMatrix diagonal(size, size);
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> A(
    dim, dim, 5);
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> B(
    dim, dim, 5);

  testproblem.five_point(A);
  A.compress(VectorOperation::insert);
  diagonal.diag(B);
  B.compress(VectorOperation::insert);

  std::vector<
    LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
                      u(n_eigenvalues,
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>(
        complete_index_set(dim), MPI_COMM_WORLD));
  std::vector<double> v(n_eigenvalues);

  LinearAlgebra::TpetraWrappers::AnasaziSolverBlockKrylovSchur
    solver<double, MemorySpace::Default>(control);

  solver.solve(A, B, u, v, v.size());

  for (int i = 0; i < v.size(), i++)
    {
      deallog << "EV" << i + 1 << ": " << v[i] << std::endl;
    }

  deallog << "OK" << std::endl;

  // now we define our
}