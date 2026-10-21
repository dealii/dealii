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

#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

// A diagonal generalized eigenproblem has known eigenvalues, independently
// of the solver interface and the shift-and-invert implementation.
struct DiagonalOperator
{
  std::vector<double> diagonal;

  void
  vmult(Vector<double> &dst, const Vector<double> &src) const
  {
    for (unsigned int i = 0; i < diagonal.size(); ++i)
      dst[i] = diagonal[i] * src[i];
  }
};

void
check(const bool legacy)
{
  const unsigned int size  = 10;
  const double       shift = 0.25;
  DiagonalOperator   system, mass, inverse;
  for (unsigned int i = 0; i < size; ++i)
    {
      mass.diagonal.push_back(i + 2.);
      system.diagonal.push_back((i + 1.) * (i + 2.));
      inverse.diagonal.push_back(1. / ((i + 1. - shift) * (i + 2.)));
    }
  SolverControl control(1000, 1e-12);
  ArpackSolver  solver(
    control,
    ArpackSolver::AdditionalData(6, ArpackSolver::largest_magnitude, true));
  solver.set_shift(shift);
  Vector<double> initial(size);
  initial = 1.;
  solver.set_initial_vector(initial);
  std::vector<std::complex<double>> eigenvalues(2);
  std::vector<Vector<double>>       eigenvectors(2, Vector<double>(size));
  if (legacy)
    solver.solve(system, mass, inverse, eigenvalues, eigenvectors);
  else
    solver.solve(mass, inverse, eigenvalues, eigenvectors);

  std::vector<double> sorted;
  for (unsigned int i = 0; i < eigenvalues.size(); ++i)
    {
      AssertThrow(std::abs(eigenvalues[i].imag()) < 1e-10, ExcInternalError());
      sorted.push_back(eigenvalues[i].real());
      Vector<double> residual(size), weighted(size);
      system.vmult(residual, eigenvectors[i]);
      mass.vmult(weighted, eigenvectors[i]);
      AssertThrow(std::abs(eigenvectors[i] * weighted - 1.) < 1e-10,
                  ExcInternalError());
      residual.add(-eigenvalues[i].real(), weighted);
      AssertThrow(residual.l2_norm() < 1e-9, ExcInternalError());
    }
  std::sort(sorted.begin(), sorted.end());
  AssertThrow(std::abs(sorted[0] - 1.) < 1e-10, ExcInternalError());
  AssertThrow(std::abs(sorted[1] - 2.) < 1e-10, ExcInternalError());
  deallog << (legacy ? "Compatibility" : "Preferred") << " interface OK"
          << std::endl;
}

int
main()
{
  initlog();
  check(false);
  check(true);
}
