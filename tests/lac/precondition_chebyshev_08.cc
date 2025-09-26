// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PreconditionChebyshev::estimate_eigenvalues() returns 1 as
// as minimum and maximum eigenvalues for fully constrained systems.

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>

#include "../tests.h"

void
test(const internal::EigenvalueAlgorithm &eigenvalue_algorithm)
{
  FullMatrix<double> matrix(2, 2);
  matrix(0, 0) = 1.0;
  matrix(1, 1) = 1.0;

  typename PreconditionChebyshev<FullMatrix<double>,
                                 Vector<double>,
                                 PreconditionIdentity>::AdditionalData ad;

  ad.preconditioner = std::make_shared<PreconditionIdentity>();
  ad.constraints.constrain_dof_to_zero(0);
  ad.constraints.constrain_dof_to_zero(1);

  ad.eigenvalue_algorithm = eigenvalue_algorithm;

  PreconditionChebyshev<FullMatrix<double>,
                        Vector<double>,
                        PreconditionIdentity>
    precon;

  precon.initialize(matrix, ad);

  Vector<double> vec(2);
  const auto     ev = precon.estimate_eigenvalues(vec);

  AssertDimension(ev.min_eigenvalue_estimate, 1.0);
  AssertDimension(ev.max_eigenvalue_estimate, 1.0);
}

int
main()
{
  initlog();

  test(internal::EigenvalueAlgorithm::lanczos);
  test(internal::EigenvalueAlgorithm::power_iteration);

  deallog << "OK!" << std::endl;
}
