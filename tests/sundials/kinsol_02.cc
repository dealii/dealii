// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"

// Solve a nonlinear system in the form accepted by Picard iteration.
//
// We solve the nonlinear problem
//
// F(u) = L u + N(u) = 0
//
// where L is a constant matrix, and N(u) is non linear.
//
// We set L = id and
//
// N_i(u) = .1*u_i^2 - i - 1
//
int
main()
{
  initlog();

  using VectorType = Vector<double>;

  // Size of the problem
  const unsigned int N = 2;

  FullMatrix<double> L(N, N);
  L(0, 0) = 1;
  L(1, 1) = 1;
  L(0, 1) = 1;

  FullMatrix<double> Linv(N, N);
  Linv.invert(L);

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/kinsol_picard.prm");
  prm.parse_input(ifile);

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [&](const VectorType &u, VectorType &F) {
    F = u;

    F[0] += .1 * u[0] * u[0] - 1;
    F[1] += .1 * u[1] * u[1] - 2;
  };

  kinsol.solve_with_jacobian =
    [&](const VectorType &rhs, VectorType &dst, double) { dst = rhs; };

  VectorType v(N);

  auto niter = kinsol.solve(v);

  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
