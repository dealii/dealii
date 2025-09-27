// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
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

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"

// Solve a nonlinear system with Newton's method in which we report
// recoverable failures when getting too far away from the solution,
// forcing KINSOL to do some backtracking before getting back into the
// region where we are willing to evaluate the residual.
//
// Specifically, we solve the nonlinear problem
//
// F(u) = atan(u)-0.5 = 0
//
// starting at u=10; This is a well-known problematic case because the
// tangent to the curve gets us far far too the left in the first
// iteration, and we use this to test the ability of KINSOL to deal
// with recoverable failures: We pretend that we can't evaluate the
// residual at that point, but the backtracking line search eventually
// gets us back to the region where things work just fine.


int
main()
{
  initlog();

  using VectorType = Vector<double>;

  // Size of the problem
  const unsigned int N = 1;

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/kinsol_newton.prm");
  prm.parse_input(ifile);

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [&](const VectorType &u, VectorType &F) {
    deallog << "Computing residual at " << u[0] << std::endl;

    if ((u[0] < -10) || (u[0] > 20))
      {
        deallog << "Reporting recoverable failure." << std::endl;
        throw RecoverableUserCallbackError();
      }


    F.reinit(u);
    F[0] = std::atan(u[0]) - 0.5;
  };

  double J_inverse;

  kinsol.setup_jacobian = [&J_inverse](const VectorType &u,
                                       const VectorType &F) {
    deallog << "Setting up Jacobian system at u=" << u[0] << std::endl;

    const double J = 1. / (1 + u[0] * u[0]);
    J_inverse      = 1. / J;
  };


  kinsol.solve_with_jacobian = [&](const VectorType &rhs,
                                   VectorType       &dst,
                                   double) { dst[0] = J_inverse * rhs[0]; };

  VectorType v(N);
  v[0] = 10;

  auto niter = kinsol.solve(v);

  deallog << "Found solution at " << std::flush;
  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
