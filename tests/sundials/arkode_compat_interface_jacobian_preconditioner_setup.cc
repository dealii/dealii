// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include "../tests.h"

// Test the backwards-compatibility proxy for
// SUNDIALS::ARKode::jacobian_preconditioner_setup.
//
// The setup callback is called by ARKStep before jacobian_preconditioner_solve
// whenever it deems the Jacobian data stale.  This test verifies the setup
// proxy routes correctly.
//
// Problem: y' = -y (implicit), y(0)=1.  J = -I.

using VectorType = Vector<double>;

int
main()
{
  initlog();
  deallog.precision(3);
  deallog << std::boolalpha;

  SUNDIALS::ARKode<VectorType>::AdditionalData data(0.0 /*t0*/,
                                                    1.0 /*tf*/,
                                                    0.01 /*dt*/);

  SUNDIALS::ARKode<VectorType> ode(data);

  bool jpu_called = false;

  ode.implicit_function = [&](const double /*t*/,
                              const VectorType &y,
                              VectorType       &f) { f[0] = -y[0]; };

  ode.jacobian_times_vector = [&](const VectorType &v,
                                  VectorType       &Jv,
                                  const double /*t*/,
                                  const VectorType & /*y*/,
                                  const VectorType & /*fy*/) { Jv[0] = -v[0]; };

  ode.solve_linearized_system =
    [&](SUNDIALS::SundialsOperator<VectorType>       &op,
        SUNDIALS::SundialsPreconditioner<VectorType> &prec,
        VectorType                                   &x,
        const VectorType                             &b,
        double                                        tol) {
      SolverControl        control(100, tol);
      SolverCG<VectorType> solver_cg(control);
      solver_cg.solve(op, x, b, prec);
    };

  // Assign via the backwards-compat proxy: store gamma for use in solve.
  ode.jacobian_preconditioner_setup = [&](const double /*t*/,
                                          const VectorType & /*y*/,
                                          const VectorType & /*fy*/,
                                          const int /*jok*/,
                                          int         &jcur,
                                          const double gamma) {
    jcur       = 1;
    jpu_called = true;
  };

  ode.jacobian_preconditioner_solve = [&](const double /*t*/,
                                          const VectorType & /*y*/,
                                          const VectorType & /*fy*/,
                                          const VectorType &r,
                                          VectorType       &z,
                                          const double      gamma,
                                          const double /*tol*/,
                                          const int /*lr*/) {
    z[0] = r[0] / (1.0 + gamma);
  };

  VectorType y(1);
  y[0] = 1.0;
  ode.solve_ode(y);

  deallog << "jacobian_preconditioner_setup_called: " << jpu_called
          << std::endl;
  deallog << "Final |y - exp(-1)| < 1e-4: "
          << (std::abs(y[0] - std::exp(-1.0)) < 1e-4) << std::endl;
}
