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

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include "../tests.h"

// Test the backwards-compatibility proxy for
// SUNDIALS::ARKode::solve_mass.
//
// Trivial ODE:  M y' = -y,  y(0)=1,  M=I  (identity mass matrix).
// solve_mass requires mass_times_vector to be set (the entire mass solver
// setup is gated by it).  For M=I the mass solve is just x = b.

using VectorType = Vector<double>;

int
main()
{
  initlog();
  deallog.precision(3);
  deallog << std::boolalpha;

  SUNDIALS::ARKode<VectorType>::AdditionalData data(0.0 /*t0*/,
                                                    1.0 /*tf*/,
                                                    0.01 /*dt*/,
                                                    0.5 /*output_period*/,
                                                    1e-8 /*min dt*/,
                                                    1e-8 /*abstol*/,
                                                    1e-6 /*reltol*/);

  SUNDIALS::ARKode<VectorType> ode(data);

  bool callback_called = false;

  ode.implicit_function = [&](const double /*t*/,
                              const VectorType &y,
                              VectorType       &f) { f[0] = -y[0]; };

  // Identity mass: Mv = v.
  ode.mass_times_vector =
    [&](const double /*t*/, const VectorType &v, VectorType &Mv) { Mv = v; };

  // Assign via the backwards-compat proxy.
  // For M=I, solving Mx=b is simply x=b.
  ode.solve_mass =
    [&](SUNDIALS::SundialsOperator<VectorType> &,
        SUNDIALS::SundialsPreconditioner<VectorType> &,
        VectorType                                   &x,
        const VectorType                             &b,
        double) {
      x               = b;
      callback_called = true;
    };

  VectorType y(1);
  y[0] = 1.0;
  ode.solve_ode(y);

  deallog << "solve_mass_called: " << callback_called << std::endl;
  deallog << "Final |y - exp(-1)| < 1e-4: "
          << (std::abs(y[0] - std::exp(-1.0)) < 1e-4) << std::endl;
}
