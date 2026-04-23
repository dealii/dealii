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
// SUNDIALS::ARKode::implicit_function.
//
// When ARKode is constructed without an explicit ARKodeStepper argument, it
// internally creates an ARKStepper and wraps its callbacks via FunctionProxy.
// Assigning to ode.implicit_function must route to the underlying ARKStepper.
// This test verifies the proxy fires correctly for a pure-implicit ODE solve.

using VectorType = Vector<double>;

int
main()
{
  initlog();
  deallog.precision(3);
  deallog << std::boolalpha;

  // Integrate  y' = -y,  y(0)=1,  implicit DIRK.  Exact: y(t)=exp(-t).
  SUNDIALS::ARKode<VectorType>::AdditionalData data(0.0 /*t0*/,
                                                    1.0 /*tf*/,
                                                    0.01 /*dt*/,
                                                    0.5 /*output_period*/,
                                                    1e-8 /*min dt*/,
                                                    1e-8 /*abstol*/,
                                                    1e-6 /*reltol*/);

  SUNDIALS::ARKode<VectorType> ode(data);

  bool implicit_callback_called = false;

  // Assign via the backwards-compat proxy (implicit-only path).
  ode.implicit_function =
    [&](const double /*t*/, const VectorType &y, VectorType &f) {
      f[0]                     = -y[0];
      implicit_callback_called = true;
    };

  VectorType y(1);
  y[0] = 1.0;
  ode.solve_ode(y);

  deallog << "implicit_callback_called: " << implicit_callback_called
          << std::endl;
  deallog << "Final |y - exp(-1)| < 1e-4: "
          << (std::abs(y[0] - std::exp(-1.0)) < 1e-4) << std::endl;
}
