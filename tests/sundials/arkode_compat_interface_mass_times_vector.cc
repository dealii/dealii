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
// SUNDIALS::ARKode::mass_times_vector.
//
// Trivial ODE:  M y' = -y,  y(0)=1,  M=I  (identity mass matrix).
// The mass matrix is the identity so SUNDIALS' built-in unpreconditioned
// SPGMR solver handles M x = b without needing solve_mass or any mass
// preconditioner callbacks.  This isolates mass_times_vector alone.

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

  bool mtv_called = false;

  ode.implicit_function = [&](const double /*t*/,
                              const VectorType &y,
                              VectorType       &f) { f[0] = -y[0]; };

  // Assign via the backwards-compat proxy.
  // Identity mass: Mv = v.
  ode.mass_times_vector =
    [&](const double /*t*/, const VectorType &v, VectorType &Mv) {
      Mv         = v;
      mtv_called = true;
    };

  VectorType y(1);
  y[0] = 1.0;
  ode.solve_ode(y);

  deallog << "mass_times_vector_called: " << mtv_called << std::endl;
  deallog << "Final |y - exp(-1)| < 1e-4: "
          << (std::abs(y[0] - std::exp(-1.0)) < 1e-4) << std::endl;
}
