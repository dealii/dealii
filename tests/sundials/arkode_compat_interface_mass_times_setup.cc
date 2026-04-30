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
// SUNDIALS::ARKode::mass_times_setup.
//
// Trivial ODE:  M y' = -y,  y(0)=1,  M=I  (identity mass matrix).
// mass_times_setup cannot be used without mass_times_vector because the
// entire mass solver setup is gated by `if (mass_times_vector)`.  We set
// mass_times_vector to the identity and flag mass_times_setup.

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

  bool mts_called = false;

  ode.implicit_function = [&](const double /*t*/,
                              const VectorType &y,
                              VectorType       &f) { f[0] = -y[0]; };

  // Identity mass: Mv = v.
  ode.mass_times_vector =
    [&](const double /*t*/, const VectorType &v, VectorType &Mv) { Mv = v; };

  // Assign via the backwards-compat proxy.
  ode.mass_times_setup = [&](const double /*t*/) { mts_called = true; };

  VectorType y(1);
  y[0] = 1.0;
  ode.solve_ode(y);

  deallog << "mass_times_setup_called: " << mts_called << std::endl;
  deallog << "Final |y - exp(-1)| < 1e-4: "
          << (std::abs(y[0] - std::exp(-1.0)) < 1e-4) << std::endl;
}
