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

#include <deal.II/sundials/arkode.h>

#include "../tests.h"

// Test implicit time stepper, no jacobian. Only implements implicit_function.

/**
 * Solve the Harmonic oscillator problem.
 *
 * u'' = -k^2 u
 * u (0) = 0
 * u'(0) = k
 *
 * write in terms of a first order ode:
 *
 * y[0]' =       y[1]
 * y[1]' = - k^2 y[0]
 *
 * That is
 *
 * y' = A y
 *
 * A = [ 0 , 1; -k^2, 0 ]
 *
 * y_0  = 0, k
 *
 * The exact solution is
 *
 * y[0](t) = sin(k t)
 * y[1](t) = k cos(k t)
 */
int
main()
{
  initlog();

  using VectorType = Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  // Use the same parameters of test 1.
  std::ifstream ifile(SOURCE_DIR "/arkode_01_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  double kappa = 1.0;

  ode.implicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = y[1];
    ydot[1] = -kappa * kappa * y[0];
  };

  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      // limit the output to every 10th step and increase the precision to make
      // the test more robust
      if (step_number % 10 == 0)
        deallog << t << ' ' << std::setprecision(7) << sol[0] << ' ' << sol[1]
                << std::endl;
    };

  Vector<double> y(2);
  y[0] = 0;
  y[1] = kappa;
  ode.solve_ode(y);
}
