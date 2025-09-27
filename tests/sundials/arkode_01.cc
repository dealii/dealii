// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
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

// Test explicit time stepper. Only implements explicit_function.

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

  // Set to true to reset input file.
  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_01_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_01_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  double kappa = 1.0;

  unsigned int n_rhs_evaluations = 0;
  ode.explicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = y[1];
    ydot[1] = -kappa * kappa * y[0];

    ++n_rhs_evaluations;
  };

  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      deallog << t << ' ' << sol[0] << ' ' << sol[1] << std::endl;
    };

  Vector<double> y(2);
  y[0] = 0;
  y[1] = kappa;

  const unsigned int n_timesteps = ode.solve_ode(y);

  deallog << "n_rhs_evaluations=" << n_rhs_evaluations << std::endl;
  deallog << "n_timesteps=" << n_timesteps << std::endl;
}
