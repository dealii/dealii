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


// Test implicit-explicit time stepper. Only jacobian_times_vector.
// Brusselator benchmark

/**
 * This test problem is called "brusselator", and is a typical benchmark for
 * ODE solvers. This problem has 3 dependent variables u, v and w, that depend
 * on the independent variable t via the IVP system
 *
 * du/dt = a - (w + 1)u + v u^2
 * dv/dt = w u - v u^2
 * dw/dt = (b - w)/eps - w u
 *
 * We integrate over the interval 0 <= t <= 10, with the initial conditions
 *
 * u(0) = 3.9, v(0) = 1.1, w(0) = 2.8,
 *
 * and parameters
 *
 * a = 1.2, b = 2.5, and eps = 10^-5
 *
 * The implicit part only contains the stiff part of the problem (the part with
 * eps in right hand side of the third equation).
 */
int
main()
{
  initlog();

  using VectorType = Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_04_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_04_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  // Parameters
  double u0 = 3.9, v0 = 1.1, w0 = 2.8, a = 1.2, b = 2.5, eps = 1e-5;
  // Explicit jacobian.
  FullMatrix<double> J(3, 3);
  J(2, 2) = -1.0 / eps;

  ode.implicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = 0;
    ydot[1] = 0;
    ydot[2] = -y[2] / eps;
  };


  ode.explicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = a - (y[2] + 1) * y[0] + y[1] * y[0] * y[0];
    ydot[1] = y[2] * y[0] - y[1] * y[0] * y[0];
    ydot[2] = b / eps - y[2] * y[0];
  };

  ode.jacobian_times_vector = [&](const VectorType &src,
                                  VectorType       &dst,
                                  double            t,
                                  const VectorType & /*y*/,
                                  const VectorType & /*fy*/) {
    J.vmult(dst, src);
  };

  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      deallog << std::setprecision(16) << t << ' ' << sol[0] << ' ' << sol[1]
              << ' ' << sol[2] << std::endl;
    };

  Vector<double> y(3);
  y[0] = u0;
  y[1] = v0;
  y[2] = w0;
  ode.solve_ode(y);
}
