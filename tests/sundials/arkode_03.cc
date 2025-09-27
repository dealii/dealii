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
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/sundials/arkode.h>

#include <arkode/arkode_arkstep.h>

#include "../tests.h"


// Test implicit-explicit time stepper, no jacobian. Use L:d:V (in serial)
// Brusselator benchmark

/**
 * This test problem is called "brusselator", and is a typical benchmark for
 * ODE solvers. This problem has 3 dependent variables u, v and w, that depend
 * on the independent variable t via the IVP system
 *
 * du/dt = a − (w + 1)u + v u^2
 * dv/dt = w u − v u^2
 * dw/dt = (b − w)/eps -w u
 *
 * We integrate over the interval 0 ≤ t ≤ 10, with the initial conditions
 *
 * u(0) = 3.9, v(0) = 1.1, w(0) = 2.8,
 *
 * and parameters
 *
 * a = 1.2, b = 2.5, and eps = 10−5
 *
 * The implicit part only contains the stiff part of the problem (the part with
 * eps in right hand side of the third equation).
 */
int
main()
{
  initlog();

  using VectorType = LinearAlgebra::distributed::Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_03_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_03_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  // Parameters
  double u0 = 3.9, v0 = 1.1, w0 = 2.8, a = 1.2, b = 2.5, eps = 1e-5;

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

  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      // limit the output to every 10th step and increase the precision to make
      // the test more robust
      if (step_number % 10 == 0)
        deallog << t << ' ' << std::setprecision(10) << sol[0] << ' ' << sol[1]
                << ' ' << sol[2] << std::endl;
    };

  // This test, for reasons I don't fully understand, generates some output
  // which varies between environments much more than the other ARKODE
  // tests. Work around it by setting a fairly stringent maximum time step.
  ode.custom_setup = [&](void *arkode_mem) {
    int ierr = ARKStepSetMinStep(arkode_mem, 1e-8);
    AssertThrow(ierr == 0, ExcInternalError());
    ierr = ARKStepSetMaxStep(arkode_mem, 1e-4);
    AssertThrow(ierr == 0, ExcInternalError());
    ierr = ARKStepSetMaxNumSteps(arkode_mem, 5000);
    AssertThrow(ierr == 0, ExcInternalError());
  };

  VectorType y(3);
  y[0] = u0;
  y[1] = v0;
  y[2] = w0;
  ode.solve_ode(y);
}
