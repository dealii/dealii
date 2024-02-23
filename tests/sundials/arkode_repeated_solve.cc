// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
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

// Test repeated calls to solve_ode with harmonic oscillator problem

using VectorType = Vector<double>;

double kappa = 1.0;

std::unique_ptr<SUNDIALS::ARKode<VectorType>>
create_solver()
{
  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/arkode_repeated_solve_in.prm");
  prm.parse_input(ifile);

  auto ode = std::make_unique<SUNDIALS::ARKode<VectorType>>(data);

  // will yield analytic solution y[0] = sin(kappa*t); y[1] = kappa*cos(kappa*t)
  ode->explicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = y[1];
    ydot[1] = -kappa * kappa * y[0];
  };

  ode->output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      deallog << std::setprecision(16) << t << ' ' << sol[0] << ' ' << sol[1]
              << std::endl;
    };
  return ode;
}

int
main()
{
  initlog();

  auto ode = create_solver();
  {
    Vector<double> y(2);
    y[0] = 0;
    y[1] = kappa;

    // solve the ode in one go
    ode->solve_ode(y);

    deallog << "First solve done" << std::endl;

    // solve again -> will call reset automatically
    y[0] = 0;
    y[1] = kappa;
    ode->solve_ode(y);
    deallog << "Second solve done" << std::endl;

    // Solve for a time greater than the originally requested final time
    ode->solve_ode_incrementally(y, 3.0);
    deallog << "Solve further than final time done" << std::endl;
  }

  // solve in two steps with new object
  ode = create_solver();

  {
    Vector<double> y(2);
    y[0] = 0;
    y[1] = kappa;

    ode->solve_ode_incrementally(y, 1.0);
    deallog << "Split solve - part 1 done" << std::endl;
  }

  {
    // N.B. initial conditions are not necessary as the internal state is reused
    Vector<double> y(2);

    ode->solve_ode_incrementally(y, 2.0);
    deallog << "Split solve - part 2 done" << std::endl;
  }

  {
    Vector<double> y(2);
    y[0] = 0;
    y[1] = kappa;
    // solving with single-argument call reruns from start
    ode->solve_ode(y);
    deallog << "Rerun with single argument call done" << std::endl;
  }

  {
    Vector<double> y0(2), y(2);
    y0[0] = 0;
    y0[1] = kappa;
    // manually reset solver
    ode->reset(0.0, 0.01, y0);
    ode->solve_ode_incrementally(y, 2.0);
  }
}
