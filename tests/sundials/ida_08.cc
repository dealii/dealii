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

#include <deal.II/sundials/ida.h>

#include "../tests.h"


// Like the _07 test, but we only provide correct initial conditions
// for Y[0]', not for Y[1]' -- this is not uncommon in cases where we
// have a first order DAE where at best we know initial conditions for
// all variables (or could compute them), but definitely not initial
// velocities Y'.

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  SUNDIALS::IDA<Vector<double>>::AdditionalData data;
  ParameterHandler                              prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/ida_08_in.prm");
  prm.parse_input(ifile);

  const double a = 1.0;
  const double p = 1.5;
  deallog << "Exponential growth factor = " << a << std::endl;

  using VectorType = Vector<double>;

  VectorType         y(2);
  VectorType         y_dot(2);
  FullMatrix<double> J(2, 2);
  FullMatrix<double> A(2, 2);
  FullMatrix<double> Jinv(2, 2);

  SUNDIALS::IDA<Vector<double>> time_stepper(data);

  time_stepper.reinit_vector = [&](VectorType &v) { v.reinit(2); };


  time_stepper.residual = [&](const double      t,
                              const VectorType &y,
                              const VectorType &y_dot,
                              VectorType       &res) {
    //  F(Y', Y, t) = [x' -a y^{1/p} ; -x^p + y]
    res    = 0;
    res[0] = y_dot[0] - a * std::pow(y[1], 1. / p);
    res[1] = -std::pow(y[0], p) + y[1];
  };

  time_stepper.setup_jacobian = [&](const double,
                                    const VectorType &y,
                                    const VectorType &,
                                    const double alpha) {
    // J = [alpha -ay^{1/p-1}/p ; -px^{p-1} 1]
    J(0, 0) = alpha;
    J(0, 1) = -a * std::pow(y[1], 1. / p - 1) / p;
    J(1, 0) = -p * std::pow(y[0], p - 1);
    J(1, 1) = 1;

    Jinv.invert(J);
  };

  time_stepper.solve_with_jacobian =
    [&](const VectorType &src, VectorType &dst, const double) {
      Jinv.vmult(dst, src);
    };

  time_stepper.output_step = [&](const double       t,
                                 const VectorType  &sol,
                                 const VectorType  &sol_dot,
                                 const unsigned int step_number) {
    deallog << t << ' ' << sol[0] << ' ' << sol[1] << ' ' << sol_dot[0] << ' '
            << sol_dot[1] << std::endl;
  };


  // Provide correct initial conditions Y(0), but incorrect initial
  // derivatives Y'(0):
  y[0] = y[1] = 1; // correct
  y_dot[0]    = 0; // wrong
  y_dot[1]    = 0; // wrong
  time_stepper.solve_dae(y, y_dot);
}
