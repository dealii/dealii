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


/**
 * Solve an ODE problem of exponential growth, written as a DAE in
 * which one could easily eliminate one variable, using a direct
 * solver for the jacobian system.
 *
 * The equation we want to solve here is
 *   x' = a y
 *   0  = x-y
 * with initial conditions
 *   x(0) = 1
 *   y(0) = 1
 *
 * That is, with Y=[x, y]:
 *   F(Y', Y, t) = [1 0 ; 0 0] Y' + [0 -a ; -1 1] Y
 *   Y(0)        = [1 1]
 *
 * The exact solution is
 *
 * x(t) = y(t) = exp(a t)
 *
 * The Jacobian to assemble is the following:
 *
 * J = dF/dY + alpha dF/dY' = [0 -a ; -1 1] + alpha [1 0 ; 0 0]
 *   = [alpha -a ; -1 1]
 */

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  SUNDIALS::IDA<Vector<double>>::AdditionalData data;
  ParameterHandler                              prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/ida_06_in.prm");
  prm.parse_input(ifile);

  const double a = 1.0;
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
    //  F(Y', Y, t) = [1 0 ; 0 0] Y' + [0 -a ; -1 1] Y
    res    = 0;
    res[0] = y_dot[0] - a * y[0];
    res[1] = -y[0] + y[1];
  };

  time_stepper.setup_jacobian = [&](const double,
                                    const VectorType &,
                                    const VectorType &,
                                    const double alpha) {
    // J = [alpha -a ; -1 1]
    J(0, 0) = alpha;
    J(0, 1) = -a;
    J(1, 0) = -1;
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


  y[0] = y[1] = 1;
  y_dot[0] = y_dot[1] = a;
  time_stepper.solve_dae(y, y_dot);
}
