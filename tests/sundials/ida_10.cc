// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/ida.h>

#include "../tests.h"


// In discretized models, one often has nodes on the boundary that
// have to satisfy boundary values. These are *algebraic* variables
// for which one would ordinarily apply boundary values that set the
// residual vector to either zero (if one solves the implicit DAE via
// a Newton method for an update vector) or to the correct boundary
// value (if one solves for a solution vector). Which one is not
// clear. Test this for a simple case where the boundary values are
// time dependent.
//
// The test case we consider looks like this:
//
// x'(t) = y(t)
// 0     = y(t) - (1+t)
//
// with initial conditions x(0)=0. (We don't provide initial
// conditions for y(0).)
//
// The exact solution is obtained by plugging y(t) = 1+t into the
// first equation and integrating, to obtain
//   x(t) = int_0^t (1+tau) dtau = t+t^2/2.

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
    //  F(Y', Y, t) = [x' -y ; y-(1+t)]
    res    = 0;
    res[0] = y_dot[0] - y[1];
    res[1] = y[1] - (1 + t);
  };

  time_stepper.setup_jacobian = [&](const double,
                                    const VectorType &y,
                                    const VectorType &,
                                    const double alpha) {
    // J = [alpha 0 ; 0 1]
    J(0, 0) = alpha;
    J(0, 1) = 0;
    J(1, 0) = 0;
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

    // Now also check that the computed solution is about correct:
    VectorType error = sol;
    error(0) -= t + t * t / 2;
    error(1) -= 1 + t;
    Assert(error.l2_norm() < 1e-9, ExcInternalError());

    // Same for Y'. Because we provide inconsistent initial conditions
    // for y(0), the computed solution is only correct for t>0. Fix
    // this by manually setting the error to zero for y at t=0.
    VectorType error_dot = sol_dot;
    error_dot(0) -= 1 + t;
    error_dot(1) -= 1;
    if (t == 0)
      error_dot(1) = 0;
    Assert(error_dot.l2_norm() < 1e-9, ExcInternalError());
  };

  time_stepper.differential_components = []() {
    IndexSet x(2);
    x.add_index(0);
    return x;
  };

  // Provide correct initial conditions x(0), but incorrect initial
  // conditions for y(0) and derivatives x'(0), y'(0):
  y[0]     = 0;  // correct
  y[1]     = 42; // wrong
  y_dot[0] = 43; // wrong
  y_dot[1] = 44; // wrong
  time_stepper.solve_dae(y, y_dot);
}
