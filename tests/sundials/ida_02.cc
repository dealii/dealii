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
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/ida.h>

#include "../tests.h"


/**
 * Solve the Harmonic oscillator problem, using an iterative solver for the
 * jacobian system (with the interface that is available only in version >
 * 4.0.0)
 *
 * u'' = -k^2 u
 * u (0) = 0
 * u'(0) = k
 *
 * written in terms of a first order ode:
 *
 * y[0]' -     y[1]  = 0
 * y[1]' + k^2 y[0]  = 0
 *
 * That is
 *
 * F(y', y, t) = y' + A y = 0
 *
 * A = [ 0 , -1; k^2, 0 ]
 *
 * y_0  = 0, k
 * y_0' = k, 0
 *
 * The exact solution is
 *
 * y[0](t) = sin(k t)
 * y[1](t) = k cos(k t)
 *
 * The Jacobian to assemble is the following:
 *
 * J = alpha I + A
 */
class HarmonicOscillator
{
public:
  HarmonicOscillator(
    double                                                        _kappa,
    const typename SUNDIALS::IDA<Vector<double>>::AdditionalData &data)
    : time_stepper(data)
    , y(2)
    , y_dot(2)
    , J(2, 2)
    , A(2, 2)
    , Jinv(2, 2)
    , kappa(_kappa)
  {
    using VectorType = Vector<double>;

    time_stepper.reinit_vector = [&](VectorType &v) { v.reinit(2); };


    time_stepper.residual = [&](const double      t,
                                const VectorType &y,
                                const VectorType &y_dot,
                                VectorType       &res) {
      res = y_dot;
      A.vmult_add(res, y);
    };

    time_stepper.setup_jacobian = [&](const double,
                                      const VectorType &,
                                      const VectorType &,
                                      const double alpha) {
      A(0, 1) = -1.0;
      A(1, 0) = kappa * kappa;

      J = A;

      J(0, 0) += alpha;
      J(1, 1) += alpha;
    };

    time_stepper.solve_with_jacobian =
      [&](const VectorType &src, VectorType &dst, const double tolerance) {
        SolverControl solver_control(1000, tolerance, false, false);
        SolverGMRES<Vector<double>> solver(solver_control);
        solver.solve(J, dst, src, PreconditionIdentity());
      };

    time_stepper.output_step = [&](const double       t,
                                   const VectorType  &sol,
                                   const VectorType  &sol_dot,
                                   const unsigned int step_number) {
      deallog << t << ' ' << sol[0] << ' ' << sol[1] << ' ' << sol_dot[0] << ' '
              << sol_dot[1] << std::endl;
    };
  }

  void
  run()
  {
    y[1]     = kappa;
    y_dot[0] = kappa;
    time_stepper.solve_dae(y, y_dot);
  }
  SUNDIALS::IDA<Vector<double>> time_stepper;

private:
  Vector<double>     y;
  Vector<double>     y_dot;
  FullMatrix<double> J;
  FullMatrix<double> A;
  FullMatrix<double> Jinv;
  double             kappa;
};


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  SUNDIALS::IDA<Vector<double>>::AdditionalData data;
  ParameterHandler                              prm;
  data.add_parameters(prm);

  // std::ofstream ofile(SOURCE_DIR "/ida_01.prm");
  // prm.print_parameters(ofile, ParameterHandler::ShortPRM);
  // ofile.close();

  std::ifstream ifile(SOURCE_DIR "/ida_01_in.prm");
  prm.parse_input(ifile);


  HarmonicOscillator ode(1.0, data);
  ode.run();
}
