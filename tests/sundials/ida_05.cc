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
 * Like the _03 test, but throw a recoverable exception if the time
 * step is longer than 0.1. We will hopefully see that the solver
 * reduces the time step and tries again.
 *
 * Solve the exponential decay problem:
 *
 * y' = -k y
 * y (0) = 1
 *
 * with k=2. That is
 *
 * F(y', y, t) = y' + k y = 0
 * y_0  = 1
 *
 * The exact solution is
 *
 * y(t) = exp(-k*t)
 *
 * The Jacobian we then need to assemble is the following:
 *
 * J = dF/dy + alpha dF/dy'
 *   = k + alpha 1
 *   = (k+alpha)
 */
class HarmonicOscillator
{
public:
  HarmonicOscillator(
    double                                                        _kappa,
    const typename SUNDIALS::IDA<Vector<double>>::AdditionalData &data)
    : time_stepper(data)
    , y(1)
    , y_dot(1)
    , kappa(_kappa)
    , last_residual_eval_time(0)
  {
    using VectorType = Vector<double>;

    time_stepper.reinit_vector = [&](VectorType &v) { v.reinit(1); };


    time_stepper.residual = [&](const double      t,
                                const VectorType &y,
                                const VectorType &y_dot,
                                VectorType       &res) {
      deallog << "Evaluating residual at t=" << t << std::endl;

      if (t > last_residual_eval_time + 0.1)
        {
          deallog << "Time step is too large: " << t - last_residual_eval_time
                  << ". Failing recoverably." << std::endl;

          throw RecoverableUserCallbackError();
        }

      res[0]                  = y_dot[0] + kappa * y[0];
      last_residual_eval_time = t;
    };

    time_stepper.setup_jacobian = [&](const double,
                                      const VectorType &,
                                      const VectorType &,
                                      const double alpha) {
      J = kappa + alpha;
    };

    time_stepper.solve_with_jacobian =
      [&](const VectorType &src, VectorType &dst, const double) {
        dst[0] = src[0] / J;
      };

    time_stepper.output_step = [&](const double       t,
                                   const VectorType  &sol,
                                   const VectorType  &sol_dot,
                                   const unsigned int step_number) {
      deallog << "Intermediate output:" << std::endl;
      deallog << "  t =" << t << std::endl;
      deallog << "  y =" << sol[0] << "  (exact: " << std::exp(-kappa * t)
              << ')' << std::endl;
      deallog << "  y'=" << sol_dot[0]
              << "  (exact: " << -kappa * std::exp(-kappa * t) << ')'
              << std::endl;
    };
  }

  void
  run()
  {
    y[0]     = 1;
    y_dot[0] = -kappa;
    try
      {
        time_stepper.solve_dae(y, y_dot);
      }
    catch (const std::exception &exc)
      {
        deallog << "Caught an irrecoverable exception in a callback:"
                << std::endl
                << exc.what() << std::endl;
      }
  }

  SUNDIALS::IDA<Vector<double>> time_stepper;

private:
  Vector<double> y;
  Vector<double> y_dot;
  double         J;
  double         kappa;
  double         last_residual_eval_time;
};


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  SUNDIALS::IDA<Vector<double>>::AdditionalData data;
  ParameterHandler                              prm;
  data.add_parameters(prm);

  // std::ofstream ofile(SOURCE_DIR "/ida_04.prm");
  // prm.print_parameters(ofile, ParameterHandler::ShortPRM);
  // ofile.close();

  std::ifstream ifile(SOURCE_DIR "/ida_04_in.prm");
  prm.parse_input(ifile);


  HarmonicOscillator ode(2.0, data);
  ode.run();
}
