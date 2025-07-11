// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>

#include "../tests.h"


// Like _08, but throw a recoverable exception when solving with the mass
// matrix.


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
      std::ofstream ofile(SOURCE_DIR "/arkode_11_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_08_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  // Explicit jacobian = stiffness matrix
  FullMatrix<double> K(2, 2);
  K(0, 0) = K(1, 1) = 0.5;
  K(1, 0) = K(0, 1) = -0.5;

  // mass matrix
  FullMatrix<double> M(2, 2);
  M(0, 0) = M(1, 1) = 2.0 / 3;
  M(1, 0) = M(0, 1) = 1.0 / 3;

  ode.implicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    K.vmult(ydot, y);
  };


  ode.explicit_function = [&](double, const VectorType &y, VectorType &ydot) {
    ydot[0] = 1;
    ydot[1] = 2;
  };

  ode.jacobian_times_vector = [&](const VectorType &v,
                                  VectorType       &Jv,
                                  double            t,
                                  const VectorType &y,
                                  const VectorType &fy) { K.vmult(Jv, v); };

  [&](SUNDIALS::SundialsOperator<VectorType> &,
      SUNDIALS::SundialsPreconditioner<VectorType> &,
      VectorType &,
      const VectorType &,
      double) {
    deallog << "Reporting recoverable failure." << std::endl;
    throw RecoverableUserCallbackError();
  };

  ode.solve_mass = [&](SUNDIALS::SundialsOperator<VectorType> &,
                       SUNDIALS::SundialsPreconditioner<VectorType> &,
                       VectorType &,
                       const VectorType &,
                       double) {
    deallog
      << "Reporting recoverable failure when solving with the mass matrix."
      << std::endl;
    throw RecoverableUserCallbackError();
  };


  ode.solve_linearized_system =
    [&](SUNDIALS::SundialsOperator<VectorType>       &op,
        SUNDIALS::SundialsPreconditioner<VectorType> &prec,
        VectorType                                   &x,
        const VectorType                             &b,
        double                                        tol) {
      ReductionControl     control;
      SolverCG<VectorType> solver_cg(control);
      solver_cg.solve(op, x, b, prec);
    };

  FullMatrix<double> M_inv(2, 2);

  ode.mass_preconditioner_solve =
    [&](double t, const VectorType &r, VectorType &z, double gamma, int lr) {
      LogStream::Prefix prefix("mass_preconditioner_solve");
      deallog << "applied" << std::endl;
      M_inv.vmult(z, r);
    };

  ode.mass_preconditioner_setup = [&](double t) {
    LogStream::Prefix prefix("mass_preconditioner_setup");
    deallog << "applied" << std::endl;
    M_inv.invert(M);
  };

  ode.mass_times_vector = [&](const double      t,
                              const VectorType &v,
                              VectorType       &Mv) { M.vmult(Mv, v); };


  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      deallog << t << ' ' << sol[0] << ' ' << sol[1] << std::endl;
    };

  Vector<double> y(2);
  y[0] = 1;
  y[1] = 0;

  try
    {
      ode.solve_ode(y);
    }
  catch (const std::exception &exc)
    {
      deallog << "Caught exception:" << std::endl << exc.what() << std::endl;
    }
}
