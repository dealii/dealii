// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>

#include "../tests.h"


// Like _08, but throw a recoverable exception when solving with the linearized
// system.


int
main()
{
  initlog();

  using VectorType = Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
  stepper_data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_10_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_10_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKStepper<VectorType> stepper(stepper_data);
  SUNDIALS::ARKode<VectorType>     ode(stepper, data);

  // Explicit jacobian = stiffness matrix
  FullMatrix<double> K(2, 2);
  K(0, 0) = K(1, 1) = 0.5;
  K(1, 0) = K(0, 1) = -0.5;

  // mass matrix
  FullMatrix<double> M(2, 2);
  M(0, 0) = M(1, 1) = 2.0 / 3;
  M(1, 0) = M(0, 1) = 1.0 / 3;

  stepper.implicit_function =
    [&](double, const VectorType &y, VectorType &ydot) { K.vmult(ydot, y); };


  stepper.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) {
      ydot[0] = 1;
      ydot[1] = 2;
    };

  bool jacobian_times_vector_called   = false;
  bool solve_linearized_system_called = false;

  stepper.jacobian_times_vector = [&](const VectorType &v,
                                      VectorType       &Jv,
                                      double            t,
                                      const VectorType &y,
                                      const VectorType &fy) { K.vmult(Jv, v); };

  [&](SUNDIALS::SundialsOperator<VectorType> &,
      SUNDIALS::SundialsPreconditioner<VectorType> &,
      VectorType &,
      const VectorType &,
      double) {
    // This should not be called at all since solve_linearized_system() always
    // throws an exception, so we will check at the end that
    // jacobian_times_vector_called remains false
    jacobian_times_vector_called = true;
    throw RecoverableUserCallbackError();
  };

  stepper.solve_linearized_system =
    [&](SUNDIALS::SundialsOperator<VectorType> &,
        SUNDIALS::SundialsPreconditioner<VectorType> &,
        VectorType &,
        const VectorType &,
        double) {
      solve_linearized_system_called = true;
      throw RecoverableUserCallbackError();
    };


  stepper.solve_mass = [&](SUNDIALS::SundialsOperator<VectorType>       &op,
                           SUNDIALS::SundialsPreconditioner<VectorType> &prec,
                           VectorType                                   &x,
                           const VectorType                             &b,
                           double                                        tol) {
    SolverControl        control(100, tol);
    SolverCG<VectorType> solver_cg(control);
    solver_cg.solve(op, x, b, prec);
  };

  FullMatrix<double> M_inv(2, 2);

  bool mass_preconditioner_solve_called = false;
  bool mass_preconditioner_setup_called = false;

  stepper.mass_preconditioner_solve =
    [&](double t, const VectorType &r, VectorType &z, double gamma, int lr) {
      mass_preconditioner_solve_called = true;
      M_inv.vmult(z, r);
    };

  stepper.mass_preconditioner_setup = [&](double t) {
    mass_preconditioner_setup_called = true;
    M_inv.invert(M);
  };

  stepper.mass_times_vector = [&](const double      t,
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
  deallog << std::boolalpha;
  deallog << "jacobian_times_vector_called : " << jacobian_times_vector_called
          << std::endl;
  deallog << "solve_linearized_system_called : "
          << solve_linearized_system_called << std::endl;
  deallog << "mass_preconditioner_setup_called : "
          << mass_preconditioner_setup_called << std::endl;
  deallog << "mass_preconditioner_solve_called : "
          << mass_preconditioner_solve_called << std::endl;
}
