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


// Like _08, but optionally report a recoverable error at a few, well-separated
// attempts to solve with the mass matrix before performing the actual solve.
// This exercises the recoverable-error code path while still allowing the
// integrator to reduce the step size, retry, and complete the time
// integration.
//
// The problem is run twice: first without any artificial exceptions
// (sub-test 1, the reference run) and then with the recoverable exceptions
// enabled (sub-test 2). Both runs write to the same output stream, and at the
// end the number of time steps taken in the two cases is compared: the
// recoverable failures force additional step-size reductions and retries, so
// the second run is expected to take at least as many steps as the first.

using VectorType = Vector<double>;

// Run the problem once. If @p throw_recoverable_errors is true, report a
// recoverable error at a few, well-separated mass-matrix solve attempts.
// Returns the number of time steps taken (0 if an exception was caught).
unsigned int
run_problem(const bool throw_recoverable_errors)
{
  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
  stepper_data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/arkode_11_in.prm");
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

  stepper.jacobian_times_vector = [&](const VectorType &v,
                                      VectorType       &Jv,
                                      double            t,
                                      const VectorType &y,
                                      const VectorType &fy) { K.vmult(Jv, v); };

  // Occasionally report a recoverable failure when solving with the mass
  // matrix to exercise ARKODE's recovery path.
  //
  // SUNDIALS "recovers" from a recoverable mass-solve failure by reducing the
  // step size and re-attempting the whole step (which calls solve_mass again).
  // The failures must therefore be sparse: they are triggered at a few fixed,
  // well-separated solve counts so that no single step accumulates more than
  // MAXNCF (=10) convergence failures (which would abort the integration).
  unsigned int n_mass_solve_runs = 0;

  stepper.solve_mass = [&](SUNDIALS::SundialsOperator<VectorType>       &op,
                           SUNDIALS::SundialsPreconditioner<VectorType> &prec,
                           VectorType                                   &x,
                           const VectorType                             &b,
                           double                                        tol) {
    ++n_mass_solve_runs;

    if (throw_recoverable_errors &&
        (n_mass_solve_runs == 30 || n_mass_solve_runs == 55 ||
         n_mass_solve_runs == 85))
      {
        deallog << "Reporting recoverable failure when solving with the "
                   "mass matrix."
                << std::endl;
        throw RecoverableUserCallbackError();
      }

    SolverControl        control(100, tol);
    SolverCG<VectorType> solver_cg(control);
    solver_cg.solve(op, x, b, prec);
  };


  stepper.solve_linearized_system =
    [&](SUNDIALS::SundialsOperator<VectorType>       &op,
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
      M_inv.vmult(z, r);
      mass_preconditioner_solve_called = true;
    };

  stepper.mass_preconditioner_setup = [&](double t) {
    M_inv.invert(M);
    mass_preconditioner_setup_called = true;
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
      const unsigned int n_timesteps = ode.solve_ode(y);

      deallog << std::boolalpha;
      deallog << "mass_preconditioner_setup_called: "
              << mass_preconditioner_setup_called << std::endl;
      deallog << "mass_preconditioner_solve_called: "
              << mass_preconditioner_solve_called << std::endl;

      return n_timesteps;
    }
  catch (const std::exception &exc)
    {
      deallog << "Caught exception:" << std::endl << exc.what() << std::endl;
      return 0;
    }
}


int
main()
{
  initlog();

  // Optionally regenerate the input parameter file.
  if (false)
    {
      ParameterHandler                             prm;
      SUNDIALS::ARKode<VectorType>::AdditionalData data;
      data.add_parameters(prm);

      SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
      stepper_data.add_parameters(prm);

      std::ofstream ofile(SOURCE_DIR "/arkode_11_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  // Sub-test 1: run the problem without any artificial exceptions. This is the
  // reference run.
  deallog << "Sub-test 1: no artificial recoverable errors" << std::endl;
  const unsigned int n_timesteps_reference =
    run_problem(/* throw_recoverable_errors = */ false);

  // Sub-test 2: run the same problem but report recoverable errors at a few
  // mass-matrix solves, forcing ARKODE to reduce the step size and retry.
  deallog << "Sub-test 2: with artificial recoverable errors" << std::endl;
  const unsigned int n_timesteps_with_errors =
    run_problem(/* throw_recoverable_errors = */ true);

  // Compare the number of time steps taken in the two cases. The recoverable
  // failures force additional step-size reductions and retries, so the second
  // run takes more steps than the reference run.
  deallog << std::boolalpha;
  deallog << "Sub-test 2 required more timesteps: "
          << (n_timesteps_with_errors - n_timesteps_reference > 2) << std::endl;
}
