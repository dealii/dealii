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
// calls to solve_linearized_system or jacobian_times_vector before performing
// the actual operation. This exercises the recoverable-error code path in
// ARKODE's linear-solver interface: each throw forces a step-size reduction and
// a retry, so the integration still completes.
//
// The problem is run three times:
//   Sub-test 1: no artificial exceptions (reference run).
//   Sub-test 2: recoverable exceptions from solve_linearized_system.
//   Sub-test 3: recoverable exceptions thrown directly from
//               jacobian_times_vector.
//               With SundialsOperator::vmult now re-throwing the positive
//               return from the ATimes callback as
//               RecoverableUserCallbackError, the exception propagates through
//               the linear solver out of solve_linearized_system, where it is
//               mapped to a recognized recoverable error code.
//
// At the end, the step counts from sub-tests 2 and 3 are compared with the
// reference: forced step-size reductions must larger number of timesteps.

using VectorType = Vector<double>;

enum class ThrowMode
{
  None,
  FromLinearSolve,
  FromJacTimesVec
};

// Run the problem once. Returns the number of time steps taken
// (0 if an exception was caught).
unsigned int
run_problem(const ThrowMode throw_mode)
{
  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
  stepper_data.add_parameters(prm);

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

  unsigned int n_jac_times_vec_runs = 0;

  stepper.jacobian_times_vector = [&](const VectorType &v,
                                      VectorType       &Jv,
                                      double            t,
                                      const VectorType &y,
                                      const VectorType &fy) {
    ++n_jac_times_vec_runs;

    if (throw_mode == ThrowMode::FromJacTimesVec &&
        (n_jac_times_vec_runs == 3 || n_jac_times_vec_runs == 12))
      {
        deallog << "Reporting recoverable failure from jacobian_times_vector."
                << std::endl;

        throw RecoverableUserCallbackError();
      }

    K.vmult(Jv, v);
  };

  unsigned int n_linear_solve_runs = 0;

  stepper.solve_linearized_system =
    [&](SUNDIALS::SundialsOperator<VectorType>       &op,
        SUNDIALS::SundialsPreconditioner<VectorType> &prec,
        VectorType                                   &x,
        const VectorType                             &b,
        double                                        tol) {
      ++n_linear_solve_runs;

      if (throw_mode == ThrowMode::FromLinearSolve &&
          (n_linear_solve_runs == 3 || n_linear_solve_runs == 12))
        {
          deallog << "Reporting recoverable failure from "
                     "solve_linearized_system."
                  << std::endl;
          throw RecoverableUserCallbackError();
        }

      SolverControl        control(100, tol);
      SolverCG<VectorType> solver_cg(control);
      solver_cg.solve(op, x, b, prec);
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

      std::ofstream ofile(SOURCE_DIR "/arkode_10_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  // Sub-test 1: run the problem without any artificial exceptions.
  deallog << "Sub-test 1: no artificial recoverable errors" << std::endl;
  const unsigned int n_steps_ref = run_problem(ThrowMode::None);

  // Sub-test 2: run the same problem but report recoverable errors from
  // solve_linearized_system at a few well-separated call counts, forcing
  // ARKODE to reduce the step size and retry.
  deallog << "Sub-test 2: recoverable errors from solve_linearized_system"
          << std::endl;
  const unsigned int n_steps_linear_solve =
    run_problem(ThrowMode::FromLinearSolve);

  // Sub-test 3: same, but the exception is thrown from jacobian_times_vector
  // and propagates through the linear solver.
  deallog << "Sub-test 3: recoverable errors from jacobian_times_vector"
          << std::endl;
  const unsigned int n_steps_jac_times_vec =
    run_problem(ThrowMode::FromJacTimesVec);

  deallog << std::boolalpha;
  deallog << "Sub-test 2 required more timesteps: "
          << (n_steps_linear_solve - n_steps_ref > 0) << std::endl;
  deallog << "Sub-test 3 required more timesteps: "
          << (n_steps_jac_times_vec - n_steps_ref > 0) << std::endl;
}
