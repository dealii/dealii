// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include <arkode/arkode_arkstep.h>

#include <array>

#include "../tests.h"


// Test ARKStepper::AdditionalData parameters for controlling the accuracy
// order and Butcher tables of ARKStep.
//
// Based on the Brusselator benchmark (arkode_03). Two sub-tests are run with
// different stepper configurations:
//   1. Explicitly set order of accuracy = 3 via AdditionalData::order.
//   2. Select order-2 paired Butcher tables (ARKODE_ARK2_DIRK_3_1_2 /
//      ARKODE_ARK2_ERK_3_1_2) via AdditionalData::implicit_butcher_table and
//      AdditionalData::explicit_butcher_table.
//
// Both runs solve the same stiff Brusselator IVP. We compare final states and
// solver work statistics to make method differences visible.

using VectorType = LinearAlgebra::distributed::Vector<double>;

struct RunResult
{
  std::array<double, 3> final_state;

  long int n_steps                    = -1;
  long int n_step_attempts            = -1;
  long int n_error_test_failures      = -1;
  long int n_explicit_rhs_evaluations = -1;
  long int n_implicit_rhs_evaluations = -1;
  long int n_nonlinear_iterations     = -1;
  long int n_nonlinear_conv_failures  = -1;
};


static void
print_work_stats(const RunResult &result)
{
  deallog << "work stats:" << std::endl;
  deallog << "  steps = " << result.n_steps << std::endl;
  deallog << "  step attempts = " << result.n_step_attempts << std::endl;
  deallog << "  error test failures = " << result.n_error_test_failures
          << std::endl;
  deallog << "  explicit rhs evals = " << result.n_explicit_rhs_evaluations
          << std::endl;
  deallog << "  implicit rhs evals = " << result.n_implicit_rhs_evaluations
          << std::endl;
  deallog << "  nonlinear iterations = " << result.n_nonlinear_iterations
          << std::endl;
  deallog << "  nonlinear convergence failures = "
          << result.n_nonlinear_conv_failures << std::endl;

  if (result.n_step_attempts > 0 && result.n_steps >= 0)
    {
      const double rejection_rate =
        static_cast<double>(result.n_step_attempts - result.n_steps) /
        static_cast<double>(result.n_step_attempts);
      deallog << "  rejection rate = " << std::setprecision(6) << rejection_rate
              << std::endl;
    }
}

// Set up ARKode AdditionalData for the Brusselator benchmark.
static SUNDIALS::ARKode<VectorType>::AdditionalData
make_arkode_data(const double absolute_tolerance,
                 const double relative_tolerance)
{
  return {0.0 /*initial_time*/,
          3.0 /*final_time*/,
          1e-6 /*initial_step_size*/,
          0.1 /*output_period*/,
          1e-7 /*minimum_step_size*/,
          absolute_tolerance,
          relative_tolerance};
}


// Run the Brusselator benchmark with the given stepper.
static RunResult
run(const SUNDIALS::ARKStepper<VectorType>::AdditionalData &stepper_data,
    const double                                            absolute_tolerance,
    const double                                            relative_tolerance)
{
  RunResult result;

  SUNDIALS::ARKStepper<VectorType> stepper(stepper_data);

  const double u0 = 3.9, v0 = 1.1, w0 = 2.8;
  const double a = 1.2, b = 2.5, eps = 1e-5;

  stepper.implicit_function =
    [&](double, const VectorType &y, VectorType &ydot) {
      ydot[0] = 0;
      ydot[1] = 0;
      ydot[2] = -y[2] / eps;
    };

  stepper.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) {
      ydot[0] = a - (y[2] + 1) * y[0] + y[1] * y[0] * y[0];
      ydot[1] = y[2] * y[0] - y[1] * y[0] * y[0];
      ydot[2] = b / eps - y[2] * y[0];
    };

  const auto arkode_data =
    make_arkode_data(absolute_tolerance, relative_tolerance);
  SUNDIALS::ARKode<VectorType> ode(stepper, arkode_data);

  void *arkode_mem = nullptr;

  ode.custom_setup = [&](void *arkode_mem_ptr) {
    arkode_mem = arkode_mem_ptr;

    int ierr = ARKStepSetMaxNumSteps(arkode_mem_ptr, 20000);
    AssertThrow(ierr == 0, ExcInternalError());
  };

  VectorType y(3);
  y[0] = u0;
  y[1] = v0;
  y[2] = w0;

  ode.solve_ode(y);

  AssertThrow(arkode_mem != nullptr, ExcInternalError());

  result.final_state = {{y[0], y[1], y[2]}};

  int ierr = ARKStepGetNumSteps(arkode_mem, &result.n_steps);
  AssertThrow(ierr == 0, ExcInternalError());

#if DEAL_II_SUNDIALS_VERSION_GTE(5, 2, 0)
  ierr = ARKStepGetNumStepAttempts(arkode_mem, &result.n_step_attempts);
  AssertThrow(ierr == 0, ExcInternalError());

  ierr = ARKStepGetNumErrTestFails(arkode_mem, &result.n_error_test_failures);
  AssertThrow(ierr == 0, ExcInternalError());

  ierr = ARKStepGetNumRhsEvals(arkode_mem,
                               &result.n_explicit_rhs_evaluations,
                               &result.n_implicit_rhs_evaluations);
  AssertThrow(ierr == 0, ExcInternalError());

  ierr =
    ARKStepGetNumNonlinSolvIters(arkode_mem, &result.n_nonlinear_iterations);
  AssertThrow(ierr == 0, ExcInternalError());

  ierr = ARKStepGetNumNonlinSolvConvFails(arkode_mem,
                                          &result.n_nonlinear_conv_failures);
  AssertThrow(ierr == 0, ExcInternalError());
#endif

  return result;
}


int
main()
{
  initlog();

  auto make_stepper_data = []() {
    SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
    stepper_data.implicit_function_is_linear           = true;
    stepper_data.implicit_function_is_time_independent = true;
    stepper_data.maximum_non_linear_iterations         = 50;
    stepper_data.anderson_acceleration_subspace        = 10;
    return stepper_data;
  };

  // Sub-test 1: set the overall order of accuracy to 3.
  {
    deallog << "=== Order 3 ===" << std::endl;

    auto stepper_data  = make_stepper_data();
    stepper_data.order = 3;

    const auto order_3_result = run(stepper_data, 1e-4, 1e-4);
    deallog << std::fixed << std::setprecision(6)
            << "final = " << order_3_result.final_state[0] << ' '
            << order_3_result.final_state[1] << ' '
            << order_3_result.final_state[2] << std::endl;
    print_work_stats(order_3_result);
  }

  // Sub-test 2: select order-2 paired ARK Butcher tables explicitly.
  // ARKODE_ARK2_DIRK_3_1_2 (implicit) and ARKODE_ARK2_ERK_3_1_2 (explicit)
  // form a second-order additive Runge-Kutta pair.
  {
    deallog << "=== Butcher tables order 2 ===" << std::endl;

    auto stepper_data                   = make_stepper_data();
    stepper_data.implicit_butcher_table = "ARKODE_ARK2_DIRK_3_1_2";
    stepper_data.explicit_butcher_table = "ARKODE_ARK2_ERK_3_1_2";

    const auto butcher_2_result = run(stepper_data, 1e-4, 1e-4);
    deallog << std::fixed << std::setprecision(6)
            << "final = " << butcher_2_result.final_state[0] << ' '
            << butcher_2_result.final_state[1] << ' '
            << butcher_2_result.final_state[2] << std::endl;
    print_work_stats(butcher_2_result);
  }
}
