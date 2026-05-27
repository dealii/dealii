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
#include <deal.II/sundials/arkode_exception.h>
#include <deal.II/sundials/arkode_stepper.h>

#include <arkode/arkode_arkstep.h>

#if DEAL_II_SUNDIALS_VERSION_LT(6, 4, 0)
#  include <arkode/arkode_butcher_dirk.h>
#  include <arkode/arkode_butcher_erk.h>
#endif

#include <array>
#include <cmath>

#include "../tests.h"


// Test ARKStepper::AdditionalData parameters for controlling the accuracy
// order and Butcher tables of ARKStep.
//
// Two sub-tests solve a 2D Prothero-Robinson IVP
//
//   y0' = lambda * y0 + (-lambda * sin(t) + cos(t))
//   y1' = lambda * y1 + (-lambda * cos(t) - sin(t))
//
// with exact solution y0(t) = sin(t), y1(t) = cos(t). The stiff part
// f_I(y) = lambda * y  (lambda = -10) is treated implicitly with the
// smooth t-dependent forcing treated explicitly, so ARKStep runs as a
// full ImEx additive Runge-Kutta method. A jacobian_times_vector callback
// (J_I * v = lambda * v) is supplied so Newton converges in one linear
// iteration per stage.
//
// Sub-tests:
//   1. AdditionalData::order = 3 (default ARK324L2SA pair, 4 stages).
//   2. Order-5 ARK548L2SA pair (8 stages) via Butcher table names
//      (requires SUNDIALS >= 6.4).
//
// At tight tolerances the 5th-order method takes far fewer steps. Although
// each step is more expensive (more stages), the total number of function
// evaluations is significantly lower, demonstrating the benefit of high-order
// ImEx methods on smooth problems.
//
// The test behaves differently for SUNDIALS >=7.3.0 since the default method
// pairs for most of the integration orders were changed (affects sub-test 1).
// Additionally, the default time step adaptivity parameters were also changed
// in SUNDIALS 7.3.0 - this modification affects both sub-tests.

using VectorType = LinearAlgebra::distributed::Vector<double>;

struct RunResult
{
  std::array<double, 2> final_state;

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
  const auto fsteps = result.n_step_attempts < 1000 ? 10 : 100;
  deallog << "work stats:" << std::endl;
  deallog << "  steps = " << fsteps * (result.n_steps / fsteps) << std::endl;
  deallog << "  step attempts = " << fsteps * (result.n_step_attempts / fsteps)
          << std::endl;
  deallog << "  error test failures = "
          << fsteps * (result.n_error_test_failures / fsteps) << std::endl;
  deallog << "  explicit rhs evals = "
          << 100 * (result.n_explicit_rhs_evaluations / 100) << std::endl;
  deallog << "  implicit rhs evals = "
          << 100 * (result.n_implicit_rhs_evaluations / 100) << std::endl;
  deallog << "  nonlinear iterations = "
          << 100 * (result.n_nonlinear_iterations / 100) << std::endl;
  deallog << "  nonlinear convergence failures = "
          << 100 * (result.n_nonlinear_conv_failures / 100) << std::endl;
}

// Set up ARKode AdditionalData for the Prothero-Robinson benchmark.
static SUNDIALS::ARKode<VectorType>::AdditionalData
make_arkode_data(const double absolute_tolerance,
                 const double relative_tolerance)
{
  return {0.0 /*initial_time*/,
          5.0 /*final_time*/,
          1e-3 /*initial_step_size*/,
          5.0 /*output_period*/,
          1e-12 /*minimum_step_size*/,
          absolute_tolerance,
          relative_tolerance};
}


// Solve the 2D Prothero-Robinson IVP
//   y0' = lambda * y0 + (-lambda * sin(t) + cos(t))
//   y1' = lambda * y1 + (-lambda * cos(t) - sin(t))
// with exact solution y0(t) = sin(t), y1(t) = cos(t).
// The stiff diagonal part f_I(y) = lambda * y is treated implicitly;
// the smooth t-dependent forcing is treated explicitly.
// A jacobian_times_vector callback (J_I * v = lambda * v) enables Newton.
static RunResult
run(const SUNDIALS::ARKStepper<VectorType>::AdditionalData &stepper_data,
    const double                                            absolute_tolerance,
    const double                                            relative_tolerance,
    const std::function<void(void *)>                      &extra_setup = {})
{
  RunResult result;

  const double lambda = -10.0;

  SUNDIALS::ARKStepper<VectorType> stepper(stepper_data);

  stepper.implicit_function =
    [&](double, const VectorType &y, VectorType &ydot) {
      ydot[0] = lambda * y[0];
      ydot[1] = lambda * y[1];
    };

  stepper.explicit_function =
    [&](double t, const VectorType & /*y*/, VectorType &ydot) {
      ydot[0] = -lambda * std::sin(t) + std::cos(t);
      ydot[1] = -lambda * std::cos(t) - std::sin(t);
    };

  // J_I = lambda * I, so J_I * v = lambda * v. With this exact Jacobian,
  // Newton converges in a single linear iteration per stage.
  stepper.jacobian_times_vector = [&](const VectorType &v,
                                      VectorType       &Jv,
                                      const double /*t*/,
                                      const VectorType & /*y*/,
                                      const VectorType & /*fy*/) {
    Jv[0] = lambda * v[0];
    Jv[1] = lambda * v[1];
  };

  const auto arkode_data =
    make_arkode_data(absolute_tolerance, relative_tolerance);
  SUNDIALS::ARKode<VectorType> ode(stepper, arkode_data);

  void *arkode_mem = nullptr;

  ode.custom_setup = [&](void *arkode_mem_ptr) {
    arkode_mem = arkode_mem_ptr;

    int status = ARKStepSetMaxNumSteps(arkode_mem_ptr, 100000);
    AssertARKode(status);

    if (extra_setup)
      extra_setup(arkode_mem_ptr);
  };

  VectorType y(2);
  y[0] = std::sin(0.0); // = 0
  y[1] = std::cos(0.0); // = 1

  ode.solve_ode(y);

  AssertThrow(arkode_mem != nullptr, ExcInternalError());

  result.final_state = {{y[0], y[1]}};

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

  const double tol = 1e-8;

  // Sub-test 1: default order-3 ARK pair (ARK324L2SA, 4 stages).
  {
    deallog << "=== Order 3 ===" << std::endl;

    SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
    stepper_data.order = 3;

    const auto result = run(stepper_data, tol, tol);
    deallog << std::fixed << std::setprecision(6)
            << "final = " << result.final_state[0] << ' '
            << result.final_state[1] << std::endl;
    print_work_stats(result);
  }

  // Sub-test 2: order-5 ARK pair by name.
  // ARKODE_ARK548L2SA_DIRK_8_4_5 / ARKODE_ARK548L2SA_ERK_8_4_5 is an
  // 8-stage 5th-order additive Runge-Kutta pair. At tight tolerances it
  // takes far fewer steps than the 3rd-order default; the total function
  // evaluation count is lower despite the higher per-step stage cost.
  {
    deallog << "=== Butcher tables order 5 ===" << std::endl;

    SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
    std::function<void(void *)>                      extra_setup;

#if DEAL_II_SUNDIALS_VERSION_GTE(6, 4, 0)
    stepper_data.implicit_butcher_table = "ARKODE_ARK548L2SA_DIRK_8_4_5";
    stepper_data.explicit_butcher_table = "ARKODE_ARK548L2SA_ERK_8_4_5";
#else
    // ARKStepSetTableName is available in SUNDIALS 6.4.0 and later, we
    // fallback to ARKStepSetTableNum.
    extra_setup = [](void *arkode_mem_ptr) {
      int status;
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      status = ARKStepSetTableNum(arkode_mem_ptr,
                                  ARKODE_ARK548L2SA_DIRK_8_4_5,
                                  ARKODE_ARK548L2SA_ERK_8_4_5);
#  else
      status = ARKStepSetTableNum(arkode_mem_ptr,
                                  ARK548L2SA_DIRK_8_4_5,
                                  ARK548L2SA_ERK_8_4_5);
#  endif
      AssertARKode(status);
    };
#endif

    const auto result = run(stepper_data, tol, tol, extra_setup);
    deallog << std::fixed << std::setprecision(6)
            << "final = " << result.final_state[0] << ' '
            << result.final_state[1] << std::endl;
    print_work_stats(result);
  }
}
