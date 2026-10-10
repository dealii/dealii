// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_tests_sundials_arkode_restart_common_h
#define dealii_tests_sundials_arkode_restart_common_h

// Common driver for the arkode_restart_* tests of
// ARKode::solver_should_restart().
//
// ARKode::solve_ode() must return the total number of steps taken, also when
// the solver is restarted. A restart re-creates the ARKODE memory, which resets
// its internal step counter, and the old objects have to be released before
// the new ones are created. The ODE y' = -y is integrated with and without
// restarts. A fixed step size is used so that restarting does not change the
// step sequence: the step count must then be identical in both runs.

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include <cmath>
#include <memory>
#include <string>
#include <type_traits>

#include "../tests.h"

using VectorType = Vector<double>;


struct RestartResult
{
  unsigned int n_steps    = 0;
  unsigned int n_restarts = 0;
  double       max_error  = 0.;
};


inline void
rhs(double, const VectorType &y, VectorType &ydot)
{
  ydot = y;
  ydot *= -1.;
}


// Set a fixed step size on the ARKODE memory of a stepper of type Stepper.
template <typename Stepper>
void
set_fixed_step(void *arkode_mem, const double step_size)
{
#if DEAL_II_SUNDIALS_VERSION_GTE(7, 1, 0)
  const int status = ARKodeSetFixedStep(arkode_mem, step_size);
#else
  const int status = std::is_same_v<Stepper, SUNDIALS::ERKStepper<VectorType>> ?
                       ERKStepSetFixedStep(arkode_mem, step_size) :
                       ARKStepSetFixedStep(arkode_mem, step_size);
#endif
  AssertThrow(status == 0, ExcInternalError());
}


template <typename Stepper>
RestartResult
integrate(SUNDIALS::ARKode<VectorType> &ode,
          const double                  step_size,
          const bool                    with_restarts)
{
  RestartResult result;

  // Called after every restart as well.
  ode.custom_setup = [&](void *arkode_mem) {
    set_fixed_step<Stepper>(arkode_mem, step_size);
  };

  // Restart each time a multiple of 2 is reached.
  double next_restart_time = 2.0;
  if (with_restarts)
    ode.solver_should_restart = [&](const double t, VectorType &) -> bool {
      if (t >= next_restart_time)
        {
          next_restart_time += 2.0;
          ++result.n_restarts;
          return true;
        }
      return false;
    };

  ode.output_step =
    [&](const double t, const VectorType &y, const unsigned int) {
      for (unsigned int i = 0; i < y.size(); ++i)
        result.max_error =
          std::max(result.max_error, std::abs(y[i] - std::exp(-t)));
    };

  VectorType y(2);
  y = 1.;

  result.n_steps = ode.solve_ode(y);
  return result;
}


// Run with and without restarts using the stepper that `make_stepper` returns,
// as a `std::unique_ptr`, and log the comparison under `name`. The step size is
// fixed for all steppers.
template <typename MakeStepper>
void
run_restart_test_with(const std::string &name,
                      MakeStepper        make_stepper,
                      const double       step_size = 0.125)
{
  using Stepper = typename decltype(make_stepper())::element_type;

  const auto run = [&](const bool with_restarts) {
    // Powers of two keep all times exactly representable.
    SUNDIALS::ARKode<VectorType>::AdditionalData data;
    data.initial_time       = 0.0;
    data.final_time         = 8.0;
    data.initial_step_size  = step_size;
    data.output_period      = 0.25;
    data.absolute_tolerance = 1e-8;
    data.relative_tolerance = 1e-6;

    const auto                   stepper = make_stepper();
    SUNDIALS::ARKode<VectorType> ode(*stepper, data);

    return integrate<Stepper>(ode, step_size, with_restarts);
  };

  const RestartResult reference = run(false);
  const RestartResult restarted = run(true);

  deallog << name << ": restarts " << restarted.n_restarts << std::endl;
  deallog << name << ": same number of steps as without restarts: "
          << (restarted.n_steps == reference.n_steps ? "yes" : "no")
          << std::endl;
  deallog << name
          << ": error is small: " << (restarted.max_error < 1e-2 ? "yes" : "no")
          << std::endl;
}


// Same as above for a default-constructible StepperType that is set up by
// `configure`.
template <typename StepperType, typename Configure>
void
run_restart_test(const std::string &name, Configure configure)
{
  run_restart_test_with(name, [&] {
    auto stepper = std::make_unique<StepperType>();
    configure(*stepper);
    return stepper;
  });
}


// Same as above for a stepper that only needs the explicit right hand side.
template <typename StepperType>
void
run_restart_test(const std::string &name)
{
  run_restart_test<StepperType>(name, [](StepperType &stepper) {
    stepper.explicit_function = rhs;
  });
}

#endif
