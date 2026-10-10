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

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>

#include "../tests.h"

// Test that ARKode::solve_ode() returns the total number of steps taken, also
// when the solver is restarted via solver_should_restart(). A restart
// re-creates the ARKODE memory, which resets its internal step counter. Use a
// fixed step size so that restarting does not change the step sequence: the
// step count must then be identical with and without restarts.

unsigned int
run(const bool with_restarts)
{
  using VectorType = Vector<double>;

  SUNDIALS::ARKode<VectorType>::AdditionalData data;

  // Powers of two keep all times exactly representable: 128 steps in total,
  // 4 steps between outputs.
  const double step_size = 0.125;
  data.initial_time      = 0.0;
  data.final_time        = 8.0;
  data.initial_step_size = step_size;
  data.output_period     = 0.25;

  SUNDIALS::ARKStepper<VectorType> stepper;
  SUNDIALS::ARKode<VectorType>     ode(stepper, data);

  const double kappa = 1.0;

  stepper.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) {
      ydot[0] = y[1];
      ydot[1] = -kappa * kappa * y[0];
    };

  // custom_setup() is also called after every restart.
  ode.custom_setup = [&](void *arkode_mem) {
#if DEAL_II_SUNDIALS_VERSION_GTE(7, 1, 0)
    const int status = ARKodeSetFixedStep(arkode_mem, step_size);
#else
    const int status = ARKStepSetFixedStep(arkode_mem, step_size);
#endif
    AssertThrow(status == 0, ExcInternalError());
  };

  unsigned int n_restarts        = 0;
  double       next_restart_time = 2.0;
  if (with_restarts)
    {
      // Restart each time a multiple of 2 is reached.
      ode.solver_should_restart = [&](const double t, VectorType &) -> bool {
        if (t >= next_restart_time)
          {
            next_restart_time += 2.0;
            ++n_restarts;
            return true;
          }
        return false;
      };
    }

  Vector<double> y(2);
  y[0]                       = 0;
  y[1]                       = kappa;
  const unsigned int n_steps = ode.solve_ode(y);

  deallog << "Restarts: " << n_restarts << std::endl;
  return n_steps;
}


int
main()
{
  initlog();

  const unsigned int n_steps_no_restart   = run(false);
  const unsigned int n_steps_with_restart = run(true);

  deallog << "Steps without restarts: " << n_steps_no_restart << std::endl;
  deallog << "Steps with restarts:    " << n_steps_with_restart << std::endl;

  if (n_steps_no_restart == n_steps_with_restart)
    deallog << "OK" << std::endl;
  else
    deallog << "Number of steps differs!" << std::endl;
}
