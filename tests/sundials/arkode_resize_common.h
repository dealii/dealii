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

#ifndef dealii_tests_sundials_arkode_resize_common_h
#define dealii_tests_sundials_arkode_resize_common_h

// Common driver for the arkode_resize_* tests of
// ARKode::solver_should_resize().
//
// The ODE y_i' = -y_i is solved with all components being identical. Doubling
// the number of components at t = 1 does not change the weighted error norm,
// so a resize that preserves the adaptivity state must reproduce the step
// sequence of the run without resize. ARKODE restarts from the output time at
// a resize, which can cost one extra step.

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include <cmath>
#include <memory>
#include <string>

#include "../tests.h"

using VectorType = Vector<double>;


struct ResizeResult
{
  unsigned int n_steps    = 0;
  unsigned int n_resizes  = 0;
  unsigned int final_size = 0;
  double       max_error  = 0.;
};


inline void
rhs(double, const VectorType &y, VectorType &ydot)
{
  ydot = y;
  ydot *= -1.;
}


inline ResizeResult
integrate(SUNDIALS::ARKode<VectorType> &ode, const bool with_resize)
{
  ResizeResult result;

  if (with_resize)
    ode.solver_should_resize = [&](const double t, VectorType &y) -> double {
      if (result.n_resizes > 0 || t < 1.0)
        return 0.;

      const unsigned int old_size = y.size();
      y.grow_or_shrink(2 * old_size);
      for (unsigned int i = 0; i < old_size; ++i)
        y[old_size + i] = y[i];
      ++result.n_resizes;
      return 1.;
    };

  ode.output_step =
    [&](const double t, const VectorType &y, const unsigned int) {
      for (unsigned int i = 0; i < y.size(); ++i)
        result.max_error =
          std::max(result.max_error, std::abs(y[i] - std::exp(-t)));
    };

  VectorType y(2);
  y = 1.;

  result.n_steps    = ode.solve_ode(y);
  result.final_size = y.size();
  return result;
}


// Integrate with and without resize using a StepperType that is set up by
// `configure` and log the comparison under `name`.
template <typename StepperType, typename Configure>
void
run_resize_test(const std::string &name, Configure configure)
{
  const auto run = [&](const bool with_resize) {
    SUNDIALS::ARKode<VectorType>::AdditionalData data;
    data.initial_time       = 0.0;
    data.final_time         = 2.0;
    data.initial_step_size  = 1e-4;
    data.output_period      = 0.5;
    data.absolute_tolerance = 1e-8;
    data.relative_tolerance = 1e-6;

    StepperType stepper;
    configure(stepper);
    SUNDIALS::ARKode<VectorType> ode(stepper, data);

    return integrate(ode, with_resize);
  };

  const ResizeResult reference = run(false);
  const ResizeResult resized   = run(true);

  deallog << name << ": resizes " << resized.n_resizes << ", final size "
          << resized.final_size << std::endl;
  deallog << name << ": step count is within one of the run without resize: "
          << (resized.n_steps >= reference.n_steps &&
                  resized.n_steps <= reference.n_steps + 1 ?
                "yes" :
                "no")
          << std::endl;
  deallog << name
          << ": error is small: " << (resized.max_error < 1e-4 ? "yes" : "no")
          << std::endl;
}


// Same as above for a stepper that only needs the explicit right hand side.
template <typename StepperType>
void
run_resize_test(const std::string &name)
{
  run_resize_test<StepperType>(name, [](StepperType &stepper) {
    stepper.explicit_function = rhs;
  });
}

#endif
