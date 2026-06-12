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

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_exception.h>
#include <deal.II/sundials/arkode_stepper.h>

#include <arkode/arkode_erkstep.h>

#if DEAL_II_SUNDIALS_VERSION_LT(6, 4, 0)
#  include <arkode/arkode_butcher_erk.h>
#endif

#include <cmath>

#include "../tests.h"


// Test ERKStepper functionality on the harmonic oscillator problem (cf.
// arkode_01). Three sub-tests are run:
//   1. Default ERKStepper settings.
//   2. ERKStepper with the order of accuracy set to 3 via
//      AdditionalData::order. Order 3 is used here because it is the only
//      order whose default Butcher table (BOGACKI_SHAMPINE_4_2_3) is
//      unchanged across all supported SUNDIALS versions (5.8.0 through 7.7.0).
//      All other orders map to different default tables in SUNDIALS 7.x.
//   3. ERKStepper with an explicit Butcher table selected by name via
//      AdditionalData::explicit_butcher_table (ARKODE_DORMAND_PRINCE_7_4_5,
//      a classical order-5 explicit Runge-Kutta method).
//
// The harmonic oscillator problem:
//
//   u'' = -k^2 u,  u(0) = 0,  u'(0) = k
//
// written as a first-order system:
//
//   y[0]' =        y[1]
//   y[1]' = -k^2 * y[0]
//
// Exact solution: y[0](t) = sin(k t), y[1](t) = k cos(k t).

using VectorType = Vector<double>;

// Common ARKode outer settings: integrate [0, 2*pi] outputting every 0.05
// time units, matching the tolerance and step-size settings of arkode_01.
static SUNDIALS::ARKode<VectorType>::AdditionalData
make_arkode_data()
{
  return {0.0 /*initial_time*/,
          6.28 /*final_time*/,
          0.01 /*initial_step_size*/,
          0.05 /*output_period*/,
          1e-6 /*minimum_step_size*/,
          1e-9 /*absolute_tolerance*/,
          1e-8 /*relative_tolerance*/};
}


static void
run(const SUNDIALS::ERKStepper<VectorType>::AdditionalData &data,
    const std::function<void(void *)>                      &extra_setup = {})
{
  SUNDIALS::ERKStepper<VectorType> stepper(data);

  const double kappa = 1.0;

  stepper.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) {
      ydot[0] = y[1];
      ydot[1] = -kappa * kappa * y[0];
    };

  SUNDIALS::ARKode<VectorType> ode(stepper, make_arkode_data());

  if (extra_setup)
    ode.custom_setup = extra_setup;

  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int /*step*/) {
      // Exact solution of the harmonic oscillator: y[0] = sin(k t),
      // y[1] = k cos(k t).
      VectorType exact(2);
      exact[0] = std::sin(kappa * t);
      exact[1] = kappa * std::cos(kappa * t);

      VectorType diff(exact);
      diff -= sol;

      deallog << t << ' ' << std::fixed << std::setprecision(4) << sol[0] << ' '
              << sol[1] << ' ' << exact[0] << ' ' << exact[1] << ' '
              << diff.l2_norm() << std::endl;
    };

  VectorType y(2);
  y[0] = 0;
  y[1] = kappa;

  ode.solve_ode(y);
}


int
main()
{
  initlog();

  // Sub-test 1: default ERKStepper settings.
  {
    deallog << "=== Default ===" << std::endl;
    SUNDIALS::ERKStepper<VectorType>::AdditionalData data;
    run(data);
  }

  // Sub-test 2: explicitly request order 3.
  // BOGACKI_SHAMPINE_4_2_3 is the default table for order 3 in all supported
  // SUNDIALS versions (5.8.0 through 7.7.0), making this sub-test stable
  // across versions without pinning a specific table by name.
  {
    deallog << "=== Order 3 ===" << std::endl;
    SUNDIALS::ERKStepper<VectorType>::AdditionalData data;
    data.order = 3;
    run(data);
  }

  // Sub-test 3: select a specific ERK Butcher table by name (SUNDIALS >= 6.4)
  // or by number (fallback).
  // ARKODE_DORMAND_PRINCE_7_4_5 is the classical Dormand-Prince 4(5) method.
  {
    deallog << "=== Butcher table ARKODE_DORMAND_PRINCE_7_4_5 ===" << std::endl;
    SUNDIALS::ERKStepper<VectorType>::AdditionalData data;
    std::function<void(void *)>                      extra_setup;

#if DEAL_II_SUNDIALS_VERSION_GTE(6, 4, 0)
    data.explicit_butcher_table = "ARKODE_DORMAND_PRINCE_7_4_5";
#else
    // ERKStepSetTableName is available in SUNDIALS 6.4.0 and later, we
    // fallback to ERKStepSetTableNum.
    extra_setup = [](void *arkode_mem_ptr) {
      int status;
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      status = ERKStepSetTableNum(arkode_mem_ptr, ARKODE_DORMAND_PRINCE_7_4_5);
#  else
      status = ERKStepSetTableNum(arkode_mem_ptr, DORMAND_PRINCE_7_4_5);
#  endif
      AssertARKode(status);
    };
#endif

    run(data, extra_setup);
  }
}
