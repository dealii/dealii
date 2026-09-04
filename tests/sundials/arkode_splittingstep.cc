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

#include <deal.II/base/config.h>

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include <cmath>
#include <memory>

#include "../tests.h"

// Reproduces the SUNDIALS example ark_analytic_partitioned.c using deal.II
// wrappers for operator-splitting (SplittingStepper).
//
// Problem (from Estep et al., SIAM J. Numer. Anal. 46 (2008)):
//   y' = f1(t,y) + f2(t,y),   y(0) = 1,   t in [0, 1]
//   f1(t,y) = -lambda * y      (linear partition, lambda = 2)
//   f2(t,y) =  y^2             (nonlinear partition)
//
// Exact solution:
//   y(t) = lambda * y(0) / (y(0) - (y(0) - lambda) * exp(lambda * t))
//         = 2 / (1 + exp(2*t))
//
// The two partitions are integrated by independent sub-steppers:
//   Partition 1 (linear):    ERKStepper, adaptive timestep
//   Partition 2 (nonlinear): ARKStepper, adaptive timestep
// The outer SplittingStep uses a fixed step dt = 0.01.
//
// SplittingStep does not support adaptive outer time-stepping; the fixed
// outer step is set via SplittingStepper::AdditionalData::step_size.
// Sub-steppers may use adaptive or fixed stepping at their own discretion.
//
// Two sub-tests:
//   1. Lie-Trotter splitting (1st order, SUNDIALS default).
//   2. Strang splitting (2nd order, "ARKODE_SPLITTING_STRANG_2_2_2").

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

using VectorType = Vector<double>;

static constexpr double lambda        = 2.0;
static constexpr double t_start       = 0.0;
static constexpr double t_end         = 1.0;
static constexpr double dt            = 0.01;
static constexpr double output_period = 0.25;

static double
exact_solution(double t)
{
  return lambda / (1.0 - (1.0 - lambda) * std::exp(lambda * t));
}

static void
run(const SUNDIALS::SplittingStepper<VectorType>::AdditionalData &split_data,
    const std::string                                            &label)
{
  deallog << "=== " << label << " ===" << std::endl;

  // Partition 1: linear, f1(t,y) = -lambda*y, integrated with ERKStepper
  auto linear_stepper = std::make_shared<SUNDIALS::ERKStepper<VectorType>>();
  linear_stepper->explicit_function =
    [](const double /*t*/, const VectorType &y, VectorType &ydot) {
      ydot[0] = -lambda * y[0];
    };

  // Partition 2: nonlinear, f2(t,y) = y^2, integrated with ARKStepper
  auto nonlinear_stepper = std::make_shared<SUNDIALS::ARKStepper<VectorType>>();
  nonlinear_stepper->implicit_function =
    [](const double /*t*/, const VectorType &y, VectorType &ydot) {
      ydot[0] = y[0] * y[0];
    };

  // Outer operator-splitting stepper; step_size in split_data sets the fixed
  // outer step via ARKodeSetFixedStep() inside SplittingStepper::reinit().
  SUNDIALS::SplittingStepper<VectorType> splitting({linear_stepper,
                                                    nonlinear_stepper},
                                                   split_data);

  // ARKode driver
  SUNDIALS::ARKode<VectorType>::AdditionalData arkode_data(
    t_start,
    t_end,
    dt, // initial_step_size (overridden by ARKodeSetFixedStep)
    output_period,
    dt * 1e-3, // minimum_step_size
    1e-8,      // absolute_tolerance
    1e-6);     // relative_tolerance

  SUNDIALS::ARKode<VectorType> ode(splitting, arkode_data);

  ode.output_step =
    [](const double t, const VectorType &y, const unsigned int /*step*/) {
      const double y_ex = exact_solution(t);
      deallog << std::fixed << std::setprecision(4) << "t=" << t
              << "  y=" << y[0] << "  exact=" << y_ex
              << "  err=" << std::abs(y[0] - y_ex) << std::endl;
    };

  VectorType y(1);
  y[0] = 1.0;
  ode.solve_ode(y);

  const double y_ex  = exact_solution(t_end);
  const double error = std::abs(y[0] - y_ex);
  deallog << "Final error: " << error << std::endl;
}

int
main()
{
  initlog();
  deallog.precision(5);

  // Sub-test 1: default (Lie-Trotter, 1st order)
  run(SUNDIALS::SplittingStepper<VectorType>::AdditionalData(dt),
      "Lie-Trotter");

  // Sub-test 2: Strang splitting (2nd order, two partitions)
  run(SUNDIALS::SplittingStepper<VectorType>::AdditionalData(
        dt, "ARKODE_SPLITTING_STRANG_2_2_2"),
      "Strang");
}

#else

int
main()
{
  return 0;
}

#endif
