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

#include "../tests.h"

// Test LSRKStepperSSP (Strong-Stability-Preserving) on the nonlinear
// Burgers equation semidiscretized with upwind differences.
//
// SSP methods preserve the strong-stability (monotonicity) of the spatial
// discretization in time, making them suitable for advection-dominated or
// hyperbolic problems. Unlike STS methods, SSP methods use a fixed number
// of stages and do not require a dominant eigenvalue estimate.
//
// Problem (inviscid Burgers, periodic boundary, one full period [0, 1]):
//   u_t + (u^2/2)_x = 0,  x in [0, 1) periodic,  t in [0, 0.4]
//   u(x, 0) = 0.5 + sin(2*pi*x)
//
// Semidiscretized with N cells, upwind flux (Godunov):
//   y_i' = -(F_{i+1/2} - F_{i-1/2}) / dx
//   F_{i+1/2} = max(y_i, 0)^2/2 + min(y_{i+1}, 0)^2/2  (upwind)
//
// For smooth initial data the exact solution before shock formation is
// u(x,t) = 0.5 + sin(2*pi*(x - (0.5 + sin(2*pi*x))*t)) (implicitly defined).
//
// Two sub-tests:
//   1. Default SSP method (SUNDIALS default: SSP(10,4)).
//   2. SSP_S_3 with 9 stages (num_stages = 9), which is a 3rd-order method
//      with 9 stages (perfect-square stage count required by SSP_S_3).

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

using VectorType = Vector<double>;

static constexpr unsigned int N  = 20;
static constexpr double       dx = 1.0 / N;

static SUNDIALS::ARKode<VectorType>::AdditionalData
make_arkode_data()
{
  return {0.0 /*initial_time*/,
          0.3 /*final_time*/,
          1e-3 /*initial_step_size*/,
          0.05 /*output_period*/,
          1e-8 /*minimum_step_size*/,
          1e-8 /*absolute_tolerance*/,
          1e-6 /*relative_tolerance*/};
}

// Godunov (upwind) flux for Burgers: F(u) = u^2/2
static double
godunov_flux(const double ul, const double ur)
{
  // Roe flux: f(ul, ur) = max(ul,0)^2/2 + min(ur,0)^2/2
  const double fl = (ul > 0.0 ? 0.5 * ul * ul : 0.0);
  const double fr = (ur < 0.0 ? 0.5 * ur * ur : 0.0);
  return fl + fr;
}

static void
run(const SUNDIALS::LSRKStepperSSP<VectorType>::AdditionalData &data)
{
  SUNDIALS::LSRKStepperSSP<VectorType> stepper(data);

  stepper.explicit_function =
    [](double, const VectorType &y, VectorType &ydot) {
      const double inv_dx = 1.0 / dx;
      for (unsigned int i = 0; i < N; ++i)
        {
          const double ul    = y[i];
          const double ur    = y[(i + 1) % N];
          const double ul_m1 = y[(i + N - 1) % N];
          ydot[i] = -inv_dx * (godunov_flux(ul, ur) - godunov_flux(ul_m1, ul));
        }
    };

  SUNDIALS::ARKode<VectorType> ode(stepper, make_arkode_data());

  ode.output_step =
    [](const double t, const VectorType &sol, const unsigned int /*step*/) {
      // Print minimum and maximum to a concise summary line.
      double vmin = sol[0], vmax = sol[0];
      for (unsigned int i = 1; i < N; ++i)
        {
          vmin = std::min(vmin, sol[i]);
          vmax = std::max(vmax, sol[i]);
        }
      deallog << t << ' ' << vmin << ' ' << vmax << std::endl;
    };

  VectorType y(N);
  for (unsigned int i = 0; i < N; ++i)
    y[i] = 0.5 + std::sin(2.0 * M_PI * (i + 0.5) * dx);

  ode.solve_ode(y);
}

int
main()
{
  initlog();

  // Sub-test 1: default SSP method (SSP(10,4)).
  {
    deallog << "=== Default SSP ===" << std::endl;
    SUNDIALS::LSRKStepperSSP<VectorType>::AdditionalData data;
    run(data);
  }

  // Sub-test 2: SSP_S_3 with 9 stages.
  // The variable-stage SSP-S-3 method requires the number of stages to be a
  // perfect square >= 4. With s=9 stages it achieves 3rd-order accuracy.
  {
    deallog << "=== SSP_S_3 with 9 stages ===" << std::endl;
    SUNDIALS::LSRKStepperSSP<VectorType>::AdditionalData data(
      "ARKODE_LSRK_SSP_S_3",
      /*num_stages=*/9);
    run(data);
  }
}

#else

int
main()
{
  return 0;
}

#endif
