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
#include <iomanip>

#include "../tests.h"

// Test LSRKStepperSSP (Strong-Stability-Preserving) on the linear advection
// equation semidiscretized with first-order upwind differences.
//
// SSP methods preserve the strong-stability (monotonicity / TVD) of the
// spatial discretization in time, making them the natural choice for
// advection-dominated or hyperbolic problems. Unlike STS methods, SSP methods
// use a fixed number of stages and do not require a dominant eigenvalue
// estimate.
//
// Problem (linear advection, periodic boundary on [0, 1)):
//   u_t + a*u_x = 0,  a = 1,  x in [0, 1) periodic,  t in [0, 0.3]
//   u(x, 0) = sin(2*pi*x)
//
// Semidiscretized at the N nodes x_j = j*dx, dx = 1/N, with first-order
// upwind differencing (a > 0):
//   y_j' = -a * (y_j - y_{j-1}) / dx          (indices taken modulo N)
//
// This is the canonical hyperbolic test problem for SSP time integrators: the
// upwind operator is strong-stability-preserving (total-variation diminishing)
// and the SSP Runge-Kutta method maintains that property in time.
//
// Because the initial datum sin(2*pi*x) is a single Fourier mode, it is an
// eigenvector of the circulant upwind operator, whose eigenvalue for this mode
// is  lambda = -(a/dx) * (1 - exp(-i*k*dx)),  k = 2*pi. The semidiscretized ODE
// system therefore has the closed-form solution
//   y_j(t) = exp(alpha*t) * sin(k*x_j + beta*t)
// with
//   alpha = -(a/dx) * (1 - cos(k*dx))   (numerical diffusion, amplitude decay)
//   beta  = -(a/dx) * sin(k*dx)         (numerical propagation, phase shift).
// This exact semidiscrete solution is used below to compute the
// time-integration error of the SSP method.
//
// Two sub-tests:
//   1. Default SSP method (SUNDIALS default: SSP(10,4)).
//   2. SSP_S_3 with 9 stages (num_stages = 9), which is a 3rd-order method
//      with 9 stages (perfect-square stage count required by SSP_S_3).

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

using VectorType = Vector<double>;

static constexpr unsigned int N  = 20;
static constexpr double       dx = 1.0 / N;
static constexpr double       a  = 1.0; // advection speed (a > 0)

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

static void
run(const SUNDIALS::LSRKStepperSSP<VectorType>::AdditionalData &data)
{
  SUNDIALS::LSRKStepperSSP<VectorType> stepper(data);

  // First-order upwind discretization of -a*u_x for a > 0.
  stepper.explicit_function =
    [](double, const VectorType &y, VectorType &ydot) {
      const double inv_dx = 1.0 / dx;
      for (unsigned int i = 0; i < N; ++i)
        {
          const double yi   = y[i];
          const double yim1 = y[(i + N - 1) % N];
          ydot[i]           = -a * inv_dx * (yi - yim1);
        }
    };

  SUNDIALS::ARKode<VectorType> ode(stepper, make_arkode_data());

  // Output t, the value at x = 1/4 (the initial peak of sin(2*pi*x)), the
  // analytical value there, and the discrete L2 error norm of the whole
  // solution. The exact semidiscrete solution is
  //   y_j(t) = exp(alpha*t) * sin(k*x_j + beta*t),  k = 2*pi,
  // so the error norm measures the time-integration error of the SSP method.
  ode.output_step =
    [](const double t, const VectorType &sol, const unsigned int /*step*/) {
      const double k     = 2.0 * M_PI;
      const double alpha = -a / dx * (1.0 - std::cos(k * dx));
      const double beta  = -a / dx * std::sin(k * dx);
      VectorType   exact(N);
      for (unsigned int i = 0; i < N; ++i)
        exact[i] = std::exp(alpha * t) * std::sin(k * i * dx + beta * t);
      VectorType diff(exact);
      diff -= sol;
      deallog << std::fixed << std::setprecision(4) << t << ' ' << sol[N / 4]
              << ' ' << exact[N / 4] << ' ' << diff.l2_norm() << std::endl;
    };

  VectorType y(N);
  for (unsigned int i = 0; i < N; ++i)
    y[i] = std::sin(2.0 * M_PI * i * dx);

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
