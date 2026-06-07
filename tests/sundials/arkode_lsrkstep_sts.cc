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
#include <complex>
#include <iomanip>

#include "../tests.h"

// Test LSRKStepperSTS (Stabilized explicit Runge-Kutta) on the
// semidiscretized 1D heat equation.
//
// LSRK-STS methods (RKC-2, RKL-2) are designed for problems whose Jacobian
// eigenvalues lie on or near the negative real axis, such as spatial
// discretizations of parabolic PDEs. The stability region grows as O(s^2)
// with the number of polynomial stages s, allowing much larger steps than
// standard explicit Runge-Kutta methods for diffusion-dominated problems.
//
// Problem:
//   u_t = u_xx,  x in (0,1), t in (0, 0.1]
//   u(0,t) = u(1,t) = 0
//   u(x,0) = sin(pi*x)
//
// Semidiscretized with N = 10 interior points using central differences
// (h = 1/11):  y_i' = (y_{i-1} - 2*y_i + y_{i+1}) / h^2
//
// Largest eigenvalue magnitude: ~4/h^2 = 484.
// Exact PDE solution:  u(x,t) = sin(pi*x) * exp(-pi^2 * t)
//
// Because the initial data sin(pi*x_i) is exactly the first eigenvector of the
// discrete Laplacian, the semidiscretized ODE system also has a closed-form
// solution y_i(t) = sin(pi*x_i) * exp(lambda*t) with the discrete eigenvalue
// lambda = -(4/h^2) * sin^2(pi*h/2). This semidiscrete solution is used below
// to compute the time-integration error.
//
// The dominant eigenvalue is supplied exactly via
// LSRKStepperSTS::dominant_eigenvalue_function so that SUNDIALS can
// determine the required number of polynomial stages at each step.
//
// Two sub-tests:
//   1. Default STS method (SUNDIALS default: RKC-2), dominant eigenvalue
//      recomputed at SUNDIALS default frequency (dom_eig_frequency = -1,
//      which leaves the SUNDIALS default of 25 steps in place).
//   2. Runge-Kutta-Legendre order-2 ("ARKODE_LSRK_RKL_2"), eigenvalue
//      recomputed every 5 steps (dom_eig_frequency = 5).

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

using VectorType = Vector<double>;

static constexpr unsigned int N = 10;
static constexpr double       h = 1.0 / (N + 1);

static SUNDIALS::ARKode<VectorType>::AdditionalData
make_arkode_data()
{
  return {0.0 /*initial_time*/,
          0.1 /*final_time*/,
          1e-3 /*initial_step_size*/,
          0.01 /*output_period*/,
          1e-8 /*minimum_step_size*/,
          1e-8 /*absolute_tolerance*/,
          1e-6 /*relative_tolerance*/};
}

static void
run(const SUNDIALS::LSRKStepperSTS<VectorType>::AdditionalData &data)
{
  SUNDIALS::LSRKStepperSTS<VectorType> stepper(data);

  stepper.explicit_function =
    [](double, const VectorType &y, VectorType &ydot) {
      const double inv_h2 = 1.0 / (h * h);
      ydot[0]             = inv_h2 * (-2.0 * y[0] + y[1]);
      for (unsigned int i = 1; i < N - 1; ++i)
        ydot[i] = inv_h2 * (y[i - 1] - 2.0 * y[i] + y[i + 1]);
      ydot[N - 1] = inv_h2 * (y[N - 2] - 2.0 * y[N - 1]);
    };

  // Provide the dominant eigenvalue of the discrete Laplacian analytically.
  // The largest-magnitude eigenvalue is -4/h^2 (purely real).
  stepper.dominant_eigenvalue_function =
    [](double /*t*/, const VectorType & /*y*/, const VectorType & /*f*/) {
      return std::complex<double>(-4.0 / (h * h), 0.0);
    };

  SUNDIALS::ARKode<VectorType> ode(stepper, make_arkode_data());

  // Output t, the midpoint value y[N/2] (i = 5, x = 6/11), the analytical
  // value there, and the discrete L2 error norm.
  //
  // The initial condition sin(pi*x_i) is exactly the first eigenvector of the
  // discrete Laplacian, so the semidiscrete ODE system y' = A*y has the exact
  // solution y_i(t) = sin(pi*x_i) * exp(lambda*t), where the eigenvalue is
  // lambda = -(4/h^2) * sin^2(pi*h/2). The error norm therefore measures the
  // time-integration error of the LSRK-STS method alone.
  ode.output_step =
    [](const double t, const VectorType &sol, const unsigned int /*step*/) {
      const double s      = std::sin(M_PI * h / 2.0);
      const double lambda = -4.0 / (h * h) * s * s;
      VectorType   exact(N);
      for (unsigned int i = 0; i < N; ++i)
        exact[i] = std::sin(M_PI * (i + 1) * h) * std::exp(lambda * t);
      VectorType diff(exact);
      diff -= sol;
      deallog << std::fixed << std::setprecision(4) << t << ' ' << sol[N / 2]
              << ' ' << exact[N / 2] << ' ' << diff.l2_norm() << std::endl;
    };

  VectorType y(N);
  for (unsigned int i = 0; i < N; ++i)
    y[i] = std::sin(M_PI * (i + 1) * h);

  ode.solve_ode(y);
}

int
main()
{
  initlog();

  // Sub-test 1: default STS method (RKC-2), eigenvalue recomputed every step.
  {
    deallog << "=== Default STS ===" << std::endl;
    SUNDIALS::LSRKStepperSTS<VectorType>::AdditionalData data;
    run(data);
  }

  // Sub-test 2: RKL-2 method, eigenvalue recomputed every 5 steps.
  // RKL-2 has a larger stability region than RKC-2 for the same stage count.
  {
    deallog << "=== RKL-2, dom_eig_frequency=5 ===" << std::endl;
    SUNDIALS::LSRKStepperSTS<VectorType>::AdditionalData data(
      "ARKODE_LSRK_RKL_2",
      /*dom_eig_frequency=*/5);
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
