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

// Like arkode_lsrkstep_sts, but demonstrates the *default* dominant-eigenvalue
// estimator for LSRKStepperSTS.
//
// LSRK-STS methods (RKC-2, RKL-2) choose the number of polynomial stages
// adaptively from the spectral radius of the Jacobian, so they need an
// estimate of the dominant (largest-magnitude) eigenvalue. In
// arkode_lsrkstep_sts this estimate is supplied analytically through
// LSRKStepperSTS::dominant_eigenvalue_function.
//
// Here that callback is deliberately left unset. In that case deal.II creates
// a SUNDIALS built-in power-iteration dominant-eigenvalue estimator
// (SUNDomEigEstimator) and attaches it automatically. The estimator derives
// the dominant eigenvalue internally from explicit_function() without any
// further user input. This automatic default is only available with SUNDIALS
// 7.5.0 or newer, which is why this test is guarded accordingly.
//
// Problem (identical to arkode_lsrkstep_sts):
//   u_t = u_xx,  x in (0,1), t in (0, 0.1]
//   u(0,t) = u(1,t) = 0
//   u(x,0) = sin(pi*x)
//
// Semidiscretized with N = 10 interior points using central differences
// (h = 1/11):  y_i' = (y_{i-1} - 2*y_i + y_{i+1}) / h^2
//
// Largest eigenvalue magnitude: ~4/h^2 = 484.
//
// Because the initial data sin(pi*x_i) is exactly the first eigenvector of the
// discrete Laplacian, the semidiscretized ODE system also has a closed-form
// solution y_i(t) = sin(pi*x_i) * exp(lambda*t) with the discrete eigenvalue
// lambda = -(4/h^2) * sin^2(pi*h/2). This semidiscrete solution is used below
// to compute the time-integration error. The fact that the error stays small
// demonstrates that the built-in estimator returns a sufficiently accurate
// dominant eigenvalue for SUNDIALS to pick a stable number of stages.
//
// Two sub-tests in arkode_lsrkstep_sts demonstrate that the built-in estimator
// reproduces the analytical-eigenvalue results. Here we additionally observe
// the effect of tuning the estimator:
//
// The number of *time steps* taken is governed by ARKODE's error-based time
// step adaptivity, not by the eigenvalue estimate, which only controls how
// many polynomial stages each step needs in order to be stable. As long as the
// estimate is accurate enough to keep the method stable, the number of time
// steps and the achieved accuracy are therefore identical regardless of how
// the estimator is tuned. What does change with the tuning is the *amount of
// work*: every power iteration (including the one-time warmup iterations)
// costs an evaluation of the right-hand side. Tuning the estimator more
// aggressively (more warmup iterations) increases the right-hand-side
// evaluation count while leaving the time step count and the solution accuracy
// intact.

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 5, 0)

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

static unsigned int
run(const SUNDIALS::LSRKStepperSTS<VectorType>::AdditionalData &data)
{
  // Number of explicit_function() (right-hand-side) evaluations performed
  // during a single run(). Reset at the start of each run().
  unsigned int n_rhs_evals = 0;

  SUNDIALS::LSRKStepperSTS<VectorType> stepper(data);

  stepper.explicit_function =
    [&n_rhs_evals](double, const VectorType &y, VectorType &ydot) {
      ++n_rhs_evals;
      const double inv_h2 = 1.0 / (h * h);
      ydot[0]             = inv_h2 * (-2.0 * y[0] + y[1]);
      for (unsigned int i = 1; i < N - 1; ++i)
        ydot[i] = inv_h2 * (y[i - 1] - 2.0 * y[i] + y[i + 1]);
      ydot[N - 1] = inv_h2 * (y[N - 2] - 2.0 * y[N - 1]);
    };

  // Note: dominant_eigenvalue_function is intentionally NOT set. deal.II
  // therefore creates and attaches the SUNDIALS built-in power-iteration
  // dominant-eigenvalue estimator, which derives the dominant eigenvalue
  // internally from explicit_function() above.

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

  return n_rhs_evals;
}

int
main()
{
  initlog();

  unsigned int n_rhs_evals_default = 0;
  unsigned int n_rhs_evals_warmup  = 0;

  // Sub-test 1: default STS method (RKC-2) with the built-in dominant-
  // eigenvalue estimator at its default settings (no eigenvalue callback and
  // no estimator tuning supplied). This is the reference run.
  {
    deallog << "=== Default STS, default estimator ===" << std::endl;
    SUNDIALS::LSRKStepperSTS<VectorType>::AdditionalData data;
    n_rhs_evals_default = run(data);
  }

  // Sub-test 2: the same RKC-2 method and the same built-in estimator, but
  // with a large number of one-time warmup power iterations. These warmups
  // refine the initial eigenvalue estimate at the cost of extra right-hand-
  // side evaluations performed once at startup. Compared with sub-test 1, the
  // solution accuracy and the number of time steps are unchanged, while the
  // number of right-hand-side evaluations is larger -- this is the visible
  // effect of tuning the estimator.
  {
    deallog << "=== Default STS, estimator with extra warmups ===" << std::endl;
    SUNDIALS::LSRKStepperSTS<VectorType>::AdditionalData data(
      /*method_name=*/"",
      /*dom_eig_frequency=*/-1,
      /*max_num_stages=*/0,
      /*dom_eig_estimator_max_iters=*/-1,
      /*dom_eig_estimator_rel_tol=*/-1.,
      /*dom_eig_estimator_num_warmups=*/200);
    n_rhs_evals_warmup = run(data);
  }

  deallog << std::boolalpha
          << "Estimator with extra warmups results in more RHS evaluations: "
          << (n_rhs_evals_warmup > n_rhs_evals_default) << std::endl;
  deallog
    << std::boolalpha
    << "Additional number of RHS evaluations less than the number of warmups: "
    << (n_rhs_evals_warmup - n_rhs_evals_default < 200) << std::endl;
}

#else

int
main()
{
  return 0;
}

#endif
