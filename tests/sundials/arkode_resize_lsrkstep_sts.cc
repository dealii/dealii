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

#include <complex>

#include "arkode_resize_common.h"

// Test ARKode::solver_should_resize() with LSRKStepperSTS. With SUNDIALS 7.5
// and later this also covers the rebuild of the built-in dominant eigenvalue
// estimator.

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

int
main()
{
  initlog();

  using Stepper = SUNDIALS::LSRKStepperSTS<VectorType>;

  run_resize_test<Stepper>("LSRKStepperSTS", [](Stepper &stepper) {
    stepper.explicit_function = rhs;
#  if !DEAL_II_SUNDIALS_VERSION_GTE(7, 5, 0)
    stepper.dominant_eigenvalue_function =
      [](double, const VectorType &, const VectorType &) {
        return std::complex<double>(-1., 0.);
      };
#  endif
  });
}

#else

int
main()
{
  return 0;
}

#endif
