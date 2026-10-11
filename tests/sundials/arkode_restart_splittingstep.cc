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

#include "arkode_restart_common.h"

// Test ARKode::solver_should_restart() with SplittingStepper. The ODE y' = -y
// is split into two parts y' = -y/2 that are integrated by adaptive
// ERKStepper sub-steppers. The outer step is fixed, and ARKode uses its initial
// step size as that step.

#if DEAL_II_SUNDIALS_VERSION_GTE(7, 2, 0)

int
main()
{
  initlog();

  const double step_size = 0.125;

  run_restart_test_with(
    "SplittingStepper",
    [&] {
      using Sub = SUNDIALS::ERKStepper<VectorType>;

      std::vector<std::shared_ptr<SUNDIALS::ARKodeStepper<VectorType>>> subs;
      for (unsigned int i = 0; i < 2; ++i)
        {
          auto sub               = std::make_shared<Sub>();
          sub->explicit_function = [](double,
                                      const VectorType &y,
                                      VectorType       &ydot) {
            ydot = y;
            ydot *= -0.5;
          };
          sub->custom_setup = [](void *arkode_mem) {
            const int status = ARKodeSStolerances(arkode_mem, 1e-8, 1e-10);
            AssertThrow(status == 0, ExcInternalError());
          };
          subs.push_back(sub);
        }

      return std::make_unique<SUNDIALS::SplittingStepper<VectorType>>(
        subs,
        SUNDIALS::SplittingStepper<VectorType>::AdditionalData(
          step_size, "ARKODE_SPLITTING_STRANG_2_2_2"));
    },
    step_size);
}

#else

int
main()
{
  return 0;
}

#endif
