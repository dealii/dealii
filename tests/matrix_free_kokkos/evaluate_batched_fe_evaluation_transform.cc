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


// Test Portable::internal::batched::FEEvaluationImplTransformToCollocation:
// evaluate() and integrate() of values and gradients for a batch of cells
// against the plain tensor product applied cell by cell, see
// evaluate_batched_fe_evaluation_common.h. The (degree, n_q_points_1d)
// pairs are the ones Portable::FEEvaluation dispatches to this kernel.

#include "evaluate_batched_fe_evaluation_common.h"


int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  Kokkos::initialize();
  {
    constexpr BatchedKernel kernel = BatchedKernel::transform_to_collocation;

    run<1, 2, 3, kernel>();
    run<1, 2, 4, kernel>();
    run<1, 3, 5, kernel>();

    run<2, 2, 3, kernel>();
    run<2, 2, 4, kernel>();
    run<2, 3, 5, kernel>();

    run<3, 2, 3, kernel>();
    run<3, 2, 4, kernel>();
    run<3, 3, 5, kernel>();
  }
  Kokkos::finalize();

  return 0;
}
