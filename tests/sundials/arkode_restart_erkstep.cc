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

// Test ARKode::solver_should_restart() with ERKStepper.

int
main()
{
  initlog();

  run_restart_test<SUNDIALS::ERKStepper<VectorType>>("ERKStepper");
}
