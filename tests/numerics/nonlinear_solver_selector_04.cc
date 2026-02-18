// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Tests that nonlinear.h can be included and that an object of
// NonlinearSolverSelector can be created, even when none of SUNDIALS,
// TRILINOS, or PETSC are enabled. Test added to trigger the bug, that this was
// not possibles.

#include <deal.II/base/function.h>

#include <deal.II/numerics/nonlinear.h>

#include "../tests.h"


int
main()
{
  initlog();

  NonlinearSolverSelector solver;

  deallog << "OK" << std::endl;
}
