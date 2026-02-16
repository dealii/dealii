// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Tests the class NonlinearSolverSelector using an example based on
// the test nonlinear_solver_selector_03. Here we use the nonlinear
// solver NOX instead of KINSOL with MPI.

#define SOLVER NonlinearSolverSelector<LA::MPI::Vector>::AdditionalData::nox

#include "nonlinear_solver_selector_03.cc"
