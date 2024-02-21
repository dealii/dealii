// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests the class NonlinearSolverSelector using an example based on
// the test nonlinear_solver_selector_01. Here we use the nonlinear
// solver NOX instead of KINSOL.

#define SOLVER NonlinearSolverSelector<Vector<double>>::AdditionalData::nox

#include "nonlinear_solver_selector_01.cc"
