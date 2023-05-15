// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Tests the class NonlinearSolverSelector using an example based on
// the test nonlinear_solver_selector_03. Here we use the nonlinear
// solver NOX instead of KINSOL with MPI.

#define SOLVER NonlinearSolverSelector<LA::MPI::Vector>::AdditionalData::nox

#include "nonlinear_solver_selector_03.cc"
