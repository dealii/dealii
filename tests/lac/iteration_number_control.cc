// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check that IterationNumber control works as expected, i.e., it terminates
// with success when reaching the maximum number of iterations


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


void
check()
{
  // Check 1: should still terminate after one iteration when using a matrix
  // with only one eigenvalue
  const unsigned int size = 20;
  FullMatrix<double> A(size, size);
  for (unsigned int i = 0; i < size; ++i)
    A(i, i) = 2;

  Vector<double> in(size), out(size);
  in = 1.;

  {
    IterationNumberControl control(5, 1e-12, false, true);
    SolverCG<>             solver(control);
    solver.solve(A, out, in, PreconditionIdentity());
    AssertThrow(control.last_step() == 1, ExcInternalError());
  }
  for (unsigned int i = 0; i < size; ++i)
    AssertThrow(std::abs(out(i) - 0.5) < 1e-12, ExcInternalError());

  // Check 2: should only do 5 iterations but the solution should not be exact
  for (unsigned int i = 0; i < size; ++i)
    A(i, i) = 1 + i;

  out = 0;
  {
    IterationNumberControl control(5, 1e-12, false, true);
    SolverCG<>             solver(control);
    solver.solve(A, out, in, PreconditionIdentity());
    AssertThrow(control.last_step() == 5, ExcInternalError());
  }
  bool solved_exactly = true;
  for (unsigned int i = 0; i < size; ++i)
    if (std::abs(out(i) - 1. / (1 + size)) > 1e-8)
      solved_exactly = false;
  AssertThrow(solved_exactly == false, ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  //  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);

  check();
}
