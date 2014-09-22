// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check that IterationNumber control works as expected, i.e., it terminates
// with success when reaching the maximum number of iterations


#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>


void check()
{
  // Check 1: should still terminate after one iteration when using a matrix
  // with only one eigenvalue
  const unsigned int size = 20;
  FullMatrix<double> A(size, size);
  for (unsigned int i=0; i<size; ++i)
    A(i,i) = 2;

  Vector<double> in(size), out(size);
  in = 1.;
  
  {
    IterationNumberControl control(5);
    SolverCG<> solver(control);
    solver.solve(A, out, in, PreconditionIdentity());
    Assert (control.last_step() == 1, ExcInternalError());
  }
  for (unsigned int i=0; i<size; ++i)
    Assert(std::abs(out(i)-0.5) < 1e-12, ExcInternalError());

  // Check 2: should only do 5 iterations but the solution should not be exact
  for (unsigned int i=0; i<size; ++i)
    A(i,i) = 1+i;

  out = 0;
  {
    IterationNumberControl control(5);
    SolverCG<> solver(control);
    solver.solve(A, out, in, PreconditionIdentity());
    Assert (control.last_step() == 5, ExcInternalError());
  }
  bool solved_exactly = true;
  for (unsigned int i=0; i<size; ++i)
    if (std::abs(out(i)-1./(1+size)) > 1e-8)
      solved_exactly = false;
  Assert(solved_exactly == false, ExcInternalError());

  deallog << "OK" << std::endl;
}

int main()
{
  std::ofstream logfile("output");
//  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check();
}

