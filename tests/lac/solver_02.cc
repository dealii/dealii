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



// test lucky breakdown in GMRES (and others)

#include "../tests.h"
#include "testmatrix.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/base/point.h>

template<class SOLVER>
void test()
{
  const unsigned int size = 3;
  SparsityPattern sparsity(size, size, 1);
  sparsity.compress();
  SparseMatrix<double> mat;
  mat.reinit(sparsity);
  mat = IdentityMatrix(size);

  Vector<double> rhs;
  Vector<double> solvec;
  solvec.reinit(size);

  rhs.reinit(size);
  rhs(size-1)=1.0;

  SolverControl solvctrl(1000, 1e-12, true);
  SOLVER solver(solvctrl);

  PreconditionIdentity precond;
  solver.solve(mat, solvec, rhs, precond);
  solvec.print(deallog);
}

int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<SolverGMRES<Vector<double> > >();
  test<SolverCG<Vector<double> > >();
  test<SolverFGMRES<Vector<double> > >();
}

