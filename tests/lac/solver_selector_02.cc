// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
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


// Test a bug in the SolverSelector when using a custom SolverControl. At one
// point the SolverControl got "sliced".

#include "../tests.h"
#include "testmatrix.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_selector.h>

#include <fstream>

class MySolverControl:
  public SolverControl
{
public:
  virtual State check(const unsigned int step,const double)
  {
    deallog << "MySolverControl " << step << std::endl;
    return SolverControl::check(step, 0);
  }
};




template <class MATRIX, class VECTOR>
void
check(const MATRIX &A, const VECTOR &f)
{
  std::vector<std::string> names;
  names.push_back("cg");
  names.push_back("bicgstab");
  names.push_back("gmres");
  names.push_back("fgmres");

  MySolverControl mycont;
  SolverSelector<VECTOR> solver;
  PreconditionSSOR<SparseMatrix<double> > pre;
  pre.initialize(A);

  VECTOR u;
  u.reinit(f);

  std::vector<std::string>::const_iterator name;

  solver.set_control(mycont);
  for (name = names.begin(); name != names.end(); ++name)
    {
      solver.select(*name);
      u = 0.;
      solver.solve(A, u, f, pre);
    }
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  unsigned int size=37;
  unsigned int dim = (size-1)*(size-1);

  deallog << "Size " << size << " Unknowns " << dim << std::endl;

  // Make matrix
  FDMatrix testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double>  A(structure);
  testproblem.five_point(A);
  Vector<double>  f(dim);
  f = 1.;

  check(A, f);
}
