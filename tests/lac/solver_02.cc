//----------------------------------------------------------------------
//    $Id: solver.cc 23710 2011-05-17 04:50:10Z bangerth $
//
//    Copyright (C) 1998 - 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


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
  std::ofstream logfile("solver_02/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<SolverGMRES<Vector<double> > >();
  test<SolverCG<Vector<double> > >();
//  test<SolverFGMRES<Vector<double> > >();
    
}

