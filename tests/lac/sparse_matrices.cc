//----------------------------  testmatrix.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  testmatrix.cc  ---------------------------


//TODO: [GK] Produce some useful output!

#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/sparse_matrix_ez.h>
#include <lac/vector.h>
#include <lac/solver_richardson.h>
#include <lac/precondition.h>

#include <lac/sparse_matrix_ez.templates.h>

#include <fstream>

template<class MATRIX>
void
check_vmult_quadratic(const MATRIX& A,
		      const char* prefix)
{
  deallog.push(prefix);
  
  Vector<double> u(A.n());
  Vector<double> f(A.m());
  GrowingVectorMemory<> mem;

  SolverControl control(20, 1.e-3, false);
  SolverRichardson<> rich(control, mem, .01);
  PreconditionIdentity prec;

  u = 0.;
  f = 1.;

  try
    {
      rich.solve(A, u, f, prec);
    }
  catch (...)
    {
    }
  deallog << "Transpose" << std::endl;
  try
    {
      rich.Tsolve(A, u, f, prec);
    }
  catch (...)
    {
    }
  deallog.pop();
}

int main()
{
  std::ofstream logfile("sparse_matrices.output");
  logfile.setf(std::ios::fixed);
  logfile.precision(2);
  deallog.attach(logfile);

				   // Switch between regression test
				   // and benchmark
#ifdef DEBUG  
  deallog.depth_console(0);
  const unsigned int size = 10;
#else
  deallog.depth_console(1000);
  deallog.log_execution_time(true);
  deallog.log_time_differences(true);
  const unsigned int size = 500;
#endif
  
  FDMatrix testproblem (size, size);
  unsigned int dim = (size-1)*(size-1);

  deallog << "Structure" << std::endl;
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double>  A(structure);
  deallog << "Assemble" << std::endl;
  testproblem.five_point(A);
  check_vmult_quadratic(A, "5-SparseMatrix<double>");

  SparseMatrixEZ<double> E(dim,dim,5,1);
  deallog << "Assemble" << std::endl;
  testproblem.five_point(E);
  check_vmult_quadratic(E, "5-SparseMatrixEZ<double>");

  A.clear();
  deallog << "Structure" << std::endl;
  structure.reinit(dim, dim, 9);
  testproblem.nine_point_structure(structure);
  structure.compress();
  A.reinit(structure);
  deallog << "Assemble" << std::endl;
  testproblem.nine_point(A);
  check_vmult_quadratic(A, "9-SparseMatrix<double>");

  E.clear();
  E.reinit(dim,dim,9,2);
  deallog << "Assemble" << std::endl;
  testproblem.nine_point(E);
  check_vmult_quadratic(E, "9-SparseMatrixEZ<double>");
  
}
