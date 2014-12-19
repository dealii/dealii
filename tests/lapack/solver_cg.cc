// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// Check eigenvalue capabilities of SolverCG

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../lac/testmatrix.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.templates.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_solve( SOLVER &solver, const MATRIX &A,
             VECTOR &u, VECTOR &f, const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A,u,f,P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_Tsolve(SOLVER &solver, const MATRIX &A,
             VECTOR &u, VECTOR &f, const PRECONDITION &P)
{
  u = 0.;
  f = 1.;
  try
    {
      solver.Tsolve(A,u,f,P);
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }
}

int main()
{
  std::ofstream logfile("output");
//  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  GrowingVectorMemory<> mem;
  SolverControl control(100, 1.e-3);
  SolverControl verbose_control(100, 1.e-3, true);
  SolverCG<>::AdditionalData data(false, true, true, true);
  SolverCG<> cg(control, mem, data);

  for (unsigned int size=4; size <= 30; size *= 3)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.five_point(A);

      PreconditionIdentity prec_no;
      PreconditionRichardson prec_richardson;
      prec_richardson.initialize(0.6);
      PreconditionSSOR<> prec_ssor;
      prec_ssor.initialize(A, 1.2);

      std::vector<unsigned int> permutation(dim);
      std::vector<unsigned int> inverse_permutation(dim);

      Vector<double>  f(dim);
      Vector<double>  u(dim);
      Vector<double> res(dim);

      f = 1.;
      u = 1.;

      try
        {
          deallog.push("no-fail");
          control.set_max_steps(10);
          check_solve(cg,A,u,f,prec_no);
          control.set_max_steps(100);
          deallog.pop();

          deallog.push("no");
          check_solve(cg,A,u,f,prec_no);
          deallog.pop();

          deallog.push("rich");
          check_solve(cg,A,u,f,prec_richardson);
          deallog.pop();

          deallog.push("ssor");
          check_solve(cg,A,u,f,prec_ssor);
          deallog.pop();
        }
      catch (std::exception &e)
        {
          std::cerr << "Exception: " << e.what() << std::endl;
        }
    }
}

