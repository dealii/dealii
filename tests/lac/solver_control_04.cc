// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// test SolverControl's get_history_data function
// This test is adapted from tests/trilinos/solver_03.cc


#include "../tests.h"
#include <deal.II/base/utilities.h>
#include "../testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/precondition.h>

template<typename MatrixType, typename VectorType, class PRECONDITION>
void
check_solve (SolverControl       &solver_control,
             const MatrixType    &A,
             VectorType          &u,
             VectorType          &f,
             const PRECONDITION  &P,
             const bool           expected_result)
{
  SolverCG<VectorType> solver(solver_control);

  u = 0.;
  f = 1.;
  bool success = false;
  try
    {
      solver.solve(A,u,f,P);
      deallog << "Success. " << std::endl;
      success = true;
    }
  catch (std::exception &e)
    {
      deallog << "Failure. " << std::endl;
    }

  deallog << "Solver history data: " << solver_control.get_history_data()
          << std::endl;
  Assert(success == expected_result, ExcMessage("Incorrect result."));
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  {
    const unsigned int size = 32;
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
    Vector<double>  u(dim);
    f = 1.;

    PreconditionJacobi<> preconditioner;
    preconditioner.initialize(A);

    deallog.push("Abs tol");
    {
      // Expects success
      SolverControl solver_control(2000, 1.e-3);
      solver_control.enable_history_data();

      check_solve (solver_control, A,u,f, preconditioner, true);
    }
    deallog.pop();
    deallog.push("Iterations");
    {
      // Expects failure
      SolverControl solver_control(20, 1.e-3);
      solver_control.enable_history_data();

      check_solve (solver_control, A,u,f, preconditioner, false);
    }
    deallog.pop();
  }
}
