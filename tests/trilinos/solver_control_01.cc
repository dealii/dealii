// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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


// test SolverControl with Trilinos solver
// This test is adapted from tests/trilinos/solver_03.cc


#include "../tests.h"
#include <deal.II/base/utilities.h>
#include "../testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <typeinfo>

template<typename MatrixType, typename VectorType, class PRECONDITION>
void
check_solve (SolverControl       &solver_control,
             const MatrixType    &A,
             VectorType          &u,
             VectorType          &f,
             const PRECONDITION  &P,
             const bool           expected_result)
{
  TrilinosWrappers::SolverCG solver(solver_control);

  u = 0.;
  f = 1.;
  bool success = false;
  try
    {
      solver.solve(A,u,f,P);
      deallog << "Success. ";
      success = true;
    }
  catch (std::exception &e)
    {
      deallog << "Failure. ";
    }

  deallog << "Solver stopped after " << solver_control.last_step()
          << " iterations" << std::endl;
  Assert(success == expected_result, ExcMessage("Incorrect result."));
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());


  {
    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix testproblem(size, size);
    DynamicSparsityPattern csp (dim, dim);
    testproblem.five_point_structure(csp);
    TrilinosWrappers::SparseMatrix  A;
    A.reinit(csp);
    testproblem.five_point(A);

    TrilinosWrappers::Vector  f(dim);
    TrilinosWrappers::Vector  u(dim);
    f = 1.;
    A.compress (VectorOperation::insert);
    f.compress (VectorOperation::insert);
    u.compress (VectorOperation::insert);

    TrilinosWrappers::PreconditionJacobi preconditioner;
    preconditioner.initialize(A);

    deallog.push("Abs tol");
    {
      // Expects success
      SolverControl solver_control(2000, 1.e-3);
      check_solve (solver_control, A,u,f, preconditioner, true);
    }
    deallog.pop();
    deallog.push("Iterations");
    {
      // Expects failure
      SolverControl solver_control(20, 1.e-3);
      check_solve (solver_control, A,u,f, preconditioner, false);
    }
    deallog.pop();
  }
}
