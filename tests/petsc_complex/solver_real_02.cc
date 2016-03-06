// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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


// test the PETSc Chebychev solver

// Note: This is (almost) a clone of the tests/petsc/solver_02.cc

#include "../tests.h"
#include "../lac/testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <typeinfo>

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_solve( SOLVER &solver, const MATRIX &A,
             VECTOR &u, VECTOR &f, const PRECONDITION &P)
{
  deallog << "Solver type: " << typeid(solver).name() << std::endl;

  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A,u,f,P);
    }
  catch (std::exception &e)
    {
      deallog << e.what() << std::endl;
      // just like for Richardson: we end up
      // here normally for this solver
    }

  deallog << "Solver stopped after " << solver.control().last_step()
          << " iterations" << std::endl;
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  {
    SolverControl control(100, 1.e-3);

    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix testproblem(size, size);
    PETScWrappers::SparseMatrix  A(dim, dim, 5);
    testproblem.five_point(A);

    PETScWrappers::Vector  f(dim);
    PETScWrappers::Vector  u(dim);
    f = 1.;
    A.compress (VectorOperation::insert);
    f.compress (VectorOperation::insert);
    u.compress (VectorOperation::insert);

    PETScWrappers::SolverChebychev solver(control);
    PETScWrappers::PreconditionJacobi preconditioner(A);
    check_solve (solver, A,u,f, preconditioner);
  }

}

