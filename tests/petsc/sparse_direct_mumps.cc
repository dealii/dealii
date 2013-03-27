//----------------------------  sparse_direct_mumps.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_direct_mumps.cc  ---------------------------

// test the PETSc SparseDirectMumps solver


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
check_solve( SOLVER& solver, const MATRIX& A,
	     VECTOR& u, VECTOR& f, const PRECONDITION& P)
{
  deallog << "Solver type: " << typeid(solver).name() << std::endl;

  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A,u,f,P);
    }
  catch (std::exception& e)
    {
      deallog << e.what() << std::endl;
      abort ();
    }

  deallog << "Solver stopped after " << solver.control().last_step()
          << " iterations" << std::endl;
}


int main(int argc, char **argv)
{
  std::ofstream logfile("sparse_direct_mumps/output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  {
    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

                                     // Make matrix
    FDMatrix testproblem(size, size);
    PETScWrappers::SparseMatrix  A(dim, dim, 5);
    testproblem.five_point(A);

    PETScWrappers::Vector  f(dim);
    PETScWrappers::Vector  u(dim);
    u = 0.;
    f = 1.;
    A.compress ();
    f.compress ();
    u.compress ();

    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn);
//    solver.set_symmetric_mode(true);
    solver.solve(A,u,f);

    PETScWrappers::Vector  tmp(dim);
    deallog << "residual = " << A.residual (tmp, u, f)
	    << std::endl;
  }
  
}

