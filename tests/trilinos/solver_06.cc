//----------------------------  trilinos_solver_06.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_solver_06.cc  ---------------------------

// test the Trilinos Bicgstab solver


#include "../tests.h" 
#include <base/utilities.h>
#include "../lac/testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <base/logstream.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_solver.h>
#include <lac/trilinos_precondition.h>
#include <lac/vector_memory.h>
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
  std::ofstream logfile("solver_06/output");
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);

  
  {
    SolverControl control(100, 1.e-3);

    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;
      
                                     // Make matrix
    FDMatrix testproblem(size, size);
    TrilinosWrappers::SparseMatrix  A(dim, dim, 5U);
    testproblem.five_point(A);

    TrilinosWrappers::Vector  f(dim);
    TrilinosWrappers::Vector  u(dim);
    f = 1.;
    A.compress ();
    f.compress ();
    u.compress ();

    TrilinosWrappers::SolverBicgstab solver(control);
    TrilinosWrappers::PreconditionJacobi preconditioner;
    preconditioner.initialize(A);
    check_solve (solver, A,u,f, preconditioner);
  }
}

