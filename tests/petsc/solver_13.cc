//----------------------------  petsc_solver_13.cc  ---------------------------
//    $Id: solver_11.cc 20442 2010-01-25 10:46:00Z heister $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_solver_13.cc  ---------------------------

// test the PETScWrapper::Precondition*::vmult


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

  solver.solve(A,u,f,P);

  deallog << "Solver stopped after " << solver.control().last_step()
          << " iterations" << std::endl;
}


int main(int argc, char **argv)
{
  std::ofstream logfile("solver_13/output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  PetscInitialize(&argc,&argv,0,0);
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
    u = 0.;
    f = 1.;
    A.compress ();
    f.compress ();
    u.compress ();
    PETScWrappers::Vector t=u;

    PETScWrappers::SolverPreOnly solver(control);
    
    PETScWrappers::PreconditionJacobi preconditioner(A);
    
    check_solve (solver, A,u,f, preconditioner);

    deallog << u.l2_norm() << std::endl;
    
    check_solve (solver, A,t,u, preconditioner);

    deallog << t.l2_norm() << std::endl;

    u = 0;
    preconditioner.vmult(u, f);
    deallog << u.l2_norm() << std::endl;
    preconditioner.vmult(t, u);
    deallog << t.l2_norm() << std::endl;    
  }
  PetscFinalize ();
}

