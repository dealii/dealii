//----------------------------  petsc_deal_solver_05.cc  ---------------------------
//    petsc_deal_solver_05.cc,v 1.33 2003/05/30 19:19:16 guido Exp
//    Version: 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_deal_solver_05.cc  ---------------------------

// test the QMRS solver


#include "../tests.h"
#include "../lac/testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <base/logstream.h>
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_vector.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_richardson.h>
#include <lac/solver_qmrs.h>
#include <lac/precondition.h>
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
  PetscInitialize(&argc,&argv,0,0);
  {
    std::ofstream logfile("petsc_deal_solver_05.output");
    logfile.precision(4);
    deallog.attach(logfile);
    deallog.depth_console(0);
  
    SolverControl control(100, 1.e-3);

    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;
      
                                     // Make matrix
    FDMatrix testproblem(size, size);
    petsc_wrappers::SparseMatrix  A(dim, dim, 5);
    testproblem.five_point(A);

    petsc_wrappers::Vector  f(dim);
    petsc_wrappers::Vector  u(dim);
    f = 1.;
    A.compress ();
    f.compress ();
    u.compress ();

    GrowingVectorMemory<petsc_wrappers::Vector> mem;
    SolverQMRS<petsc_wrappers::Vector> solver(control,mem);
    PreconditionIdentity preconditioner;
    check_solve (solver, A,u,f, preconditioner);
  }
  PetscFinalize ();
}

