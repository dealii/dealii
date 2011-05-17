//----------------------------  petsc_deal_solver_06.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_deal_solver_06.cc  ---------------------------

// test the Richardson solver using the PETSc matrix and vector classes


#include "../tests.h"
#include "../lac/testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/precondition.h>
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
    }

  deallog << "Solver stopped after " << solver.control().last_step()
          << " iterations" << std::endl;
}


int main(int argc, char **argv)
{
  std::ofstream logfile("deal_solver_06/output");
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  PetscInitialize(&argc,&argv,0,0);
  {
    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);
    SolverControl control(10000, 1.e-3);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

                                     // Make matrix
    FDMatrix testproblem(size, size);
    PETScWrappers::SparseMatrix  A(dim, dim, 5);
    testproblem.five_point(A);

    PETScWrappers::Vector  f(dim);
    PETScWrappers::Vector  u(dim);
    f = 1.;
    A.compress ();
    f.compress ();
    u.compress ();

    GrowingVectorMemory<PETScWrappers::Vector> mem;
    SolverRichardson<PETScWrappers::Vector> solver(control,mem);
    PreconditionIdentity preconditioner;
    check_solve (solver, A,u,f, preconditioner);
  }
  GrowingVectorMemory<PETScWrappers::Vector>::release_unused_memory ();
  PetscFinalize ();
}

