//----------------------------------------------------------------------
//    $Id: solver.cc 21548 2010-07-20 19:58:37Z kanschat $
//
//    Copyright (C) 1998 - 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// test all solver with a start value that is already a solution (i.e. 0
// iterations). This caused a memory leak in FGMRES.


#include "../tests.h"
#include "testmatrix.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/solver_minres.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_richardson.h>
#include <lac/solver_qmrs.h>
#include <lac/precondition.h>

template<class SOLVER, class MATRIX, class VECTOR>
void
check_solve( const MATRIX& A,
	     VECTOR& u, VECTOR& f)
{
  GrowingVectorMemory<> mem;
  SolverControl control(100, 1.e-3);
  SOLVER solver(control, mem);
  PreconditionIdentity prec_no;
  u = 0.;
  f = 0.;
  
  try 
    {
      solver.solve(A, u, f, prec_no);
    }
  catch (std::exception& e)
    {
      deallog << e.what() << std::endl;
    }  
}

int main()
{
  std::ofstream logfile("solver_leak/output");
//  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
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
      
      Vector<double>  f(dim);
      Vector<double>  u(dim);
      Vector<double> res(dim);

      deallog.push("alreadydone");      
      check_solve<SolverCG<> >(A,u,f);
      check_solve<SolverGMRES<> >(A,u,f);
//      check_solve<SolverFGMRES<> >(A,u,f);
      check_solve<SolverBicgstab<> >(A,u,f);
      check_solve<SolverQMRS<> >(A,u,f);

      deallog.pop();
    }

}

