//----------------------------  solver.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver.cc  ---------------------------


#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_richardson.h>
#include <lac/solver_qmrs.h>
#include <lac/precondition.h>

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_solve( SOLVER& solver, const MATRIX& A,
	     VECTOR& u, VECTOR& f, const PRECONDITION& P)
{
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
}

template<class SOLVER, class MATRIX, class VECTOR, class PRECONDITION>
void
check_Tsolve(SOLVER& solver, const MATRIX& A,
	     VECTOR& u, VECTOR& f, const PRECONDITION& P)
{
  u = 0.;
  f = 1.;
  try 
    {
      solver.Tsolve(A,u,f,P);
    }
  catch (std::exception& e)
    {
      deallog << e.what() << std::endl;
    }  
}

int main()
{
  std::ofstream logfile("solver.output");
//  logfile.setf(std::ios::fixed);
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  GrowingVectorMemory<> mem;
  SolverControl control(100, 1.e-3);
  SolverControl verbose_control(100, 1.e-3, true);
  SolverCG<> cg(control, mem);
  SolverGMRES<> gmres(control, mem, 8);
  SolverGMRES<>::AdditionalData(8, true);
  SolverGMRES<> gmresright(control, mem, 8);  
  SolverBicgstab<> bicgstab(control, mem);
  SolverRichardson<> rich(control, mem);
  SolverQMRS<> qmrs(control, mem);

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
      PreconditionSOR<> prec_sor;
      prec_sor.initialize(A, 1.2);
      PreconditionSSOR<> prec_ssor;
      prec_ssor.initialize(A, 1.2);

      std::vector<unsigned int> permutation(dim);
      std::vector<unsigned int> inverse_permutation(dim);

				       // Create a permutation: Blocks
				       // backwards and every second
				       // block backwards
      unsigned int k = 0;
      for (unsigned int i=0;i<size-1;++i)
	for (unsigned int j=0;j<size-1;++j)
	  {
	    if (i % 2)
	      permutation[k++] = (size-i-2) * (size-1) + j;
	    else
	      permutation[k++] = (size-i-2) * (size-1) + size-j-2;
	  }
      

      for (unsigned int i=0;i<dim;++i)
	inverse_permutation[permutation[i]] = i;

      PreconditionPSOR<> prec_psor;
      prec_psor.initialize(A, permutation, inverse_permutation, 1.2);
      
      Vector<double>  f(dim);
      Vector<double>  u(dim);
      Vector<double> res(dim);

      f = 1.;
      u = 1.;
      
      A.residual(res,u,f);
      A.SOR(res);
      res.add(1.,u);
      A.SOR_step(u,f);
      res.add(-1.,u);
    
      deallog << "SOR-diff:" << res*res << std::endl;

      try
	{
	  deallog.push("no-fail");

	  control.set_max_steps(10);
	  check_solve(cg,A,u,f,prec_no);
	  check_solve(bicgstab,A,u,f,prec_no);
	  check_solve(gmres,A,u,f,prec_no);
	  check_solve(gmresright,A,u,f,prec_no);
	  check_solve(qmrs,A,u,f,prec_no);
	  control.set_max_steps(100);
	  
	  deallog.pop();
	  
	  deallog.push("no");
	  
	  check_solve(cg,A,u,f,prec_no);
	  check_solve(bicgstab,A,u,f,prec_no);
	  check_solve(gmres,A,u,f,prec_no);
	  check_solve(gmresright,A,u,f,prec_no);
	  check_solve(qmrs,A,u,f,prec_no);
	  
	  deallog.pop();
	  
	  deallog.push("ssor");
	  
	  check_Tsolve(rich,A,u,f,prec_ssor);
	  check_solve(rich,A,u,f,prec_ssor);
	  check_solve(cg,A,u,f,prec_ssor);
	  check_solve(bicgstab,A,u,f,prec_ssor);
	  check_solve(gmres,A,u,f,prec_ssor);
	  check_solve(gmresright,A,u,f,prec_ssor);
	  check_solve(qmrs,A,u,f,prec_ssor);
	  
	  deallog.pop();
	  
	  deallog.push("sor");
	  
	  check_Tsolve(rich,A,u,f,prec_sor);
	  check_solve(rich,A,u,f,prec_sor);
	  check_solve(cg,A,u,f,prec_sor);
	  check_solve(bicgstab,A,u,f,prec_sor);
	  check_solve(gmres,A,u,f,prec_sor);
	  check_solve(gmresright,A,u,f,prec_sor);
	  
	  deallog.pop();
	  
	  deallog.push("psor");
	  
	  check_Tsolve(rich,A,u,f,prec_psor);
	  check_solve(rich,A,u,f,prec_psor);
	  check_solve(cg,A,u,f,prec_psor);
	  check_solve(bicgstab,A,u,f,prec_psor);
	  check_solve(gmres,A,u,f,prec_psor);
	  check_solve(gmresright,A,u,f,prec_psor);
	  
	  deallog.pop();
	}
      catch (std::exception& e)
	{
	  std::cerr << e.what() << std::endl;
	}
    };
}

