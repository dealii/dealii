//----------------------------  solver.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver.cc  ---------------------------


#include <cmath>
#include <fstream>
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
check_method( SOLVER& solver, const MATRIX& A,
	     VECTOR& u, VECTOR& f, const PRECONDITION& P)
{
  u = 0.;
  f = 1.;
  solver.solve(A,u,f,P);
}

int main()
{
  ofstream logfile("solver.output");
  logfile.setf(ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  GrowingVectorMemory<> mem;
  SolverControl control(100, 1.e-3);
  SolverControl verbose_control(100, 1.e-3, true);
  SolverCG<> cg(control, mem);
  SolverGMRES<> gmres(control, mem,20);
  SolverBicgstab<> bicgstab(control, mem);
  SolverRichardson<> rich(control, mem);
  SolverQMRS<> qmrs(control, mem);

  for (unsigned int size=4; size <= 40; size *= 3)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;
      
				       // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.build_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.laplacian(A);

      PreconditionIdentity prec_no;
      PreconditionSOR<> prec_sor;
      prec_sor.initialize(A, 1.2);
      PreconditionSSOR<> prec_ssor;
      prec_ssor.initialize(A, 1.2);
      
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
      
      deallog.push("no");

      check_method(cg,A,u,f,prec_no);
      check_method(bicgstab,A,u,f,prec_no);
      check_method(gmres,A,u,f,prec_no);
      check_method(qmrs,A,u,f,prec_no);

      deallog.pop();
      
      deallog.push("ssor");

      check_method(rich,A,u,f,prec_ssor);
      check_method(cg,A,u,f,prec_ssor);
      check_method(bicgstab,A,u,f,prec_ssor);
      check_method(gmres,A,u,f,prec_ssor);
      check_method(qmrs,A,u,f,prec_ssor);

      deallog.pop();

      deallog.push("sor");

      check_method(rich,A,u,f,prec_sor);
      check_method(cg,A,u,f,prec_sor);
      check_method(bicgstab,A,u,f,prec_sor);
      check_method(gmres,A,u,f,prec_sor);
      check_method(qmrs,A,u,f,prec_sor);

      deallog.pop();
    };
};

