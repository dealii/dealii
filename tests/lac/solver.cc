// $Id$
//
// Test program for linear solvers.

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
  SolverControl control(100, 1.e-5);
  SolverControl verbose_control(100, 1.e-5, true);
  SolverCG<> cg(control, mem);
  SolverGMRES<> gmres(control, mem,20);
  SolverBicgstab<> bicgstab(control, mem);
  SolverRichardson<> rich(control, mem);
  SolverQMRS<> qmrs(control, mem);

  for (unsigned int size=4; size <= 40; size *= 3)
    {
      unsigned int dim = (size-1)*(size-1);

      deallog << "Size " << size << " Unknowns " << dim << endl;
      
				       // Make matrix
      FDMatrix testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.build_structure(structure);
      structure.compress();
      SparseMatrix<double>  A(structure);
      testproblem.laplacian(A);

      PreconditionIdentity
	prec_no;
      PreconditionRelaxation<>
	prec_sor(A, &SparseMatrix<double>::template precondition_SOR<double>, 1.2);
      PreconditionRelaxation<>
	prec_ssor(A, &SparseMatrix<double>::template precondition_SSOR<double>, 1.2);
      
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
    
      deallog << "SOR-diff:" << res*res << endl;
      
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

