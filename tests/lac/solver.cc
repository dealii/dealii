// $Id$
//
// Test program for linear solvers.

#include <cmath>
#include <fstream>
#include <iomanip>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparsematrix.h>
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

main()
{
  ofstream logfile("solver.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  PrimitiveVectorMemory<Vector<double>  > mem;
  SolverControl control(50, 1.e-5);
  SolverCG<SparseMatrix<float> , Vector<double>  > cg(control, mem);
  SolverGMRES<SparseMatrix<float> , Vector<double>  > gmres(control, mem,20);
  SolverBicgstab<SparseMatrix<float> , Vector<double>  > bicgstab(control, mem);
  SolverRichardson<SparseMatrix<float> , Vector<double>  > rich(control, mem);
  SolverQMRS<SparseMatrix<float> , Vector<double>  > qmrs(control, mem);

  for (unsigned int size=4; size <= 40; size *= 3)
    {
      deallog << "Size " << size << endl;
      
      unsigned int dim = (size-1)*(size-1);

				       // Make matrix
      FDMatrix testproblem(size, size);
      SparseMatrixStruct structure(dim, dim, 5);
      testproblem.build_structure(structure);
      structure.compress();
      SparseMatrix<float>  A(structure);
      testproblem.laplacian(A);

      PreconditionIdentity<Vector<double>  >
	prec_no;
      PreconditionRelaxation<SparseMatrix<float> , Vector<double> >
	prec_sor(A, &SparseMatrix<float> ::template precondition_SOR<double>, 1.2);
      PreconditionRelaxation<SparseMatrix<float> , Vector<double> >
	prec_ssor(A, &SparseMatrix<float> ::template precondition_SSOR<double>, 1.2);
      
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
      {
	check_method(cg,A,u,f,prec_no);
	check_method(bicgstab,A,u,f,prec_no);
	check_method(gmres,A,u,f,prec_no);
	check_method(qmrs,A,u,f,prec_no);
      };
      deallog.pop();
      
      deallog.push("ssor");      
      {
	check_method(rich,A,u,f,prec_ssor);
	check_method(cg,A,u,f,prec_ssor);
	check_method(bicgstab,A,u,f,prec_ssor);
	check_method(gmres,A,u,f,prec_ssor);
	check_method(qmrs,A,u,f,prec_ssor);
      };
      deallog.pop();

      deallog.push("sor");      
      {
	check_method(rich,A,u,f,prec_sor);
	check_method(cg,A,u,f,prec_sor);
	check_method(bicgstab,A,u,f,prec_sor);
	check_method(gmres,A,u,f,prec_sor);
	check_method(qmrs,A,u,f,prec_sor);
      };
      deallog.pop();
    };
};

