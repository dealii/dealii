// $Id$
//
// Test program for linear solvers.

#include <cmath>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparsematrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/solver_pcg.h>
#include <lac/solver_pgmres.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_richardson.h>
#include <lac/precondition.h>

#define MATRIX SparseMatrix<float> 
#define VECTOR Vector<double> 

main()
{
  PrimitiveVectorMemory<VECTOR > mem;
  SolverControl control(100, 1.e-5, true);
  SolverPCG<MATRIX, VECTOR > cg(control, mem);
  SolverGMRES<MATRIX, VECTOR > gmres(control, mem,20);
  SolverBicgstab<MATRIX, VECTOR > bicgstab(control, mem);
  SolverRichardson<MATRIX, VECTOR > rich(control, mem);

  for (unsigned int size=10; size <= 10; size *= 10)
    {
      deallog << "Size " << size << endl;
      
      unsigned int dim = (size+1)*(size+1);

				       // Make matrix
      FDMatrix testproblem(size, size);
      SparseMatrixStruct structure(dim, dim, 5);
      testproblem.build_structure(structure);
      structure.compress();
      MATRIX A(structure);
      testproblem.laplacian(A);

      PreconditionIdentity<VECTOR >
	prec_no;
      PreconditionRelaxation<MATRIX, VECTOR>
	prec_ssor(A, &MATRIX::template precondition_SSOR<double>, 1.2);
      
      VECTOR f(dim);
      VECTOR u(dim);
      
      deallog.push("no");

      f = 1.;
      u = 0.;
      cg.solve(A,u,f,prec_no);

      f = 1.;
      u = 0.;
      bicgstab.solve(A,u,f,prec_no);

      f = 1.;
      u = 0.;
      gmres.solve(A,u,f,prec_no);

      deallog.pop();
      deallog.push("ssor");      

      f = 1.;
      u = 0.;
      rich.solve(A,u,f,prec_ssor);

      f = 1.;
      u = 0.;
      cg.solve(A,u,f,prec_ssor);

      f = 1.;
      u = 0.;
      bicgstab.solve(A,u,f,prec_ssor);

      f = 1.;
      u = 0.;
      gmres.solve(A,u,f,prec_ssor);

      deallog.pop();
    }
}
