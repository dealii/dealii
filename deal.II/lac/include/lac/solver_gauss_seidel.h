/*----------------------------   solver_gauss_seidel.h     ---------------------------*/
/*      $Id$             */
#ifndef __solver_gauss_seidel_H
#define __solver_gauss_seidel_H
/*----------------------------   solver_gauss_seidel.h     ---------------------------*/

#include <base/trace.h>
#include <lac/solver_control.h>
#include <lac/solver.h>

/**
 * Implementation of the Gauss-Seidel method. The stopping criterion
 * is the norm of the defect, i.e. if x is the iteration vector, b the
 * rhs and A the matrix, the $res = \| A x - b \|$. 
 */
template<class Matrix, class Vector>
class SolverGaussSeidel : public Solver<Matrix, Vector> {
 public:
  /**
   * Constructor.  Takes a SolverControl and some
   VectorMemory. VectorMemory is not used. 
   */
  SolverGaussSeidel (SolverControl &cn, VectorMemory<Vector> &mem) :
    Solver<Matrix,Vector> (cn,mem)
    {};
  
  /**
   * Solver method. Just specify the vectors A, x and b and let it
   * run. x should contain a reasonable starting value.  
   */
    virtual ReturnState solve (const Matrix &A,
			       Vector       &x,
			       const Vector &b);

  protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();
    
				     /**
				      * Within the iteration loop, the
				      * square of the residual vector is
				      * stored in this variable. The
				      * function #criterion# uses this
				      * variable to compute the convergence
				      * value, which in this class is the
				      * norm of the residual vector and thus
				      * the square root of the #res2# value.
				      */
    double res2;
};

/*----------------- Implementation of the Gauss-Seidel Method ------------------------*/

template<class Matrix,class Vector> 
SolverGaussSeidel<Matrix,Vector>::ReturnState 
SolverGaussSeidel<Matrix,Vector>::solve (const Matrix &A,
					 Vector       &x,
					 const Vector &b) {

  deallog.push("lac/solver_gauss_seidel.h");
  deallog.push("SolverGaussSeidel::solve");

  unsigned int i;
  unsigned int n = b.size();
  double res;

  SolverControl::State conv=SolverControl::iterate;

  Vector defect(n);

  deallog.push("Main loop");

  // Main loop
  for(int iter=0; conv==SolverControl::iterate; iter++)
    {
      // Calculate defect, i.e. defect = A x - b
      A.vmult(defect,x);      
      defect.add(-1,b);    
      
      // Calculate residual
      res2 = defect * defect;
     
      // Apply Gauss-Seidel preconditioner
      for (i=0;i<n;i++)
	{
	  defect(i) = defect(i) / A(i,i);
	}

      // Correct defect
      x.add(-1,defect);

      // Check residual
      res = criterion();
      conv = (control().check(iter, res));
      if (conv != SolverControl::iterate)
	break;
    }
  
  // Output
  
  deallog.pop();

  if (conv == SolverControl::failure)
    {
      TRACEMSG("*** exceeded maximum number of iterations");
      deallog.pop();
      deallog.pop();
      return exceeded;
    }
  else
    {
      TRACEMSG("success");
      deallog.pop();
      deallog.pop();
      return success;
    }
} 

template<class Matrix,class Vector>
inline double
SolverGaussSeidel<Matrix,Vector>::criterion()
{
  return sqrt(res2);
}

/*----------------------------   solver_gauss_seidel.h     ---------------------------*/
/* end of #ifndef __solver_gauss_seidel_H */
#endif
/*----------------------------   solver_gauss_seidel.h     ---------------------------*/
