/*----------------------------   solver_richardson.h     ---------------------------*/
/*   $Id$              */
/*            Ralf Hartmann, University of Heidelberg                               */
#ifndef __solver_richardson_H
#define __solver_richardson_H
/*----------------------------   solver_richardson.h     ---------------------------*/



#include <lac/solver.h>
#include <lac/solver_control.h>



/**
 * Implementation of the richardson iteration method. The stopping criterion
 * is the norm of the residual.
 *
 * @author Ralf Hartmann
 */
template<class Matrix, class Vector>
class SolverRichardson : public Solver<Matrix, Vector>
{
  public:
				     /**
				      * Constructor.
				      */
    SolverRichardson (SolverControl &cn, VectorMemory<Vector> &mem) :
		    Solver<Matrix,Vector> (cn,mem), omega(1.)
      {};

				     /**
				      * Solve $Ax=b$ for $x$.
				      */
    template<class Preconditioner>
    ReturnState solve (const Matrix &A,
		       Vector       &x,
		       const Vector &b,
		       const Preconditioner& precondition);

				     /**
				      * Set the damping-coefficient.
				      * Default is 1., i.e. no damping.
				      */
    void set_omega(double om=1.);
    
  protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();
    
				     /**
				      * Temporary vectors, allocated through
				      * the #VectorMemory# object at the start
				      * of the actual solution process and
				      * deallocated at the end.
				      */
    Vector *Vr, *Vd;

				     /**
				      * Damping-coefficient.
				      */
    double omega;  ;
    
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




/*----------------- Implementation of the Richardson Method ------------------*/


template<class Matrix,class Vector>
template<class Preconditioner>
Solver<Matrix,Vector>::ReturnState 
SolverRichardson<Matrix,Vector>::solve (const Matrix &A,
					Vector       &x,
					const Vector &b,
					const Preconditioner& precondition)
{
  SolverControl::State conv=SolverControl::iterate;

				   // Memory allocation
  Vr  = memory.alloc(); Vector& r  = *Vr; r.reinit(x);
  Vd  = memory.alloc(); Vector& d  = *Vd; d.reinit(x);

  deallog.push("Richardson");

				   // Main loop
  for(int iter=0; conv==SolverControl::iterate; iter++)
    {
      A.residual(r,x,b);

      res2 = r*r;
      conv = control().check (iter, criterion());
      if (conv != SolverControl::iterate)
	break;

      precondition(d,r);
      x.add(omega,d);
    }

				   // Deallocate Memory
  memory.free(Vr);
  memory.free(Vd);

  deallog.pop();
				   // Output
  if (conv == SolverControl::failure)
    return exceeded;
  else
    return success;
}


template<class Matrix,class Vector>
inline double
SolverRichardson<Matrix,Vector>::criterion()
{
  return sqrt(res2);
}

template<class Matrix,class Vector>
inline void
SolverRichardson<Matrix,Vector>::set_omega(double om)
{
  omega=om;
}


/*------------------   solver_richardson.h     ----------------------*/
/* end of #ifndef __solver_richardson_H */
#endif
/*------------------   solver_richardson.h     ----------------------*/
