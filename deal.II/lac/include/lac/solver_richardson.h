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
 * The use of the #AdditionalData# struct is described in the #solver#
 * base class.
 *
 * @author Ralf Hartmann
 */
template<class Matrix, class Vector>
class SolverRichardson : public Solver<Matrix, Vector>
{
  public:
				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver.
				      */
    struct AdditionalData 
    {
	AdditionalData(double omega=1):
			omega(omega) {};
	
					 /**
					  * Damping parameter.
					  */
	double omega;
    };
	
				     /**
				      * Constructor.
				      */
    SolverRichardson (SolverControl &cn,
		      VectorMemory<Vector> &mem,
		      const AdditionalData &data=AdditionalData());

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
				      * Damping parameter.
				      */
    AdditionalData additional_data;
    
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


template<class Matrix, class Vector>
SolverRichardson<Matrix,Vector>::SolverRichardson(SolverControl &cn,
						  VectorMemory<Vector> &mem,
						  const AdditionalData &data):
		Solver<Matrix,Vector> (cn,mem),
		additional_data(data)  {};


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
      x.add(additional_data.omega,d);
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
  additional_data.omega=om;
}


/*------------------   solver_richardson.h     ----------------------*/
/* end of #ifndef __solver_richardson_H */
#endif
/*------------------   solver_richardson.h     ----------------------*/
