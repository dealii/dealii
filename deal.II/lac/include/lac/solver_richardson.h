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
 * Like all other solver classes, this class has a local structure called
 * #AdditionalData# which is used to pass additional parameters to the
 * solver, like damping parameters or the number of temporary vectors. We
 * use this additional structure instead of passing these values directly
 * to the constructor because this makes the use of the #SolverSelector# and
 * other classes much easier and guarantees that these will continue to
 * work even if number or type of the additional parameters for a certain
 * solver changes.
 *
 * For the Richardson method, the additional data is the damping parameter,
 * which is the only content of the #AdditionalData# structure. By default,
 * the constructor of the structure sets it to one.
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
					 /**
					  * Constructor. By default,
					  * set the damping parameter
					  * to one.
					  */
	AdditionalData(double omega=1):
			omega(omega)
	  {};
	
					 /**
					  * Relaxation parameter.
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
				      * norm of the residual is
				      * stored in this variable. The
				      * function #criterion# uses this
				      * variable to compute the convergence
				      * value, which in this class is the
				      * norm of the residual vector and thus
				      * the square root of the #res2# value.
				      */
    double res;
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
      res=A.residual(r,x,b);

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
  return res;
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
