//----------------------------  solver_richardson.h  ---------------------------
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
//----------------------------  solver_richardson.h  ---------------------------
#ifndef __deal2__solver_richardson_h
#define __deal2__solver_richardson_h


#include <lac/solver.h>
#include <lac/solver_control.h>


/**
 * Implementation of the richardson iteration method. The stopping criterion
 * is the norm of the residual.
 *
 * Like all other solver classes, this class has a local structure called
 * @p{AdditionalData} which is used to pass additional parameters to the
 * solver, like damping parameters or the number of temporary vectors. We
 * use this additional structure instead of passing these values directly
 * to the constructor because this makes the use of the @p{SolverSelector} and
 * other classes much easier and guarantees that these will continue to
 * work even if number or type of the additional parameters for a certain
 * solver changes.
 *
 * For the Richardson method, the additional data is the damping parameter,
 * which is the only content of the @p{AdditionalData} structure. By default,
 * the constructor of the structure sets it to one.
 *
 * @author Ralf Hartmann
 */
template <class VECTOR = Vector<double> >
class SolverRichardson : private Solver<VECTOR>
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
		      VectorMemory<VECTOR> &mem,
		      const AdditionalData &data=AdditionalData());

				     /**
				      * Solve $Ax=b$ for $x$.
				      */
    template<class MATRIX, class PRECONDITIONER>
    typename Solver<VECTOR>::ReturnState solve (const MATRIX &A,
						VECTOR       &x,
						const VECTOR &b,
						const PRECONDITIONER& precondition);

				     /**
				      * Solve $A^Tx=b$ for $x$.
				      */
    template<class MATRIX, class PRECONDITIONER>
    typename Solver<VECTOR>::ReturnState Tsolve (const MATRIX &A,
						 VECTOR       &x,
						 const VECTOR &b,
						 const PRECONDITIONER& precondition);

				     /**
				      * Set the damping-coefficient.
				      * Default is 1., i.e. no damping.
				      */
    void set_omega(double om=1.);

				     /**
				      * Interface for derived class.
				      * This function gets the current
				      * iteration vector, the residual
				      * and the update vector in each
				      * step. It can be used for a
				      * graphical output of the
				      * convergence history.
				      */
    virtual void print_vectors(const unsigned int step,
			       const VECTOR& x,
			       const VECTOR& r,
			       const VECTOR& d) const;
    
  protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();
    
				     /**
				      * Temporary vectors, allocated through
				      * the @p{VectorMemory} object at the start
				      * of the actual solution process and
				      * deallocated at the end.
				      */
    VECTOR *Vr;
    VECTOR *Vd;

				     /**
				      * Damping parameter.
				      */
    AdditionalData additional_data;
    
				     /**
				      * Within the iteration loop, the
				      * norm of the residual is
				      * stored in this variable. The
				      * function @p{criterion} uses this
				      * variable to compute the convergence
				      * value, which in this class is the
				      * norm of the residual vector and thus
				      * the square root of the @p{res2} value.
				      */
    double res;
};


/*----------------- Implementation of the Richardson Method ------------------*/


template<class VECTOR>
SolverRichardson<VECTOR>::SolverRichardson(SolverControl &cn,
					   VectorMemory<VECTOR> &mem,
					   const AdditionalData &data):
		Solver<VECTOR> (cn,mem),
		additional_data(data)  {};


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
typename Solver<VECTOR>::ReturnState 
SolverRichardson<VECTOR>::solve (const MATRIX &A,
				 VECTOR       &x,
				 const VECTOR &b,
				 const PRECONDITIONER& precondition)
{
  SolverControl::State conv=SolverControl::iterate;

				   // Memory allocation
  Vr  = memory.alloc(); VECTOR& r  = *Vr; r.reinit(x);
  Vd  = memory.alloc(); VECTOR& d  = *Vd; d.reinit(x);

  deallog.push("Richardson");

				   // Main loop
  for(int iter=0; conv==SolverControl::iterate; iter++)
    {
      res=A.residual(r,x,b);

      conv = control().check (iter, criterion());
      if (conv != SolverControl::iterate)
	break;

      precondition.vmult(d,r);
      x.add(additional_data.omega,d);
      print_vectors(iter,x,r,d);
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


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
typename Solver<VECTOR>::ReturnState 
SolverRichardson<VECTOR>::Tsolve (const MATRIX &A,
				  VECTOR       &x,
				  const VECTOR &b,
				  const PRECONDITIONER& precondition)
{
  SolverControl::State conv=SolverControl::iterate;

				   // Memory allocation
  Vr  = memory.alloc(); VECTOR& r  = *Vr; r.reinit(x);
  Vd  = memory.alloc(); VECTOR& d  = *Vd; d.reinit(x);

  deallog.push("Richardson");

				   // Main loop
  for(int iter=0; conv==SolverControl::iterate; iter++)
    {
				       // Do not use Tresidual,
				       // but do it in 2 steps
      A.Tvmult(r,x);
      r.sadd(-1.,1.,b);
      res=sqrt(r*r);

      conv = control().check (iter, criterion());
      if (conv != SolverControl::iterate)
	break;

      precondition.Tvmult(d,r);
      x.add(additional_data.omega,d);
      print_vectors(iter,x,r,d);
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


template<class VECTOR>
void
SolverRichardson<VECTOR>::print_vectors(const unsigned int,
					const VECTOR&,
					const VECTOR&,
					const VECTOR&) const
{}



template<class VECTOR>
inline double
SolverRichardson<VECTOR>::criterion()
{
  return res;
}


template<class VECTOR>
inline void
SolverRichardson<VECTOR>::set_omega(double om)
{
  additional_data.omega=om;
}


#endif
