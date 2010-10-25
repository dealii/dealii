//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__solver_relaxation_h
#define __deal2__solver_relaxation_h


#include <base/config.h>
#include <base/logstream.h>
#include <lac/solver.h>
#include <lac/solver_control.h>
#include <base/subscriptor.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of an iterative solver based on relaxation
 * methods. The stopping criterion is the norm of the residual.
 *
 * For the requirements on matrices and vectors in order to work with
 * this class, see the documentation of the Solver base class.
 *
 * Like all other solver classes, this class has a local structure called
 * @p AdditionalData which is used to pass additional parameters to the
 * solver, like damping parameters or the number of temporary vectors. We
 * use this additional structure instead of passing these values directly
 * to the constructor because this makes the use of the @p SolverSelector and
 * other classes much easier and guarantees that these will continue to
 * work even if number or type of the additional parameters for a certain
 * solver changes. AdditionalData of this class currently does not
 * contain any data.
 *
 * @ingroup Solvers
 * @author Guido Kanschat
 * @date 2010
 */
template <class VECTOR = Vector<double> >
class SolverRelaxation : public Solver<VECTOR>
{
  public:
				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver. There is no data in
				      * here for relaxation methods.
				      */
    struct AdditionalData {};
	
				     /**
				      * Constructor.
				      */
    SolverRelaxation (SolverControl        &cn,
		      const AdditionalData &data=AdditionalData());

				     /**
				      * Virtual destructor.
				      */
    virtual ~SolverRelaxation ();
    
				     /**
				      * Solve the system $Ax = b$
				      * using the relaxation method
				      * $x_{k+1} = R(x_k,b)$. The
				      * amtrix <i>A</i> itself is only
				      * used to compute the residual.
				      */
    template<class MATRIX, class RELAXATION>
    void
    solve (const MATRIX& A,
	   VECTOR& x,
	   const VECTOR& b,
	   const RELAXATION& R);
};

//----------------------------------------------------------------------//

template <class VECTOR>
SolverRelaxation<VECTOR>::SolverRelaxation(SolverControl &cn,
					   const AdditionalData &)
		:
		Solver<VECTOR> (cn)
{}



template <class VECTOR>
SolverRelaxation<VECTOR>::~SolverRelaxation()
{}


template <class VECTOR>
template <class MATRIX, class RELAXATION>
void
SolverRelaxation<VECTOR>::solve (
  const MATRIX& A,
  VECTOR& x,
  const VECTOR& b,
  const RELAXATION &R)
{
  GrowingVectorMemory<VECTOR> mem;
  SolverControl::State conv=SolverControl::iterate;

				   // Memory allocation
  typename VectorMemory<VECTOR>::Pointer Vr(mem); VECTOR& r  = *Vr; r.reinit(x);
  typename VectorMemory<VECTOR>::Pointer Vd(mem); VECTOR& d  = *Vd; d.reinit(x);
  
  deallog.push("Relaxation");
  
  try 
    {
				       // Main loop
      for(int iter=0; conv==SolverControl::iterate; iter++)
	{
					   // Compute residual
	  A.vmult(r,x);
	  r.sadd(-1.,1.,b);

					   // The required norm of the
					   // (preconditioned)
					   // residual is computed in
					   // criterion() and stored
					   // in res.
	  conv = this->control().check (iter, r.l2_norm());
	  if (conv != SolverControl::iterate)
	    break;
	  R.step(x,b);
	}
    }
  catch (...)
    {
      deallog.pop();
      throw;
    }
  deallog.pop();

				   // in case of failure: throw
				   // exception
  if (this->control().last_check() != SolverControl::success)
    throw SolverControl::NoConvergence (this->control().last_step(),
					this->control().last_value());
				   // otherwise exit as normal
}


DEAL_II_NAMESPACE_CLOSE

#endif
