//----------------------------  eigen.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  eigen.h  ---------------------------
#ifndef __deal2__eigen_h
#define __deal2__eigen_h


#include <lac/shifted_matrix.h>
#include <lac/solver.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

/**
 * Power method (von Mises).
 *
 * This method determines the largest eigenvalue of a matrix by
 * applying increasing powers of this matrix to a vector. If there is
 * an eigenvalue $l$ with dominant absolute value, the iteration vectors
 * will become aligned to its eigenspace and $Ax = lx$.
 *
 * A shift parameter allows to shift the spectrum, so it is possible
 * to compute the smallest eigenvalue, too.
 *
 * Convergence of this method is known to be slow.
 *
 * @author Guido Kanschat, 2000
 */
template <class VECTOR = Vector<double> >
class EigenPower : private Solver<VECTOR>
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
					  * Shift parameter. This
					  * parameter allows to shift
					  * the spectrum to compute a
					  * different eigenvalue.
					  */
	double shift;
					 /**
					  * Constructor. Set the shift parameter.
					  */
	AdditionalData (const double shift = 0.):
			shift(shift)
	  {}
	
    };

				     /**
				      * Constructor.
				      */
    EigenPower (SolverControl &cn,
		VectorMemory<VECTOR> &mem,
		const AdditionalData &data=AdditionalData());

				     /**
				      * Virtual destructor.
				      */
    virtual ~EigenPower ();

				     /**
				      * Power method. @p{x} is the
				      * (not necessarily normalized)
				      * start vector for the power
				      * method. After the iteration,
				      * @p{value} is the approximated
				      * eigenvalue and @p{x} is the
				      * corresponding eigenvector,
				      * normalized with respect to the
				      * l2-norm.
				      */
    template <class MATRIX>
    void
    solve (double       &value,
	   const MATRIX &A,
	   VECTOR       &x);

  protected:
				     /**
				      * Shift parameter.
				      */
    AdditionalData additional_data;
};


/**
 * Inverse iteration (Wieland).
 *
 * This class implements an adaptive version of the inverse iteration by Wieland.
 *
 * There are two choices for the stopping criterion: by default, the
 * norm of the residual $A x - l x$ is computed. Since this might not
 * converge to zero for non-symmetric matrices with non-trivial Jordan
 * blocks, it can be replaced by checking the difference of successive
 * eigenvalues. Use @p{AdditionalData::use_residual} for switching
 * these options.
 *
 * @author Guido Kanschat, 2000
 */
template <class VECTOR = Vector<double> >
class EigenInverse : private Solver<VECTOR>
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
					  * Flag for the stopping criterion.
					  */
	bool use_residual;
					 /**
					  * Constructor.
					  */
	AdditionalData (bool use_residual = true):
			use_residual(use_residual)
	  {}
	
    };
    
  				     /**
				      * Constructor.
				      */
    EigenInverse (SolverControl &cn,
		  VectorMemory<VECTOR> &mem,
		  const AdditionalData &data=AdditionalData());


				     /**
				      * Virtual destructor.
				      */
    virtual ~EigenInverse ();

				     /**
				      * Inverse method. @p{value} is
				      * the start guess for the
				      * eigenvalue and @p{x} is the
				      * (not necessarily normalized)
				      * start vector for the power
				      * method. After the iteration,
				      * @p{value} is the approximated
				      * eigenvalue and @p{x} is the
				      * corresponding eigenvector,
				      * normalized with respect to the
				      * l2-norm.
				      */
    template <class MATRIX>
    void
    solve (double       &value,
	   const MATRIX &A,
	   VECTOR       &x);

  protected:
				     /**
				      * Shift parameter.
				      */
    AdditionalData additional_data;
};

//----------------------------------------------------------------------


template <class VECTOR>
EigenPower<VECTOR>::EigenPower (SolverControl &cn,
				VectorMemory<VECTOR> &mem,
				const AdditionalData &data):
		Solver<VECTOR>(cn, mem),
		additional_data(data)
{}


template <class VECTOR>
EigenPower<VECTOR>::~EigenPower ()
{}


template <class VECTOR>
template <class MATRIX>
void
EigenPower<VECTOR>::solve (double       &value,
			   const MATRIX &A,
			   VECTOR       &x)
{
  SolverControl::State conv=SolverControl::iterate;

  deallog.push("Power method");

  VECTOR* Vy = memory.alloc (); VECTOR& y = *Vy; y.reinit (x);
  VECTOR* Vr = memory.alloc (); VECTOR& r = *Vr; r.reinit (x);
  
  double length = x.l2_norm ();
  double old_length = 0.;
  x.scale(1./length);
  
  A.vmult (y,x);
  
				   // Main loop
  for(int iter=0; conv==SolverControl::iterate; iter++)
    {
      y.add(additional_data.shift, x);
      
				       // Compute absolute value of eigenvalue
      old_length = length;
      length = y.l2_norm ();

				       // do a little trick to compute the sign
				       // with not too much effect of round-off errors.
      double entry = 0.;
      unsigned int i = 0;
      double thresh = length/x.size();
      do 
	{
	  Assert (i<x.size(), ExcInternalError());
	  entry = y (i++);
	}
      while (fabs(entry) < thresh);

      --i;

				       // Compute unshifted eigenvalue
      value = (entry * x (i) < 0.) ? -length : length;
      value -= additional_data.shift;

				       // Update normalized eigenvector
      x.equ (1/length, y);

				       // Compute residual
      A.vmult (y,x);

				       // Check the change of the eigenvalue
				       // Brrr, this is not really a good criterion
      conv = control().check (iter, fabs(1./length-1./old_length));
    }
  
  memory.free(Vy);
  memory.free(Vr);

  deallog.pop();

				   // in case of failure: throw
				   // exception
  if (control().last_check() != SolverControl::success)
    throw SolverControl::NoConvergence (control().last_step(),
					control().last_value());
				   // otherwise exit as normal
}

//----------------------------------------------------------------------//

template <class VECTOR>
EigenInverse<VECTOR>::EigenInverse (SolverControl &cn,
					    VectorMemory<VECTOR> &mem,
					    const AdditionalData &data):
		Solver<VECTOR>(cn, mem),
		additional_data(data)
{}


template <class VECTOR>
EigenInverse<VECTOR>::~EigenInverse ()
{}


template <class VECTOR>
template <class MATRIX>
void
EigenInverse<VECTOR>::solve (double       &value,
			     const MATRIX &A,
			     VECTOR       &x)
{
  deallog.push("Wieland");

  SolverControl::State conv=SolverControl::iterate;

				   // Prepare matrix for solver
  ShiftedMatrix <MATRIX> A_s(A, -value);

				   // Define solver
  ReductionControl inner_control (100, 1.e-16, 1.e-8, false, false);
  PreconditionIdentity prec;
  SolverCG<VECTOR>
    solver(inner_control, memory);

				   // Next step for recomputing the shift
  unsigned int goal = 10;
  
				   // Auxiliary vector
  VECTOR* Vy = memory.alloc (); VECTOR& y = *Vy; y.reinit (x);
  VECTOR* Vr = memory.alloc (); VECTOR& r = *Vr; r.reinit (x);
  
  double length = x.l2_norm ();
  double old_length = 0.;
  double old_value = value;
  
  x.scale(1./length);
  
				   // Main loop
  for (unsigned int iter=0; conv==SolverControl::iterate; iter++)
    {
      solver.solve (A_s, y, x, prec);
      
				       // Compute absolute value of eigenvalue
      old_length = length;
      length = y.l2_norm ();

				       // do a little trick to compute the sign
				       // with not too much effect of round-off errors.
      double entry = 0.;
      unsigned int i = 0;
      double thresh = length/x.size();
      do 
	{
	  Assert (i<x.size(), ExcInternalError());
	  entry = y (i++);
	}
      while (fabs(entry) < thresh);

      --i;

				       // Compute unshifted eigenvalue
      value = (entry * x (i) < 0.) ? -length : length;
      value = 1./value;      
      value -= A_s.shift ();

      if (iter==goal)
	{
	  A_s.shift(-value);
	  ++goal;
	}
      
				       // Update normalized eigenvector
      x.equ (1./length, y);
				       // Compute residual
      if (additional_data.use_residual)
	{
	  y.equ (value, x);
	  double res = A.residual (r,x,y);
					   // Check the residual
	  conv = control().check (iter, res);
	} else {
	  conv = control().check (iter, fabs(1./value-1./old_value));
	}
      old_value = value;
    }

  memory.free(Vy);
  memory.free(Vr);
  
  deallog.pop();

				   // in case of failure: throw
				   // exception
  if (control().last_check() != SolverControl::success)
    throw SolverControl::NoConvergence (control().last_step(),
					control().last_value());
				   // otherwise exit as normal
}

#endif



