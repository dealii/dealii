//----------------------------  eigen.h  ---------------------------
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
//----------------------------  eigen.h  ---------------------------
#ifndef __deal2__eigen_h
#define __deal2__eigen_h


#include <lac/forward_declarations.h>
#include <lac/solver.h>
#include <lac/solver_control.h>
#include <lac/vector_memory.h>

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
template <class MATRIX = SparseMatrix<double>,
          class VECTOR = Vector<double> >
class EigenPower : public Solver<MATRIX,VECTOR>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver. This solver does not
				      * need additional data yet.
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
	AdditionalData (const double shift):
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
				      * Power method. @p x is the (not
				      * necessarily normalized) start
				      * vector for the power
				      * method. After the iteration,
				      * @p value is the approximated
				      * eigenvalue and @p x is the
				      * corresponding eigenvector,
				      * normalized with respect to the l2-norm.
				      */
    typename Solver<MATRIX,VECTOR>::ReturnState
    solve (double       &value,
	   const MATRIX &A,
	   VECTOR       &x);

  protected:
				     /**
				      * Shift parameter.
				      */
    AdditionalData additional_data;
};

//----------------------------------------------------------------------//

template <class MATRIX, class VECTOR>
EigenPower<MATRIX, VECTOR>::EigenPower (SolverControl &cn,
					VectorMemory<VECTOR> &mem,
					const AdditionalData &data):
		Solver<MATRIX, VECTOR>(cn, mem),
		additional_data(data)
{}


template <class MATRIX, class VECTOR>
EigenPower<MATRIX, VECTOR>::~EigenPower ()
{}


template <class MATRIX, class VECTOR>
typename Solver<MATRIX,VECTOR>::ReturnState
EigenPower<MATRIX, VECTOR>::solve (double       &value,
				   const MATRIX &A,
				   VECTOR       &x)
{
  SolverControl::State conv=SolverControl::iterate;

  deallog.push("Power method");

  VECTOR* Vy = memory.alloc (); VECTOR& y = *Vy; y.reinit (x);
  
  double length = x.l2_norm ();
  double old_length = 0.;
  x.scale(1./length);
  
  
				   // Main loop
  for(int iter=0; conv==SolverControl::iterate; iter++)
    {
      A.vmult (y,x);
      y.add(additional_data.shift, x);
      
				       // Compute absolute value of eigenvalue
      old_length = length;
      length = y.l2_norm ();

				       // do a little trick to compute the sign
				       // with not too much round-off errors.
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

      conv = control().check (iter, fabs(1.-length/old_length));
    }
  
  memory.free(Vy);

  deallog.pop();
				   // Output
  if (conv == SolverControl::failure)
    return exceeded;
  else
    return success;

}


#endif
