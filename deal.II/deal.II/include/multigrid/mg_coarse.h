//------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------------------------------------------------------
#ifndef __deal2__mg_coarse_h
#define __deal2__mg_coarse_h


#include <multigrid/mg_base.h>


/**
 * Coarse grid solver using LAC iterative methods.
 * This is a little wrapper, transforming a triplet of iterative
 * solver, matrix and preconditioner into a coarse grid solver.
 *
 * The type of the matrix (i.e. the template parameter @p{MATRIX})
 * should be derived from @p{Subscriptor} to allow for the use of a
 * smart pointer to it.
 *
 * @author Guido Kanschat, 1999
 */
template<class SOLVER, class MATRIX, class PRECOND, class VECTOR = Vector<double> >
class MGCoarseGridLACIteration :  public MGCoarseGrid<VECTOR>
{
  public:
				     /**
				      * Constructor.
				      * Store solver, matrix and
				      * preconditioning method for later
				      * use.
				      */
    MGCoarseGridLACIteration (SOLVER        &,
			      const MATRIX  &,
			      const PRECOND &);
    
				     /**
				      * Implementation of the abstract
				      * function.
				      * Calls the solver method with
				      * matrix, vectors and
				      * preconditioner.
				      */
    void operator() (const unsigned int   level,
		     VECTOR       &dst,
		     const VECTOR &src) const;
  private:
				     /**
				      * Reference to the solver.
				      */
    SOLVER& solver;
    
				     /**
				      * Reference to the matrix.
				      */
    const SmartPointer<const MATRIX> matrix;
    
				     /**
				      * Reference to the preconditioner.
				      */
    const PRECOND& precondition;
};


/* ------------------ Functions for MGCoarseGridLACIteration ------------ */


template<class SOLVER, class MATRIX, class PRECOND, class VECTOR>
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND, VECTOR>
::MGCoarseGridLACIteration(SOLVER& s,
			   const MATRIX  &m,
			   const PRECOND &p)
		:
		solver(s),
		matrix(&m),
		precondition(p)
{};


template<class SOLVER, class MATRIX, class PRECOND, class VECTOR>
void
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND, VECTOR>
::operator() (const unsigned int    /* level */,
	      VECTOR       &dst,
	      const VECTOR &src) const
{
  solver.solve(*matrix, dst, src, precondition);
}



#endif
