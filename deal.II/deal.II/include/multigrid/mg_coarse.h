//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_coarse_h
#define __deal2__mg_coarse_h


#include <lac/full_matrix.h>
#include <lac/matrix_lib.h>
#include <multigrid/mg_base.h>

/*!@addtogroup mg */
/*@{*/

/**
 * Coarse grid solver using LAC iterative methods.
 * This is a little wrapper, transforming a triplet of iterative
 * solver, matrix and preconditioner into a coarse grid solver.
 *
 * The type of the matrix (i.e. the template parameter @p MATRIX)
 * should be derived from @p Subscriptor to allow for the use of a
 * smart pointer to it.
 *
 * @author Guido Kanschat, 1999, Ralf Hartmann, 2002.
 */
template<class SOLVER, class VECTOR = Vector<double> >
class MGCoarseGridLACIteration :  public MGCoarseGridBase<VECTOR>
{
  public:
				     /**
				      * Default constructor.
				      */
    MGCoarseGridLACIteration ();
    
				     /**
				      * Constructor.
				      * Store solver, matrix and
				      * preconditioning method for later
				      * use.
				      */
    template<class MATRIX, class PRECOND>
    MGCoarseGridLACIteration (SOLVER        &,
			      const MATRIX  &,
			      const PRECOND &);
    
				     /**
				      * Destructor freeing the pointers.
				      */
    ~MGCoarseGridLACIteration ();
    
				     /**
				      * Initialize new data.
				      */
    template<class MATRIX, class PRECOND>
    void initialize (SOLVER        &,
		     const MATRIX  &,
		     const PRECOND &);

				     /**
				      * Clear all pointers.
				      */
    void clear ();
    
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

				     /**
				      * Sets the matrix. This gives
				      * the possibility to replace the
				      * matrix that was given to the
				      * constructor by a new matrix.
				      */
    template <class MATRIX>
    void set_matrix (const MATRIX &);
    
  private:
				     /**
				      * Reference to the solver.
				      */
    SmartPointer<SOLVER> solver;
    
				     /**
				      * Reference to the matrix.
				      */
    PointerMatrixBase<VECTOR>* matrix;
    
				     /**
				      * Reference to the preconditioner.
				      */
    PointerMatrixBase<VECTOR>* precondition;
};



//TODO: Improve QR-factorization to make this more efficient.

/**
 * Coarse grid solver by QR factorization
 *
 * This is a direct solver for possibly indefinite coarse grid
 * problems, using the @p least_squares function of the class
 * FullMatrix.
 *
 * Since the currently implemented Householder algorithm transforms
 * the right hand side immediately, the transformation has to be
 * computed for each coarse grid solution. Therefore, this coarse grid
 * solver may be inefficient for larger coarse grid systems.
 *
 * @author Guido Kanschat, 2003
 */
template<typename number, class VECTOR>
class MGCoarseGridHouseholder : public MGCoarseGridBase<VECTOR>
{
  public:
				     /**
				      * Constructor, taking the coarse
				      * grid matrix.
				      */
    MGCoarseGridHouseholder (const FullMatrix<number>* A = 0);

				     /**
				      * Initialize for a new matrix.
				      */
    void initialize (const FullMatrix<number>& A);

				     /**
				      * Release matrix pointer.
				      */
    void clear ();
    
				     /**
				      * Solution operator, defined in
				      * the base class.
				      */
    void operator() (const unsigned int   level,
		     VECTOR       &dst,
		     const VECTOR &src) const;

  private:
				     /**
				      * Pointer to the coarse grid
				      * matrix.
				      */
    SmartPointer<const FullMatrix<number> > matrix;

				     /**
				      * Matrix for QR-factorization.
				      */
    mutable FullMatrix<number> work;

				     /**
				      * Auxiliary vector.
				      */
    mutable VECTOR aux;
};

/*@}*/
  
/// @if NoDoc
/* ------------------ Functions for MGCoarseGridLACIteration ------------ */


template<class SOLVER, class VECTOR>
MGCoarseGridLACIteration<SOLVER, VECTOR>
::MGCoarseGridLACIteration()
		:
		solver(0, typeid(*this).name()),
		matrix(0),
		precondition(0)
{}


template<class SOLVER, class VECTOR>
template<class MATRIX, class PRECOND>
MGCoarseGridLACIteration<SOLVER, VECTOR>
::MGCoarseGridLACIteration(SOLVER& s,
			   const MATRIX  &m,
			   const PRECOND &p)
		:
		solver(&s, typeid(*this).name())
{
  matrix = new PointerMatrix<MATRIX, VECTOR>(&m);
  precondition = new PointerMatrix<PRECOND, VECTOR>(&p);
}


template<class SOLVER, class VECTOR>
MGCoarseGridLACIteration<SOLVER, VECTOR>
::~MGCoarseGridLACIteration()
{
  clear();
}


template<class SOLVER, class VECTOR>
template<class MATRIX, class PRECOND>
void
MGCoarseGridLACIteration<SOLVER, VECTOR>
::initialize(SOLVER& s,
	     const MATRIX  &m,
	     const PRECOND &p)
{
  solver = &s;
  if (matrix)
    delete matrix;
  matrix = new PointerMatrix<MATRIX, VECTOR>(&m);
  if (precondition)
    delete precondition;
  precondition = new PointerMatrix<PRECOND, VECTOR>(&p);
}


template<class SOLVER, class VECTOR>
void
MGCoarseGridLACIteration<SOLVER, VECTOR>
::clear()
{
  solver = 0;
  if (matrix)
    delete matrix;
  matrix = 0;
  if (precondition)
    delete precondition;
  precondition = 0;
}


template<class SOLVER, class VECTOR>
void
MGCoarseGridLACIteration<SOLVER, VECTOR>
::operator() (const unsigned int    /* level */,
	      VECTOR       &dst,
	      const VECTOR &src) const
{
  Assert(solver!=0, ExcNotInitialized());
  Assert(matrix!=0, ExcNotInitialized());
  Assert(precondition!=0, ExcNotInitialized());
  solver->solve(*matrix, dst, src, *precondition);
}


template<class SOLVER, class VECTOR>
template<class MATRIX>
void
MGCoarseGridLACIteration<SOLVER, VECTOR>
::set_matrix(const MATRIX &m)
{
  if (matrix)
    delete matrix;
  matrix = new PointerMatrix<MATRIX, VECTOR>(&m);
}

//---------------------------------------------------------------------------

template<typename number, class VECTOR>
MGCoarseGridHouseholder<number, VECTOR>::MGCoarseGridHouseholder(
  const FullMatrix<number>* A)
		:
		matrix(A, typeid(*this).name())
{}



template<typename number, class VECTOR>
void
MGCoarseGridHouseholder<number, VECTOR>::initialize(
  const FullMatrix<number>& A)
{
  matrix = &A;
}



template<typename number, class VECTOR>
void
MGCoarseGridHouseholder<number, VECTOR>::clear()
{
  matrix = 0;
}



template<typename number, class VECTOR>
void
MGCoarseGridHouseholder<number, VECTOR>::operator() (
  const unsigned int,
  VECTOR       &dst,
  const VECTOR &src) const
{
  work = *matrix;
  aux = src;
  work.least_squares(dst, aux);
}

/// @endif

#endif
