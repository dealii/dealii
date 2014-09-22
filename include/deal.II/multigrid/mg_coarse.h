// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__mg_coarse_h
#define __deal2__mg_coarse_h


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/householder.h>
#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

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
  MGCoarseGridLACIteration (SOLVER &,
                            const MATRIX &,
                            const PRECOND &);

  /**
   * Destructor freeing the pointers.
   */
  ~MGCoarseGridLACIteration ();

  /**
   * Initialize new data.
   */
  template<class MATRIX, class PRECOND>
  void initialize (SOLVER &,
                   const MATRIX &,
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
  SmartPointer<SOLVER,MGCoarseGridLACIteration<SOLVER,VECTOR> > solver;

  /**
   * Reference to the matrix.
   */
  PointerMatrixBase<VECTOR> *matrix;

  /**
   * Reference to the preconditioner.
   */
  PointerMatrixBase<VECTOR> *precondition;
};



/**
 * Coarse grid solver by QR factorization implemented in the class Householder.
 *
 * Upon initialization, the QR decomposition of the matrix is
 * computed. then, the operator() uses Householder::least_squares() to
 * compute the action of the inverse.
 *
 * @author Guido Kanschat, 2003, 2012
 */
template<typename number = double, class VECTOR = Vector<number> >
class MGCoarseGridHouseholder : public MGCoarseGridBase<VECTOR>
{
public:
  /**
   * Constructor, taking the coarse
   * grid matrix.
   */
  MGCoarseGridHouseholder (const FullMatrix<number> *A = 0);

  /**
   * Initialize for a new matrix.
   */
  void initialize (const FullMatrix<number> &A);

  void operator() (const unsigned int   level,
                   VECTOR       &dst,
                   const VECTOR &src) const;

private:
  /**
   * Matrix for QR-factorization.
   */
  Householder<number> householder;
};

/**
 * Coarse grid solver using singular value decomposition of LAPACK matrices.
 *
 * Upon initialization, the singular value decomposition of the matrix is
 * computed. then, the operator() uses
 *
 * @author Guido Kanschat, 2003, 2012
 */
template<typename number = double, class VECTOR = Vector<number> >
class MGCoarseGridSVD : public MGCoarseGridBase<VECTOR>
{
public:
  /**
   * Constructor leaving an
   * uninitialized object.
   */
  MGCoarseGridSVD ();

  /**
   * Initialize for a new
   * matrix. This resets the
   * dimensions to the
   */
  void initialize (const FullMatrix<number> &A, const double threshold = 0);

  void operator() (const unsigned int   level,
                   VECTOR       &dst,
                   const VECTOR &src) const;

  /**
   * Write the singular values to @p deallog.
   */
  void log () const;

private:

  /**
   * Matrix for singular value decomposition.
   */
  LAPACKFullMatrix<number> matrix;
};

/*@}*/

#ifndef DOXYGEN
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
::MGCoarseGridLACIteration(SOLVER &s,
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
::initialize(SOLVER &s,
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
  const FullMatrix<number> *A)
{
  if (A != 0) householder.initialize(*A);
}



template<typename number, class VECTOR>
void
MGCoarseGridHouseholder<number, VECTOR>::initialize(
  const FullMatrix<number> &A)
{
  householder.initialize(A);
}



template<typename number, class VECTOR>
void
MGCoarseGridHouseholder<number, VECTOR>::operator() (
  const unsigned int /*level*/,
  VECTOR       &dst,
  const VECTOR &src) const
{
  householder.least_squares(dst, src);
}

//---------------------------------------------------------------------------

template<typename number, class VECTOR>
inline
MGCoarseGridSVD<number, VECTOR>::MGCoarseGridSVD()
{}



template<typename number, class VECTOR>
void
MGCoarseGridSVD<number, VECTOR>::initialize(
  const FullMatrix<number> &A,
  double threshold)
{
  matrix.reinit(A.n_rows(), A.n_cols());
  matrix = A;
  matrix.compute_inverse_svd(threshold);
}


template<typename number, class VECTOR>
void
MGCoarseGridSVD<number, VECTOR>::operator() (
  const unsigned int /*level*/,
  VECTOR       &dst,
  const VECTOR &src) const
{
  matrix.vmult(dst, src);
}


template<typename number, class VECTOR>
void
MGCoarseGridSVD<number, VECTOR>::log() const
{
  const unsigned int n = std::min(matrix.n_rows(), matrix.n_cols());

  for (unsigned int i=0; i<n; ++i)
    deallog << ' ' << matrix.singular_value(i);
  deallog << std::endl;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
