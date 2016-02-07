// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2016 by the deal.II authors
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

#ifndef dealii__mg_coarse_h
#define dealii__mg_coarse_h


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/householder.h>
#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mg */
/*@{*/

/**
 * Coarse grid solver using LAC iterative methods. This is a little wrapper,
 * transforming a triplet of iterative solver, matrix and preconditioner into
 * a coarse grid solver.
 *
 * The type of the matrix (i.e. the template parameter @p MatrixType) should
 * be derived from @p Subscriptor to allow for the use of a smart pointer to
 * it.
 *
 * @author Guido Kanschat, 1999, Ralf Hartmann, 2002.
 */
template<typename SolverType, class VectorType = Vector<double> >
class MGCoarseGridLACIteration :  public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Default constructor.
   */
  MGCoarseGridLACIteration ();

  /**
   * Constructor. Store solver, matrix and preconditioning method for later
   * use.
   */
  template<typename MatrixType, typename PreconditionerType>
  MGCoarseGridLACIteration (SolverType &,
                            const MatrixType &,
                            const PreconditionerType &);

  /**
   * Destructor freeing the pointers.
   */
  ~MGCoarseGridLACIteration ();

  /**
   * Initialize new data.
   */
  template<typename MatrixType, typename PreconditionerType>
  void initialize (SolverType &,
                   const MatrixType &,
                   const PreconditionerType &);

  /**
   * Clear all pointers.
   */
  void clear ();

  /**
   * Implementation of the abstract function. Calls the solver method with
   * matrix, vectors and preconditioner.
   */
  void operator() (const unsigned int level,
                   VectorType         &dst,
                   const VectorType   &src) const;

  /**
   * Sets the matrix. This gives the possibility to replace the matrix that
   * was given to the constructor by a new matrix.
   */
  template <typename MatrixType>
  void set_matrix (const MatrixType &);

private:
  /**
   * Reference to the solver.
   */
  SmartPointer<SolverType,MGCoarseGridLACIteration<SolverType,VectorType> > solver;

  /**
   * Reference to the matrix.
   */
  PointerMatrixBase<VectorType> *matrix;

  /**
   * Reference to the preconditioner.
   */
  PointerMatrixBase<VectorType> *precondition;
};



/**
 * Coarse grid solver by QR factorization implemented in the class
 * Householder.
 *
 * Upon initialization, the QR decomposition of the matrix is computed. then,
 * the operator() uses Householder::least_squares() to compute the action of
 * the inverse.
 *
 * @author Guido Kanschat, 2003, 2012
 */
template<typename number = double, class VectorType = Vector<number> >
class MGCoarseGridHouseholder : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Constructor, taking the coarse grid matrix.
   */
  MGCoarseGridHouseholder (const FullMatrix<number> *A = 0);

  /**
   * Initialize for a new matrix.
   */
  void initialize (const FullMatrix<number> &A);

  void operator() (const unsigned int level,
                   VectorType         &dst,
                   const VectorType   &src) const;

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
template<typename number = double, class VectorType = Vector<number> >
class MGCoarseGridSVD : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Constructor leaving an uninitialized object.
   */
  MGCoarseGridSVD ();

  /**
   * Initialize for a new matrix. This resets the dimensions to the
   */
  void initialize (const FullMatrix<number> &A, const double threshold = 0);

  void operator() (const unsigned int level,
                   VectorType         &dst,
                   const VectorType   &src) const;

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


template<typename SolverType, class VectorType>
MGCoarseGridLACIteration<SolverType, VectorType>
::MGCoarseGridLACIteration()
  :
  solver(0, typeid(*this).name()),
  matrix(0),
  precondition(0)
{}


template<typename SolverType, class VectorType>
template<typename MatrixType, typename PreconditionerType>
MGCoarseGridLACIteration<SolverType, VectorType>
::MGCoarseGridLACIteration (SolverType               &s,
                            const MatrixType         &m,
                            const PreconditionerType &p)
  :
  solver(&s, typeid(*this).name())
{
  matrix = new PointerMatrix<MatrixType, VectorType>(&m);
  precondition = new PointerMatrix<PreconditionerType, VectorType>(&p);
}


template<typename SolverType, class VectorType>
MGCoarseGridLACIteration<SolverType, VectorType>
::~MGCoarseGridLACIteration()
{
  clear();
}


template<typename SolverType, class VectorType>
template<typename MatrixType, typename PreconditionerType>
void
MGCoarseGridLACIteration<SolverType, VectorType>
::initialize (SolverType               &s,
              const MatrixType         &m,
              const PreconditionerType &p)
{
  solver = &s;
  if (matrix)
    delete matrix;
  matrix = new PointerMatrix<MatrixType, VectorType>(&m);
  if (precondition)
    delete precondition;
  precondition = new PointerMatrix<PreconditionerType, VectorType>(&p);
}


template<typename SolverType, class VectorType>
void
MGCoarseGridLACIteration<SolverType, VectorType>
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


template<typename SolverType, class VectorType>
void
MGCoarseGridLACIteration<SolverType, VectorType>
::operator() (const unsigned int /* level */,
              VectorType         &dst,
              const VectorType   &src) const
{
  Assert(solver!=0, ExcNotInitialized());
  Assert(matrix!=0, ExcNotInitialized());
  Assert(precondition!=0, ExcNotInitialized());
  solver->solve(*matrix, dst, src, *precondition);
}


template<typename SolverType, class VectorType>
template<typename MatrixType>
void
MGCoarseGridLACIteration<SolverType, VectorType>
::set_matrix(const MatrixType &m)
{
  if (matrix)
    delete matrix;
  matrix = new PointerMatrix<MatrixType, VectorType>(&m);
}

//---------------------------------------------------------------------------

template<typename number, class VectorType>
MGCoarseGridHouseholder<number, VectorType>::MGCoarseGridHouseholder
(const FullMatrix<number> *A)
{
  if (A != 0) householder.initialize(*A);
}



template<typename number, class VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::initialize(const FullMatrix<number> &A)
{
  householder.initialize(A);
}



template<typename number, class VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::operator() (const unsigned int /*level*/,
                                                         VectorType         &dst,
                                                         const VectorType   &src) const
{
  householder.least_squares(dst, src);
}

//---------------------------------------------------------------------------

template<typename number, class VectorType>
inline
MGCoarseGridSVD<number, VectorType>::MGCoarseGridSVD()
{}



template<typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::initialize (const FullMatrix<number> &A,
                                                 double                    threshold)
{
  matrix.reinit(A.n_rows(), A.n_cols());
  matrix = A;
  matrix.compute_inverse_svd(threshold);
}


template<typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::operator() (
  const unsigned int /*level*/,
  VectorType         &dst,
  const VectorType   &src) const
{
  matrix.vmult(dst, src);
}


template<typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::log() const
{
  const unsigned int n = std::min(matrix.n_rows(), matrix.n_cols());

  for (unsigned int i=0; i<n; ++i)
    deallog << ' ' << matrix.singular_value(i);
  deallog << std::endl;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
