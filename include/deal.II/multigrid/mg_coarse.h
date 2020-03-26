// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_mg_coarse_h
#define dealii_mg_coarse_h


#include <deal.II/base/config.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/linear_operator.h>

#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mg */
/*@{*/

/**
 * Coarse grid solver using smoother only. This is a little wrapper,
 * transforming a smoother into a coarse grid solver.
 *
 * @author Denis Davydov, 2016.
 */
template <class VectorType = Vector<double>>
class MGCoarseGridApplySmoother : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Default constructor.
   */
  MGCoarseGridApplySmoother();

  /**
   * Constructor. Store a pointer to the smoother for later use.
   */
  MGCoarseGridApplySmoother(const MGSmootherBase<VectorType> &coarse_smooth);

  /**
   * Clear the pointer.
   */
  void
  clear();

  /**
   * Initialize new data.
   */
  void
  initialize(const MGSmootherBase<VectorType> &coarse_smooth);

  /**
   * Implementation of the abstract function.
   */
  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const override;

private:
  /**
   * Reference to the smoother.
   */
  SmartPointer<const MGSmootherBase<VectorType>,
               MGCoarseGridApplySmoother<VectorType>>
    coarse_smooth;
};


/**
 * Coarse grid solver using LAC iterative methods. This is a little wrapper,
 * transforming a triplet of iterative solver, matrix and preconditioner into
 * a coarse grid solver.
 *
 * The type of the matrix (i.e. the template parameter @p MatrixType) should
 * be derived from @p Subscriptor to allow for the use of a smart pointer to
 * it.
 *
 * @deprecated Use MGCoarseGridIterativeSolver instead.
 *
 * @author Guido Kanschat, 1999, Ralf Hartmann, 2002.
 */
template <typename SolverType, class VectorType = Vector<double>>
class DEAL_II_DEPRECATED MGCoarseGridLACIteration
  : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Default constructor.
   */
  MGCoarseGridLACIteration();

  /**
   * Constructor. Store solver, matrix and preconditioning method for later
   * use.
   */
  template <typename MatrixType, typename PreconditionerType>
  MGCoarseGridLACIteration(SolverType &,
                           const MatrixType &,
                           const PreconditionerType &);

  /**
   * Destructor freeing the pointers.
   */
  ~MGCoarseGridLACIteration();

  /**
   * Initialize new data.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  initialize(SolverType &, const MatrixType &, const PreconditionerType &);

  /**
   * Clear all pointers.
   */
  void
  clear();

  /**
   * Implementation of the abstract function. Calls the solver method with
   * matrix, vectors and preconditioner.
   */
  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const;

  /**
   * Set the matrix. This gives the possibility to replace the matrix that
   * was given to the constructor by a new matrix.
   */
  template <typename MatrixType>
  void
  set_matrix(const MatrixType &);

private:
  /**
   * Reference to the solver.
   */
  SmartPointer<SolverType, MGCoarseGridLACIteration<SolverType, VectorType>>
    solver;

  /**
   * LinearOperator wrapping a reference to the matrix.
   */
  LinearOperator<VectorType> matrix;

  /**
   * LinearOperator wrapping a reference to the preconditioner.
   */
  LinearOperator<VectorType> precondition;
};



/**
 * Coarse grid multigrid operator for an iterative solver.
 *
 * This class provides a wrapper for a deal.II iterative solver with a given
 * matrix and preconditioner as a coarse grid operator.
 */
template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
class MGCoarseGridIterativeSolver : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Default constructor.
   */
  MGCoarseGridIterativeSolver();

  /**
   * Constructor. Only a reference to these objects is stored, so
   * their lifetime needs to exceed the usage in this class.
   */
  MGCoarseGridIterativeSolver(SolverType &              solver,
                              const MatrixType &        matrix,
                              const PreconditionerType &precondition);

  /**
   * Initialize with new data, see the corresponding constructor for more
   * details.
   */
  void
  initialize(SolverType &              solver,
             const MatrixType &        matrix,
             const PreconditionerType &precondition);

  /**
   * Clear all pointers.
   */
  void
  clear();

  /**
   * Implementation of the abstract function. Calls the solve method of the
   * given solver with matrix, vectors, and preconditioner.
   */
  virtual void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const override;

private:
  /**
   * Reference to the solver.
   */
  SmartPointer<SolverType,
               MGCoarseGridIterativeSolver<VectorType,
                                           SolverType,
                                           MatrixType,
                                           PreconditionerType>>
    solver;

  /**
   * Reference to the matrix.
   */
  SmartPointer<const MatrixType,
               MGCoarseGridIterativeSolver<VectorType,
                                           SolverType,
                                           MatrixType,
                                           PreconditionerType>>
    matrix;

  /**
   * Reference to the preconditioner.
   */
  SmartPointer<const PreconditionerType,
               MGCoarseGridIterativeSolver<VectorType,
                                           SolverType,
                                           MatrixType,
                                           PreconditionerType>>
    preconditioner;
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
template <typename number = double, class VectorType = Vector<number>>
class MGCoarseGridHouseholder : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Constructor, taking the coarse grid matrix.
   */
  MGCoarseGridHouseholder(const FullMatrix<number> *A = nullptr);

  /**
   * Initialize for a new matrix.
   */
  void
  initialize(const FullMatrix<number> &A);

  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const override;

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
template <typename number = double, class VectorType = Vector<number>>
class MGCoarseGridSVD : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Constructor leaving an uninitialized object.
   */
  MGCoarseGridSVD() = default;

  /**
   * Initialize for a new matrix. This resets the dimensions to the
   */
  void
  initialize(const FullMatrix<number> &A, const double threshold = 0);

  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const;

  /**
   * Write the singular values to @p deallog.
   */
  void
  log() const;

private:
  /**
   * Matrix for singular value decomposition.
   */
  LAPACKFullMatrix<number> matrix;
};

/*@}*/

#ifndef DOXYGEN
/* ------------------ Functions for MGCoarseGridApplySmoother -----------*/
template <class VectorType>
MGCoarseGridApplySmoother<VectorType>::MGCoarseGridApplySmoother()
  : coarse_smooth(nullptr)
{}

template <class VectorType>
MGCoarseGridApplySmoother<VectorType>::MGCoarseGridApplySmoother(
  const MGSmootherBase<VectorType> &coarse_smooth)
  : coarse_smooth(nullptr)
{
  initialize(coarse_smooth);
}


template <class VectorType>
void
MGCoarseGridApplySmoother<VectorType>::initialize(
  const MGSmootherBase<VectorType> &coarse_smooth_)
{
  coarse_smooth =
    SmartPointer<const MGSmootherBase<VectorType>,
                 MGCoarseGridApplySmoother<VectorType>>(&coarse_smooth_,
                                                        typeid(*this).name());
}


template <class VectorType>
void
MGCoarseGridApplySmoother<VectorType>::clear()
{
  coarse_smooth = nullptr;
}


template <class VectorType>
void
MGCoarseGridApplySmoother<VectorType>::operator()(const unsigned int level,
                                                  VectorType &       dst,
                                                  const VectorType & src) const
{
  coarse_smooth->smooth(level, dst, src);
}

/* ------------------ Functions for MGCoarseGridLACIteration ------------ */


template <typename SolverType, class VectorType>
MGCoarseGridLACIteration<SolverType, VectorType>::MGCoarseGridLACIteration()
  : solver(0, typeid(*this).name())
  , matrix(0)
  , precondition(0)
{}


template <typename SolverType, class VectorType>
template <typename MatrixType, typename PreconditionerType>
MGCoarseGridLACIteration<SolverType, VectorType>::MGCoarseGridLACIteration(
  SolverType &              s,
  const MatrixType &        m,
  const PreconditionerType &p)
  : solver(&s, typeid(*this).name())
{
  // Workaround: Unfortunately, not every "m" object has a rich enough
  // interface to populate reinit_(domain|range)_vector. Thus, supply an
  // empty LinearOperator exemplar.
  matrix       = linear_operator<VectorType>(LinearOperator<VectorType>(), m);
  precondition = linear_operator<VectorType>(matrix, p);
}


template <typename SolverType, class VectorType>
MGCoarseGridLACIteration<SolverType, VectorType>::~MGCoarseGridLACIteration()
{
  clear();
}


template <typename SolverType, class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
MGCoarseGridLACIteration<SolverType, VectorType>::initialize(
  SolverType &              s,
  const MatrixType &        m,
  const PreconditionerType &p)
{
  solver = &s;
  // Workaround: Unfortunately, not every "m" object has a rich enough
  // interface to populate reinit_(domain|range)_vector. Thus, supply an
  // empty LinearOperator exemplar.
  matrix       = linear_operator<VectorType>(LinearOperator<VectorType>(), m);
  precondition = linear_operator<VectorType>(matrix, p);
}


template <typename SolverType, class VectorType>
void
MGCoarseGridLACIteration<SolverType, VectorType>::clear()
{
  solver       = nullptr;
  matrix       = LinearOperator<VectorType>();
  precondition = LinearOperator<VectorType>();
}


template <typename SolverType, class VectorType>
void
MGCoarseGridLACIteration<SolverType, VectorType>::
operator()(const unsigned int /* level */,
           VectorType &      dst,
           const VectorType &src) const
{
  Assert(solver != nullptr, ExcNotInitialized());
  solver->solve(matrix, dst, src, precondition);
}


template <typename SolverType, class VectorType>
template <typename MatrixType>
void
MGCoarseGridLACIteration<SolverType, VectorType>::set_matrix(
  const MatrixType &m)
{
  // Workaround: Unfortunately, not every "m" object has a rich enough
  // interface to populate reinit_(domain|range)_vector. Thus, supply an
  // empty LinearOperator exemplar.
  matrix = linear_operator<VectorType>(LinearOperator<VectorType>(), m);
}



/* ------------------ Functions for MGCoarseGridIterativeSolver ------------ */

template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::MGCoarseGridIterativeSolver()
  : solver(0, typeid(*this).name())
  , matrix(0, typeid(*this).name())
  , preconditioner(0, typeid(*this).name())
{}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::
  MGCoarseGridIterativeSolver(SolverType &              solver,
                              const MatrixType &        matrix,
                              const PreconditionerType &preconditioner)
  : solver(&solver, typeid(*this).name())
  , matrix(&matrix, typeid(*this).name())
  , preconditioner(&preconditioner, typeid(*this).name())
{}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
MGCoarseGridIterativeSolver<
  VectorType,
  SolverType,
  MatrixType,
  PreconditionerType>::initialize(SolverType &              solver_,
                                  const MatrixType &        matrix_,
                                  const PreconditionerType &preconditioner_)
{
  solver         = &solver_;
  matrix         = &matrix_;
  preconditioner = &preconditioner_;
}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::clear()
{
  solver         = 0;
  matrix         = 0;
  preconditioner = 0;
}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::
operator()(const unsigned int /*level*/,
           VectorType &      dst,
           const VectorType &src) const
{
  Assert(solver != nullptr, ExcNotInitialized());
  Assert(matrix != nullptr, ExcNotInitialized());
  Assert(preconditioner != nullptr, ExcNotInitialized());
  solver->solve(*matrix, dst, src, *preconditioner);
}



/* ------------------ Functions for MGCoarseGridHouseholder ------------ */

template <typename number, class VectorType>
MGCoarseGridHouseholder<number, VectorType>::MGCoarseGridHouseholder(
  const FullMatrix<number> *A)
{
  if (A != nullptr)
    householder.initialize(*A);
}



template <typename number, class VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::initialize(
  const FullMatrix<number> &A)
{
  householder.initialize(A);
}



template <typename number, class VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::
operator()(const unsigned int /*level*/,
           VectorType &      dst,
           const VectorType &src) const
{
  householder.least_squares(dst, src);
}

//---------------------------------------------------------------------------



template <typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::initialize(const FullMatrix<number> &A,
                                                double threshold)
{
  matrix.reinit(A.n_rows(), A.n_cols());
  matrix = A;
  matrix.compute_inverse_svd(threshold);
}


template <typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::operator()(const unsigned int /*level*/,
                                                VectorType &      dst,
                                                const VectorType &src) const
{
  matrix.vmult(dst, src);
}


template <typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::log() const
{
  const unsigned int n = std::min(matrix.n_rows(), matrix.n_cols());

  for (unsigned int i = 0; i < n; ++i)
    deallog << ' ' << matrix.singular_value(i);
  deallog << std::endl;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
