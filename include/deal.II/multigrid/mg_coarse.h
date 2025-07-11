// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_coarse_h
#define dealii_mg_coarse_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/linear_operator.h>

#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup mg
 * @{
 */

/**
 * Coarse grid solver using smoother only. This is a little wrapper,
 * transforming a smoother into a coarse grid solver.
 */
template <typename VectorType = Vector<double>>
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
             VectorType        &dst,
             const VectorType  &src) const override;

private:
  /**
   * Reference to the smoother.
   */
  ObserverPointer<const MGSmootherBase<VectorType>,
                  MGCoarseGridApplySmoother<VectorType>>
    coarse_smooth;
};



/**
 * Coarse grid solver using an operator, e.g., preconditioner, only.
 * This is a little wrapper, transforming an operator into a coarse grid solver.
 */
template <class VectorType, class MatrixType>
class MGCoarseGridApplyOperator : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * Default constructor.
   */
  MGCoarseGridApplyOperator();

  /**
   * Constructor. Store a pointer to the operator for later use.
   */
  MGCoarseGridApplyOperator(const MatrixType &matrix);

  /**
   * Clear the pointer.
   */
  void
  clear();

  /**
   * Initialize new data.
   */
  void
  initialize(const MatrixType &matrix);

  /**
   * Implementation of the abstract function.
   */
  virtual void
  operator()(const unsigned int level,
             VectorType        &dst,
             const VectorType  &src) const override;

private:
  /**
   * Reference to the operator.
   */
  ObserverPointer<const MatrixType,
                  MGCoarseGridApplyOperator<VectorType, MatrixType>>
    matrix;
};



/**
 * Coarse grid multigrid operator for an iterative solver.
 *
 * This class provides a wrapper for a deal.II iterative solver with a given
 * matrix and preconditioner as a coarse grid operator.
 */
template <typename VectorType,
          typename SolverType,
          typename MatrixType,
          typename PreconditionerType>
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
  MGCoarseGridIterativeSolver(SolverType               &solver,
                              const MatrixType         &matrix,
                              const PreconditionerType &precondition);

  /**
   * Initialize with new data, see the corresponding constructor for more
   * details.
   */
  void
  initialize(SolverType               &solver,
             const MatrixType         &matrix,
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
             VectorType        &dst,
             const VectorType  &src) const override;

private:
  /**
   * Reference to the solver.
   */
  ObserverPointer<SolverType,
                  MGCoarseGridIterativeSolver<VectorType,
                                              SolverType,
                                              MatrixType,
                                              PreconditionerType>>
    solver;

  /**
   * Reference to the matrix.
   */
  ObserverPointer<const MatrixType,
                  MGCoarseGridIterativeSolver<VectorType,
                                              SolverType,
                                              MatrixType,
                                              PreconditionerType>>
    matrix;

  /**
   * Reference to the preconditioner.
   */
  ObserverPointer<const PreconditionerType,
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
 */
template <typename number = double, typename VectorType = Vector<number>>
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
             VectorType        &dst,
             const VectorType  &src) const override;

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
 */
template <typename number = double, typename VectorType = Vector<number>>
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
             VectorType        &dst,
             const VectorType  &src) const;

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

/** @} */

#ifndef DOXYGEN
/* ------------------ Functions for MGCoarseGridApplySmoother -----------*/
template <typename VectorType>
MGCoarseGridApplySmoother<VectorType>::MGCoarseGridApplySmoother()
  : coarse_smooth(nullptr)
{}

template <typename VectorType>
MGCoarseGridApplySmoother<VectorType>::MGCoarseGridApplySmoother(
  const MGSmootherBase<VectorType> &coarse_smooth)
  : coarse_smooth(nullptr)
{
  initialize(coarse_smooth);
}


template <typename VectorType>
void
MGCoarseGridApplySmoother<VectorType>::initialize(
  const MGSmootherBase<VectorType> &coarse_smooth_)
{
  coarse_smooth = ObserverPointer<const MGSmootherBase<VectorType>,
                                  MGCoarseGridApplySmoother<VectorType>>(
    &coarse_smooth_, typeid(*this).name());
}


template <typename VectorType>
void
MGCoarseGridApplySmoother<VectorType>::clear()
{
  coarse_smooth = nullptr;
}


template <typename VectorType>
void
MGCoarseGridApplySmoother<VectorType>::operator()(const unsigned int level,
                                                  VectorType        &dst,
                                                  const VectorType  &src) const
{
  coarse_smooth->apply(level, dst, src);
}

/* ------------------ Functions for MGCoarseGridApplyOperator ----------*/

template <class VectorType, class PreconditionerType>
MGCoarseGridApplyOperator<VectorType,
                          PreconditionerType>::MGCoarseGridApplyOperator()
  : matrix(0, typeid(*this).name())
{}



template <class VectorType, class PreconditionerType>
MGCoarseGridApplyOperator<VectorType, PreconditionerType>::
  MGCoarseGridApplyOperator(const PreconditionerType &matrix)
  : matrix(&matrix, typeid(*this).name())
{}



template <class VectorType, class PreconditionerType>
void
MGCoarseGridApplyOperator<VectorType, PreconditionerType>::initialize(
  const PreconditionerType &matrix_)
{
  matrix = &matrix_;
}



template <class VectorType, class PreconditionerType>
void
MGCoarseGridApplyOperator<VectorType, PreconditionerType>::clear()
{
  matrix = 0;
}


template <class VectorType, class PreconditionerType>
void
MGCoarseGridApplyOperator<VectorType, PreconditionerType>::operator()(
  const unsigned int /*level*/,
  VectorType       &dst,
  const VectorType &src) const
{
  Assert(matrix != nullptr, ExcNotInitialized());

  matrix->vmult(dst, src);
}

/* ------------------ Functions for MGCoarseGridIterativeSolver ------------ */

template <typename VectorType,
          typename SolverType,
          typename MatrixType,
          typename PreconditionerType>
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::MGCoarseGridIterativeSolver()
  : solver(0, typeid(*this).name())
  , matrix(0, typeid(*this).name())
  , preconditioner(0, typeid(*this).name())
{}



template <typename VectorType,
          typename SolverType,
          typename MatrixType,
          typename PreconditionerType>
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::
  MGCoarseGridIterativeSolver(SolverType               &solver,
                              const MatrixType         &matrix,
                              const PreconditionerType &preconditioner)
  : solver(&solver, typeid(*this).name())
  , matrix(&matrix, typeid(*this).name())
  , preconditioner(&preconditioner, typeid(*this).name())
{}



template <typename VectorType,
          typename SolverType,
          typename MatrixType,
          typename PreconditionerType>
void
MGCoarseGridIterativeSolver<
  VectorType,
  SolverType,
  MatrixType,
  PreconditionerType>::initialize(SolverType               &solver_,
                                  const MatrixType         &matrix_,
                                  const PreconditionerType &preconditioner_)
{
  solver         = &solver_;
  matrix         = &matrix_;
  preconditioner = &preconditioner_;
}



template <typename VectorType,
          typename SolverType,
          typename MatrixType,
          typename PreconditionerType>
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



template <typename VectorType,
          typename SolverType,
          typename MatrixType,
          typename PreconditionerType>
void
MGCoarseGridIterativeSolver<
  VectorType,
  SolverType,
  MatrixType,
  PreconditionerType>::operator()(const unsigned int /*level*/,
                                  VectorType       &dst,
                                  const VectorType &src) const
{
  Assert(solver != nullptr, ExcNotInitialized());
  Assert(matrix != nullptr, ExcNotInitialized());
  Assert(preconditioner != nullptr, ExcNotInitialized());

  dst = 0;

  if constexpr (std::is_same_v<VectorType, typename SolverType::vector_type>)
    {
      solver->solve(*matrix, dst, src, *preconditioner);
    }
  else
    {
      typename SolverType::vector_type src_;
      typename SolverType::vector_type dst_;

      src_ = src;
      dst_ = dst;

      solver->solve(*matrix, dst_, src_, *preconditioner);

      dst = dst_;
    }
}



/* ------------------ Functions for MGCoarseGridHouseholder ------------ */

template <typename number, typename VectorType>
MGCoarseGridHouseholder<number, VectorType>::MGCoarseGridHouseholder(
  const FullMatrix<number> *A)
{
  if (A != nullptr)
    householder.initialize(*A);
}



template <typename number, typename VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::initialize(
  const FullMatrix<number> &A)
{
  householder.initialize(A);
}



template <typename number, typename VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::operator()(
  const unsigned int /*level*/,
  VectorType       &dst,
  const VectorType &src) const
{
  householder.least_squares(dst, src);
}

//---------------------------------------------------------------------------



template <typename number, typename VectorType>
void
MGCoarseGridSVD<number, VectorType>::initialize(const FullMatrix<number> &A,
                                                double threshold)
{
  matrix.reinit(A.n_rows(), A.n_cols());
  matrix = A;
  matrix.compute_inverse_svd(threshold);
}


template <typename number, typename VectorType>
void
MGCoarseGridSVD<number, VectorType>::operator()(const unsigned int /*level*/,
                                                VectorType       &dst,
                                                const VectorType &src) const
{
  matrix.vmult(dst, src);
}


template <typename number, typename VectorType>
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
