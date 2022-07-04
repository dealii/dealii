// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2021 by the deal.II authors
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
  coarse_smooth->apply(level, dst, src);
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



namespace internal
{
  namespace MGCoarseGridIterativeSolver
  {
    template <
      class VectorType,
      class SolverType,
      class MatrixType,
      class PreconditionerType,
      typename std::enable_if<
        std::is_same<VectorType, typename SolverType::vector_type>::value,
        VectorType>::type * = nullptr>
    void
    solve(SolverType &              solver,
          const MatrixType &        matrix,
          const PreconditionerType &preconditioner,
          VectorType &              dst,
          const VectorType &        src)
    {
      solver.solve(matrix, dst, src, preconditioner);
    }

    template <
      class VectorType,
      class SolverType,
      class MatrixType,
      class PreconditionerType,
      typename std::enable_if<
        !std::is_same<VectorType, typename SolverType::vector_type>::value,
        VectorType>::type * = nullptr>
    void
    solve(SolverType &              solver,
          const MatrixType &        matrix,
          const PreconditionerType &preconditioner,
          VectorType &              dst,
          const VectorType &        src)
    {
      typename SolverType::vector_type src_;
      typename SolverType::vector_type dst_;

      src_ = src;
      dst_ = dst;

      solver.solve(matrix, dst_, src_, preconditioner);

      dst = dst_;
    }
  } // namespace MGCoarseGridIterativeSolver
} // namespace internal



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
                       MGCoarseGridIterativeSolver<
                         VectorType,
                         SolverType,
                         MatrixType,
                         PreconditionerType>::operator()(const unsigned int /*level*/,
                                  VectorType &      dst,
                                  const VectorType &src) const
{
  Assert(solver != nullptr, ExcNotInitialized());
  Assert(matrix != nullptr, ExcNotInitialized());
  Assert(preconditioner != nullptr, ExcNotInitialized());

  dst = 0;
  internal::MGCoarseGridIterativeSolver::solve(
    *solver, *matrix, *preconditioner, dst, src);
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
MGCoarseGridHouseholder<number, VectorType>::operator()(
  const unsigned int /*level*/,
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
