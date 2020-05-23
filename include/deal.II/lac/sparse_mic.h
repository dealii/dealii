// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2019 by the deal.II authors
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

#ifndef dealii_sparse_mic_h
#define dealii_sparse_mic_h

#include <deal.II/base/config.h>

#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Implementation of the Modified Incomplete Cholesky (MIC(0)) preconditioner
 * for symmetric matrices. This class conforms to the state and usage
 * specification in SparseLUDecomposition.
 *
 *
 * <h3>The decomposition</h3>
 *
 * Let a symmetric, positive-definite, sparse matrix $A$ be in the form $A = D
 * - L - L^T$, where $D$ is the diagonal part of $A$ and $-L$ is a strictly
 * lower triangular matrix. The MIC(0) decomposition of the matrix $A$ is
 * defined by $B = (X-L)X^{-1}(X-L^T)$, where $X$ is a diagonal matrix defined
 * by the condition $\text{rowsum}(A) = \text{rowsum}(B)$.
 *
 * @author Stephen "Cheffo" Kolaroff, 2002, unified interface: Ralf Hartmann
 * 2003; extension for full compatibility with LinearOperator class: Jean-Paul
 * Pelteret, 2015.
 */
template <typename number>
class SparseMIC : public SparseLUDecomposition<number>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Constructor. Does nothing, so you have to call @p decompose sometimes
   * afterwards.
   */
  SparseMIC();

  /**
   * Destructor.
   */
  virtual ~SparseMIC() override;

  /**
   * Deletes all member variables. Leaves the class in the state that it had
   * directly after calling the constructor
   */
  virtual void
  clear() override;

  /**
   * Make the @p AdditionalData type in the base class accessible to this
   * class as well.
   */
  using AdditionalData = typename SparseLUDecomposition<number>::AdditionalData;

  /**
   * Perform the incomplete LU factorization of the given matrix.
   *
   * This function needs to be called before an object of this class is used
   * as preconditioner.
   *
   * For more details about possible parameters, see the class documentation
   * of SparseLUDecomposition and the documentation of the @p
   * SparseLUDecomposition::AdditionalData class.
   *
   * According to the @p parameters, this function creates a new
   * SparsityPattern or keeps the previous sparsity or takes the sparsity
   * given by the user to @p data. Then, this function performs the MIC
   * decomposition.
   *
   * After this function is called the preconditioner is ready to be used.
   */
  template <typename somenumber>
  void
  initialize(const SparseMatrix<somenumber> &matrix,
             const AdditionalData &          parameters = AdditionalData());

  /**
   * Apply the incomplete decomposition, i.e. do one forward-backward step
   * $dst=(LU)^{-1}src$.
   *
   * Call @p initialize before calling this function.
   */
  template <typename somenumber>
  void
  vmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;

  /**
   * Apply the transpose of the incomplete decomposition, i.e. do one forward-
   * backward step $dst=(LU)^{-1}src$.
   *
   * Call @p initialize before calling this function.
   *
   * @note This function has not yet been implemented
   *
   */
  template <typename somenumber>
  void
  Tvmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const override;

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException0(ExcStrengthenDiagonalTooSmall);
  /**
   * Exception
   */
  DeclException1(ExcInvalidStrengthening,
                 double,
                 << "The strengthening parameter " << arg1
                 << " is not greater or equal than zero!");
  /**
   * Exception
   */
  DeclException2(ExcDecompositionNotStable,
                 int,
                 double,
                 << "The diagonal element (" << arg1 << "," << arg1 << ") is "
                 << arg2 << ", but must be positive");

  //@}
private:
  /**
   * Values of the computed diagonal.
   */
  std::vector<number> diag;

  /**
   * Inverses of the diagonal: precomputed for faster vmult.
   */
  std::vector<number> inv_diag;

  /**
   * Values of the computed "inner sums", i.e. per-row sums of the elements
   * laying on the right side of the diagonal.
   */
  std::vector<number> inner_sums;

  /**
   * Compute the row-th "inner sum".
   */
  number
  get_rowsum(const size_type row) const;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_
