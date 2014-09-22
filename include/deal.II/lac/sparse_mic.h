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

#ifndef __deal2__sparse_mic_h
#define __deal2__sparse_mic_h

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_decomposition.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Modified incomplete Cholesky (MIC(0)) preconditioner.  This class
 * conforms to the state and usage specification in
 * SparseLUDecomposition.
 *
 *
 * <h3>The decomposition</h3>
 *
 * Let a sparse matrix $A$ be in the form $A = - L - U + D$, where $-L$ and
 * $-U$ are strictly lower and upper triangular matrices. The MIC(0)
 * decomposition of the matrix $A$ is defined by $B = (X-L)X^(-1)(X-U)$,
 * where $X$ is a diagonal matrix, defined by the condition $\text{rowsum}(A) =
 * \text{rowsum}(B)$.
 *
 * @author Stephen "Cheffo" Kolaroff, 2002, unified interface: Ralf
 * Hartmann 2003.
 */
template <typename number>
class SparseMIC : public SparseLUDecomposition<number>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor. Does nothing, so
   * you have to call @p decompose
   * sometimes afterwards.
   */
  SparseMIC ();

  /**
   * @deprecated This method is deprecated, and left for backward
   * compatibility. It will be removed in later versions.  Instead,
   * pass the sparsity pattern that you want used for the
   * decomposition in the AdditionalData structure.
   */
  SparseMIC (const SparsityPattern &sparsity) DEAL_II_DEPRECATED;

  /**
   * Destructor.
   */
  virtual ~SparseMIC();

  /**
   * Deletes all member
   * variables. Leaves the class in
   * the state that it had directly
   * after calling the constructor
   */
  virtual void clear();

  /**
   * Make the @p AdditionalData
   * type in the base class
   * accessible to this class as
   * well.
   */
  typedef
  typename SparseLUDecomposition<number>::AdditionalData
  AdditionalData;

  /**
   * @deprecated This method is deprecated, and
   * left for backward
   * compatibility. It will be
   * removed in later versions.
   */
  void reinit (const SparsityPattern &sparsity) DEAL_II_DEPRECATED;

  /**
   * Perform the incomplete LU
   * factorization of the given
   * matrix.
   *
   * This function needs to be
   * called before an object of
   * this class is used as
   * preconditioner.
   *
   * For more details about
   * possible parameters, see the
   * class documentation of
   * SparseLUDecomposition and the
   * documentation of the
   * @p SparseLUDecomposition::AdditionalData
   * class.
   *
   * According to the
   * @p parameters, this function
   * creates a new SparsityPattern
   * or keeps the previous sparsity
   * or takes the sparsity given by
   * the user to @p data. Then,
   * this function performs the MIC
   * decomposition.
   *
   * After this function is called
   * the preconditioner is ready to
   * be used.
   */
  template <typename somenumber>
  void initialize (const SparseMatrix<somenumber> &matrix,
                   const AdditionalData &parameters = AdditionalData());

  /**
   * @deprecated This method is deprecated, and left for backward
   * compability. It will be removed in later versions. Use
   * initialize() instead.
   */
  template <typename somenumber>
  void decompose (const SparseMatrix<somenumber> &matrix,
                  const double                   strengthen_diagonal=0.) DEAL_II_DEPRECATED;

  /**
   * Apply the incomplete decomposition,
   * i.e. do one forward-backward step
   * $dst=(LU)^{-1}src$.
   *
   * Call @p initialize before
   * calling this function.
   */
  template <typename somenumber>
  void vmult (Vector<somenumber>       &dst,
              const Vector<somenumber> &src) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /** @addtogroup Exceptions
   * @{ */

  /**
   * Exception
   */
  DeclException0 (ExcStrengthenDiagonalTooSmall);
  /**
   * Exception
   */
  DeclException1 (ExcInvalidStrengthening,
                  double,
                  << "The strengthening parameter " << arg1
                  << " is not greater or equal than zero!");
  /**
   * Exception
   */
  DeclException2(ExcDecompositionNotStable, int, double,
                 << "The diagonal element (" <<arg1<<","<<arg1<<") is "
                 << arg2 <<", but must be positive");

  //@}
private:
  /**
   * Values of the computed
   * diagonal.
   */
  std::vector<number> diag;

  /**
   * Inverses of the the diagonal:
   * precomputed for faster vmult.
   */
  std::vector<number> inv_diag;

  /**
   * Values of the computed "inner
   * sums", i.e. per-row sums of
   * the elements laying on the
   * right side of the diagonal.
   */
  std::vector<number> inner_sums;

  /**
   * Compute the row-th "inner
   * sum".
   */
  number get_rowsum (const size_type row) const;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif  // __deal2__
