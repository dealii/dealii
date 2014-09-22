// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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


#ifndef __deal2__sparse_ilu_h
#define __deal2__sparse_ilu_h


#include <deal.II/base/config.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * This class computes an Incomplete LU (ILU) decomposition of a sparse matrix,
 * using either the same sparsity pattern or a different one.
 * By incomplete we mean that unlike the exact decomposition, the incomplete
 * one is also computed using sparse factors, and entries
 * in the decomposition that do not fit into the given sparsity structure
 * are discarded.
 *
 * The algorithm used by this class is essentially a copy of the
 * algorithm given in the book Y. Saad: "Iterative methods for sparse
 * linear systems", second edition, in section 10.3.2.
 *
 *
 * <h3>Usage and state management</h3>
 *
 * Refer to SparseLUDecomposition documentation for suggested
 * usage and state management. This class is used in the @ref
 * step_22 "step-22" tutorial program.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @author Wolfgang Bangerth, 2008, 2009; unified interface: Ralf Hartmann
 */
template <typename number>
class SparseILU : public SparseLUDecomposition<number>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor. Does nothing.
   *
   * Call the @p initialize
   * function before using this
   * object as preconditioner.
   */
  SparseILU ();

  /**
   * @deprecated This method is deprecated, and
   * left for backward
   * compatibility. It will be removed
   * in later versions.
   * Instead, pass the sparsity pattern that you want used for
   * the decomposition in the AdditionalData structure.
   */
  SparseILU (const SparsityPattern &sparsity) DEAL_II_DEPRECATED;

  /**
   * Make
   * SparseLUDecomposition::AdditionalData
   * accessible to this class as
   * well.
   */
  typedef
  typename SparseLUDecomposition<number>::AdditionalData
  AdditionalData;

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
   * this function performs the LU
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
   * compability. It will be removed in later versions.
   */
  template <typename somenumber>
  void decompose (const SparseMatrix<somenumber> &matrix,
                  const double                    strengthen_diagonal=0.) DEAL_II_DEPRECATED;

  /**
   * @deprecated This method is deprecated, and
   * left for backward
   * compatibility. It will be
   * removed in later versions.
   */
  template <typename somenumber>
  void apply_decomposition (Vector<somenumber>       &dst,
                            const Vector<somenumber> &src) const DEAL_II_DEPRECATED;

  /**
   * Apply the incomplete decomposition,
   * i.e. do one forward-backward step
   * $dst=(LU)^{-1}src$.
   *
   * The initialize() function
   * needs to be called before.
   */
  template <typename somenumber>
  void vmult (Vector<somenumber>       &dst,
              const Vector<somenumber> &src) const;


  /**
   * Apply the transpose of the
   * incomplete decomposition,
   * i.e. do one forward-backward step
   * $dst=(LU)^{-T}src$.
   *
   * The initialize() function
   * needs to be called before.
   */
  template <typename somenumber>
  void Tvmult (Vector<somenumber>       &dst,
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
  DeclException1 (ExcInvalidStrengthening,
                  double,
                  << "The strengthening parameter " << arg1
                  << " is not greater or equal than zero!");
  //@}
};

/*@}*/
//---------------------------------------------------------------------------

template <typename number>
template <typename somenumber>
inline
void
SparseILU<number>::apply_decomposition (Vector<somenumber>       &dst,
                                        const Vector<somenumber> &src) const
{
  vmult (dst, src);
}




DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__sparse_ilu_h
