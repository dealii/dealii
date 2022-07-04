// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#ifndef dealii_trilinos_linear_operator_h
#define dealii_trilinos_linear_operator_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/block_linear_operator.h>
#  include <deal.II/lac/linear_operator.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // Forward declarations:
#  ifndef DOXYGEN
  class SparseMatrix;
  class PreconditionBase;
  class BlockSparseMatrix;

  namespace internal
  {
    namespace LinearOperatorImplementation
    {
      class TrilinosPayload;
    }

    namespace BlockLinearOperatorImplementation
    {
      template <typename PayloadBlockType>
      class TrilinosBlockPayload;
    }
  } // namespace internal
#  endif

  /**
   * @name Creation of a LinearOperator
   */
  //@{


  /**
   * @relatesalso LinearOperator
   *
   * A function that encapsulates generic @p matrix objects, based on an
   * @p operator_exemplar, that act on a compatible Vector type into a
   * LinearOperator.
   *
   * This function is the equivalent of the dealii::linear_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <typename Range, typename Domain = Range, typename Matrix>
  inline LinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>
  linear_operator(const TrilinosWrappers::SparseMatrix &operator_exemplar,
                  const Matrix &                        matrix)
  {
    using OperatorExemplar = TrilinosWrappers::SparseMatrix;
    using Payload =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    return dealii::
      linear_operator<Range, Domain, Payload, OperatorExemplar, Matrix>(
        operator_exemplar, matrix);
  }


  /**
   * @relatesalso LinearOperator
   *
   * A function that encapsulates generic @p matrix objects that act on a
   * compatible Vector type into a LinearOperator.
   *
   * This function is the equivalent of the dealii::linear_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <typename Range, typename Domain = Range>
  inline LinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>
  linear_operator(const TrilinosWrappers::SparseMatrix &matrix)
  {
    using Matrix = TrilinosWrappers::SparseMatrix;
    using Payload =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    return dealii::linear_operator<Range, Domain, Payload, Matrix, Matrix>(
      matrix, matrix);
  }


  /**
   * @relatesalso LinearOperator
   *
   * A function that encapsulates generic @p matrix objects, based on an
   * @p operator_exemplar, that act on a compatible Vector type into a
   * LinearOperator.
   *
   * This function is the equivalent of the dealii::linear_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <typename Range, typename Domain, typename Matrix>
  inline LinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>
  linear_operator(
    const LinearOperator<
      Range,
      Domain,
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>
      &           operator_exemplar,
    const Matrix &matrix)
  {
    using Payload =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using OperatorExemplar = LinearOperator<Range, Domain, Payload>;
    return dealii::
      linear_operator<Range, Domain, Payload, OperatorExemplar, Matrix>(
        operator_exemplar, matrix);
  }


  //@}
  /**
   * @name Creation of a BlockLinearOperator
   */
  //@{


  /**
   * @relatesalso BlockLinearOperator
   *
   * A function that encapsulates a @p block_matrix into a BlockLinearOperator.
   *
   * This function is the equivalent of the dealii::block_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <typename Range, typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_operator(const TrilinosWrappers::BlockSparseMatrix &block_matrix)
  {
    using BlockMatrix = TrilinosWrappers::BlockSparseMatrix;
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::block_operator<Range, Domain, BlockPayload, BlockMatrix>(
      block_matrix);
  }


  /**
   * @relatesalso BlockLinearOperator
   *
   * A variant of above function that builds up a block diagonal linear operator
   * from an array @p ops of diagonal elements (off-diagonal blocks are assumed
   * to be 0).
   *
   * This function is the equivalent of the dealii::block_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <std::size_t m,
            std::size_t n,
            typename Range,
            typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_operator(
    const std::array<
      std::array<
        LinearOperator<typename Range::BlockType,
                       typename Domain::BlockType,
                       TrilinosWrappers::internal::
                         LinearOperatorImplementation::TrilinosPayload>,
        n>,
      m> &ops)
  {
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::block_operator<m, n, Range, Domain, BlockPayload>(ops);
  }


  /**
   * @relatesalso BlockLinearOperator
   *
   * This function extracts the diagonal blocks of @p block_matrix (either a
   * block matrix type or a BlockLinearOperator) and creates a
   * BlockLinearOperator with the diagonal. Off-diagonal elements are
   * initialized as null_operator (with correct reinit_range_vector and
   * reinit_domain_vector methods).
   *
   * This function is the equivalent of the dealii::block_diagonal_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <typename Range, typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_diagonal_operator(
    const TrilinosWrappers::BlockSparseMatrix &block_matrix)
  {
    using BlockMatrix = TrilinosWrappers::BlockSparseMatrix;
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::
      block_diagonal_operator<Range, Domain, BlockPayload, BlockMatrix>(
        block_matrix);
  }


  /**
   * @relatesalso BlockLinearOperator
   *
   * A variant of above function that builds up a block diagonal linear operator
   * from an array @p ops of diagonal elements (off-diagonal blocks are assumed
   * to be 0).
   *
   * This function is the equivalent of the dealii::block_diagonal_operator, but
   * ensures full compatibility with Trilinos operations by preselecting the
   * appropriate template parameters.
   *
   *
   * @ingroup TrilinosWrappers
   */
  template <std::size_t m, typename Range, typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_diagonal_operator(
    const std::array<
      LinearOperator<typename Range::BlockType,
                     typename Domain::BlockType,
                     TrilinosWrappers::internal::LinearOperatorImplementation::
                       TrilinosPayload>,
      m> &ops)
  {
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::block_diagonal_operator<m, Range, Domain, BlockPayload>(ops);
  }

  //@}

} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
#endif
