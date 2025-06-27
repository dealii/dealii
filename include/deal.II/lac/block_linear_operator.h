// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_block_linear_operator_h
#define dealii_block_linear_operator_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/linear_operator.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
namespace internal
{
  namespace BlockLinearOperatorImplementation
  {
    template <typename PayloadBlockType =
                internal::LinearOperatorImplementation::EmptyPayload>
    class EmptyBlockPayload;
  }
} // namespace internal

template <typename Number>
class BlockVector;

template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
class BlockLinearOperator;
#endif

template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(const BlockMatrixType &matrix);

template <std::size_t m,
          std::size_t n,
          typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(
  const std::array<std::array<LinearOperator<typename Range::BlockType,
                                             typename Domain::BlockType,
                                             typename BlockPayload::BlockType>,
                              n>,
                   m> &);

template <std::size_t m,
          typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const std::array<LinearOperator<typename Range::BlockType,
                                  typename Domain::BlockType,
                                  typename BlockPayload::BlockType>,
                   m> &);

template <std::size_t m,
          typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const LinearOperator<typename Range::BlockType,
                       typename Domain::BlockType,
                       typename BlockPayload::BlockType> &op);



/**
 * A class to store the concept of a block linear operator.
 *
 * This class increases the interface of LinearOperator (which encapsulates
 * the  @p Matrix interface) by three additional functions:
 * @code
 *   std::function<unsigned int()> n_block_rows;
 *   std::function<unsigned int()> n_block_cols;
 *   std::function<BlockType(unsigned int, unsigned int)> block;
 * @endcode
 * that describe the underlying block structure (of an otherwise opaque)
 * linear operator.
 *
 * Objects of type BlockLinearOperator can be created similarly to
 * LinearOperator with a wrapper function:
 * @code
 * dealii::BlockSparseMatrix<double> A;
 * const auto block_op_a = block_operator(A);
 * @endcode
 *
 * Alternatively, there are several helper functions available for creating
 * instances from multiple independent matrices of possibly different types.
 * Here is an example of a block diagonal matrix created from a FullMatrix and
 * a SparseMatrixEZ:
 *
 * @code
 * FullMatrix<double> top_left(2, 2);
 * top_left(0, 0) = 2.0;
 * top_left(0, 1) = -1.0;
 * top_left(1, 0) = -1.0;
 * top_left(1, 1) = 2.0;
 *
 * SparseMatrixEZ<double> bottom_right(4, 4, 4);
 * for (std::size_t row_n = 0; row_n < 4; ++row_n)
 *   {
 *     bottom_right.add(row_n, row_n, 1.0);
 *     if (row_n < 3)
 *       bottom_right.add(row_n, row_n + 1, -1.0);
 *   }
 *
 * auto top_left_op = linear_operator(top_left);
 * auto bottom_right_op = linear_operator(bottom_right);
 * std::array<decltype(top_left_op), 2> operators {{top_left_op,
 *                                                  bottom_right_op}};
 * auto block_op = block_diagonal_operator (operators);
 *
 * std::vector<BlockVector<double>::size_type> block_sizes {2, 4};
 * BlockVector<double> src(block_sizes);
 * src = 2.0;
 * BlockVector<double> dst(block_sizes);
 * block_op.vmult(dst, src); // now equal to 2, 2, 0, 0, 0, 2
 * @endcode
 *
 *
 * A BlockLinearOperator can be sliced to a LinearOperator at any time. This
 * removes all information about the underlying block structure (because above
 * <code>std::function</code> objects are no longer available) - the linear
 * operator interface, however, remains intact.
 *
 * @note This class makes heavy use of <code>std::function</code> objects and
 * lambda functions. This flexibility comes with a run-time penalty. Only use
 * this object to encapsulate object with medium to large individual block
 * sizes, and small block structure (as a rule of thumb, matrix blocks greater
 * than $1000\times1000$).
 *
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename BlockPayload>
class BlockLinearOperator
  : public LinearOperator<Range, Domain, typename BlockPayload::BlockType>
{
public:
  using BlockType = LinearOperator<typename Range::BlockType,
                                   typename Domain::BlockType,
                                   typename BlockPayload::BlockType>;

  /**
   * Create an empty BlockLinearOperator object.
   *
   * All<code>std::function</code> member objects of this class and its base
   * class LinearOperator are initialized with default variants that throw an
   * exception upon invocation.
   */
  BlockLinearOperator(const BlockPayload &payload)
    : LinearOperator<Range, Domain, typename BlockPayload::BlockType>(
        typename BlockPayload::BlockType(payload, payload))
  {
    n_block_rows = []() -> unsigned int {
      Assert(
        false,
        ExcMessage(
          "Uninitialized BlockLinearOperator<Range, Domain>::n_block_rows called"));
      return 0;
    };

    n_block_cols = []() -> unsigned int {
      Assert(
        false,
        ExcMessage(
          "Uninitialized BlockLinearOperator<Range, Domain>::n_block_cols called"));
      return 0;
    };

    block = [](unsigned int, unsigned int) -> BlockType {
      Assert(
        false,
        ExcMessage(
          "Uninitialized BlockLinearOperator<Range, Domain>::block called"));
      return BlockType();
    };
  }

  /**
   * Default copy constructor.
   */
  BlockLinearOperator(
    const BlockLinearOperator<Range, Domain, BlockPayload> &) = default;

  /**
   * Templated copy constructor that creates a BlockLinearOperator object from
   * an object @p op for which the conversion function
   * <code>block_operator</code> is defined.
   */
  template <typename Op>
  BlockLinearOperator(const Op &op)
  {
    *this = block_operator<Range, Domain, BlockPayload, Op>(op);
  }

  /**
   * Create a BlockLinearOperator from a two-dimensional array @p ops of
   * LinearOperator. This constructor calls the corresponding block_operator()
   * specialization.
   */
  template <std::size_t m, std::size_t n>
  BlockLinearOperator(const std::array<std::array<BlockType, n>, m> &ops)
  {
    *this = block_operator<m, n, Range, Domain, BlockPayload>(ops);
  }

  /**
   * Create a block-diagonal BlockLinearOperator from a one-dimensional array
   * @p ops of LinearOperator. This constructor calls the corresponding
   * block_operator() specialization.
   */
  template <std::size_t m>
  BlockLinearOperator(const std::array<BlockType, m> &ops)
  {
    *this = block_diagonal_operator<m, Range, Domain, BlockPayload>(ops);
  }

  /**
   * Default copy assignment operator.
   */
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const BlockLinearOperator<Range, Domain, BlockPayload> &) = default;

  /**
   * Templated copy assignment operator for an object @p op for which the
   * conversion function <code>block_operator</code> is defined.
   */
  template <typename Op>
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const Op &op)
  {
    *this = block_operator<Range, Domain, BlockPayload, Op>(op);
    return *this;
  }

  /**
   * Copy assignment from a two-dimensional array @p ops of LinearOperator.
   * This assignment operator calls the corresponding block_operator()
   * specialization.
   */
  template <std::size_t m, std::size_t n>
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const std::array<std::array<BlockType, n>, m> &ops)
  {
    *this = block_operator<m, n, Range, Domain, BlockPayload>(ops);
    return *this;
  }

  /**
   * Copy assignment from a one-dimensional array @p ops of LinearOperator
   * that creates a block-diagonal BlockLinearOperator. This assignment
   * operator calls the corresponding block_operator() specialization.
   */
  template <std::size_t m>
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const std::array<BlockType, m> &ops)
  {
    *this = block_diagonal_operator<m, Range, Domain, BlockPayload>(ops);
    return *this;
  }

  /**
   * Return the number of blocks in a column (i.e, the number of "block rows",
   * or the number $m$, if interpreted as a $m\times n$ block system).
   */
  std::function<unsigned int()> n_block_rows;

  /**
   * Return the number of blocks in a row (i.e, the number of "block columns",
   * or the number $n$, if interpreted as a $m\times n$ block system).
   */
  std::function<unsigned int()> n_block_cols;

  /**
   * Access the block with the given coordinates. This
   * <code>std::function</code> object returns a LinearOperator representing
   * the $(i,j)$-th block of the BlockLinearOperator.
   */
  std::function<BlockType(unsigned int, unsigned int)> block;
};


namespace internal
{
  namespace BlockLinearOperatorImplementation
  {
    // A helper function to apply a given vmult, or Tvmult to a vector with
    // intermediate storage, similar to the corresponding helper
    // function for LinearOperator. Here, two operators are used.
    // The first one takes care of the first "column" and typically doesn't add.
    // On the other hand, the second operator is normally an adding one.
    template <typename Function1,
              typename Function2,
              typename Range,
              typename Domain>
    void
    apply_with_intermediate_storage(const Function1 &first_op,
                                    const Function2 &loop_op,
                                    Range           &v,
                                    const Domain    &u,
                                    bool             add)
    {
      GrowingVectorMemory<Range> vector_memory;

      typename VectorMemory<Range>::Pointer tmp(vector_memory);
      tmp->reinit(v, /*bool omit_zeroing_entries =*/true);

      const unsigned int n = u.n_blocks();
      const unsigned int m = v.n_blocks();

      for (unsigned int i = 0; i < m; ++i)
        {
          first_op(*tmp, u, i, 0);
          for (unsigned int j = 1; j < n; ++j)
            loop_op(*tmp, u, i, j);
        }

      if (add)
        v += *tmp;
      else
        v = *tmp;
    }

    // Populate the LinearOperator interfaces with the help of the
    // BlockLinearOperator functions
    template <typename Range, typename Domain, typename BlockPayload>
    inline void
    populate_linear_operator_functions(
      dealii::BlockLinearOperator<Range, Domain, BlockPayload> &op)
    {
      op.reinit_range_vector = [=](Range &v, bool omit_zeroing_entries) {
        const unsigned int m = op.n_block_rows();

        // Reinitialize the block vector to m blocks:
        v.reinit(m);

        // And reinitialize every individual block with reinit_range_vectors:
        for (unsigned int i = 0; i < m; ++i)
          op.block(i, 0).reinit_range_vector(v.block(i), omit_zeroing_entries);

        v.collect_sizes();
      };

      op.reinit_domain_vector = [=](Domain &v, bool omit_zeroing_entries) {
        const unsigned int n = op.n_block_cols();

        // Reinitialize the block vector to n blocks:
        v.reinit(n);

        // And reinitialize every individual block with reinit_domain_vectors:
        for (unsigned int i = 0; i < n; ++i)
          op.block(0, i).reinit_domain_vector(v.block(i), omit_zeroing_entries);

        v.collect_sizes();
      };

      op.vmult = [&op](Range &v, const Domain &u) {
        const unsigned int m = op.n_block_rows();
        const unsigned int n = op.n_block_cols();
        Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
        Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range             &v,
                                        const Domain      &u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(i, j).vmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range             &v,
                                       const Domain      &u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(i, j).vmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, false);
          }
        else
          {
            for (unsigned int i = 0; i < m; ++i)
              {
                op.block(i, 0).vmult(v.block(i), u.block(0));
                for (unsigned int j = 1; j < n; ++j)
                  op.block(i, j).vmult_add(v.block(i), u.block(j));
              }
          }
      };

      op.vmult_add = [&op](Range &v, const Domain &u) {
        const unsigned int m = op.n_block_rows();
        const unsigned int n = op.n_block_cols();
        Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
        Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range             &v,
                                        const Domain      &u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(i, j).vmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range             &v,
                                       const Domain      &u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(i, j).vmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, true);
          }
        else
          {
            for (unsigned int i = 0; i < m; ++i)
              for (unsigned int j = 0; j < n; ++j)
                op.block(i, j).vmult_add(v.block(i), u.block(j));
          }
      };

      op.Tvmult = [&op](Domain &v, const Range &u) {
        const unsigned int n = op.n_block_cols();
        const unsigned int m = op.n_block_rows();
        Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
        Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range             &v,
                                        const Domain      &u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(j, i).Tvmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range             &v,
                                       const Domain      &u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(j, i).Tvmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, false);
          }
        else
          {
            for (unsigned int i = 0; i < n; ++i)
              {
                op.block(0, i).Tvmult(v.block(i), u.block(0));
                for (unsigned int j = 1; j < m; ++j)
                  op.block(j, i).Tvmult_add(v.block(i), u.block(j));
              }
          }
      };

      op.Tvmult_add = [&op](Domain &v, const Range &u) {
        const unsigned int n = op.n_block_cols();
        const unsigned int m = op.n_block_rows();
        Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
        Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range             &v,
                                        const Domain      &u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(j, i).Tvmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range             &v,
                                       const Domain      &u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(j, i).Tvmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, true);
          }
        else
          {
            for (unsigned int i = 0; i < n; ++i)
              for (unsigned int j = 0; j < m; ++j)
                op.block(j, i).Tvmult_add(v.block(i), u.block(j));
          }
      };
    }



    /**
     * A dummy class for BlockLinearOperators that do not require any
     * extensions to facilitate the operations of the block matrix or its
     * subblocks.
     *
     * This is the Payload class typically associated with deal.II's native
     * BlockSparseMatrix. To use either TrilinosWrappers::BlockSparseMatrix or
     * PETScWrappers::BlockSparseMatrix one must initialize a
     * BlockLinearOperator with their associated BlockPayload.
     *
     *
     * @ingroup LAOperators
     */
    template <typename PayloadBlockType>
    class EmptyBlockPayload
    {
    public:
      /**
       * Type of payload held by each subblock
       */
      using BlockType = PayloadBlockType;

      /**
       * Default constructor
       *
       * Since this class does not do anything in particular and needs no
       * special configuration, we have only one generic constructor that can
       * be called under any conditions.
       */
      template <typename... Args>
      EmptyBlockPayload(const Args &...)
      {}
    };

  } // namespace BlockLinearOperatorImplementation
} // namespace internal



/**
 * @name Creation of a BlockLinearOperator
 */
/** @{ */

/**
 * @relatesalso BlockLinearOperator
 *
 * A function that encapsulates a @p block_matrix into a BlockLinearOperator.
 *
 * All changes made on the block structure and individual blocks of @p
 * block_matrix after the creation of the BlockLinearOperator object are
 * reflected by the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename Domain,
          typename BlockPayload,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(const BlockMatrixType &block_matrix)
{
  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  BlockLinearOperator<Range, Domain, BlockPayload> return_op{
    BlockPayload(block_matrix, block_matrix)};

  return_op.n_block_rows = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_rows();
  };

  return_op.n_block_cols = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_cols();
  };

  return_op.block = [&block_matrix](unsigned int i,
                                    unsigned int j) -> BlockType {
    if constexpr (running_in_debug_mode())
      {
        const unsigned int m = block_matrix.n_block_rows();
        const unsigned int n = block_matrix.n_block_cols();
        AssertIndexRange(i, m);
        AssertIndexRange(j, n);
      }

    return BlockType(block_matrix.block(i, j));
  };

  populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relatesalso BlockLinearOperator
 *
 * A variant of above function that encapsulates a given collection @p ops of
 * LinearOperators into a block structure. Here, it is assumed that Range and
 * Domain are block vectors, i.e., derived from
 * @ref BlockVectorBase.
 * The individual linear operators in @p ops must act on the underlying vector
 * type of the block vectors, i.e., on Domain::BlockType yielding a result in
 * Range::BlockType.
 *
 * The list @p ops is best passed as an initializer list. Consider for example
 * a linear operator block (acting on Vector<double>)
 * @code
 *  op_a00 | op_a01
 *         |
 *  ---------------
 *         |
 *  op_a10 | op_a11
 * @endcode
 * The corresponding block_operator invocation takes the form
 * @code
 * block_operator<2, 2, BlockVector<double>>({op_a00, op_a01, op_a10, op_a11});
 * @endcode
 *
 * @ingroup LAOperators
 */
template <std::size_t m,
          std::size_t n,
          typename Range,
          typename Domain,
          typename BlockPayload>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(
  const std::array<std::array<LinearOperator<typename Range::BlockType,
                                             typename Domain::BlockType,
                                             typename BlockPayload::BlockType>,
                              n>,
                   m> &ops)
{
  static_assert(m > 0 && n > 0,
                "a blocked LinearOperator must consist of at least one block");

  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  // TODO: Create block payload so that this can be initialized correctly
  BlockLinearOperator<Range, Domain, BlockPayload> return_op{BlockPayload()};

  return_op.n_block_rows = []() -> unsigned int { return m; };

  return_op.n_block_cols = []() -> unsigned int { return n; };

  return_op.block = [ops](unsigned int i, unsigned int j) -> BlockType {
    AssertIndexRange(i, m);
    AssertIndexRange(j, n);

    return ops[i][j];
  };

  populate_linear_operator_functions(return_op);
  return return_op;
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
 * All changes made on the individual diagonal blocks of @p block_matrix after
 * the creation of the BlockLinearOperator object are reflected by the
 * operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(const BlockMatrixType &block_matrix)
{
  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  BlockLinearOperator<Range, Domain, BlockPayload> return_op{
    BlockPayload(block_matrix, block_matrix)};

  return_op.n_block_rows = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_rows();
  };

  return_op.n_block_cols = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_cols();
  };

  return_op.block = [&block_matrix](unsigned int i,
                                    unsigned int j) -> BlockType {
    if constexpr (running_in_debug_mode())
      {
        const unsigned int m = block_matrix.n_block_rows();
        const unsigned int n = block_matrix.n_block_cols();
        Assert(m == n, ExcDimensionMismatch(m, n));
        AssertIndexRange(i, m);
        AssertIndexRange(j, n);
      }
    if (i == j)
      return BlockType(block_matrix.block(i, j));
    else
      return null_operator(BlockType(block_matrix.block(i, j)));
  };

  populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relatesalso BlockLinearOperator
 *
 * A variant of above function that builds up a block diagonal linear operator
 * from an array @p ops of diagonal elements (off-diagonal blocks are assumed
 * to be 0).
 *
 * The list @p ops is best passed as an initializer list. Consider for example
 * a linear operator block (acting on Vector<double>) <code>diag(op_a0, op_a1,
 * ..., op_am)</code>. The corresponding block_operator invocation takes the
 * form
 * @code
 * block_diagonal_operator<m, BlockVector<double>>({op_00, op_a1, ..., op_am});
 * @endcode
 *
 * @ingroup LAOperators
 */
template <std::size_t m, typename Range, typename Domain, typename BlockPayload>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const std::array<LinearOperator<typename Range::BlockType,
                                  typename Domain::BlockType,
                                  typename BlockPayload::BlockType>,
                   m> &ops)
{
  static_assert(
    m > 0, "a blockdiagonal LinearOperator must consist of at least one block");

  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  std::array<std::array<BlockType, m>, m> new_ops;

  // This is a bit tricky. We have to make sure that the off-diagonal
  // elements of return_op.ops are populated correctly. They must be
  // null_operators, but with correct reinit_domain_vector and
  // reinit_range_vector functions.
  for (unsigned int i = 0; i < m; ++i)
    for (unsigned int j = 0; j < m; ++j)
      if (i == j)
        {
          // diagonal elements are easy:
          new_ops[i][j] = ops[i];
        }
      else
        {
          // create a null-operator...
          new_ops[i][j] = null_operator(ops[i]);
          // ... and fix up reinit_domain_vector:
          new_ops[i][j].reinit_domain_vector = ops[j].reinit_domain_vector;
        }

  return block_operator<m, m, Range, Domain>(new_ops);
}



/**
 * @relatesalso BlockLinearOperator
 *
 * A variant of above function that only takes a single LinearOperator
 * argument @p op and creates a blockdiagonal linear operator with @p m copies
 * of it.
 *
 * @ingroup LAOperators
 */
template <std::size_t m, typename Range, typename Domain, typename BlockPayload>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const LinearOperator<typename Range::BlockType,
                       typename Domain::BlockType,
                       typename BlockPayload::BlockType> &op)
{
  static_assert(m > 0,
                "a blockdiagonal LinearOperator must consist of at least "
                "one block");

  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;
  std::array<BlockType, m> new_ops;
  new_ops.fill(op);

  return block_diagonal_operator(new_ops);
}



/** @} */
/**
 * @name Manipulation of a BlockLinearOperator
 */
/** @{ */

/**
 * @relatesalso LinearOperator
 *
 * This function implements forward substitution to invert a lower block
 * triangular matrix. As arguments, it takes a BlockLinearOperator @p
 * block_operator representing a block lower triangular matrix, as well as a
 * BlockLinearOperator @p diagonal_inverse representing inverses of diagonal
 * blocks of @p block_operator.
 *
 * Let us assume we have a linear system with the following block structure:
 *
 * @code
 * A00 x0 + ...                   = y0
 * A01 x0 + A11 x1 + ...          = y1
 * ...        ...
 * A0n x0 + A1n x1 + ... + Ann xn = yn
 * @endcode
 *
 * First of all, <code>x0 = A00^-1 y0</code>. Then, we can use x0 to recover
 * x1:
 * @code
 *    x1 = A11^-1 ( y1 - A01 x0 )
 * @endcode
 * and therefore:
 * @code
 *    xn = Ann^-1 ( yn - A0n x0 - ... - A(n-1)n x(n-1) )
 * @endcode
 *
 * @note We are not using all blocks of the BlockLinearOperator arguments:
 * Just the lower triangular block matrix of @p block_operator is used as well
 * as the diagonal of @p diagonal_inverse.
 *
 * @ingroup LAOperators
 */
template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
LinearOperator<Domain, Range, typename BlockPayload::BlockType>
block_forward_substitution(
  const BlockLinearOperator<Range, Domain, BlockPayload> &block_operator,
  const BlockLinearOperator<Domain, Range, BlockPayload> &diagonal_inverse)
{
  LinearOperator<Range, Range, typename BlockPayload::BlockType> return_op{
    typename BlockPayload::BlockType(diagonal_inverse)};

  return_op.reinit_range_vector  = diagonal_inverse.reinit_range_vector;
  return_op.reinit_domain_vector = diagonal_inverse.reinit_domain_vector;

  return_op.vmult = [block_operator, diagonal_inverse](Range       &v,
                                                       const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    if (m == 0)
      return;

    diagonal_inverse.block(0, 0).vmult(v.block(0), u.block(0));
    for (unsigned int i = 1; i < m; ++i)
      {
        auto &dst = v.block(i);
        dst       = u.block(i);
        dst *= -1.;
        for (unsigned int j = 0; j < i; ++j)
          block_operator.block(i, j).vmult_add(dst, v.block(j));
        dst *= -1.;
        diagonal_inverse.block(i, i).vmult(dst,
                                           dst); // uses intermediate storage
      }
  };

  return_op.vmult_add = [block_operator, diagonal_inverse](Range       &v,
                                                           const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    if (m == 0)
      return;

    GrowingVectorMemory<typename Range::BlockType>            vector_memory;
    typename VectorMemory<typename Range::BlockType>::Pointer tmp(
      vector_memory);

    diagonal_inverse.block(0, 0).vmult_add(v.block(0), u.block(0));

    for (unsigned int i = 1; i < m; ++i)
      {
        diagonal_inverse.block(i, i).reinit_range_vector(
          *tmp, /*bool omit_zeroing_entries=*/true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j = 0; j < i; ++j)
          block_operator.block(i, j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        diagonal_inverse.block(i, i).vmult_add(v.block(i), *tmp);
      }
  };

  return return_op;
}



/**
 * @relatesalso BlockLinearOperator
 *
 * This function implements back substitution to invert an upper block
 * triangular matrix. As arguments, it takes a BlockLinearOperator @p
 * block_operator representing an upper block triangular matrix, as well as a
 * BlockLinearOperator @p diagonal_inverse representing inverses of diagonal
 * blocks of @p block_operator.
 *
 * Let us assume we have a linear system with the following block structure:
 *
 * @code
 * A00 x0 + A01 x1 + ... + A0n xn = yn
 *          A11 x1 + ...          = y1
 *                          ...     ..
 *                         Ann xn = yn
 * @endcode
 *
 * First of all, <code>xn = Ann^-1 yn</code>. Then, we can use xn to recover
 * x(n-1):
 * @code
 *    x(n-1) = A(n-1)(n-1)^-1 ( y(n-1) - A(n-1)n x(n-1) )
 * @endcode
 * and therefore:
 * @code
 *    x0 = A00^-1 ( y0 - A0n xn - ... - A01 x1 )
 * @endcode
 *
 * @note We are not using all blocks of the BlockLinearOperator arguments:
 * Just the upper triangular block matrix of @p block_operator is used as well
 * as the diagonal of @p diagonal_inverse.
 *
 * @ingroup LAOperators
 */
template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
LinearOperator<Domain, Range, typename BlockPayload::BlockType>
block_back_substitution(
  const BlockLinearOperator<Range, Domain, BlockPayload> &block_operator,
  const BlockLinearOperator<Domain, Range, BlockPayload> &diagonal_inverse)
{
  LinearOperator<Range, Range, typename BlockPayload::BlockType> return_op{
    typename BlockPayload::BlockType(diagonal_inverse)};

  return_op.reinit_range_vector  = diagonal_inverse.reinit_range_vector;
  return_op.reinit_domain_vector = diagonal_inverse.reinit_domain_vector;

  return_op.vmult = [block_operator, diagonal_inverse](Range       &v,
                                                       const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    if (m == 0)
      return;

    diagonal_inverse.block(m - 1, m - 1).vmult(v.block(m - 1), u.block(m - 1));

    for (int i = m - 2; i >= 0; --i)
      {
        auto &dst = v.block(i);
        dst       = u.block(i);
        dst *= -1.;
        for (unsigned int j = i + 1; j < m; ++j)
          block_operator.block(i, j).vmult_add(dst, v.block(j));
        dst *= -1.;
        diagonal_inverse.block(i, i).vmult(dst,
                                           dst); // uses intermediate storage
      }
  };

  return_op.vmult_add = [block_operator, diagonal_inverse](Range       &v,
                                                           const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));
    GrowingVectorMemory<typename Range::BlockType>            vector_memory;
    typename VectorMemory<typename Range::BlockType>::Pointer tmp(
      vector_memory);

    if (m == 0)
      return;

    diagonal_inverse.block(m - 1, m - 1)
      .vmult_add(v.block(m - 1), u.block(m - 1));

    for (int i = m - 2; i >= 0; --i)
      {
        diagonal_inverse.block(i, i).reinit_range_vector(
          *tmp, /*bool omit_zeroing_entries=*/true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j = i + 1; j < m; ++j)
          block_operator.block(i, j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        diagonal_inverse.block(i, i).vmult_add(v.block(i), *tmp);
      }
  };

  return return_op;
}

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
