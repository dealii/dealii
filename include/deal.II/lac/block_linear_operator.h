// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

#ifndef dealii__block_linear_operator_h
#define dealii__block_linear_operator_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/linear_operator.h>

#ifdef DEAL_II_WITH_CXX11

DEAL_II_NAMESPACE_OPEN

// Forward declarations:

template <typename Number> class BlockVector;

template <typename Range = BlockVector<double>,
          typename Domain = Range>
class BlockLinearOperator;

template <typename Range = BlockVector<double>,
          typename Domain = Range,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain>
block_operator(const BlockMatrixType &matrix);

template <size_t m, size_t n,
          typename Range = BlockVector<double>,
          typename Domain = Range>
BlockLinearOperator<Range, Domain>
block_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, m> &);

template <size_t m,
          typename Range = BlockVector<double>,
          typename Domain = Range>
BlockLinearOperator<Range, Domain>
block_diagonal_operator(const std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, m> &);

template <size_t m,
          typename Range = BlockVector<double>,
          typename Domain = Range>
BlockLinearOperator<Range, Domain>
block_diagonal_operator(const LinearOperator<typename Range::BlockType, typename Domain::BlockType> &op);

// This is a workaround for a bug in <=gcc-4.7 that does not like partial
// template default values in combination with local lambda expressions [1]
//
// [1] https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53624
//
// Forward declare functions with partial template defaults:

template <typename Range = BlockVector<double>,
          typename Domain = Range,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain>
block_diagonal_operator(const BlockMatrixType &block_matrix);

template <typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Domain, Range>
block_forward_substitution(const BlockLinearOperator<Range, Domain> &,
                           const BlockLinearOperator<Domain, Range> &);

template <typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Domain, Range>
block_back_substitution(const BlockLinearOperator<Range, Domain> &,
                        const BlockLinearOperator<Domain, Range> &);

// end of workaround



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
 * A BlockLinearOperator can be sliced to a LinearOperator at any time.
 * This removes all information about the underlying block structure
 * (beacuse above <code>std::function</code> objects are no longer
 * available) - the linear operator interface, however, remains intact.
 *
 * @note This class makes heavy use of <code>std::function</code> objects
 * and lambda functions. This flexibiliy comes with a run-time penalty.
 * Only use this object to encapsulate object with medium to large
 * individual block sizes, and small block structure (as a rule of thumb,
 * matrix blocks greater than $1000\times1000$).
 *
 * @note This class is only available if deal.II was configured with C++11
 * support, i.e., if <code>DEAL_II_WITH_CXX11</code> is enabled during cmake
 * configure.
 *
 * @author Matthias Maier, 2015
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
class BlockLinearOperator : public LinearOperator<Range, Domain>
{
public:

  typedef LinearOperator<typename Range::BlockType, typename Domain::BlockType> BlockType;

  /**
   * Create an empty BlockLinearOperator object.
   *
   * All<code>std::function</code> member objects of this class and its
   * base class LinearOperator are initialized with default variants that
   * throw an exception upon invocation.
   */
  BlockLinearOperator()
    : LinearOperator<Range, Domain>()
  {

    n_block_rows = []() -> unsigned int
    {
      Assert(false, ExcMessage("Uninitialized BlockLinearOperator<Range, Domain>::n_block_rows called"));
      return 0;
    };

    n_block_cols = []() -> unsigned int
    {
      Assert(false, ExcMessage("Uninitialized BlockLinearOperator<Range, Domain>::n_block_cols called"));
      return 0;
    };

    block = [](unsigned int, unsigned int) -> BlockType
    {
      Assert(false, ExcMessage("Uninitialized BlockLinearOperator<Range, Domain>::block called"));
      return BlockType();
    };
  }

  /**
   * Default copy constructor.
   */
  BlockLinearOperator(const BlockLinearOperator<Range, Domain> &) =
    default;

  /**
   * Templated copy constructor that creates a BlockLinearOperator object
   * from an object @p op for which the conversion function
   * <code>block_operator</code> is defined.
   */
  template<typename Op>
  BlockLinearOperator(const Op &op)
  {
    *this = block_operator<Range, Domain, Op>(op);
  }

  /**
   * Create a BlockLinearOperator from a two-dimensional array @p ops of
   * LinearOperator. This constructor calls the corresponding
   * block_operator() specialization.
   */
  template<size_t m, size_t n>
  BlockLinearOperator(const std::array<std::array<BlockType, n>, m> &ops)
  {
    *this = block_operator<m, n, Range, Domain>(ops);
  }

  /**
   * Create a block-diagonal BlockLinearOperator from a one-dimensional
   * array @p ops of LinearOperator. This constructor calls the
   * corresponding block_operator() specialization.
   */
  template<size_t m>
  BlockLinearOperator(const std::array<BlockType, m> &ops)
  {
    *this = block_diagonal_operator<m, Range, Domain>(ops);
  }

  /**
   * Default copy assignment operator.
   */
  BlockLinearOperator<Range, Domain> &
  operator=(const BlockLinearOperator<Range, Domain> &) = default;

  /**
   * Templated copy assignment operator for an object @p op for which the
   * conversion function <code>block_operator</code> is defined.
   */
  template <typename Op>
  BlockLinearOperator<Range, Domain> &operator=(const Op &op)
  {
    *this = block_operator<Range, Domain, Op>(op);
    return *this;
  }

  /**
   * Copy assignment from a two-dimensional array @p ops of LinearOperator.
   * This assignment operator calls the corresponding block_operator()
   * specialization.
   */
  template <size_t m, size_t n>
  BlockLinearOperator<Range, Domain> &
  operator=(const std::array<std::array<BlockType, n>, m> &ops)
  {
    *this = block_operator<m, n, Range, Domain>(ops);
    return *this;
  }

  /**
   * Copy assignment from a one-dimensional array @p ops of LinearOperator
   * that creates a block-diagonal BlockLinearOperator.
   * This assignment operator calls the corresponding block_operator()
   * specialization.
   */
  template <size_t m>
  BlockLinearOperator<Range, Domain> &
  operator=(const std::array<BlockType, m> &ops)
  {
    *this = block_diagonal_operator<m, Range, Domain>(ops);
    return *this;
  }

  /**
   * Return the number of blocks in a column (i.e, the number of "block
   * rows", or the number $m$, if interpreted as a $m\times n$ block
   * system).
   */
  std::function<unsigned int()> n_block_rows;

  /**
   * Return the number of blocks in a row (i.e, the number of "block
   * columns", or the number $n$, if interpreted as a $m\times n$ block
   * system).
   */
  std::function<unsigned int()> n_block_cols;

  /**
   * Access the block with the given coordinates. This
   * <code>std::function</code> object returns a LinearOperator
   * representing the $(i,j)$-th block of the BlockLinearOperator.
   */
  std::function<BlockType(unsigned int, unsigned int)> block;

  //@}
};



namespace internal
{
  namespace BlockLinearOperator
  {
    // Populate the LinearOperator interfaces with the help of the
    // BlockLinearOperator functions
    template <typename Range, typename Domain>
    inline void
    populate_linear_operator_functions(
      dealii::BlockLinearOperator<Range, Domain> &op)
    {
      op.reinit_range_vector = [=](Range &v, bool omit_zeroing_entries)
      {
        const unsigned int m = op.n_block_rows();

        // Reinitialize the block vector to m blocks:
        v.reinit(m);

        // And reinitialize every individual block with reinit_range_vectors:
        for (unsigned int i = 0; i < m; ++i)
          op.block(i, 0).reinit_range_vector(v.block(i), omit_zeroing_entries);

        v.collect_sizes();
      };

      op.reinit_domain_vector = [=](Domain &v, bool omit_zeroing_entries)
      {
        const unsigned int n = op.n_block_cols();

        // Reinitialize the block vector to n blocks:
        v.reinit(n);

        // And reinitialize every individual block with reinit_domain_vectors:
        for (unsigned int i = 0; i < n; ++i)
          op.block(0, i).reinit_domain_vector(v.block(i), omit_zeroing_entries);

        v.collect_sizes();
      };

      op.vmult = [=](Range &v, const Domain &u)
      {
        const unsigned int m = op.n_block_rows();
        const unsigned int n = op.n_block_cols();
        Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
        Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

        for (unsigned int i = 0; i < m; ++i)
          {
            op.block(i, 0).vmult(v.block(i), u.block(0));
            for (unsigned int j = 1; j < n; ++j)
              op.block(i, j).vmult_add(v.block(i), u.block(j));
          }
      };

      op.vmult_add = [=](Range &v, const Domain &u)
      {
        const unsigned int m = op.n_block_rows();
        const unsigned int n = op.n_block_cols();
        Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
        Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

        for (unsigned int i = 0; i < m; ++i)
          for (unsigned int j = 0; j < n; ++j)
            op.block(i, j).vmult_add(v.block(i), u.block(j));
      };

      op.Tvmult = [=](Domain &v, const Range &u)
      {
        const unsigned int n = op.n_block_cols();
        const unsigned int m = op.n_block_rows();
        Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
        Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

        for (unsigned int i = 0; i < n; ++i)
          {
            op.block(0, i).Tvmult(v.block(i), u.block(0));
            for (unsigned int j = 1; j < m; ++j)
              op.block(j, i).Tvmult_add(v.block(i), u.block(j));
          }
      };

      op.Tvmult_add = [=](Domain &v, const Range &u)
      {
        const unsigned int n = op.n_block_cols();
        const unsigned int m = op.n_block_rows();
        Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
        Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int j = 0; j < m; ++j)
            op.block(j, i).Tvmult_add(v.block(i), u.block(j));
      };
    }
  } /*namespace BlockLinearOperator*/
} /*namespace internal*/



/**
 * @name Creation of a BlockLinearOperator
 */
//@{

/**
 * @relates BlockLinearOperator
 *
 * A function that encapsulates a @p block_matrix into a
 * BlockLinearOperator.
 *
 * All changes made on the block structure and individual blocks of @p
 * block_matrix after the creation of the BlockLinearOperator object are
 * reflected by the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename Domain,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain>
block_operator(const BlockMatrixType &block_matrix)
{
  typedef typename BlockLinearOperator<Range, Domain>::BlockType BlockType;

  BlockLinearOperator<Range, Domain> return_op;

  return_op.n_block_rows = [&block_matrix]() -> unsigned int
  {
    return block_matrix.n_block_rows();
  };

  return_op.n_block_cols = [&block_matrix]() -> unsigned int
  {
    return block_matrix.n_block_cols();
  };

  return_op.block = [&block_matrix](unsigned int i, unsigned int j) -> BlockType
  {
#ifdef DEBUG
    const unsigned int m = block_matrix.n_block_rows();
    const unsigned int n = block_matrix.n_block_cols();
    Assert(i < m, ExcIndexRange (i, 0, m));
    Assert(j < n, ExcIndexRange (j, 0, n));
#endif

    return BlockType(block_matrix.block(i, j));
  };

  internal::BlockLinearOperator::populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relates BlockLinearOperator
 *
 * A variant of above function that encapsulates a given collection @p ops
 * of LinearOperators into a block structure. Here, it is assumed that
 * Range and Domain are blockvectors, i.e., derived from @ref
 * BlockVectorBase. The individual linear operators in @p ops must act on
 * the underlying vector type of the block vectors, i.e., on
 * Domain::BlockType yielding a result in Range::BlockType.
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
template <size_t m, size_t n, typename Range, typename Domain>
BlockLinearOperator<Range, Domain>
block_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, m> &ops)
{
  static_assert(m > 0 && n > 0,
                "a blocked LinearOperator must consist of at least one block");

  typedef typename BlockLinearOperator<Range, Domain>::BlockType BlockType;

  BlockLinearOperator<Range, Domain> return_op;

  return_op.n_block_rows = []() -> unsigned int
  {
    return m;
  };

  return_op.n_block_cols = []() -> unsigned int
  {
    return n;
  };

  return_op.block = [ops](unsigned int i, unsigned int j) -> BlockType
  {
    Assert(i < m, ExcIndexRange (i, 0, m));
    Assert(j < n, ExcIndexRange (j, 0, n));

    return ops[i][j];
  };

  internal::BlockLinearOperator::populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relates BlockLinearOperator
 *
 * This function extracts the diagonal blocks of @p block_matrix (either a
 * block matrix type or a BlockLinearOperator) and creates a
 * BlockLinearOperator with the diagonal. Off-diagonal elements are
 * initialized as null_operator (with correct reinit_range_vector and
 * reinit_domain_vector methods).
 *
 * All changes made on the individual diagonal blocks of @p block_matrix
 * after the creation of the BlockLinearOperator object are reflected by
 * the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename Domain,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain>
block_diagonal_operator(const BlockMatrixType &block_matrix)
{
  typedef typename BlockLinearOperator<Range, Domain>::BlockType BlockType;

  BlockLinearOperator<Range, Domain> return_op;

  return_op.n_block_rows = [&block_matrix]() -> unsigned int
  {
    return block_matrix.n_block_rows();
  };

  return_op.n_block_cols = [&block_matrix]() -> unsigned int
  {
    return block_matrix.n_block_cols();
  };

  return_op.block = [&block_matrix](unsigned int i, unsigned int j) -> BlockType
  {
#ifdef DEBUG
    const unsigned int m = block_matrix.n_block_rows();
    const unsigned int n = block_matrix.n_block_cols();
    Assert(m == n, ExcDimensionMismatch(m, n));
    Assert(i < m, ExcIndexRange (i, 0, m));
    Assert(j < n, ExcIndexRange (j, 0, n));
#endif
    if (i == j)
      return BlockType(block_matrix.block(i, j));
    else
      return null_operator(BlockType(block_matrix.block(i, j)));
  };

  internal::BlockLinearOperator::populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relates BlockLinearOperator
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
template <size_t m, typename Range, typename Domain>
BlockLinearOperator<Range, Domain>
block_diagonal_operator(const std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, m> &ops)
{
  static_assert(m > 0,
                "a blockdiagonal LinearOperator must consist of at least one block");

  typedef typename BlockLinearOperator<Range, Domain>::BlockType BlockType;

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

  return block_operator<m,m,Range,Domain>(new_ops);
}



/**
 * @relates BlockLinearOperator
 *
 * A variant of above function that only takes a single LinearOperator
 * argument @p op and creates a blockdiagonal linear operator with @p m copies
 * of it.
 *
 * @ingroup LAOperators
 */
template <size_t m, typename Range, typename Domain>
BlockLinearOperator<Range, Domain>
block_diagonal_operator(const LinearOperator<typename Range::BlockType, typename Domain::BlockType> &op)
{
  static_assert(m > 0,
                "a blockdiagonal LinearOperator must consist of at least "
                "one block");

  typedef typename BlockLinearOperator<Range, Domain>::BlockType BlockType;
  std::array<BlockType, m> new_ops;
  new_ops.fill(op);

  return block_diagonal_operator(new_ops);
}



//@}
/**
 * @name Manipulation of a BlockLinearOperator
 */
//@{

/**
 * @relates LinearOperator
 * @relates BlockLinearOperator
 *
 * This function implements forward substitution to invert a lower block
 * triangular matrix. As arguments, it takes a BlockLinearOperator
 * @p block_operator representing a block lower triangular matrix, as well as
 * a BlockLinearOperator @p diagonal_inverse representing inverses of
 * diagonal blocks of @p block_operator.
 *
 * Let us assume we have a linear system with the following block structure:
 *
 * A00 x0 + ...                   = y0
 * A01 x0 + A11 x1 + ...          = y1
 * ...        ...
 * A0n x0 + A1n x1 + ... + Ann xn = yn
 *
 * First of all, x0 = A00^-1 y0
 * Then, we can use x0 to recover x1:
 *    x1 = A11^-1 ( y1 - A01 x0 )
 * and therefore:
 *    xn = Ann^-1 ( yn - A0n x0 - ... - A(n-1)n x(n-1) )
 *
 * @note We are not using all blocks of the BlockLinearOperator arguments:
 * Just the lower triangular block matrix of @p block_operator is used as
 * well as the diagonal of @p diagonal_inverse.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
LinearOperator<Domain, Range>
block_forward_substitution(const BlockLinearOperator<Range, Domain> &block_operator,
                           const BlockLinearOperator<Domain, Range> &diagonal_inverse)
{
  LinearOperator<Range, Range> return_op;

  return_op.reinit_range_vector = diagonal_inverse.reinit_range_vector;
  return_op.reinit_domain_vector = diagonal_inverse.reinit_domain_vector;

  return_op.vmult = [block_operator, diagonal_inverse](Range &v, const Range &u)
  {
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
        dst = u.block(i);
        dst *= -1.;
        for (unsigned int j = 0; j < i; ++j)
          block_operator.block(i, j).vmult_add(dst, v.block(j));
        dst *= -1.;
        diagonal_inverse.block(i, i).vmult(dst, dst); // uses intermediate storage
      }
  };

  return_op.vmult_add = [block_operator, diagonal_inverse](Range &v, const Range &u)
  {
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

    static GrowingVectorMemory<typename  Range::BlockType> vector_memory;
    typename Range::BlockType *tmp = vector_memory.alloc();

    diagonal_inverse.block(0, 0).vmult_add(v.block(0), u.block(0));

    for (unsigned int i = 1; i < m; ++i)
      {
        diagonal_inverse.block(i, i).reinit_range_vector(*tmp, /*bool omit_zeroing_entries=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j = 0; j < i; ++j)
          block_operator.block(i, j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        diagonal_inverse.block(i, i).vmult_add(v.block(i),*tmp);
      }

    vector_memory.free(tmp);
  };

  return return_op;
}



/**
 * @relates LinearOperator
 * @relates BlockLinearOperator
 *
 * This function implements back substitution to invert an upper block
 * triangular matrix. As arguments, it takes a BlockLinearOperator
 * @p block_operator representing an upper block triangular matrix, as well
 * as a BlockLinearOperator @p diagonal_inverse representing inverses of
 * diagonal blocks of @p block_operator.
 *
 * Let us assume we have a linear system with the following block structure:
 *
 * A00 x0 + A01 x1 + ... + A0n xn = yn
 *          A11 x1 + ...          = y1
 *                          ...     ..
 *                         Ann xn = yn
 *
 * First of all, xn = Ann^-1 yn
 * Then, we can use xn to recover x(n-1):
 *    x(n-1) = A(n-1)(n-1)^-1 ( y(n-1) - A(n-1)n x(n-1) )
 * and therefore:
 *    x0 = A00^-1 ( y0 - A0n xn - ... - A01 x1 )
 *
 * @note We are not using all blocks of the BlockLinearOperator arguments:
 * Just the upper triangular block matrix of @p block_operator is used as
 * well as the diagonal of @p diagonal_inverse.

 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
LinearOperator<Domain, Range>
block_back_substitution(const BlockLinearOperator<Range, Domain> &block_operator,
                        const BlockLinearOperator<Domain, Range> &diagonal_inverse)
{
  LinearOperator<Range, Range> return_op;

  return_op.reinit_range_vector = diagonal_inverse.reinit_range_vector;
  return_op.reinit_domain_vector = diagonal_inverse.reinit_domain_vector;

  return_op.vmult = [block_operator, diagonal_inverse](Range &v, const Range &u)
  {
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

    diagonal_inverse.block(m-1, m-1).vmult(v.block(m-1),u.block(m-1));

    for (int i = m - 2; i >= 0; --i)
      {
        auto &dst = v.block(i);
        dst = u.block(i);
        dst *= -1.;
        for (unsigned int j = i + 1; j < m; ++j)
          block_operator.block(i, j).vmult_add(dst, v.block(j));
        dst *= -1.;
        diagonal_inverse.block(i, i).vmult(dst, dst); // uses intermediate storage
      }
  };

  return_op.vmult_add = [block_operator, diagonal_inverse](Range &v, const Range &u)
  {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));
    static GrowingVectorMemory<typename  Range::BlockType> vector_memory;
    typename  Range::BlockType *tmp = vector_memory.alloc();

    if (m == 0)
      return;

    diagonal_inverse.block(m-1, m-1).vmult_add(v.block(m-1),u.block(m-1));

    for (int i = m - 2; i >= 0; --i)
      {
        diagonal_inverse.block(i, i).reinit_range_vector(*tmp, /*bool omit_zeroing_entries=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j = i + 1; j < m; ++j)
          block_operator.block(i, j).vmult_add(*tmp,v.block(j));
        *tmp *= -1.;
        diagonal_inverse.block(i, i).vmult_add(v.block(i),*tmp);
      }

    vector_memory.free(tmp);
  };

  return return_op;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_CXX11
#endif
