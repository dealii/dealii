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

/**
 * @name Creation of LinearOperator block structures
 */
//@{

/**
 * @relates LinearOperator
 *
 * A function that encapsulates a given collection @p ops of
 * LinearOperators into a block structure. Hereby, it is assumed that Range
 * and Domain are blockvectors, i.e., derived from @ref BlockVectorBase.
 * The individual linear operators in @p ops must act on a the underlying
 * vector type of the block vectors, i.e., on Domain::BlockType yielding a
 * result in Range::BlockType.
 *
 * The list @p ops is best passed as an initializer list. Consider for
 * example a linear operator block (acting on Vector<double>)
 * @code
 *  op_a00 | op_a01
 *         |
 *  ---------------
 *         |
 *  op_a10 | op_a11
 * @endcode
 * The coresponding block_operator invocation takes the form
 * @code
 * block_operator<2, 2, BlockVector<double>>({op_a00, op_a01, op_a10, op_a11});
 * @endcode
 *
 * @ingroup LAOperators
 */
template <unsigned int m, unsigned int n,
          typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Range, Domain>
block_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, m> &ops)
{
  static_assert(m > 0 && n > 0,
                "a blocked LinearOperator must consist of at least one block");

  LinearOperator<Range, Domain> return_op;

  return_op.reinit_range_vector = [ops](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(m);

    // And reinitialize every individual block with reinit_range_vectors:
    for (unsigned int i = 0; i < m; ++i)
      ops[i][0].reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [ops](Domain &v, bool fast)
  {
    // Reinitialize the block vector to n blocks:
    v.reinit(n);

    // And reinitialize every individual block with reinit_domain_vectors:
    for (unsigned int i = 0; i < n; ++i)
      ops[0][i].reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [ops](Range &v, const Domain &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

    for (unsigned int i = 0; i < m; ++i)
      {
        ops[i][0].vmult(v.block(i), u.block(0));
        for (unsigned int j = 1; j < n; ++j)
          ops[i][j].vmult_add(v.block(i), u.block(j));
      }
  };

  return_op.vmult_add = [ops](Range &v, const Domain &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

    for (unsigned int i = 0; i < m; ++i)
      for (unsigned int j = 0; j < n; ++j)
        ops[i][j].vmult_add(v.block(i), u.block(j));
  };

  return_op.Tvmult = [ops](Domain &v, const Range &u)
  {
    Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < n; ++i)
      {
        ops[0][i].Tvmult(v.block(i), u.block(0));
        for (unsigned int j = 1; j < m; ++j)
          ops[j][i].Tvmult_add(v.block(i), u.block(j));
      }
  };

  return_op.Tvmult_add = [ops](Domain &v, const Range &u)
  {
    Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < m; ++j)
        ops[j][i].Tvmult_add(v.block(i), u.block(j));
  };

  return return_op;
}


/**
 * @relates LinearOperator
 *
 * A variant of above function that builds up a block diagonal linear
 * operator from an array @p ops of diagonal elements (off-diagonal blocks
 * are assumed to be 0).
 *
 * The list @p ops is best passed as an initializer list. Consider for
 * example a linear operator block (acting on Vector<double>)
 * <code>diag(op_a0, op_a1, ..., op_am)</code>. The coresponding
 * block_operator invocation takes the form
 * @code
 * block_diagonal_operator<m, BlockVector<double>>({op_00, op_a1, ..., op_am});
 * @endcode
 *
 * @ingroup LAOperators
 */
template <unsigned int m,
          typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Range, Domain>
block_diagonal_operator(const std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, m> &ops)
{
  static_assert(m > 0,
                "a blockdiagonal LinearOperator must consist of at least one block");

  LinearOperator<Range, Domain> return_op;

  return_op.reinit_range_vector = [ops](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(m);

    // And reinitialize every individual block with reinit_range_vectors:
    for (unsigned int i = 0; i < m; ++i)
      ops[i].reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [ops](Domain &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(m);

    // And reinitialize every individual block with reinit_domain_vectors:
    for (unsigned int i = 0; i < m; ++i)
      ops[i].reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [ops](Range &v, const Domain &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      ops[i].vmult(v.block(i), u.block(i));
  };

  return_op.vmult_add = [ops](Range &v, const Domain &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      ops[i].vmult_add(v.block(i), u.block(i));
  };

  return_op.Tvmult = [ops](Domain &v, const Range &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      ops[i].Tvmult(v.block(i), u.block(i));
  };

  return_op.Tvmult_add = [ops](Domain &v, const Range &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      ops[i].Tvmult_add(v.block(i), u.block(i));
  };

  return return_op;
}


/**
 * @relates LinearOperator
 *
 * A variant of above function that only takes a single LinearOperator
 * argument @p op and creates a blockdiagonal linear operator with @p m
 * copies of it.
 *
 * @ingroup LAOperators
 */
template <unsigned int m,
          typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Range, Domain>
block_diagonal_operator(const LinearOperator<typename Range::BlockType, typename Domain::BlockType> &op)
{

  static_assert(m > 0,
                "a blockdiagonal LinearOperator must consist of at least "
                "one block");

  LinearOperator<Range, Domain> return_op;

  return_op.reinit_range_vector = [op](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(m);

    // And reinitialize every individual block with reinit_range_vectors:
    for (unsigned int i = 0; i < m; ++i)
      op.reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [op](Domain &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(m);

    // And reinitialize every individual block with reinit_domain_vectors:
    for (unsigned int i = 0; i < m; ++i)
      op.reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [op](Range &v, const Domain &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      op.vmult(v.block(i), u.block(i));
  };

  return_op.vmult_add = [op](Range &v, const Domain &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      op.vmult_add(v.block(i), u.block(i));
  };

  return_op.Tvmult = [op](Domain &v, const Range &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      op.Tvmult(v.block(i), u.block(i));
  };

  return_op.Tvmult_add = [op](Domain &v, const Range &u)
  {
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    for (unsigned int i = 0; i < m; ++i)
      op.Tvmult_add(v.block(i), u.block(i));
  };

  return return_op;
}


/**
 * @relates LinearOperator
 *
 * A function that takes an array of arrays of LinearOperators @p block_matrix
 * and returns its associated lower triangular matrix operator
 * (diagonal is not included).
 *
 * @code
 *  a00 | a01 | a02                 |     |
 *  ---------------             ---------------
 *  a10 | a11 | a12     ->      a10 |     |
 *  ---------------             ---------------
 *  a20 | a21 | a22             a20 | a21 |
 * @endcode
 *
 * @ingroup LAOperators
 */

// This is a workaround for a bug in <=gcc-4.7 that does not like partial
// template default values in function definitions in combination with
// local lambda expressions [1] in the function body. As a workaround
// declare the function with all default types and parameters first such
// that the function definition is without default types and parameters.
//
// [1] https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53624

template <unsigned int n,
          typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Range, Domain>
lower_triangular_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, n> &);


template <unsigned int n, typename Range, typename Domain>
LinearOperator<Range, Domain>
lower_triangular_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, n> &block_matrix)
{
  LinearOperator<Range, Domain> return_op;

  return_op.reinit_range_vector = [block_matrix](Range &v, bool fast)
  {
    v.reinit(n);

    for (unsigned int i = 0; i < n; ++i)
      block_matrix[i][i].reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [block_matrix](Domain &v, bool fast)
  {
    v.reinit(n);

    for (unsigned int i = 0; i < n; ++i)
      block_matrix[i][i].reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [block_matrix](Range &v, const Domain &u)
  {
    v.block(0) = 0;
    for (unsigned int i = 1; i < n; ++i)
      {
        block_matrix[i][0].vmult(v.block(i), u.block(0));
        for (unsigned int j = 1; j < i; ++j)
          block_matrix[i][j].vmult_add(v.block(i), u.block(j));
      }
  };

  return_op.vmult_add = [block_matrix](Range &v, const Domain &u)
  {
    for (unsigned int i = 1; i < n; ++i)
      {
        block_matrix[i][0].vmult_add(v.block(i), u.block(0));
        for (unsigned int j = 1; j < i; ++j)
          block_matrix[i][j].vmult_add(v.block(i), u.block(j));
      }
  };

  return_op.Tvmult = [block_matrix](Domain &v, const Range &u)
  {
    for (unsigned int i = 0; i < n; ++i)
      {
        v.block(i) = 0;
        for (unsigned int j = i+1; j < n; ++j)
          block_matrix[j][i].Tvmult_add(v.block(i), u.block(j));
      }
  };

  return_op.Tvmult_add = [block_matrix](Domain &v, const Range &u)
  {
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = i+1; j < n; ++j)
        block_matrix[j][i].Tvmult_add(v.block(i), u.block(j));
  };

  return return_op;
}

/**
 * @relates LinearOperator
 *
 * This function is a specification of the above function that
 * allows to work with block matrices @p block_matrix .
 *
 * @ingroup LAOperators
 */

template <unsigned int n,
          typename Range = BlockVector<double>,
          typename Domain = Range,
          typename BlockMatrix>
LinearOperator<Range, Domain>
lower_triangular_operator(const BlockMatrix &);


template <unsigned int n, typename Range, typename Domain, typename BlockMatrix>
LinearOperator<Range, Domain>
lower_triangular_operator(const BlockMatrix &block_matrix)
{
  Assert(block_matrix.n_block_rows() == block_matrix.n_block_cols(),
         ExcDimensionMismatch(block_matrix.n_block_rows(),block_matrix.n_block_cols()) );

  std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, n> M;
  for (unsigned int i = 0; i<n; ++i)
    {
      for (unsigned int j = 0; j<i; ++j)
        M[i][j] = linear_operator<typename Range::BlockType, typename Domain::BlockType>(block_matrix.block(i,j));
      M[i][i] = null_operator(linear_operator<typename Range::BlockType, typename Domain::BlockType>(block_matrix.block(i,i)).reinit_range_vector);
    }
  return lower_triangular_operator<n, Range, Domain>(M);
}

/**
 * @relates LinearOperator
 *
 * A function that takes an array of arrays of LinearOperators @p block_matrix
 * and returns its associated upper triangular matrix operator
 * (diagonal is not included).
 *
 * @code
 *  a00 | a01 | a02                 | a01 | a02
 *  ---------------             ---------------
 *  a10 | a11 | a12     ->          |     | a12
 *  ---------------             ---------------
 *  a20 | a21 | a22                 |     |
 * @endcode
 *
 * @ingroup LAOperators
 */

// This is a workaround for a bug in <=gcc-4.7 that does not like partial
// template default values in function definitions in combination with
// local lambda expressions [1] in the function body. As a workaround
// declare the function with all default types and parameters first such
// that the function definition is without default types and parameters.
//
// [1] https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53624

template <unsigned int n,
          typename Range = BlockVector<double>,
          typename Domain = Range>
LinearOperator<Range, Domain>
upper_triangular_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, n> &);


template <unsigned int n, typename Range, typename Domain>
LinearOperator<Range, Domain>
upper_triangular_operator(const std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, n> &block_matrix)
{
  LinearOperator<Range, Domain> return_op;

  return_op.reinit_range_vector = [block_matrix](Range &v, bool fast)
  {
    v.reinit(n);

    for (unsigned int i = 0; i < n; ++i)
      block_matrix[i][i].reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [block_matrix](Domain &v, bool fast)
  {
    v.reinit(n);

    for (unsigned int i = 0; i < n; ++i)
      block_matrix[i][i].reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [block_matrix](Range &v, const Domain &u)
  {
    for (unsigned int i = 0; i < n - 1; ++i)
      {
        block_matrix[i][n - 1].vmult(v.block(i), u.block(n - 1));
        for (unsigned int j = n - 2; j > i; --j)
          block_matrix[i][j].vmult_add(v.block(i), u.block(j));
      }
    v.block(n - 1) = 0;
  };

  return_op.vmult_add = [block_matrix](Range &v, const Domain &u)
  {
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = i + 1; j < n; ++j)
        block_matrix[i][j].vmult_add(v.block(i), u.block(j));
  };

  return_op.Tvmult = [block_matrix](Domain &v, const Range &u)
  {
    for (unsigned int i = 0; i < n; ++i)
      {
        v.block(i) = 0;
        for (unsigned int j = 0; j < i; ++j)
          block_matrix[j][i].Tvmult_add(v.block(i), u.block(j));
      }
  };

  return_op.Tvmult_add = [block_matrix](Domain &v, const Range &u)
  {
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < i; ++j)
        block_matrix[j][i].Tvmult_add(v.block(i), u.block(j));
  };

  return return_op;
}


/**
 * @relates LinearOperator
 *
 * This function is a specification of the above function that
 * allows to work with block matrices @p block_matrix .
 *
 * @ingroup LAOperators
 */

template <unsigned int n,
          typename Range = BlockVector<double>,
          typename Domain = Range,
          typename BlockMatrix>
LinearOperator<Range, Domain>
upper_triangular_operator(const BlockMatrix &);


template <unsigned int n, typename Range, typename Domain, typename BlockMatrix>
LinearOperator<Range, Domain>
upper_triangular_operator(const BlockMatrix &block_matrix)
{
  Assert(block_matrix.n_block_rows() == block_matrix.n_block_cols(),
         ExcDimensionMismatch(block_matrix.n_block_rows(),block_matrix.n_block_cols()) );

  std::array<std::array<LinearOperator<typename Range::BlockType, typename Domain::BlockType>, n>, n> M;
  for (unsigned int i = 0; i<n; ++i)
    {
      M[i][i] = null_operator(linear_operator<typename Range::BlockType, typename Domain::BlockType>(block_matrix.block(i,i)).reinit_range_vector);

      for (unsigned int j = i+1; j<n; ++j)
        M[i][j] = linear_operator<typename Range::BlockType, typename Domain::BlockType>(block_matrix.block(i,j));
    }
  return upper_triangular_operator<n, Range, Domain>(M);
}


/**
 * @relates LinearOperator
 *
 * Let M be a block matrix of the form Id + T made of nxn blocks and where
 * T is a lower / upper triangular (without diagonal). Then, its inverse is
 * of the form:
 * @code
 *  Id + sum_{i=1}^{n-1} (-1)^i T^i
 * @endcode
 * This formula can be used to invert all triangular matrices (diagonal
 * included of course).
 *
 * This function takes a block matrix @p block_matrix (possibly full) and a
 * linear block operator  @p inverse_diagonal made of the inverse of the
 * diagonal blocks inverses. The output is the inverse of the matrix in the
 * case of a triangular matrix and as inverse_diagonal its diagonal blocks
 * inverses. Otherwise, the result is a preconditioner.
 *
 * The parameter @p lower is a bool that allows to specify if we want to
 * use lower triangular part of @p block_matrix (true, this is the default
 * value) or to use upper triangular part of @p block_matrix (false).
 *
 * @ingroup LAOperators
 */

// workaround for a bug in <=gcc-4.7 that does not like partial template
// default values in combination with local lambda expressions [1]
// [1] https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53624
template <unsigned int n,
          typename Range = BlockVector<double>,
          typename Domain = Range,
          typename BlockMatrix>
LinearOperator<Range, Domain>
block_triangular_inverse(const BlockMatrix &,
                         const LinearOperator<Range, Domain> &,
                         bool lower = true);


template <unsigned int n, typename Range, typename Domain, typename BlockMatrix>
LinearOperator<Range, Domain>
block_triangular_inverse(const BlockMatrix &block_matrix,
                         const LinearOperator<Range, Domain> &inverse_diagonal,
                         bool lower)
{
  Assert(block_matrix.n_block_rows() == n,
         ExcDimensionMismatch(block_matrix.n_block_rows(), n));
  Assert(block_matrix.n_block_rows() == block_matrix.n_block_cols(),
         ExcDimensionMismatch(block_matrix.n_block_rows(),
                              block_matrix.n_block_cols()));

  LinearOperator<Range, Domain> op_a;

  if (lower)
    {
      op_a = lower_triangular_operator<n, Range, Domain, BlockMatrix>(block_matrix);
    }
  else
    {
      op_a = upper_triangular_operator<n, Range, Domain, BlockMatrix>(block_matrix);
    }

  auto id     = identity_operator(op_a.reinit_range_vector);
  auto result = identity_operator(op_a.reinit_range_vector);

  // Notice that the following formula is recursive. We are evaluating:
  // Id - T + T^2 - T^3 ... (- T)^n
  for (unsigned int i = 0; i < n - 1; ++i)
    result = id - inverse_diagonal * op_a * result;

  return result * inverse_diagonal;
}

/**
 * @relates LinearOperator
 *
 * This function implement a forward substitution argument to invert a lower
 * block triangular matrix.
 * It takes as argement a block matrix @p block_matrix and an array of LinearOperators
 * @p inverse_diagonal representing inverses of digonal blocks of @p block_matrix.
 *
 * Let us assume we have a linear system where each coefficient of the system is a
 * matrix:
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
 * Caveat: Tvmult and Tvmult_add have not been implemented, yet. This may lead to mistakes.
 * @ingroup LAOperators
 */

template <typename BlockMatrix,
          unsigned int n  = BlockMatrix::n ,
          typename Range  = typename BlockMatrix::BlockType>
LinearOperator<Range, Range>
block_forward_substitution(const BlockMatrix &block_matrix,
                           const std::array<LinearOperator<typename Range::BlockType, typename Range::BlockType>, n> &inverse_diagonal)
{
  LinearOperator<Range, Range> return_op;

  return_op.reinit_range_vector = [inverse_diagonal](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(n);

    // And reinitialize every individual block with reinit_range_vectors:
    for (unsigned int i = 0; i < n; ++i)
      inverse_diagonal[i].reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [inverse_diagonal](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(n);

    // And reinitialize every individual block with reinit_domain_vectors:
    for (unsigned int i = 0; i < n; ++i)
      inverse_diagonal[i].reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [&block_matrix, inverse_diagonal](Range &v, const Range &u)
  {
    static GrowingVectorMemory<typename  Range::BlockType> vector_memory;
    typename  Range::BlockType *tmp = vector_memory.alloc();

    inverse_diagonal[0].vmult(v.block(0), u.block(0));
    for (unsigned int i=1; i<n; ++i)
      {
        inverse_diagonal[i].reinit_range_vector(*tmp, /*bool fast=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j=0; j<i; ++j)
          block_matrix.block(i,j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        inverse_diagonal[i].vmult(v.block(i),*tmp);
      }

    vector_memory.free(tmp);
  };

  return_op.vmult_add = [&block_matrix, inverse_diagonal](Range &v, const Range &u)
  {
    static GrowingVectorMemory<typename  Range::BlockType> vector_memory;
    typename  Range::BlockType *tmp = vector_memory.alloc();

    inverse_diagonal[0].vmult_add(v.block(0), u.block(0));
    for (unsigned int i=1; i<n; ++i)
      {
        inverse_diagonal[i].reinit_range_vector(*tmp, /*bool fast=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j=0; j<i; ++j)
          block_matrix.block(i,j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        inverse_diagonal[i].vmult_add(v.block(i),*tmp);
      }

    vector_memory.free(tmp);
  };

  return return_op;
}

/**
 * @relates LinearOperator
 *
 * This function implement a back substitution argument to invert an upper
 * block triangular matrix.
 * It takes as argement a block matrix @p block_matrix and an array of LinearOperators
 * @p inverse_diagonal representing inverses of digonal blocks of @p block_matrix.
 *
 * Let us assume we have a linear system where each coefficient of the system is a
 * matrix:

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
 * Caveat: Tvmult and Tvmult_add have not been implemented, yet. This may lead to mistakes.
 * @ingroup LAOperators
 */

template <typename BlockMatrix,
          unsigned int n  = BlockMatrix::n ,
          typename Range  = typename BlockMatrix::BlockType>
LinearOperator<Range, Range>
block_back_substitution(const BlockMatrix &block_matrix,
                        const std::array<LinearOperator<typename Range::BlockType, typename Range::BlockType>, n> &inverse_diagonal)
{
  LinearOperator<Range, Range> return_op;

  return_op.reinit_range_vector = [inverse_diagonal](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(n);

    // And reinitialize every individual block with reinit_range_vectors:
    for (unsigned int i = 0; i < n; ++i)
      inverse_diagonal[i].reinit_range_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.reinit_domain_vector = [inverse_diagonal](Range &v, bool fast)
  {
    // Reinitialize the block vector to m blocks:
    v.reinit(n);

    // And reinitialize every individual block with reinit_domain_vectors:
    for (unsigned int i = 0; i < n; ++i)
      inverse_diagonal[i].reinit_domain_vector(v.block(i), fast);

    v.collect_sizes();
  };

  return_op.vmult = [&block_matrix, inverse_diagonal](Range &v, const Range &u)
  {
    static GrowingVectorMemory<typename  Range::BlockType> vector_memory;
    typename  Range::BlockType *tmp = vector_memory.alloc();

    inverse_diagonal[n-1].vmult(v.block(n-1),u.block(n-1));
    for (int i=n-2; i>=0; --i)
      {
        inverse_diagonal[i].reinit_range_vector(*tmp, /*bool fast=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (int j=i+1; j<n; ++j)
          block_matrix.block(i,j).vmult_add(*tmp,v.block(j));
        *tmp *= -1.;
        inverse_diagonal[i].vmult(v.block(i),*tmp);
      }

    vector_memory.free(tmp);
  };

  return_op.vmult_add = [&block_matrix, inverse_diagonal](Range &v, const Range &u)
  {
    static GrowingVectorMemory<typename  Range::BlockType> vector_memory;
    typename  Range::BlockType *tmp = vector_memory.alloc();

    inverse_diagonal[n-1].vmult_add(v.block(n-1),u.block(n-1));
    for (int i=n-2; i>=0; --i)
      {
        inverse_diagonal[i].reinit_range_vector(*tmp, /*bool fast=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (int j=i+1; j<n; ++j)
          block_matrix.block(i,j).vmult_add(*tmp,v.block(j));
        *tmp *= -1.;
        inverse_diagonal[i].vmult_add(v.block(i),*tmp);
      }

    vector_memory.free(tmp);
  };

  return return_op;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_CXX11
#endif
