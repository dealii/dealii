// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_householder_h
#define dealii_householder_h


#include <deal.II/base/config.h>

#include <deal.II/lac/full_matrix.h>

#include <cmath>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
#endif

/**
 * @addtogroup Matrix2
 * @{
 */


/**
 * QR-decomposition of a full matrix.
 *
 * This class computes the QR-decomposition of given matrix by the
 * Householder algorithm. Then, the function least_squares() can be
 * used to compute the vector $x$ minimizing $\|Ax-b\|$ for a given
 * vector $b$. The QR decomposition of $A$ is useful for this purpose
 * because the minimizer is given by the equation
 * $x=(A^TA)^{-1}A^Tb=(R^TQ^TQR)^{-1}R^TQ^Tb$ which is easy to compute
 * because $Q$ is an orthogonal matrix, and consequently
 * $Q^TQ=I$. Thus,
 * $x=(R^TR)^{-1}R^TQ^Tb=R^{-1}R^{-T}R^TQ^Tb=R^{-1}Q^Tb$. Furthermore,
 * $R$ is triangular, so applying $R^{-1}$ to a vector only involves a
 * backward or forward solve.
 *
 *
 * <h3>Implementation details</h3>
 *
 * The class does not in fact store the $Q$ and $R$ factors explicitly
 * as matrices. It does store $R$, but the $Q$ factor is stored as the
 * product of Householder reflections of the form $Q_i = I-v_i v_i^T$
 * where the vectors $v_i$ are so that they can be stored in the
 * lower-triangular part of an underlying matrix object, whereas $R$
 * is stored in the upper triangular part.
 *
 * The $v_i$ vectors and the $R$ matrix now are in conflict because they
 * both want to use the diagonal entry of the matrix, but we can only
 * store one in these positions, of course. Consequently, the entries
 * $(v_i)_i$ are stored separately in the `diagonal` member variable.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on
 * @ref Instantiations
 * in the manual).
 */
template <typename number>
class Householder
{
public:
  /**
   * Declare type of container size type.
   */
  using size_type = types::global_dof_index;

  /**
   * Create an empty object.
   */
  Householder() = default;

  /**
   * Create an object holding the QR-decomposition of the matrix $A$.
   */
  template <typename number2>
  Householder(const FullMatrix<number2> &A);

  /**
   * Compute the QR-decomposition of the given matrix $A$.
   *
   * This overwrites any previously computed QR decomposition.
   */
  template <typename number2>
  void
  initialize(const FullMatrix<number2> &A);

  /**
   * Solve the least-squares problem for the right hand side <tt>src</tt>. The
   * returned scalar value is the Euclidean norm of the approximation error.
   *
   * @arg @c dst contains the solution of the least squares problem on return.
   *
   * @arg @c src contains the right hand side <i>b</i> of the least squares
   * problem. It will be changed during the algorithm and is unusable on
   * return.
   */
  template <typename number2>
  number2
  least_squares(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * This function does the same as the previous one, but for BlockVectors.
   */
  template <typename number2>
  double
  least_squares(BlockVector<number2>       &dst,
                const BlockVector<number2> &src) const;

  /**
   * A wrapper to least_squares(), implementing the standard MatrixType
   * interface.
   */
  template <typename VectorType>
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * A wrapper to least_squares() that implements multiplication with
   * the transpose matrix.
   */
  template <typename VectorType>
  void
  Tvmult(VectorType &dst, const VectorType &src) const;


private:
  /**
   * Storage for the diagonal elements of the orthogonal
   * transformation. See the class documentation for more information.
   */
  std::vector<number> diagonal;

  /**
   * Storage that is internally used for the Householder transformation.
   */
  FullMatrix<number> storage;
};

/** @} */

#ifndef DOXYGEN
/*-------------------------Inline functions -------------------------------*/

// QR-transformation cf. Stoer 1 4.8.2 (p. 191)

template <typename number>
template <typename number2>
void
Householder<number>::initialize(const FullMatrix<number2> &M)
{
  const size_type m = M.n_rows(), n = M.n_cols();
  storage.reinit(m, n);
  storage.fill(M);
  Assert(!storage.empty(), typename FullMatrix<number2>::ExcEmptyMatrix());
  diagonal.resize(m);

  // m > n, src.n() = m
  Assert(storage.n_cols() <= storage.n_rows(),
         ExcDimensionMismatch(storage.n_cols(), storage.n_rows()));

  for (size_type j = 0; j < n; ++j)
    {
      number2   sigma = 0;
      size_type i;
      // sigma = ||v||^2
      for (i = j; i < m; ++i)
        sigma += storage(i, j) * storage(i, j);
      // We are ready if the column is
      // empty. Are we?
      if (std::abs(sigma) < 1.e-15)
        return;

      number2 s;
      if constexpr (numbers::NumberTraits<number2>::is_complex)
        s = storage(j, j).real() < 0 ? std::sqrt(sigma) : -std::sqrt(sigma);
      else
        s = storage(j, j) < 0 ? std::sqrt(sigma) : -std::sqrt(sigma);
      //
      number2 beta = std::sqrt(1. / (sigma - s * storage(j, j)));

      // Make column j the Householder
      // vector, store first entry in
      // diagonal
      diagonal[j]   = beta * (storage(j, j) - s);
      storage(j, j) = s;

      for (i = j + 1; i < m; ++i)
        storage(i, j) *= beta;


      // For all subsequent columns do
      // the Householder reflection
      for (size_type k = j + 1; k < n; ++k)
        {
          number2 sum = diagonal[j] * storage(j, k);
          for (i = j + 1; i < m; ++i)
            sum += storage(i, j) * storage(i, k);

          storage(j, k) -= sum * this->diagonal[j];
          for (i = j + 1; i < m; ++i)
            storage(i, k) -= sum * storage(i, j);
        }
    }
}



template <typename number>
template <typename number2>
Householder<number>::Householder(const FullMatrix<number2> &M)
{
  initialize(M);
}



template <typename number>
template <typename number2>
number2
Householder<number>::least_squares(Vector<number2>       &dst,
                                   const Vector<number2> &src) const
{
  Assert(!storage.empty(), typename FullMatrix<number2>::ExcEmptyMatrix());
  AssertDimension(dst.size(), storage.n());
  AssertDimension(src.size(), storage.m());

  const size_type m = storage.m(), n = storage.n();

  Vector<number2> aux(src);
  // m > n, m = src.n, n = dst.n

  // Multiply Q_n ... Q_2 Q_1 src
  // Where Q_i = I - v_i v_i^T
  for (size_type j = 0; j < n; ++j)
    {
      // sum = v_i^T dst
      number2 sum = diagonal[j] * aux(j);
      for (size_type i = j + 1; i < m; ++i)
        sum += static_cast<number2>(storage(i, j)) * aux(i);
      // dst -= v * sum
      aux(j) -= sum * diagonal[j];
      for (size_type i = j + 1; i < m; ++i)
        aux(i) -= sum * static_cast<number2>(storage(i, j));
    }
  // Compute norm of residual
  number2 sum = 0.;
  for (size_type i = n; i < m; ++i)
    sum += aux(i) * aux(i);
  AssertIsFinite(sum);

  // Compute solution
  storage.backward(dst, aux);

  return std::sqrt(sum);
}



template <typename number>
template <typename number2>
double
Householder<number>::least_squares(BlockVector<number2>       &dst,
                                   const BlockVector<number2> &src) const
{
  Assert(!storage.empty(), typename FullMatrix<number2>::ExcEmptyMatrix());
  AssertDimension(dst.size(), storage.n());
  AssertDimension(src.size(), storage.m());

  const size_type m = storage.m(), n = storage.n();

  BlockVector<number2> aux;
  aux.reinit(src, true);
  aux = src;
  // m > n, m = src.n, n = dst.n

  // Multiply Q_n ... Q_2 Q_1 src
  // Where Q_i = I-v_i v_i^T
  for (size_type j = 0; j < n; ++j)
    {
      // sum = v_i^T dst
      number2 sum = diagonal[j] * aux(j);
      for (size_type i = j + 1; i < m; ++i)
        sum += storage(i, j) * aux(i);
      // dst -= v * sum
      aux(j) -= sum * diagonal[j];
      for (size_type i = j + 1; i < m; ++i)
        aux(i) -= sum * storage(i, j);
    }
  // Compute norm of residual
  number2 sum = 0.;
  for (size_type i = n; i < m; ++i)
    sum += *aux(i) * aux(i);
  AssertIsFinite(sum);

  // backward works for Vectors only, so copy them before
  Vector<number2> v_dst, v_aux;
  v_dst = dst;
  v_aux = aux;
  // Compute solution
  storage.backward(v_dst, v_aux);
  // copy the result back to the BlockVector
  dst = v_dst;

  return std::sqrt(sum);
}


template <typename number>
template <typename VectorType>
void
Householder<number>::vmult(VectorType &dst, const VectorType &src) const
{
  least_squares(dst, src);
}


template <typename number>
template <typename VectorType>
void
Householder<number>::Tvmult(VectorType &, const VectorType &) const
{
  DEAL_II_NOT_IMPLEMENTED();
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
