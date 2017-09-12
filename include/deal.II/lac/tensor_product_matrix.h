// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii_tensor_product_matrix_h
#define dealii_tensor_product_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

DEAL_II_NAMESPACE_OPEN

template <typename> class Vector;
template <typename> class FullMatrix;

/**
 * This is a special matrix class defined as the tensor product (or Kronecker
 * product) of 1D matrices of the type
 * @f{align*}{
 * L &= A \otimes M + M \otimes A
 * @f}
 * in 2D and
 * @f{align*}{
 * L &= A \otimes M \otimes M + M \otimes A \otimes M + M \otimes M \otimes A
 * @f}
 * in 3D. The typical application setting is a discretization of the Laplacian
 * $L$ on a Cartesian (axis-aligned) geometry, where it can be exactly
 * represented by the Kronecker or tensor product of a 1D mass matrix $M$ and
 * a 1D Laplace matrix $A$ in each dimension. The dimension of the resulting
 * class is the product of the one-dimensional matrices.
 *
 * This class implements two basic operations, namely the usual multiplication
 * by a vector and the inverse. For both operations, fast tensorial techniques
 * can be applied that implement the operator evaluation in
 * $\text{size}(M)^{d+1}$ arithmetic operations, considerably less than
 * $\text{size}(M)^{2d}$ for the naive forward transformation and
 * $\text{size}(M)^{3d}$ for setting up the inverse of $L$.
 *
 * Interestingly, the exact inverse of the matrix $L$ can be found through
 * tensor products due to an article by <a
 * href="http://dl.acm.org/citation.cfm?id=2716130">R. E. Lynch, J. R. Rice,
 * D. H. Thomas, Direct solution of partial difference equations by tensor
 * product methods, Numerische Mathematik 6, 185-199</a> from 1964,
 * @f{align*}{
 * L^{-1} &= S \otimes S (\Lambda \otimes I + I \otimes \Lambda)^{-1}
 * S^\mathrm T \otimes S^\mathrm T,
 * @f}
 * where $S$ is the matrix of eigenvectors to the generalized eigenvalue problem
 * @f{align*}{
 * A s  &= \lambda M s,
 * @f}
 * and $\Lambda$ is the diagonal matrix representing the generalized
 * eigenvalues $\lambda$. Note that the vectors $s$ are such that they
 * simultaneously diagonalize $A$ and $M$, $S^{\mathrm T} A S = \Lambda$ and
 * $S^{\mathrm T} B S = I$. This method of matrix inversion is called fast
 * diagonalization method.
 *
 * This class requires LAPACK support.
 *
 * Note that this class allows for two modes of usage. The first is a use case
 * with run time constants for the matrix dimensions that is achieved by
 * setting the optional template parameter for the size to -1. The second mode
 * of usage that is faster allows to set the template parameter as a compile
 * time constant, giving significantly faster code in particular for small
 * sizes of the matrix.
 *
 * @note This class uses a temporary array for storing intermediate results
 * that is a class member. A mutex is used to protect access to this array and
 * ensure correct results. If several threads run parallel instances of this
 * class, it is recommended that each threads holds its own matrix version.
 *
 * @tparam dim Dimension of the problem. Currently, 1D, 2D, and 3D codes are
 * implemented.
 *
 * @tparam Number Type of the underlying array elements. Note that the
 * underlying LAPACK implementation supports only float and double numbers, so
 * only these two types are currently supported.
 *
 * @tparam size Compile-time array lengths. By default at -1, which means that
 * the run-time info stored in the matrices passed to the reinit()
 * function is used.
 *
 * @author Martin Kronbichler, 2017
 */
template <int dim, typename Number, int size = -1>
class TensorProductMatrixSymmetricSum
{
public:
  /**
   * Constructor.
   */
  TensorProductMatrixSymmetricSum();

  /**
   * Constructor that is equivalent to the previous constructor and
   * immediately calling reinit().
   */
  TensorProductMatrixSymmetricSum(const FullMatrix<Number> &mass_matrix,
                                  const FullMatrix<Number> &derivative_matrix);

  /**
   * Initializes the matrix to the given mass matrix $M$ and derivative matrix
   * $A$. Note that the current implementation requires $M$ to be symmetric
   * and positive definite and $A$ to be symmetric and invertible but not
   * necessarily positive defininte.
   */
  void reinit (const FullMatrix<Number> &mass_matrix,
               const FullMatrix<Number> &derivative_matrix);

  /**
   * Returns the number of rows of this matrix, given by the dim-th power of
   * the size of the 1D matrices passed to the constructor.
   */
  unsigned int m() const;

  /**
   * Returns the number of columns of this matrix, given by the dim-th power
   * of the size of the 1D matrices passed to the constructor.
   */
  unsigned int n() const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class.
   */
  void vmult (Vector<Number> &dst,
              const Vector<Number> &src) const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class. Same as the other
   * vmult() function, but operating on plain pointers rather than a vector
   * (no check of array bounds possible).
   */
  void vmult (Number *dst,
              const Number *src) const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class.
   */
  void apply_inverse (Vector<Number> &dst,
                      const Vector<Number> &src) const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class. Same as the other
   * apply_inverse() function, but operating on plain pointers rather than a
   * vector (no check of array bounds possible).
   */
  void apply_inverse (Number *dst,
                      const Number *src) const;

private:
  /**
   * A copy of the @p mass_matrix object passed to the reinit() method.
   */
  FullMatrix<Number> mass_matrix;

  /**
   * A copy of the @p derivative_matrix object passed to the reinit() method.
   */
  FullMatrix<Number> derivative_matrix;

  /**
   * A vector containing the generalized eigenvalues of A s = lambda B s.
   */
  AlignedVector<Number> eigenvalues;

  /**
   * The matrix containing the generalized eigenvectors.
   */
  Table<2,Number> eigenvectors;

  /**
   * An array for temporary data.
   */
  mutable AlignedVector<Number> tmp_array;

  /**
   * A mutex that guards access to the array @p tmp_array.
   */
  mutable Threads::Mutex mutex;
};


/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN


template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,Number,size>
::TensorProductMatrixSymmetricSum() = default;



template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,Number,size>
::TensorProductMatrixSymmetricSum(const FullMatrix<Number> &mass_matrix,
                                  const FullMatrix<Number> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::reinit(const FullMatrix<Number> &mass_matrix,
         const FullMatrix<Number> &derivative_matrix)
{
  Assert(size == -1 ||
         (size > 0 && static_cast<unsigned int>(size) == mass_matrix.m()),
         ExcDimensionMismatch(size, mass_matrix.m()));
  AssertDimension(mass_matrix.m(), mass_matrix.n());
  AssertDimension(mass_matrix.m(), derivative_matrix.m());
  AssertDimension(mass_matrix.m(), derivative_matrix.n());

  this->mass_matrix = mass_matrix;
  this->derivative_matrix = derivative_matrix;

  std::vector<Vector<Number> > eigenvecs(mass_matrix.m());
  LAPACKFullMatrix<Number> mass_copy(mass_matrix.m(), mass_matrix.n());
  LAPACKFullMatrix<Number> deriv_copy(derivative_matrix.m(), derivative_matrix.n());
  mass_copy = mass_matrix;
  deriv_copy = derivative_matrix;

  deriv_copy.compute_generalized_eigenvalues_symmetric(mass_copy, eigenvecs);
  AssertDimension(eigenvecs.size(), mass_matrix.m());
  eigenvectors.reinit(mass_matrix.m(), mass_matrix.m());
  for (unsigned int i=0; i<mass_matrix.m(); ++i)
    for (unsigned int j=0; j<mass_matrix.n(); ++j)
      eigenvectors(i,j) = eigenvecs[j][i];

  eigenvalues.resize(mass_matrix.m());
  for (unsigned int i=0; i<mass_matrix.m(); ++i)
    eigenvalues[i] = deriv_copy.eigenvalue(i).real();
}



template <int dim, typename Number, int size>
inline
unsigned int
TensorProductMatrixSymmetricSum<dim,Number,size>::m() const
{
  return Utilities::fixed_power<dim>(mass_matrix.m());
}



template <int dim, typename Number, int size>
inline
unsigned int
TensorProductMatrixSymmetricSum<dim,Number,size>::n() const
{
  return Utilities::fixed_power<dim>(mass_matrix.n());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::vmult(Vector<Number> &dst,
        const Vector<Number> &src) const
{
  AssertDimension(dst.size(), Utilities::fixed_power<dim>(eigenvalues.size()));
  AssertDimension(src.size(), Utilities::fixed_power<dim>(eigenvalues.size()));
  vmult(dst.begin(), src.begin());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::apply_inverse(Vector<Number> &dst,
                const Vector<Number> &src) const
{
  AssertDimension(dst.size(), Utilities::fixed_power<dim>(eigenvalues.size()));
  AssertDimension(src.size(), Utilities::fixed_power<dim>(eigenvalues.size()));
  apply_inverse(dst.begin(), src.begin());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::vmult(Number *dst,
        const Number *src) const
{
  Threads::Mutex::ScopedLock lock(this->mutex);
  const unsigned int n = Utilities::fixed_power<dim>(size > 0 ? size : eigenvalues.size());
  tmp_array.resize_fast(n*2);
  const int kernel_size = size > 0 ? size-1 : -1;
  internal::EvaluatorTensorProduct<internal::evaluate_general,dim,kernel_size,kernel_size+1,Number>
  eval(AlignedVector<Number>(), AlignedVector<Number>(),
       AlignedVector<Number>(), mass_matrix.m()-1, mass_matrix.m());
  const Number *A = &derivative_matrix(0,0);
  const Number *M = &mass_matrix(0,0);
  Number *t = tmp_array.begin();
  if (dim == 1)
    eval.template apply<0, true, false>(A, src, dst);
  else if (dim == 2)
    {
      eval.template apply<0, true, false>(M, src, t);
      eval.template apply<1, true, false>(A, t, dst);
      eval.template apply<0, true, false>(A, src, t);
      eval.template apply<1, true, true> (M, t, dst);
    }
  else if (dim == 3)
    {
      eval.template apply<0, true, false>(M, src, t+n);
      eval.template apply<1, true, false>(M, t+n, t);
      eval.template apply<2, true, false>(A, t, dst);
      eval.template apply<1, true, false>(A, t+n, t);
      eval.template apply<0, true, false>(A, src, t+n);
      eval.template apply<1, true, true> (M, t+n, t);
      eval.template apply<2, true, true> (M, t, dst);
    }
  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::apply_inverse(Number *dst,
                const Number *src) const
{
  Threads::Mutex::ScopedLock lock(this->mutex);
  const unsigned int n = size > 0 ? size : eigenvalues.size();
  tmp_array.resize_fast(Utilities::fixed_power<dim>(n));
  const int kernel_size = size > 0 ? size-1 : -1;
  internal::EvaluatorTensorProduct<internal::evaluate_general,dim,kernel_size,kernel_size+1,Number>
  eval(AlignedVector<Number>(), AlignedVector<Number>(),
       AlignedVector<Number>(), mass_matrix.m()-1, mass_matrix.m());
  const Number *S = &eigenvectors(0,0);
  Number *t = tmp_array.begin();

  switch (dim)
    {
    case 1:
      eval.template apply<0, true, false> (S, src, t);
      for (unsigned int i=0; i<n; ++i)
        t[i] /= eigenvalues[i];
      eval.template apply<0, false, false> (S, t, dst);
      break;

    case 2:
      eval.template apply<0, true, false> (S, src, t);
      eval.template apply<1, true, false> (S, t, dst);
      for (unsigned int i=0, c=0; i<n; ++i)
        for (unsigned int j=0; j<n; ++j, ++c)
          dst[c] /= (eigenvalues[i] + eigenvalues[j]);
      eval.template apply<1, false, false> (S, dst, t);
      eval.template apply<0, false, false> (S, t, dst);
      break;

    case 3:
      eval.template apply<0, true, false> (S, src, t);
      eval.template apply<1, true, false> (S, t, dst);
      eval.template apply<2, true, false> (S, dst, t);
      for (unsigned int i=0, c=0; i<n; ++i)
        for (unsigned int j=0; j<n; ++j)
          for (unsigned int k=0; k<n; ++k, ++c)
            t[c] /= (eigenvalues[i] + eigenvalues[j] + eigenvalues[k]);
      eval.template apply<2, false, false> (S, t, dst);
      eval.template apply<1, false, false> (S, dst, t);
      eval.template apply<0, false, false> (S, t, dst);
      break;

    default:
      Assert(false, ExcNotImplemented());
    }
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif
