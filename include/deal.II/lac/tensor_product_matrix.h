// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_tensor_product_matrix_h
#define dealii_tensor_product_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <deal.II/matrix_free/tensor_product_kernels.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename>
class Vector;
template <typename>
class FullMatrix;
#endif

/**
 * This is an abstract base class used for a special matrix class, namely the
 * TensorProductMatrixSymmetricSum.
 *
 * First, the base class acts like a container storing 1D mass matrices and
 * 1D derivative matrices as well as the generalized eigenvalues and
 * eigenvectors for each tensor direction. For a detailed definition of these
 * matrices and corresponding generalized eigenproblems we refer to the main
 * documentation of TensorProductMatrixSymmetricSum.
 *
 * @note This base class has no functionality to calculate eigenvalues and
 * eigenvectors for mass and derivative matrices given. The responsibility of
 * initializing the data members completely lies with the derived class.
 *
 * Second, it implements the matrix-vector product with the tensor product
 * matrix (vmult()) and its inverse (apply_inverse()) as described in the
 * main documentation of TensorProductMatrixSymmetricSum.
 *
 * @note This class uses a temporary array for storing intermediate results
 * that is a class member. A mutex is used to protect access to this array and
 * ensure correct results. If several threads run parallel instances of this
 * class, it is recommended that each threads holds its own matrix version.
 *
 * @tparam dim Dimension of the problem. Currently, 1D, 2D, and 3D codes are
 * implemented.
 *
 * @tparam Number Arithmetic type of the underlying array elements.
 *
 * @tparam n_rows_1d Compile-time number of rows of 1D matrices (only
 * valid if the number of rows and columns coincide for each
 * dimension). By default at -1, which means that the number of rows
 * is determined at run-time by means of the matrices passed to the
 * reinit() function.
 *
 * @author Martin Kronbichler and Julius Witte, 2017
 */
template <int dim, typename Number, int n_rows_1d = -1>
class TensorProductMatrixSymmetricSumBase
{
public:
  /**
   * Type of matrix entries. This alias is analogous to <tt>value_type</tt>
   * in the standard library containers.
   */
  using value_type = Number;

  /**
   * The static number of rows of the 1D matrices. For more details,
   * see the description of the template parameter <tt>n_rows_1d</tt>.
   */
  static constexpr int n_rows_1d_static = n_rows_1d;

  /**
   * Return the number of rows of the tensor product matrix
   * resulting from the Kronecker product of 1D matrices, which is described
   * in the main documentation of TensorProductMatrixSymmetricSum.
   */
  unsigned int
  m() const;

  /**
   * Return the number of columns of the tensor product matrix
   * resulting from the Kronecker product of 1D matrices, which is described
   * in the main documentation of TensorProductMatrixSymmetricSum.
   */
  unsigned int
  n() const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of TensorProductMatrixSymmetricSum.
   * This function is operating on ArrayView to allow checks of
   * array bounds with respect to @p dst and @p src.
   */
  void
  vmult(const ArrayView<Number> &dst, const ArrayView<const Number> &src) const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of TensorProductMatrixSymmetricSum.
   * This function is operating on ArrayView to allow checks of
   * array bounds with respect to @p dst and @p src.
   */
  void
  apply_inverse(const ArrayView<Number> &      dst,
                const ArrayView<const Number> &src) const;

protected:
  /**
   * Default constructor.
   */
  TensorProductMatrixSymmetricSumBase() = default;

  /**
   * An array containing a mass matrix for each tensor direction.
   */
  std::array<Table<2, Number>, dim> mass_matrix;

  /**
   * An array containing a derivative matrix for each tensor direction.
   */
  std::array<Table<2, Number>, dim> derivative_matrix;

  /**
   * An array storing the generalized eigenvalues
   * for each tensor direction.
   */
  std::array<AlignedVector<Number>, dim> eigenvalues;

  /**
   * An array storing the generalized eigenvectors
   * for each tensor direction.
   */
  std::array<Table<2, Number>, dim> eigenvectors;

private:
  /**
   * An array for temporary data.
   */
  mutable AlignedVector<Number> tmp_array;

  /**
   * A mutex that guards access to the array @p tmp_array.
   */
  mutable Threads::Mutex mutex;
};



/**
 * This is a special matrix class defined as the tensor product (or Kronecker
 * product) of 1D matrices of the type
 * @f{align*}{
 * L &= A_1 \otimes M_0 + M_1 \otimes A_0
 * @f}
 * in 2D and
 * @f{align*}{
 * L &= A_2 \otimes M_1 \otimes M_0 + M_2 \otimes A_1 \otimes M_0 + M_2 \otimes
 * M_1 \otimes A_0
 * @f}
 * in 3D. The typical application setting is a discretization of the Laplacian
 * $L$ on a Cartesian (axis-aligned) geometry, where it can be exactly
 * represented by the Kronecker or tensor product of a 1D mass matrix $M$ and
 * a 1D Laplace matrix $A$ in each tensor direction (due to symmetry $M$ and $A$
 * are the same in each dimension). The dimension of the resulting class is the
 * product of the one-dimensional matrices.
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
 * L^{-1} &= S_1 \otimes S_0 (\Lambda_1 \otimes I + I \otimes \Lambda_0)^{-1}
 * S_1^\mathrm T \otimes S_0^\mathrm T,
 * @f}
 * where $S_d$ is the matrix of eigenvectors to the generalized eigenvalue
 * problem in the given tensor direction $d$:
 * @f{align*}{
 * A_d s  &= \lambda M_d s, d = 0, \quad \ldots,\mathrm{dim},
 * @f}
 * and $\Lambda_d$ is the diagonal matrix representing the generalized
 * eigenvalues $\lambda$. Note that the vectors $s$ are such that they
 * simultaneously diagonalize $A_d$ and $M_d$, i.e. $S_d^{\mathrm T} A_d S_d =
 * \Lambda_d$ and $S_d^{\mathrm T} M_d S_d = I$. This method of matrix inversion
 * is called fast diagonalization method.
 *
 * This class requires LAPACK support.
 *
 * Note that this class allows for two modes of usage. The first is a use case
 * with run time constants for the matrix dimensions that is achieved by
 * setting the optional template parameter <tt>n_rows_1d</tt> to -1. The second
 * mode of usage that is faster allows to set the template parameter as a
 * compile time constant, giving significantly faster code in particular for
 * small sizes of the matrix.
 *
 * @tparam dim Dimension of the problem. Currently, 1D, 2D, and 3D codes are
 * implemented.
 *
 * @tparam Number Arithmetic type of the underlying array elements. Note that the
 * underlying LAPACK implementation supports only float and double numbers, so
 * only these two types are currently supported by the generic class.
 * Nevertheless, a template specialization for the vectorized types
 * VectorizedArray<float> and VectorizedArray<double> exists. This is necessary
 * to perform LAPACK calculations for each vectorization lane, i.e. for the
 * supported float and double numbers.
 *
 * @tparam n_rows_1d Compile-time number of rows of 1D matrices (only
 * valid if the number of rows and columns coincide for each
 * dimension). By default at -1, which means that the number of rows
 * is determined at run-time by means of the matrices passed to the
 * reinit() function.
 *
 * @author Martin Kronbichler and Julius Witte, 2017
 */
template <int dim, typename Number, int n_rows_1d = -1>
class TensorProductMatrixSymmetricSum
  : public TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>
{
public:
  /**
   * Default constructor.
   */
  TensorProductMatrixSymmetricSum() = default;

  /**
   * Constructor that is equivalent to the empty constructor and
   * immediately calling
   * reinit(const std::array<Table<2,Number>, dim>&,const
   * std::array<Table<2,Number>, dim>&).
   */
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, Number>, dim> &mass_matrix,
    const std::array<Table<2, Number>, dim> &derivative_matrix);

  /**
   * Constructor that is equivalent to the empty constructor and
   * immediately calling
   * reinit(const std::array<FullMatrix<Number>,dim>&,const
   * std::array<FullMatrix<Number>,dim>&).
   */
  TensorProductMatrixSymmetricSum(
    const std::array<FullMatrix<Number>, dim> &mass_matrix,
    const std::array<FullMatrix<Number>, dim> &derivative_matrix);

  /**
   * Constructor that is equivalent to the empty constructor and
   * immediately calling reinit(const Table<2,Number>&,const Table<2,Number>&).
   */
  TensorProductMatrixSymmetricSum(const Table<2, Number> &mass_matrix,
                                  const Table<2, Number> &derivative_matrix);

  /**
   * Initializes the tensor product matrix by copying the arrays of 1D mass
   * matrices @p mass_matrix and 1D derivative matrices @p derivative_matrix into its
   * base class counterparts, respectively, and by assembling the regarding
   * generalized eigenvalues and eigenvectors in
   * TensorProductMatrixSymmetricSumBase::eigenvalues
   * and TensorProductMatrixSymmetricSumBase::eigenvectors, respectively.
   * Note that the current implementation requires each $M_{d}$ to be symmetric
   * and positive definite and every $A_{d}$ to be symmetric and invertible but
   * not necessarily positive definite.
   */
  void
  reinit(const std::array<Table<2, Number>, dim> &mass_matrix,
         const std::array<Table<2, Number>, dim> &derivative_matrix);

  /**
   * This function is equivalent to the previous reinit() except that
   * the 1D matrices in @p mass_matrix and @p derivative_matrix are
   * passed in terms of a FullMatrix, respectively.
   */
  void
  reinit(const std::array<FullMatrix<Number>, dim> &mass_matrix,
         const std::array<FullMatrix<Number>, dim> &derivative_matrix);

  /**
   * This function is equivalent to the first reinit() except that
   * we consider the same 1D mass matrix @p mass_matrix and the same 1D
   * derivative matrix @p derivative_matrix for each tensor direction.
   */
  void
  reinit(const Table<2, Number> &mass_matrix,
         const Table<2, Number> &derivative_matrix);

private:
  /**
   * A generic implementation of all reinit() functions based on
   * perfect forwarding, that allows to pass lvalue as well
   * as rvalue arguments.
   * @tparam MatrixArray Has to be convertible to the underlying
   * type of TensorProductMatrixSymmetricSumBase::mass_matrix and
   * TensorProductMatrixSymmetricSumBase::derivative_matrix.
   */
  template <typename MatrixArray>
  void
  reinit_impl(MatrixArray &&mass_matrix, MatrixArray &&derivative_matrix);
};



/**
 * This is the template specialization for VectorizedArray<Number>
 * being the arithmetic template. For a detailed description see
 * the main documentation of the generic
 * TensorProductMatrixSymmetricSum class.
 *
 * @author Martin Kronbichler and Julius Witte, 2017
 */
template <int dim, typename Number, int n_rows_1d>
class TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>
  : public TensorProductMatrixSymmetricSumBase<dim,
                                               VectorizedArray<Number>,
                                               n_rows_1d>
{
public:
  /**
   * Default constructor.
   */
  TensorProductMatrixSymmetricSum() = default;

  /**
   * Constructor that is equivalent to the empty constructor and
   * immediately calling
   * reinit(const std::array<Table<2,VectorizedArray<Number> >, dim>&,const
   * std::array<Table<2,VectorizedArray<Number> >, dim>&).
   */
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
    const std::array<Table<2, VectorizedArray<Number>>, dim>
      &derivative_matrix);

  /**
   * Constructor that is equivalent to the empty constructor and
   * immediately calling
   * reinit(const Table<2,VectorizedArray<Number> >&,const
   * Table<2,VectorizedArray<Number> >&).
   */
  TensorProductMatrixSymmetricSum(
    const Table<2, VectorizedArray<Number>> &mass_matrix,
    const Table<2, VectorizedArray<Number>> &derivative_matrix);

  /**
   * Initializes the tensor product matrix by copying the arrays of 1D mass
   * matrices @p mass_matrix and 1D derivative matrices @p derivative_matrix into its
   * base class counterparts, respectively, and by assembling the regarding
   * generalized eigenvalues and eigenvectors in
   * TensorProductMatrixSymmetricSumBase::eigenvalues
   * and TensorProductMatrixSymmetricSumBase::eigenvectors, respectively.
   * Note that the current implementation requires each $M_{d}$ to be symmetric
   * and positive definite and every $A_{d}$ to be symmetric and invertible but
   * not necessarily positive definite.
   */
  void
  reinit(const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
         const std::array<Table<2, VectorizedArray<Number>>, dim>
           &derivative_matrix);

  /**
   * This function is equivalent to the previous reinit() except that
   * we consider the same 1D mass matrix @p mass_matrix and the same 1D
   * derivative matrix @p derivative_matrix for each tensor direction.
   */
  void
  reinit(const Table<2, VectorizedArray<Number>> &mass_matrix,
         const Table<2, VectorizedArray<Number>> &derivative_matrix);

private:
  /**
   * A generic implementation of all reinit() functions based on
   * perfect forwarding, that allows to pass lvalue as well
   * as rvalue arguments.
   * @tparam MatrixArray Has to be convertible to the underlying
   * type of TensorProductMatrixSymmetricSumBase::mass_matrix and
   * TensorProductMatrixSymmetricSumBase::derivative_matrix.
   */
  template <typename MatrixArray>
  void
  reinit_impl(MatrixArray &&mass_matrix, MatrixArray &&derivative_matrix);
};


/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

namespace internal
{
  namespace TensorProductMatrix
  {
    /**
     * Compute generalized eigenvalues and eigenvectors of the real
     * generalized symmetric eigenproblem $A v = \lambda M v$. Since we are
     * operating on plain pointers we require the size of the matrices
     * beforehand. Note that the data arrays for the eigenvalues and
     * eigenvectors have to be initialized to a proper size, too. (no check of
     * array bounds possible)
     */
    template <typename Number>
    void
    spectral_assembly(const Number *     mass_matrix,
                      const Number *     derivative_matrix,
                      const unsigned int n_rows,
                      const unsigned int n_cols,
                      Number *           eigenvalues,
                      Number *           eigenvectors)
    {
      Assert(n_rows == n_cols, ExcNotImplemented());

      auto &&transpose_fill_nm = [](Number *           out,
                                    const Number *     in,
                                    const unsigned int n,
                                    const unsigned int m) {
        for (unsigned int mm = 0; mm < m; ++mm)
          for (unsigned int nn = 0; nn < n; ++nn)
            out[mm + nn * m] = *(in++);
      };

      std::vector<dealii::Vector<Number>> eigenvecs(n_rows);
      LAPACKFullMatrix<Number>            mass_copy(n_rows, n_cols);
      LAPACKFullMatrix<Number>            deriv_copy(n_rows, n_cols);

      transpose_fill_nm(&(mass_copy(0, 0)), mass_matrix, n_rows, n_cols);
      transpose_fill_nm(&(deriv_copy(0, 0)), derivative_matrix, n_rows, n_cols);

      deriv_copy.compute_generalized_eigenvalues_symmetric(mass_copy,
                                                           eigenvecs);
      AssertDimension(eigenvecs.size(), n_rows);
      for (unsigned int i = 0; i < n_rows; ++i)
        for (unsigned int j = 0; j < n_cols; ++j, ++eigenvectors)
          *eigenvectors = eigenvecs[j][i];

      for (unsigned int i = 0; i < n_rows; ++i, ++eigenvalues)
        *eigenvalues = deriv_copy.eigenvalue(i).real();
    }
  } // namespace TensorProductMatrix
} // namespace internal


template <int dim, typename Number, int n_rows_1d>
inline unsigned int
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::m() const
{
  unsigned int m = mass_matrix[0].n_rows();
  for (unsigned int d = 1; d < dim; ++d)
    m *= mass_matrix[d].n_rows();
  return m;
}



template <int dim, typename Number, int n_rows_1d>
inline unsigned int
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::n() const
{
  unsigned int n = mass_matrix[0].n_cols();
  for (unsigned int d = 1; d < dim; ++d)
    n *= mass_matrix[d].n_cols();
  return n;
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::vmult(
  const ArrayView<Number> &      dst_view,
  const ArrayView<const Number> &src_view) const
{
  AssertDimension(dst_view.size(), this->m());
  AssertDimension(src_view.size(), this->n());
  std::lock_guard<std::mutex> lock(this->mutex);
  const unsigned int          n = Utilities::fixed_power<dim>(
    n_rows_1d > 0 ? n_rows_1d : eigenvalues[0].size());
  tmp_array.resize_fast(n * 2);
  constexpr int kernel_size = n_rows_1d > 0 ? n_rows_1d : 0;
  internal::EvaluatorTensorProduct<internal::evaluate_general,
                                   dim,
                                   kernel_size,
                                   kernel_size,
                                   Number>
                eval(AlignedVector<Number>{},
         AlignedVector<Number>{},
         AlignedVector<Number>{},
         mass_matrix[0].n_rows(),
         mass_matrix[0].n_rows());
  Number *      t   = tmp_array.begin();
  const Number *src = src_view.begin();
  Number *      dst = dst_view.data();

  if (dim == 1)
    {
      const Number *A = &derivative_matrix[0](0, 0);
      eval.template apply<0, false, false>(A, src, dst);
    }

  else if (dim == 2)
    {
      const Number *A0 = &derivative_matrix[0](0, 0);
      const Number *M0 = &mass_matrix[0](0, 0);
      const Number *A1 = &derivative_matrix[1](0, 0);
      const Number *M1 = &mass_matrix[1](0, 0);
      eval.template apply<0, false, false>(M0, src, t);
      eval.template apply<1, false, false>(A1, t, dst);
      eval.template apply<0, false, false>(A0, src, t);
      eval.template apply<1, false, true>(M1, t, dst);
    }

  else if (dim == 3)
    {
      const Number *A0 = &derivative_matrix[0](0, 0);
      const Number *M0 = &mass_matrix[0](0, 0);
      const Number *A1 = &derivative_matrix[1](0, 0);
      const Number *M1 = &mass_matrix[1](0, 0);
      const Number *A2 = &derivative_matrix[2](0, 0);
      const Number *M2 = &mass_matrix[2](0, 0);
      eval.template apply<0, false, false>(M0, src, t + n);
      eval.template apply<1, false, false>(M1, t + n, t);
      eval.template apply<2, false, false>(A2, t, dst);
      eval.template apply<1, false, false>(A1, t + n, t);
      eval.template apply<0, false, false>(A0, src, t + n);
      eval.template apply<1, false, true>(M1, t + n, t);
      eval.template apply<2, false, true>(M2, t, dst);
    }

  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSumBase<dim, Number, n_rows_1d>::apply_inverse(
  const ArrayView<Number> &      dst_view,
  const ArrayView<const Number> &src_view) const
{
  AssertDimension(dst_view.size(), this->n());
  AssertDimension(src_view.size(), this->m());
  std::lock_guard<std::mutex> lock(this->mutex);
  const unsigned int n = n_rows_1d > 0 ? n_rows_1d : eigenvalues[0].size();
  tmp_array.resize_fast(Utilities::fixed_power<dim>(n));
  constexpr int kernel_size = n_rows_1d > 0 ? n_rows_1d : 0;
  internal::EvaluatorTensorProduct<internal::evaluate_general,
                                   dim,
                                   kernel_size,
                                   kernel_size,
                                   Number>
                eval(AlignedVector<Number>(),
         AlignedVector<Number>(),
         AlignedVector<Number>(),
         mass_matrix[0].n_rows(),
         mass_matrix[0].n_rows());
  Number *      t   = tmp_array.begin();
  const Number *src = src_view.data();
  Number *      dst = dst_view.data();

  // NOTE: dof_to_quad has to be interpreted as 'dof to eigenvalue index'
  //       --> apply<.,true,.> (S,src,dst) calculates dst = S^T * src,
  //       --> apply<.,false,.> (S,src,dst) calculates dst = S * src,
  //       while the eigenvectors are stored column-wise in S, i.e.
  //       rows correspond to dofs whereas columns to eigenvalue indices!
  if (dim == 1)
    {
      const Number *S = &eigenvectors[0](0, 0);
      eval.template apply<0, true, false>(S, src, t);
      for (unsigned int i = 0; i < n; ++i)
        t[i] /= eigenvalues[0][i];
      eval.template apply<0, false, false>(S, t, dst);
    }

  else if (dim == 2)
    {
      const Number *S0 = &(eigenvectors[0](0, 0));
      const Number *S1 = &(eigenvectors[1](0, 0));
      eval.template apply<0, true, false>(S0, src, t);
      eval.template apply<1, true, false>(S1, t, dst);
      for (unsigned int i1 = 0, c = 0; i1 < n; ++i1)
        for (unsigned int i0 = 0; i0 < n; ++i0, ++c)
          dst[c] /= (eigenvalues[1][i1] + eigenvalues[0][i0]);
      eval.template apply<0, false, false>(S0, dst, t);
      eval.template apply<1, false, false>(S1, t, dst);
    }

  else if (dim == 3)
    {
      const Number *S0 = &eigenvectors[0](0, 0);
      const Number *S1 = &eigenvectors[1](0, 0);
      const Number *S2 = &eigenvectors[2](0, 0);
      eval.template apply<0, true, false>(S0, src, t);
      eval.template apply<1, true, false>(S1, t, dst);
      eval.template apply<2, true, false>(S2, dst, t);
      for (unsigned int i2 = 0, c = 0; i2 < n; ++i2)
        for (unsigned int i1 = 0; i1 < n; ++i1)
          for (unsigned int i0 = 0; i0 < n; ++i0, ++c)
            t[c] /=
              (eigenvalues[2][i2] + eigenvalues[1][i1] + eigenvalues[0][i0]);
      eval.template apply<0, false, false>(S0, t, dst);
      eval.template apply<1, false, false>(S1, dst, t);
      eval.template apply<2, false, false>(S2, t, dst);
    }

  else
    Assert(false, ExcNotImplemented());
}


//---------------------- TensorProductMatrixSymmetricSum ----------------------

template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, Number>, dim> &mass_matrix,
    const std::array<Table<2, Number>, dim> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const std::array<FullMatrix<Number>, dim> &mass_matrix,
    const std::array<FullMatrix<Number>, dim> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::
  TensorProductMatrixSymmetricSum(const Table<2, Number> &mass_matrix,
                                  const Table<2, Number> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
template <typename MatrixArray>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit_impl(
  MatrixArray &&mass_matrices_,
  MatrixArray &&derivative_matrices_)
{
  auto &&mass_matrices       = std::forward<MatrixArray>(mass_matrices_);
  auto &&derivative_matrices = std::forward<MatrixArray>(derivative_matrices_);
  this->mass_matrix          = mass_matrices;
  this->derivative_matrix    = derivative_matrices;

  for (int dir = 0; dir < dim; ++dir)
    {
      Assert(n_rows_1d == -1 ||
               (n_rows_1d > 0 && static_cast<unsigned int>(n_rows_1d) ==
                                   mass_matrices[dir].n_rows()),
             ExcDimensionMismatch(n_rows_1d, mass_matrices[dir].n_rows()));
      AssertDimension(mass_matrices[dir].n_rows(), mass_matrices[dir].n_cols());
      AssertDimension(mass_matrices[dir].n_rows(),
                      derivative_matrices[dir].n_rows());
      AssertDimension(mass_matrices[dir].n_rows(),
                      derivative_matrices[dir].n_cols());

      this->eigenvectors[dir].reinit(mass_matrices[dir].n_cols(),
                                     mass_matrices[dir].n_rows());
      this->eigenvalues[dir].resize(mass_matrices[dir].n_cols());
      internal::TensorProductMatrix::spectral_assembly<Number>(
        &(mass_matrices[dir](0, 0)),
        &(derivative_matrices[dir](0, 0)),
        mass_matrices[dir].n_rows(),
        mass_matrices[dir].n_cols(),
        this->eigenvalues[dir].begin(),
        &(this->eigenvectors[dir](0, 0)));
    }
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const std::array<Table<2, Number>, dim> &mass_matrix,
  const std::array<Table<2, Number>, dim> &derivative_matrix)
{
  reinit_impl(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const std::array<FullMatrix<Number>, dim> &mass_matrix,
  const std::array<FullMatrix<Number>, dim> &derivative_matrix)
{
  std::array<Table<2, Number>, dim> mass_copy;
  std::array<Table<2, Number>, dim> deriv_copy;

  std::transform(mass_matrix.cbegin(),
                 mass_matrix.cend(),
                 mass_copy.begin(),
                 [](const FullMatrix<Number> &m) -> Table<2, Number> {
                   return m;
                 });
  std::transform(derivative_matrix.cbegin(),
                 derivative_matrix.cend(),
                 deriv_copy.begin(),
                 [](const FullMatrix<Number> &m) -> Table<2, Number> {
                   return m;
                 });

  reinit_impl(std::move(mass_copy), std::move(deriv_copy));
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, Number, n_rows_1d>::reinit(
  const Table<2, Number> &mass_matrix,
  const Table<2, Number> &derivative_matrix)
{
  std::array<Table<2, Number>, dim> mass_matrices;
  std::array<Table<2, Number>, dim> derivative_matrices;

  std::fill(mass_matrices.begin(), mass_matrices.end(), mass_matrix);
  std::fill(derivative_matrices.begin(),
            derivative_matrices.end(),
            derivative_matrix);

  reinit_impl(std::move(mass_matrices), std::move(derivative_matrices));
}



//------------- vectorized spec.: TensorProductMatrixSymmetricSum -------------

template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim,
                                       VectorizedArray<Number>,
                                       n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
    const std::array<Table<2, VectorizedArray<Number>>, dim> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline TensorProductMatrixSymmetricSum<dim,
                                       VectorizedArray<Number>,
                                       n_rows_1d>::
  TensorProductMatrixSymmetricSum(
    const Table<2, VectorizedArray<Number>> &mass_matrix,
    const Table<2, VectorizedArray<Number>> &derivative_matrix)
{
  reinit(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
template <typename MatrixArray>
inline void
TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>::
  reinit_impl(MatrixArray &&mass_matrices_, MatrixArray &&derivative_matrices_)
{
  auto &&mass_matrix       = std::forward<MatrixArray>(mass_matrices_);
  auto &&derivative_matrix = std::forward<MatrixArray>(derivative_matrices_);
  this->mass_matrix        = mass_matrix;
  this->derivative_matrix  = derivative_matrix;

  constexpr unsigned int macro_size = VectorizedArray<Number>::size();
  std::size_t            n_rows_max = (n_rows_1d > 0) ? n_rows_1d : 0;
  if (n_rows_1d == -1)
    for (unsigned int d = 0; d < dim; ++d)
      n_rows_max = std::max(n_rows_max, mass_matrix[d].n_rows());
  const std::size_t nm_flat_size_max = n_rows_max * n_rows_max * macro_size;
  const std::size_t n_flat_size_max  = n_rows_max * macro_size;

  std::vector<Number> mass_matrix_flat;
  std::vector<Number> deriv_matrix_flat;
  std::vector<Number> eigenvalues_flat;
  std::vector<Number> eigenvectors_flat;
  mass_matrix_flat.resize(nm_flat_size_max);
  deriv_matrix_flat.resize(nm_flat_size_max);
  eigenvalues_flat.resize(n_flat_size_max);
  eigenvectors_flat.resize(nm_flat_size_max);
  std::array<unsigned int, macro_size> offsets_nm;
  std::array<unsigned int, macro_size> offsets_n;
  for (int dir = 0; dir < dim; ++dir)
    {
      Assert(n_rows_1d == -1 ||
               (n_rows_1d > 0 && static_cast<unsigned int>(n_rows_1d) ==
                                   mass_matrix[dir].n_rows()),
             ExcDimensionMismatch(n_rows_1d, mass_matrix[dir].n_rows()));
      AssertDimension(mass_matrix[dir].n_rows(), mass_matrix[dir].n_cols());
      AssertDimension(mass_matrix[dir].n_rows(),
                      derivative_matrix[dir].n_rows());
      AssertDimension(mass_matrix[dir].n_rows(),
                      derivative_matrix[dir].n_cols());

      const unsigned int n_rows = mass_matrix[dir].n_rows();
      const unsigned int n_cols = mass_matrix[dir].n_cols();
      const unsigned int nm     = n_rows * n_cols;
      for (unsigned int vv = 0; vv < macro_size; ++vv)
        offsets_nm[vv] = nm * vv;

      vectorized_transpose_and_store(false,
                                     nm,
                                     &(mass_matrix[dir](0, 0)),
                                     offsets_nm.cbegin(),
                                     mass_matrix_flat.data());
      vectorized_transpose_and_store(false,
                                     nm,
                                     &(derivative_matrix[dir](0, 0)),
                                     offsets_nm.cbegin(),
                                     deriv_matrix_flat.data());

      const Number *mass_cbegin    = mass_matrix_flat.data();
      const Number *deriv_cbegin   = deriv_matrix_flat.data();
      Number *      eigenvec_begin = eigenvectors_flat.data();
      Number *      eigenval_begin = eigenvalues_flat.data();
      for (unsigned int lane = 0; lane < macro_size; ++lane)
        internal::TensorProductMatrix::spectral_assembly<Number>(
          mass_cbegin + nm * lane,
          deriv_cbegin + nm * lane,
          n_rows,
          n_cols,
          eigenval_begin + n_rows * lane,
          eigenvec_begin + nm * lane);

      this->eigenvalues[dir].resize(n_rows);
      this->eigenvectors[dir].reinit(n_rows, n_cols);
      for (unsigned int vv = 0; vv < macro_size; ++vv)
        offsets_n[vv] = n_rows * vv;
      vectorized_load_and_transpose(n_rows,
                                    eigenvalues_flat.data(),
                                    offsets_n.cbegin(),
                                    this->eigenvalues[dir].begin());
      vectorized_load_and_transpose(nm,
                                    eigenvectors_flat.data(),
                                    offsets_nm.cbegin(),
                                    &(this->eigenvectors[dir](0, 0)));
    }
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>::
  reinit(
    const std::array<Table<2, VectorizedArray<Number>>, dim> &mass_matrix,
    const std::array<Table<2, VectorizedArray<Number>>, dim> &derivative_matrix)
{
  reinit_impl(mass_matrix, derivative_matrix);
}



template <int dim, typename Number, int n_rows_1d>
inline void
TensorProductMatrixSymmetricSum<dim, VectorizedArray<Number>, n_rows_1d>::
  reinit(const Table<2, VectorizedArray<Number>> &mass_matrix,
         const Table<2, VectorizedArray<Number>> &derivative_matrix)
{
  std::array<Table<2, VectorizedArray<Number>>, dim> mass_matrices;
  std::array<Table<2, VectorizedArray<Number>>, dim> derivative_matrices;

  std::fill(mass_matrices.begin(), mass_matrices.end(), mass_matrix);
  std::fill(derivative_matrices.begin(),
            derivative_matrices.end(),
            derivative_matrix);

  reinit_impl(std::move(mass_matrices), std::move(derivative_matrices));
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif
