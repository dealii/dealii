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
class TensorProductMatrixSymmetricSumBase
{
public:
  /**
   * Returns the number of rows of this matrix, given by the dim-th power of
   * the size of the 1D matrices passed to the constructor.
   */
  unsigned int m () const;

  /**
   * Returns the number of columns of this matrix, given by the dim-th power
   * of the size of the 1D matrices passed to the constructor.
   */
  unsigned int n () const;

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
   * described in the main documentation of this class. Same as the other
   * apply_inverse() function, but operating on plain pointers rather than a
   * vector (no check of array bounds possible).
   */
  void apply_inverse (Number *dst,
                      const Number *src) const;

protected:
  /**
   * Constructor.
   */
  TensorProductMatrixSymmetricSumBase () = default ;

  /**
   * A copy of the @p mass_matrix object passed to the reinit() method.
   */
  std::array<Table<2,Number>,dim> mass_matrix;

  /**
   * A copy of the @p derivative_matrix object passed to the reinit() method.
   */
  std::array<Table<2,Number>,dim> derivative_matrix;

  /**
   * A vector containing the generalized eigenvalues of A s = lambda B s.
   */
  std::array<AlignedVector<Number>,dim> eigenvalues;

  /**
   * The matrix containing the generalized eigenvectors.
   */
  std::array<Table<2,Number>,dim> eigenvectors;

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
 * ... new TensorProductMatrixSymmetricSum using the base class as tensor product
 * container and interface to arithmetic operations for a generic Number type ...
 */
template <int dim, typename Number, int size = -1>
class TensorProductMatrixSymmetricSum
  : public TensorProductMatrixSymmetricSumBase<dim,Number,size>
{
public:
  /**
   * Constructor.
   */
  TensorProductMatrixSymmetricSum () ;

  /**
   * Constructor that is equivalent to the previous constructor and
   * immediately calling the corresponding reinit().
   */
  TensorProductMatrixSymmetricSum (const std::array<Table<2,Number>, dim> &mass_matrix,
                                   const std::array<Table<2,Number>, dim> &derivative_matrix) ;

  /**
   * Constructor that is equivalent to the first constructor and
   * immediately calling the corresponding reinit().
   */
  TensorProductMatrixSymmetricSum (const std::array<FullMatrix<Number>,dim> &mass_matrix,
                                   const std::array<FullMatrix<Number>,dim> &derivative_matrix) ;

  /**
   * Constructor that is equivalent to the first constructor and
   * immediately calling the corresponding reinit().
   */
  TensorProductMatrixSymmetricSum (const FullMatrix<Number> &mass_matrix,
                                   const FullMatrix<Number> &derivative_matrix) ;

  /**
   * Initializes the tensor product matrix to the given mass matrices $M_0,\ldots,M_{dim}$
   * and derivative matrices $A_0,\ldots,A_{dim}$.
   * Note that the current implementation requires each $M_{d}$ to be symmetric
   * and positive definite and every $A_{d}$ to be symmetric and invertible but not
   * necessarily positive defininte.
   */
  void reinit (const std::array<Table<2,Number>,dim> &mass_matrix,
               const std::array<Table<2,Number>,dim> &derivative_matrix) ;

  /**
   * Equivalent to the previous reinit() unless that the mass and derivative
   * matrices are passed by Table instead of FullMatrix.
   */
  void reinit (const std::array<FullMatrix<Number>,dim> &mass_matrix,
               const std::array<FullMatrix<Number>,dim> &derivative_matrix) ;

  /**
   * Initializes the same mass matrix $M$ and derivative matrix $A$ to the given array
   * of mass matrices and array of derivative matrices, respectively.
   * Note that the current implementation requires $M$ to be symmetric
   * and positive definite and $A$ to be symmetric and invertible but not
   * necessarily positive defininte.
   */
  void reinit (const FullMatrix<Number> &mass_matrix,
               const FullMatrix<Number> &derivative_matrix) ;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class.
   */
  void vmult (Vector<Number> &dst,
              const Vector<Number> &src) const;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class.
   */
  void apply_inverse (Vector<Number> &dst,
                      const Vector<Number> &src) const;

  /**
   * ... for compability to MappingQGeneric
   */
  using TensorProductMatrixSymmetricSumBase<dim,Number,size>::vmult ;

  /**
   * ... for compability to MappingQGeneric
   */
  using TensorProductMatrixSymmetricSumBase<dim,Number,size>::apply_inverse ;

private:
  /**
   * A generic implementation of all reinit() functions based on
   * perfect forwarding, that makes it possible to pass lvalue as well
   * as rvalue arguments. MatrixArray has to be convertible to the underlying
   * type of the bass class' members mass_matrices and derivative_matrices.
   */
  template <typename MatrixArray>
  void reinit_impl (MatrixArray &&mass_matrix,
                    MatrixArray &&derivative_matrix) ;
};


/**
 * ... same as previous class but based on a vectorized value type, namely
 * VectorizedArray<Number> ...
 */
template <int dim, typename Number, int size>
class TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
  : public TensorProductMatrixSymmetricSumBase<dim,VectorizedArray<Number>,size>
{
public:
  /**
   * Constructor.
   */
  TensorProductMatrixSymmetricSum () ;

  /**
   * Constructor that is equivalent to the previous constructor and
   * immediately calling reinit().
   */
  TensorProductMatrixSymmetricSum (const std::array<Table<2,VectorizedArray<Number> >,dim> &mass_matrix,
                                   const std::array<Table<2,VectorizedArray<Number> >,dim> &derivative_matrix) ;

  /**
   * Constructor that is equivalent to the first constructor and
   * immediately calling the corresponding reinit().
   */
  TensorProductMatrixSymmetricSum (const Table<2,VectorizedArray<Number> > &mass_matrix,
                                   const Table<2,VectorizedArray<Number> > &derivative_matrix) ;

  /**
   * Initializes the tensor product matrix to the given mass matrices $M_0,\ldots,M_{dim}$
   * and derivative matrices $A_0,\ldots,A_{dim}$.
   * Note that the current implementation requires each $M_{d}$ to be symmetric
   * and positive definite and every $A_{d}$ to be symmetric and invertible but not
   * necessarily positive defininte.
   */
  void reinit (const std::array<Table<2,VectorizedArray<Number> >,dim> &mass_matrix,
               const std::array<Table<2,VectorizedArray<Number> >,dim> &derivative_matrix) ;

  /**
   * Initializes the same mass matrix $M$ and derivative matrix $A$ to the given array
   * of mass matrices and array of derivative matrices, respectively.
   * Note that the current implementation requires $M$ to be symmetric
   * and positive definite and $A$ to be symmetric and invertible but not
   * necessarily positive defininte.
   */
  void reinit (const Table<2,VectorizedArray<Number> > &mass_matrix,
               const Table<2,VectorizedArray<Number> > &derivative_matrix) ;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class.
   */
  void vmult (AlignedVector<VectorizedArray<Number> > &dst,
              const AlignedVector<VectorizedArray<Number> > &src) const ;

  /**
   * Implements a matrix-vector product with the underlying matrix as
   * described in the main documentation of this class.
   */
  void apply_inverse (AlignedVector<VectorizedArray<Number> > &dst,
                      const AlignedVector<VectorizedArray<Number> > &src) const ;

private:
  /**
   * A generic implementation of all reinit() functions based on
   * perfect forwarding, that makes it possible to pass lvalue as well
   * as rvalue arguments. MatrixArray has to be convertible to the underlying
   * type of the bass class' members mass_matrices and derivative_matrices.
   */
  template <typename MatrixArray>
  void reinit_impl (MatrixArray &&mass_matrix,
                    MatrixArray &&derivative_matrix) ;
};


/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

namespace
{
  /**
   * Compute generalized eigenvalues and eigenvectors of the real
   * generalized symmetric eigenproblem $M v = \lambda A v$. Since we are
   * operating on plain pointers we require the size of the matrices beforehand.
   * Note that the data arrays for the eigenvalues and eigenvectors
   * have to be initialized to a proper size, too. (no check of array bounds
   * possible)
   */
  template <typename Number>
  void
  spectral_assembly (const Number *mass_matrix,
                     const Number *derivative_matrix,
                     const unsigned int n_rows,
                     const unsigned int n_cols,
                     Number *eigenvalues,
                     Number *eigenvectors)
  {
    Assert (n_rows == n_cols, ExcNotImplemented()) ;

    auto &&transpose_fill_nm
      = [](Number *out, const Number *in, const unsigned int n, const unsigned int m)
    {
      for (unsigned int mm = 0; mm < m; ++mm)
        for (unsigned int nn = 0; nn < n; ++nn)
          out[mm+nn*m] = *(in++) ;
    };

    std::vector<Vector<Number> > eigenvecs(n_rows) ;
    LAPACKFullMatrix<Number> mass_copy(n_rows, n_cols) ;
    LAPACKFullMatrix<Number> deriv_copy(n_rows, n_cols) ;

    transpose_fill_nm (&(mass_copy(0,0)), mass_matrix, n_rows, n_cols) ;
    transpose_fill_nm (&(deriv_copy(0,0)), derivative_matrix, n_rows, n_cols) ;

    deriv_copy.compute_generalized_eigenvalues_symmetric (mass_copy, eigenvecs);
    AssertDimension (eigenvecs.size(), n_rows) ;
    for (unsigned int i=0; i<n_rows; ++i)
      for (unsigned int j=0; j<n_cols; ++j, ++eigenvectors)
        *eigenvectors = eigenvecs[j][i] ;

    for (unsigned int i=0; i<n_rows; ++i, ++eigenvalues)
      *eigenvalues = deriv_copy.eigenvalue(i).real();
  }
}



template <int dim, typename Number, int size>
inline
unsigned int
TensorProductMatrixSymmetricSumBase<dim,Number,size>::m() const
{
  unsigned int m = mass_matrix[0].n_rows() ;
  for (unsigned int d = 1; d < dim; ++d)
    m *= mass_matrix[d].n_rows() ;
  return m ;
}



template <int dim, typename Number, int size>
inline
unsigned int
TensorProductMatrixSymmetricSumBase<dim,Number,size>::n() const
{
  unsigned int n = mass_matrix[0].n_cols() ;
  for (unsigned int d = 1; d < dim; ++d)
    n *= mass_matrix[d].n_cols() ;
  return n ;
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSumBase<dim,Number,size>
::vmult(Number *dst,
        const Number *src) const
{
  Threads::Mutex::ScopedLock lock(this->mutex);
  const unsigned int n = Utilities::fixed_power<dim>(size > 0 ? size : eigenvalues[0].size());
  tmp_array.resize_fast(n*2);
  constexpr int kernel_size = size > 0 ? size-1 : -1;
  internal::EvaluatorTensorProduct<internal::evaluate_general,dim,kernel_size,kernel_size+1,Number>
  eval(AlignedVector<Number> {}, AlignedVector<Number> {},
       AlignedVector<Number> {}, mass_matrix[0].n_rows()-1, mass_matrix[0].n_rows());
  Number *t = tmp_array.begin();

  if (dim == 1)
    {
      const Number *A = &derivative_matrix[0](0,0);
      eval.template apply<0, false, false> (A, src, dst);
    }

  else if (dim == 2)
    {
      const Number *A0 = &derivative_matrix[0](0,0);
      const Number *M0 = &mass_matrix[0](0,0);
      const Number *A1 = &derivative_matrix[1](0,0);
      const Number *M1 = &mass_matrix[1](0,0);
      eval.template apply<0, false, false> (M0, src, t);
      eval.template apply<1, false, false> (A1, t, dst);
      eval.template apply<0, false, false> (A0, src, t);
      eval.template apply<1, false, true>  (M1, t, dst);
    }

  else if (dim == 3)
    {
      const Number *A0 = &derivative_matrix[0](0,0);
      const Number *M0 = &mass_matrix[0](0,0);
      const Number *A1 = &derivative_matrix[1](0,0);
      const Number *M1 = &mass_matrix[1](0,0);
      const Number *A2 = &derivative_matrix[2](0,0);
      const Number *M2 = &mass_matrix[2](0,0);
      eval.template apply<0, false, false> (M0, src, t+n);
      eval.template apply<1, false, false> (M1, t+n, t);
      eval.template apply<2, false, false> (A2, t, dst);
      eval.template apply<1, false, false> (A1, t+n, t);
      eval.template apply<0, false, false> (A0, src, t+n);
      eval.template apply<1, false, true> (M1, t+n, t);
      eval.template apply<2, false, true> (M2, t, dst);
    }

  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSumBase<dim,Number,size>
::apply_inverse(Number *dst,
                const Number *src) const
{
  Threads::Mutex::ScopedLock lock(this->mutex);
  const unsigned int n = size > 0 ? size : eigenvalues[0].size();
  tmp_array.resize_fast (Utilities::fixed_power<dim>(n));
  constexpr int kernel_size = size > 0 ? size-1 : -1;
  internal::EvaluatorTensorProduct<internal::evaluate_general,dim,kernel_size,kernel_size+1,Number>
  eval(AlignedVector<Number>(), AlignedVector<Number>(),
       AlignedVector<Number>(), mass_matrix[0].n_rows()-1, mass_matrix[0].n_rows());
  Number *t = tmp_array.begin();

  // NOTE: dof_to_quad has to be interpreted as 'dof to eigenvalue index'
  //       --> apply<.,true,.> (S,src,dst) calculates dst = S^T * src,
  //       --> apply<.,false,.> (S,src,dst) calculates dst = S * src,
  //       while the eigenvectors are stored column-wise in S, i.e.
  //       rows correspond to dofs whereas columns to eigenvalue indices!
  if (dim == 1)
    {
      const Number *S = &eigenvectors[0](0,0);
      eval.template apply<0, true, false> (S, src, t);
      for (unsigned int i=0; i<n; ++i)
        t[i] /= eigenvalues[0][i];
      eval.template apply<0, false, false> (S, t, dst);
    }

  else if (dim == 2)
    {
      const Number *S0 = &(eigenvectors[0](0,0));
      const Number *S1 = &(eigenvectors[1](0,0));
      eval.template apply<0, true, false> (S0, src, t);
      eval.template apply<1, true, false> (S1, t, dst);
      for (unsigned int i1=0, c=0; i1<n; ++i1)
        for (unsigned int i0=0; i0<n; ++i0, ++c)
          dst[c] /= (eigenvalues[1][i1] + eigenvalues[0][i0]);
      eval.template apply<0, false, false> (S0, dst, t);
      eval.template apply<1, false, false> (S1, t, dst);
    }

  else if (dim == 3)
    {
      const Number *S0 = &eigenvectors[0](0,0);
      const Number *S1 = &eigenvectors[1](0,0);
      const Number *S2 = &eigenvectors[2](0,0);
      eval.template apply<0, true, false> (S0, src, t);
      eval.template apply<1, true, false> (S1, t, dst);
      eval.template apply<2, true, false> (S2, dst, t);
      for (unsigned int i2=0, c=0; i2<n; ++i2)
        for (unsigned int i1=0; i1<n; ++i1)
          for (unsigned int i0=0; i0<n; ++i0, ++c)
            t[c] /= (eigenvalues[2][i2] + eigenvalues[1][i1] + eigenvalues[0][i0]);
      eval.template apply<0, false, false> (S0, t, dst);
      eval.template apply<1, false, false> (S1, dst, t);
      eval.template apply<2, false, false> (S2, t, dst);
    }

  else
    Assert(false, ExcNotImplemented());
}


// ------------------------------   TensorProductMatrixSymmetricSum   ------------------------------

template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,Number,size>
::TensorProductMatrixSymmetricSum ()
  : TensorProductMatrixSymmetricSumBase<dim,Number,size>()
{}



template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,Number,size>
::TensorProductMatrixSymmetricSum (const std::array<Table<2,Number>, dim> &mass_matrix,
                                   const std::array<Table<2,Number>, dim> &derivative_matrix)
{
  reinit (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,Number,size>
::TensorProductMatrixSymmetricSum(const std::array<FullMatrix<Number>, dim> &mass_matrix,
                                  const std::array<FullMatrix<Number>, dim> &derivative_matrix)
{
  reinit (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,Number,size>
::TensorProductMatrixSymmetricSum (const FullMatrix<Number> &mass_matrix,
                                   const FullMatrix<Number> &derivative_matrix)
{
  reinit (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
template <typename MatrixArray>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::reinit_impl (MatrixArray &&mass_matrices_,
               MatrixArray &&derivative_matrices_)
{
  auto &&mass_matrices = std::forward<MatrixArray>(mass_matrices_) ;
  auto &&derivative_matrices = std::forward<MatrixArray>(derivative_matrices_) ;
  this->mass_matrix = mass_matrices ;
  this->derivative_matrix = derivative_matrices ;

  for (int dir = 0; dir < dim; ++dir)
    {
      Assert (size == -1 || (size > 0 && static_cast<unsigned int>(size) == mass_matrices[dir].n_rows()),
              ExcDimensionMismatch(size, mass_matrices[dir].n_rows()));
      AssertDimension (mass_matrices[dir].n_rows(), mass_matrices[dir].n_cols());
      AssertDimension (mass_matrices[dir].n_rows(), derivative_matrices[dir].n_rows());
      AssertDimension (mass_matrices[dir].n_rows(), derivative_matrices[dir].n_cols());

      this->eigenvectors[dir].reinit (mass_matrices[dir].n_cols(), mass_matrices[dir].n_rows()) ;
      this->eigenvalues[dir].resize (mass_matrices[dir].n_cols()) ;
      spectral_assembly<Number> (&(mass_matrices[dir](0,0))
                                 , &(derivative_matrices[dir](0,0))
                                 , mass_matrices[dir].n_rows()
                                 , mass_matrices[dir].n_cols()
                                 , this->eigenvalues[dir].begin()
                                 , &(this->eigenvectors[dir](0,0))) ;
    }
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::reinit (const std::array<Table<2,Number>, dim> &mass_matrix,
          const std::array<Table<2,Number>, dim> &derivative_matrix)
{
  reinit_impl (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::reinit (const std::array<FullMatrix<Number>, dim> &mass_matrix,
          const std::array<FullMatrix<Number>, dim> &derivative_matrix)
{
  std::array<Table<2,Number>,dim> mass_copy ;
  std::array<Table<2,Number>,dim> deriv_copy ;

  std::transform (mass_matrix.cbegin(), mass_matrix.cend(), mass_copy.begin(),
                  [] (const FullMatrix<Number> &m) ->Table<2,Number> {return m;}) ;
  std::transform (derivative_matrix.cbegin(), derivative_matrix.cend(), deriv_copy.begin(),
                  [] (const FullMatrix<Number> &m) ->Table<2,Number> {return m;}) ;

  reinit_impl (std::move(mass_copy), std::move(deriv_copy)) ;
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::reinit (const FullMatrix<Number> &mass_matrix,
          const FullMatrix<Number> &derivative_matrix)
{
  std::array<Table<2,Number>,dim> mass_matrices ;
  std::array<Table<2,Number>,dim> derivative_matrices ;

  std::fill (mass_matrices.begin(), mass_matrices.end(), mass_matrix) ;
  std::fill (derivative_matrices.begin(), derivative_matrices.end(), derivative_matrix) ;

  reinit_impl (std::move(mass_matrices), std::move(derivative_matrices)) ;
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::vmult (Vector<Number> &dst,
         const Vector<Number> &src) const
{
  AssertDimension(dst.size(), this->m()) ;
  AssertDimension(src.size(), this->n()) ;
  TensorProductMatrixSymmetricSumBase<dim,Number,size>::vmult (dst.begin(), src.begin());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,Number,size>
::apply_inverse (Vector<Number> &dst,
                 const Vector<Number> &src) const
{
  AssertDimension (dst.size(), this->n()) ;
  AssertDimension (src.size(), this->m()) ;
  TensorProductMatrixSymmetricSumBase<dim,Number,size>::apply_inverse (dst.begin(), src.begin());
}



// ------------------------------ vectorized spec.: TensorProductMatrixSymmetricSum   ------------------------------

template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::TensorProductMatrixSymmetricSum ()
  : TensorProductMatrixSymmetricSumBase<dim,VectorizedArray<Number>,size>()
{}



template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::TensorProductMatrixSymmetricSum (const std::array<Table<2,VectorizedArray<Number> >,dim> &mass_matrix,
                                   const std::array<Table<2,VectorizedArray<Number> >,dim> &derivative_matrix)
{
  reinit (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
inline
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::TensorProductMatrixSymmetricSum (const Table<2,VectorizedArray<Number> > &mass_matrix,
                                   const Table<2,VectorizedArray<Number> > &derivative_matrix)
{
  reinit (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
template <typename MatrixArray>
inline
void
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::reinit_impl (MatrixArray &&mass_matrices_,
               MatrixArray &&derivative_matrices_)
{
  auto &&mass_matrix = std::forward<MatrixArray>(mass_matrices_) ;
  auto &&derivative_matrix = std::forward<MatrixArray>(derivative_matrices_) ;
  this->mass_matrix = mass_matrix ;
  this->derivative_matrix = derivative_matrix ;

  constexpr unsigned int macro_size = VectorizedArray<Number>::n_array_elements ;
  const unsigned int nm_flat_size
    = (size > 0)
      ? (Utilities::fixed_int_power<size,dim>::value
         * Utilities::fixed_int_power<size,dim>::value * macro_size)
      : (Utilities::fixed_power<dim>(mass_matrix[0].n_rows())
         * Utilities::fixed_power<dim>(mass_matrix[0].n_rows()) * macro_size) ;
  const unsigned int n_flat_size
    = (size > 0)
      ? Utilities::fixed_int_power<size,dim>::value * macro_size
      : Utilities::fixed_power<dim>(mass_matrix[0].n_rows()) * macro_size ;

  std::vector<Number> mass_matrix_flat ;
  std::vector<Number> deriv_matrix_flat ;
  std::vector<Number> eigenvalues_flat ;
  std::vector<Number> eigenvectors_flat ;
  mass_matrix_flat.reserve (nm_flat_size) ;
  deriv_matrix_flat.reserve (nm_flat_size) ;
  eigenvalues_flat.reserve (n_flat_size) ;
  eigenvectors_flat.reserve (nm_flat_size) ;
  std::array<unsigned int,macro_size> offsets_nm ;
  std::array<unsigned int,macro_size> offsets_n ;
  for (int dir = 0; dir < dim; ++dir)
    {
      Assert (size == -1 ||
              (size > 0 && static_cast<unsigned int>(size) == mass_matrix[dir].n_rows()),
              ExcDimensionMismatch(size, mass_matrix[dir].n_rows()));
      AssertDimension (mass_matrix[dir].n_rows(), mass_matrix[dir].n_cols());
      AssertDimension (mass_matrix[dir].n_rows(), derivative_matrix[dir].n_rows());
      AssertDimension (mass_matrix[dir].n_rows(), derivative_matrix[dir].n_cols());

      const unsigned int n_rows = mass_matrix[dir].n_rows() ;
      const unsigned int n_cols = mass_matrix[dir].n_cols() ;
      const unsigned int nm = n_rows * n_cols ;

      mass_matrix_flat.resize (macro_size*nm) ;
      deriv_matrix_flat.resize (macro_size*nm) ;
      eigenvalues_flat.resize (macro_size*n_rows) ;
      eigenvectors_flat.resize (macro_size*nm) ;
      for (unsigned int vv=0; vv<macro_size; ++vv)
        offsets_nm[vv] = nm * vv ;
      for (unsigned int vv=0; vv<macro_size; ++vv)
        offsets_n[vv] = n_rows * vv ;

      vectorized_transpose_and_store (false, nm, &(mass_matrix[dir](0,0))
                                      , offsets_nm.cbegin(), mass_matrix_flat.data()) ;
      vectorized_transpose_and_store (false, nm, &(derivative_matrix[dir](0,0))
                                      , offsets_nm.cbegin(), deriv_matrix_flat.data()) ;

      const Number *mass_cbegin = mass_matrix_flat.data() ;
      const Number *deriv_cbegin = deriv_matrix_flat.data() ;
      Number *eigenvec_begin = eigenvectors_flat.data() ;
      Number *eigenval_begin = eigenvalues_flat.data() ;

      spectral_assembly<Number> (mass_cbegin, deriv_cbegin, n_rows, n_cols
                                 , eigenval_begin, eigenvec_begin) ;
      for (unsigned int lane = 1; lane < macro_size; ++lane)
        {
          std::advance (mass_cbegin, nm) ;
          std::advance (deriv_cbegin, nm) ;
          std::advance (eigenvec_begin, nm) ;
          std::advance (eigenval_begin, n_rows) ;
          spectral_assembly<Number> (mass_cbegin, deriv_cbegin, n_rows, n_cols
                                     , eigenval_begin, eigenvec_begin) ;
        }

      this->eigenvalues[dir].resize (n_rows) ;
      this->eigenvectors[dir].reinit (n_rows, n_cols) ;
      vectorized_load_and_transpose (n_rows, eigenvalues_flat.data()
                                     , offsets_n.cbegin(), this->eigenvalues[dir].begin()) ;
      vectorized_load_and_transpose (nm, eigenvectors_flat.data()
                                     , offsets_nm.cbegin(), &(this->eigenvectors[dir](0,0))) ;
    }
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::reinit (const std::array<Table<2,VectorizedArray<Number> >,dim> &mass_matrix,
          const std::array<Table<2,VectorizedArray<Number> >,dim> &derivative_matrix)
{
  reinit_impl (mass_matrix, derivative_matrix) ;
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::reinit (const Table<2,VectorizedArray<Number> > &mass_matrix,
          const Table<2,VectorizedArray<Number> > &derivative_matrix)
{
  std::array<Table<2,VectorizedArray<Number> >,dim> mass_matrices ;
  std::array<Table<2,VectorizedArray<Number> >,dim> derivative_matrices ;

  std::fill (mass_matrices.begin(), mass_matrices.end(), mass_matrix) ;
  std::fill (derivative_matrices.begin(), derivative_matrices.end(), derivative_matrix) ;

  reinit_impl (std::move(mass_matrices), std::move(derivative_matrices)) ;
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::vmult (AlignedVector<VectorizedArray<Number> > &dst,
         const AlignedVector<VectorizedArray<Number> > &src) const
{
  AssertDimension(dst.size(), this->m()) ;
  AssertDimension(src.size(), this->n()) ;
  TensorProductMatrixSymmetricSumBase<dim,VectorizedArray<Number>,size>::vmult (dst.begin(), src.begin());
}



template <int dim, typename Number, int size>
inline
void
TensorProductMatrixSymmetricSum<dim,VectorizedArray<Number>,size>
::apply_inverse (AlignedVector<VectorizedArray<Number> > &dst,
                 const AlignedVector<VectorizedArray<Number> > &src) const
{
  AssertDimension (dst.size(), this->n()) ;
  AssertDimension (src.size(), this->m()) ;
  TensorProductMatrixSymmetricSumBase<dim,VectorizedArray<Number>,size>::apply_inverse (dst.begin(), src.begin());
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif
